
"""Create a summary that shows the assignment of 3' end sites to annotation features"""

__date__ = "2017-02-14"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import os
import time
from argparse import ArgumentParser, RawTextHelpFormatter
import gzip
import re
import bisect

out = open("out.table","w")

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("-g",
                    "--gtf",
                    dest="gtf",
                    help="gene annotation file in uncompressed gtf format")

parser.add_argument("--b",
                    "--bed",
                    dest="bed_file",
                    help="BED-file of the poly(A) sites (gzipped or uncompressed; raw 3' ends or clustered poly(A) sites")

parser.add_argument("-u",
                    "--utr_name",
                    dest="utr_name",
                    help="file name providing the poly(A) sites in extended BED ")


syserr = sys.stderr.write
sysout = sys.stdout.write

def insert_feature( coordinates_list, sets_list, start, end, term):

    # end is in BED-format -> this means, that it is actually the first
    # pos of the next region

    # set initial values for both lists if necessary
    if len(coordinates_list) == 0:
        coordinates_list.append(0)
        sets_list.append(set())

    inserted = False

    insert_index = bisect.bisect(coordinates_list, start)
    if term not in sets_list[ insert_index - 1 ]:
        if coordinates_list[ insert_index - 1] == start:
            # add term to set of features
            sets_list[ insert_index - 1].add(term)
            # special case because start and coordinate are the same
            insert_index -= 1
            inserted = True
        else:
            # start is bigger than coordinate before
            # create a new step
            coordinates_list.insert(insert_index, start)
            sets_list.insert(insert_index, {term} | sets_list[insert_index - 1] )
            inserted = True

    else:
        # term is already present in previous entry
        # -> no insertion needed
        # but next entry upstream must not be skipped
        insert_index -= 1

    # start-pos handling finished

    # now, iterate over all coordinates that are smaller than the end
    duplicates = set()
    end_index = bisect.bisect( coordinates_list, end )
    for idx in range(insert_index + 1, end_index):
        if idx == end_index - 1:
            # last iteration
            # treat potential special case that end and already
            # present coordinate are equal
            if coordinates_list[ idx ] == end:
                # return both lists
                return( coordinates_list, sets_list )
            
        if term in sets_list[ idx ]:
            inserted = False
        else:
            sets_list[ idx ].add(term)
            inserted = True
            if sets_list[idx] == sets_list[idx-1]:
                duplicates.add( idx )
            if (idx+1) < len(sets_list):
                if sets_list[idx] == sets_list[ idx+1 ]:
                    duplicates.add(idx+1)

    # handle potential duplicate entries
    for idx in sorted(list(duplicates))[::-1]:
        del coordinates_list[idx]
        del sets_list[idx]

    end_index = bisect.bisect( coordinates_list, end )

    # handle current region end
    # insert entry to coordinates_list/sets_list if necessary
    if inserted:
        # term was inserted for the last step
        # insert new entry 
        coordinates_list.insert(end_index, end)
        sets_list.insert(end_index, sets_list[end_index - 1].difference({term}) )

    # # test -->
    # syserr("[TMP] Intermediate: coords: %s\n" % str(coordinates_list))
    # syserr("[TMP] Intermediate: sets: %s\n\n" % str(sets_list))
    # # <-- test
            
    return( coordinates_list, sets_list)

def insert_feature_curr_transcript( coordinates_list, value_list, start, end, term ):
    '''same function as above except for the fact that one region/step
    is only allowed to take one value with the 3' utr having the highest prio
    '''

    # initialize lists if necessary
    if len(coordinates_list) == 0:
        return( [start,end], [term,''])

    insert_index = bisect.bisect(coordinates_list, start)

    if term != "utr":
        # now, only insertions at the beginning and at the end
        # of the coordinates_list are allowed (due to the fact,
        # that a transcript is processed based on sorted exons)
        if insert_index != 0 and coordinates_list[insert_index - 1] != start:
            syserr("[ERROR] The current region is not inserted properly in the transcript annotation\n")
            syserr("[INFO] start: %i, end: %i, term: %s\n" %(start,end, term))
            syserr("[INFO] coods so far: %s\n" % str(coordinates_list))
            sys.exit(2)
        if insert_index == 0:
            coordinates_list.insert(0, start)
            value_list.insert(0,term)
            return( coordinates_list, value_list)
        else:
            if value_list[-1] != "":
                # this is a consistency check only
                syserr("[ERROR] transcript map was not build properly! The last feature does not end\n")
                syserr("[INFO] coords: %s, vals: %s\n" % (str(coordinates_list), str(value_list)))
                sys.exit(2)
            value_list[-1] = term
            value_list.insert(insert_index,'')
            coordinates_list.append(end)
            return( coordinates_list, value_list)

    else:
        # handle utr entries
        # which might occur anywhere

        inserted = False
        prev_entry = None
        to_delete = []

        # first part: insert start of utr-region
        if insert_index == 0:
            coordinates_list.insert(0, start)
            value_list.insert(0, term)
        else:
            # utr region starts within already defined regions
            if coordinates_list[ insert_index - 1 ] == start:
                if value_list[ insert_index -1 ] != term:
                    prev_entry = value_list[ insert_index -1 ]
                    value_list[ insert_index - 1 ] = term
                    insert_index -= 1
                    inserted = True
            else:
                if value_list[ insert_index -1 ] != term:
                    coordinates_list.insert(insert_index, start)
                    value_list.insert(insert_index, term)
                    inserted = True
                    prev_entry = value_list[insert_index - 1]
                
        # now, check for coordinate steps smaller than the end
        end_index = bisect.bisect( coordinates_list, end )
        
        for idx in range(insert_index + 1, end_index):
            if idx == end_index - 1:
                # last iteration
                # treat potential special case that end and already
                # present coordinate are equal
                if coordinates_list[ idx ] == end:
                    # return both lists
                    return( coordinates_list, value_list )
            
                if term == value_list[ idx ]:
                    inserted = False
                else:
                    if value_list[ idx ] == '':
                        if not inserted:
                            # in the step before, nothing was inserted
                            # hence, it was utr-annotated already
                            # since it stays utr, this border can be deleted afterwards
                            to_delete.append(idx)
                        else:
                            # entering this case would mean that the utr region
                            # spans multiple features
                            # this is not possible
                            syserr("[ERROR] utr spans different regions\n")
                            syserr("[ERROR] first region: %i - %i, %s\n"
                                   % (coordinates_list[ idx - 1], coordinates_list[idx], value_list[ idx -1]))
                            syserr("[ERROR] next region: %i - %i, %s\n"
                                   % (coordinates_list[ idx], coordinates_list[idx + 1], value_list[ idx]))
                            syserr("[ERROR] UTR: %i - %i, %s\n"
                                   % (start, end, term))
                            sys.exit(2)
                        
                    inserted = True

        # handle end of region
        # last, potentially insert a new end part
        if inserted:
            if len(to_delete) > 0:
                # utr regions extends previous annotation
                syserr("[Warning] utr region extends previous region: %i - %i\n" % (coordinates_list[to_delete[0]-1], coordinates_list[to_delete[0]]))
                sys.exit(2)
            # term was inserted in the last step
            # insert new entry
            coordinates_list.insert(end_index, end)
            value_list.insert(end_index, prev_entry )
        
        return( coordinates_list, value_list)        

def read_gtf(gtf, three_utr_name):
    '''iterate through gtf file and return a dict that contains for every 
    genomic position the feature at this position
    '''

    # The following features are used:
    # exon
    # intron
    # terminal exon
    # 3' UTR
    # intergenic

    feature_dict = {}
    previous_exon = None
    tmp_transcript_dict = {}
    utr_entries = {}
    gtf_order = None

    with open(gtf, 'r') as in_gtf:
        for line in in_gtf:
            if line.startswith("#"):
                continue
            F = line.rstrip().split("\t")

            if F[2] == "exon":
                # get transcript id
                mo = re.match('.+transcript_id\s\"([^\"]+)', F[8])
                assert mo
                tr_id = mo.groups()[0]
                chrom = F[0]
                strand = F[6]
                # use chr:strand as key
                key = chrom + ":" + strand
                # change coordinates to BED format
                start = int(F[3]) - 1
                end = int(F[4])

                if key not in feature_dict:
                    feature_dict[ key ] = {"coords": [], "sets": []}

                if key not in tmp_transcript_dict:
                    tmp_transcript_dict[ key ] = {"coords": [], "vals": []}
                
                if previous_exon is not None:
                    # process previous exon

                    # define the ordering of the gtf file if not yet done (three_to_five or increasing)
                    if gtf_order is None:
                        if tr_id == previous_exon[0] and strand == "-":
                            if start > previous_exon[3]:
                                gtf_order = "increasing"
                            else:
                                gtf_order = "three_to_five"

                            # this is the first transcript with at least two exons
                            # check if the previous exon needs to be set as terminal
                            if gtf_order == "increasing":
                                (tmp_transcript_dict[key]["coords"],
                                 tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                     tmp_transcript_dict[key]["vals"],
                                                                                                     previous_exon[3],
                                                                                                     previous_exon[4],
                                                                                                     "term_exon"
                                 )
                                # append a note that the exon was processed already
                                previous_exon.append(1)

                    prev_tr_id  = previous_exon[0]
                    prev_chr    = previous_exon[1]
                    prev_strand = previous_exon[2]
                    prev_start  = previous_exon[3]
                    prev_end    = previous_exon[4]
                    prev_key    = previous_exon[5]
                    # for the case of increasing ordered gtfs
                    # it is saved whether the last exon was the first
                    # for the corresponding transcript
                    if len(previous_exon) == 7:
                        first_exon = True
                    else:
                        first_exon = False
                        
                    if prev_tr_id == tr_id:
                        # same transcript as previous exon
                        if strand == "+":
                            # store previous exon as normal exon
                            (tmp_transcript_dict[key]["coords"],
                             tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                 tmp_transcript_dict[key]["vals"],
                                                                                                 prev_start,
                                                                                                 prev_end,
                                                                                                 "exon"
                             )
                            # store the space inbetween as intron
                            (tmp_transcript_dict[key]["coords"],
                             tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                 tmp_transcript_dict[key]["vals"],
                                                                                                 prev_end,
                                                                                                 start,
                                                                                                "intron"
                            )

                        else:
                            # on negative strand it's important to respect
                            # the gtf-order
                            if gtf_order is "three_to_five":
                                # similar to plus-strand case
                                # only store previous exon and intron
                                (tmp_transcript_dict[key]["coords"],
                                 tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                     tmp_transcript_dict[key]["vals"],
                                                                                                     prev_start,
                                                                                                     prev_end,
                                                                                                     "exon"
                                )
                                # store the space inbetween as intron
                                (tmp_transcript_dict[key]["coords"],
                                 tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                     tmp_transcript_dict[key]["vals"],
                                                                                                     end,
                                                                                                     prev_start,
                                                                                                     "intron"
                                 )

                            elif gtf_order is "increasing":
                                # store the previous exon if it was not
                                # the first of the transcript;
                                # anyway, the intron region
                                if not first_exon:
                                    (tmp_transcript_dict[key]["coords"],
                                     tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                         tmp_transcript_dict[key]["vals"],
                                                                                                         prev_start,
                                                                                                         prev_end,
                                                                                                         "exon"
                                     )
                                # store the space inbetween as intron
                                (tmp_transcript_dict[key]["coords"],
                                 tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                     tmp_transcript_dict[key]["vals"],
                                                                                                     prev_end,
                                                                                                     start,
                                                                                                     "intron"
                                 )

                        # store information for new exon
                        previous_exon = [tr_id, chrom, strand, start, end, key]
            
                    else:
                        # previous and current exon belong to different transcript
                        
                        # depending on the strand of the previous transcript
                        # and the order of the gtf the annotation changes
                        if (prev_strand == "+" or
                            ( prev_strand == "-" and gtf_order == "three_to_five" ) or
                            gtf_order is None):
                            # test -->
                            if prev_key not in tmp_transcript_dict:
                                syserr("[TMP] %s not found in hash for entry: %s\n" % (prev_key, previous_exon))
                            # annotate the previous exon as terminal exon
                            (tmp_transcript_dict[prev_key]["coords"],
                             tmp_transcript_dict[prev_key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[prev_key]["coords"],
                                                                                                      tmp_transcript_dict[prev_key]["vals"],
                                                                                                      prev_start,
                                                                                                      prev_end,
                                                                                                      "term_exon"
                             )
                                
                        elif prev_strand == "-" and gtf_order == "increasing":
                            # store previous exon as normal exon
                            (tmp_transcript_dict[prev_key]["coords"],
                             tmp_transcript_dict[prev_key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[prev_key]["coords"],
                                                                                                      tmp_transcript_dict[prev_key]["vals"],
                                                                                                      prev_start,
                                                                                                      prev_end,
                                                                                                      "exon"
                             )
                        else:
                            # consistency check:
                            # this part should never be accessed
                            syserr("[ERROR] Previous exon was not properly recognized and annotated\n")
                            syserr("[ERROR] Previous exon: %s\n" % "\t".join(previous_exon))
                            sys.exit(2)

                        # previous transcript has been processed
                        if len(tmp_transcript_dict) <= 2:
                            # include utr annotation to tmp_transcript_dict
                            if len(utr_entries) > 0:
                                for utr_entry in utr_entries[ prev_tr_id ]:
                                    u_start = utr_entry[0]
                                    u_end = utr_entry[1]
                                    (tmp_transcript_dict[prev_key]["coords"],
                                     tmp_transcript_dict[prev_key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[prev_key]["coords"],
                                                                                                              tmp_transcript_dict[prev_key]["vals"],
                                                                                                              u_start,
                                                                                                              u_end,
                                                                                                              "utr"
                                     )

                                # clean utr_entries for next transcript
                                utr_entries = {}


                            for idx in range(len(tmp_transcript_dict[prev_key]["coords"])-1):
                                (feature_dict[ prev_key ]["coords"],
                                 feature_dict[ prev_key ]["sets"]) = insert_feature( feature_dict[ prev_key ]["coords"],
                                                                                     feature_dict[ prev_key ]["sets"],
                                                                                     tmp_transcript_dict[prev_key]["coords"][idx],
                                                                                     tmp_transcript_dict[prev_key]["coords"][idx+1],
                                                                                     tmp_transcript_dict[prev_key]["vals"][idx]
                                 )

                            # # test -->
                            # syserr("[TMP] coords: %s\n" % str(feature_dict[prev_key]["coords"]))
                            # syserr("[TMP] sets: %s\n" % str(feature_dict[prev_key]["sets"]))
                            # # <-- test

                        else:
                            # consistency check
                            syserr("[ERROR] intermediate transcript information belongs to unambiguous number of keys (%i)\n"
                                   % len(tmp_transcript_dict))
                            sys.exit(2)

                        tmp_transcript_dict = {}
                        tmp_transcript_dict[ key ] = {"coords":[], "vals": []}

                        # store information for new exon
                        # count it as terminal in the case of increased gtf order
                        # and negative strand
                        if strand == "-" and gtf_order == "increasing":
                            # count this first exon as terminal exon
                            # and mark that it has been considered already
                            (tmp_transcript_dict[key]["coords"],
                             tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                                      tmp_transcript_dict[key]["vals"],
                                                                                                      start,
                                                                                                      end,
                                                                                                      "term_exon"
                             )
                            previous_exon = [tr_id, chrom, strand, start, end, key, 1]
                        else:
                            # only save the current exon
                            previous_exon = [tr_id, chrom, strand, start, end, key]

                else:
                    # no previous exon
                    # store the current information:
                    # [transcript_id, chr, strand, start, end, key(i.e. "chrom:strand")]
                    previous_exon = [tr_id, chrom, strand, start, end, key]

                
            #######################################################
            #######################################################

            if F[2] == three_utr_name:
                                    
                # 3' UTR has highest priority

                # get transcript id
                mo = re.match('.+transcript_id\s\"([^\"]+)', F[8])
                assert mo
                tr_id = mo.groups()[0]
                chrom = F[0]
                strand = F[6]
                # use chr:strand as key
                key = chrom + ":" + strand
                # change coordinates to BED format
                start = int(F[3]) - 1
                end = int(F[4])

                if key not in feature_dict:
                    feature_dict[ key ] = {"coords": [], "sets": []}

                if key not in tmp_transcript_dict:
                    tmp_transcript_dict[ key ] = {}

                if previous_exon is not None:

                    prev_tr_id  = previous_exon[0]
                    prev_chr    = previous_exon[1]
                    prev_strand = previous_exon[2]
                    prev_start  = previous_exon[3]
                    prev_end    = previous_exon[4]
                    prev_key    = previous_exon[5]
                    
                    # compare previous and current transcript
                    if prev_tr_id == tr_id:
                        # 3' UTR belongs to the current transcript
                        # store utr infos

                        if tr_id not in utr_entries:
                            utr_entries[tr_id] = []
                        utr_entries[tr_id].append([ start,end])

                    else :
                        # current utr section marks a new transcript
                        
                        # all information from the previous transcript was considered
                        if len(tmp_transcript_dict) > 0:
                            # test -->
                            syserr("[TMP] utr marks new transcript\n")
                            # <-- test

                            if len(tmp_transcript_dict) <= 2:
                                # include utr annotation to tmp_transcript_dict
                                if len( utr_entries ) > 0:
                                    for utr_entry in utr_entries[ prev_tr_id ]:
                                        u_start = utr_entry[0]
                                        u_end = utr_entry[1]
                                        (tmp_transcript_dict[prev_key]["coords"],
                                         tmp_transcript_dict[prev_key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[prev_key]["coords"],
                                                                                                                  tmp_transcript_dict[prev_key]["vals"],
                                                                                                                  u_start,
                                                                                                                  u_end,
                                                                                                                  "utr"
                                         )

                                    # clean utr_entries for next transcript
                                    utr_entries = {}

                                for idx in range(len(tmp_transcript_dict[prev_key]["coords"])-1):
                                    (feature_dict[ prev_key ]["coords"],
                                     feature_dict[ prev_key ]["sets"]) = insert_feature( feature_dict[ prev_key ]["coords"],
                                                                                         feature_dict[ prev_key ]["sets"],
                                                                                         tmp_transcript_dict[prev_key]["coords"][idx],
                                                                                         tmp_transcript_dict[prev_key]["coords"][idx+1],
                                                                                         tmp_transcript_dict[prev_key]["vals"][idx]
                                     )

                            else:
                                # consistency check
                                syserr("[ERROR] intermediate transcript information belongs to unambiguous number of keys (%i)\n"
                                       % len( tmp_transcript_dict))
                                sys.exit(2)

                            tmp_transcript_dict = {}
                            tmp_transcript_dict[ key ] = {"coords":[], "vals": []}
                            
                            if tr_id not in utr_entries:
                                utr_entries[tr_id] = []
                            utr_entries[tr_id].append([ start,end])

                else:
                    # no previous exon is available
                    # simply store the utr information
                    if tr_id not in utr_entries:
                            utr_entries[tr_id] = []
                    utr_entries[tr_id].append([ start,end])

    ###
    # process last exon
    ###
    if previous_exon is not None:
        prev_tr_id  = previous_exon[0]
        prev_chr    = previous_exon[1]
        prev_strand = previous_exon[2]
        prev_start  = previous_exon[3]
        prev_end    = previous_exon[4]
        prev_key    = previous_exon[5]
        # for the case of increasing ordered gtfs
        if len(previous_exon) == 7:
            first_exon = True
        else:
            first_exon = False
        if (strand == "+" or
            ( strand == "-" and gtf_order == "three_to_five" ) or
            gtf_order is None):
            # annotate the exon as terminal exon
            (tmp_transcript_dict[key]["coords"],
             tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                 tmp_transcript_dict[key]["vals"],
                                                                                 start,
                                                                                 end,
                                                                                 "term_exon"
             )
                                
        elif strand == "-" and gtf_order == "increasing":
            # store exon as normal exon
            (tmp_transcript_dict[key]["coords"],
             tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                 tmp_transcript_dict[key]["vals"],
                                                                                 start,
                                                                                 end,
                                                                                 "exon"
             )

    ###
    # process last transcript
    ###
    if len( utr_entries ) > 0:
        for utr_entry in utr_entries[ tr_id ]:
            u_start = utr_entry[0]
            u_end = utr_entry[1]
            (tmp_transcript_dict[key]["coords"],
             tmp_transcript_dict[key]["vals"]) = insert_feature_curr_transcript( tmp_transcript_dict[key]["coords"],
                                                                                      tmp_transcript_dict[key]["vals"],
                                                                                      u_start,
                                                                                      u_end,
                                                                                      "utr"
             )

    for idx in range(len(tmp_transcript_dict[key]["coords"])-1):
        (feature_dict[ key ]["coords"],
         feature_dict[ key ]["sets"]) = insert_feature( feature_dict[ key ]["coords"],
                                                        feature_dict[ key ]["sets"],
                                                        tmp_transcript_dict[key]["coords"][idx],
                                                        tmp_transcript_dict[key]["coords"][idx+1],
                                                        tmp_transcript_dict[key]["vals"][idx]
         )

    # # test -->
    # syserr("[TMP] coords: %s\n" % str(feature_dict[prev_key]["coords"]))
    # syserr("[TMP] sets: %s\n" % str(feature_dict[prev_key]["sets"]))
    # # <-- test
    
    return(feature_dict)

def get_overlapping_features( feat_dict, start, end):

    coords_list = feat_dict[ 'coords' ]
    feature_sets_list = feat_dict[ 'sets' ]

    return_set = set()

    start_idx = bisect.bisect( coords_list, start) - 1
    end_idx = bisect.bisect( coords_list, end)

    if start_idx == end_idx:
        end_idx += 1

    for region_idx in range( start_idx, end_idx ):
        if region_idx == end_idx - 1:
            # only count this features if
            # the corresponding coordinate is << then
            # the end coordinate (BED format: end coord is not included in the feature)
            if end ==  coords_list[ region_idx ]:
                continue
        return_set = return_set.union( feature_sets_list[ region_idx ] )

    return( return_set )


def main(options):

    # read the gtf_file
    # gives back region annotated as:
    # intron
    # exon
    # terminal exon
    # 3' UTR
    feature_dict = read_gtf( options.gtf, options.utr_name )

    # total number of sites
    total_cnt = 0
    # dict to count the overlap of sites with speicifc features types
    feature_cnt = {"intergenic": 0}

    # read input file
    if options.bed_file.endswith(".gz"):
        input_bed = gzip.open( options.bed_file, "rt")
    else:
        input_bed = open( options.bed_file, "r")
    for line in input_bed:
        if line.startswith("#"):
            continue
        total_cnt += 1
        line_list = line.rstrip().split("\t")
        curr_key = line_list[0] + ":" + line_list[5]
        id_list = line_list[3].split(":")
        curr_start = int(line_list[1])
        curr_end = int(line_list[2])

        if curr_key not in feature_dict:
            syserr("[INFO] No feature found for %s in gtf. Site %s is counted as intergenic\n"
                   % (curr_key, line_list[0] + ":" + line_list[1] + ":" + line_list[2] + ":" + line_list[5]) )
            feature_cnt["intergenic"] += 1
            out.write(curr_key + "\t" + "intergenic")
	    continue

        feature_set = get_overlapping_features( feature_dict[ curr_key ], curr_start, curr_end)
        
        # # test -->
        # if "utr" in feature_set:
        #     syserr("[tmp] %s overlaps with 3UTR\n" % line_list[3])
        # # <-- test
        
        number_of_features = len(feature_set)
        if (number_of_features == 0 or
            ( number_of_features == 1 and "" in feature_set)):
            feature_cnt["intergenic"] += 1
	    out.write(line_list[0] + ":" + line_list[1] + "-" + line_list[2] + "\t" + "intergenic" + "\n")
        else:
            for feat in feature_set:
                if feat == "":
                    continue
                if feat not in feature_cnt:
                    feature_cnt[ feat ] = 0
                feature_cnt[ feat ] += 1
		out.write(line_list[0] + ":" + line_list[1] + "-" + line_list[2] + "\t" + feat + "\n")
    input_bed.close()
    
    # output overview table
    sysout("feature\t#_sites\tfraction_of_total\n")
    for feat in sorted(feature_cnt):
        sysout("%s\t%.2f\t%.2f\n"
               %( feat, float( feature_cnt[ feat ] ), float( feature_cnt[ feat ] ) / total_cnt ) )

if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception as e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)

        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)

out.close()
