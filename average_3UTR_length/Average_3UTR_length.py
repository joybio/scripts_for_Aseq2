#!/bin/miniconda3/bin/python

"""Create a summary of average 3'UTR length in Control and Treat"""

__date__ = "2019-10-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

#imports
import re
import os
import optparse
#import bisect
#sort  package
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-c','--Control',
                dest='control',
                help='enrichment peak file in bed format')
parser.add_option('-t','--Treat',
                dest='treat',
                help='enrichment peak file in bed format')


(options,args) = parser.parse_args()

with open("distance_control.csv","w") as out:
        with open(options.control,'r') as data:
                for i in data:
                        i = i.strip().split('\t')
                        n=3
			stop_codon = int(i[2])
			TPM_rep1 = 0
			TPM_rep2 = 0
			UTR3_rep1 = 0
			UTR3_rep2 = 0
			if i[1] == '+':
				while n < len(i):
					a = i[n].split(',')
					distance_TPM_rep1 = (int(a[0]) - stop_codon) * float(a[1])
					distance_TPM_rep2 = (int(a[0]) - stop_codon) * float(a[2])
					TPM_rep1 += float(a[1])
					TPM_rep2 += float(a[2])
					UTR3_rep1 += distance_TPM_rep1
					UTR3_rep2 += distance_TPM_rep2
					n += 1
			else:
				while n < len(i):
					a = i[n].split(',')
					distance_TPM_rep1 = (stop_codon - int(a[0])) * float(a[1])
					distance_TPM_rep2 = (stop_codon - int(a[0])) * float(a[2])
					TPM_rep1 += float(a[1])
					TPM_rep2 += float(a[2])
					UTR3_rep1 += distance_TPM_rep1
					UTR3_rep2 += distance_TPM_rep2
					n += 1
			UTR3_average_rep1 = int(UTR3_rep1/TPM_rep1)
			UTR3_average_rep2 = int(UTR3_rep2/TPM_rep2)
			out.write(i[0] + '\t' + i[1] + '\t' + str(UTR3_average_rep1) + '\t' + str(UTR3_average_rep2) + '\n')
	data.close()
out.close()

with open("distance_treat.csv","w") as out:
        with open(options.treat,'r') as data:
                for i in data:
                        i = i.strip().split('\t')
                        n=3     
                        stop_codon = int(i[2])
			TPM_rep1 = 0
                        TPM_rep2 = 0
                        UTR3_rep1 = 0
                        UTR3_rep2 = 0
                        if i[1] == '+':
                                while n < len(i):
                                        a = i[n].split(',')
                                        distance_TPM_rep1 = (int(a[0]) - stop_codon) * float(a[1])
                                        distance_TPM_rep2 = (int(a[0]) - stop_codon) * float(a[2])
                                        TPM_rep1 += float(a[1])
                                        TPM_rep2 += float(a[2])
                                        UTR3_rep1 += distance_TPM_rep1
                                        UTR3_rep2 += distance_TPM_rep2
                                        n += 1
                        else:
                                while n < len(i):
                                        a = i[n].split(',')
                                        distance_TPM_rep1 = (stop_codon - int(a[0])) * float(a[1])
                                        distance_TPM_rep2 = (stop_codon - int(a[0])) * float(a[2])
                                        TPM_rep1 += float(a[1])
                                        TPM_rep2 += float(a[2])
                                        UTR3_rep1 += distance_TPM_rep1
                                        UTR3_rep2 += distance_TPM_rep2
                                        n += 1
                        UTR3_average_rep1 = int(UTR3_rep1/TPM_rep1)
                        UTR3_average_rep2 = int(UTR3_rep2/TPM_rep2)
                        out.write(i[0] + '\t' + i[1] + '\t' + str(UTR3_average_rep1) + '\t' + str(UTR3_average_rep2) + '\n')
	data.close()
out.close()

control_dict_rep1 = {}
control_dict_rep2 = {}
with open('distance_control.csv','r') as f:
	for i in f:
		i = i.strip().split('\t')
		if i[0] in control_dict_rep1.keys():
			control_dict_rep1[i[0]] += ',' + i[2]
			control_dict_rep2[i[0]] += ',' + i[3]
		else:
			control_dict_rep1[i[0]] = i[2]
			control_dict_rep2[i[0]] = i[3]
f.close()

treat_dict_rep1 = {}
treat_dict_rep2 = {}
out = open('distance_treat2.csv','w')
with open('distance_treat.csv','r') as f:
	for i in f:
		i = i.strip().split('\t')
		if i[0] in treat_dict_rep1.keys():
			treat_dict_rep1[i[0]] += ',' + i[2]
			treat_dict_rep2[i[0]] += ',' + i[3]
		else:
			treat_dict_rep1[i[0]] = i[2]
			treat_dict_rep2[i[0]] = i[3]
	for i in treat_dict_rep1.keys():
		out.write(i + '\t' + treat_dict_rep1[i] + '\t' + treat_dict_rep2[i] + '\n')
f.close()
out.close()
out = open('final_distance.csv','w')
with open('distance_treat2.csv','r') as fh:
	for i in fh:
		i = i.strip().split('\t')
		if i[0] in control_dict_rep1.keys():
			out.write(i[0] + "\t" + control_dict_rep1[i[0]] + "\t" + control_dict_rep2[i[0]] + '\t' + i[1] + '\t' + i[2] + '\n')
fh.close()
