#!/bin/miniconda3/bin/python
"""Create a summary of AAUAAA percentage and 1-nt viaration percentage"""

__date__ = "2019-10-2"
__author__ = "Junbo Yang"
__email__ = "yang_junbo_hi@126.com"
__license__ = "PKU.jia.group"

import os
import re
import optparse
import bisect
#sort  package
from optparse import OptionParser

parser = OptionParser('Usage: %prog ')
parser.add_option('-f','--fasta',
                dest='fa',
                help='fasta format')

parser.add_option('-o','--out',
                dest='out',
                help='percentage file')
(options,args) = parser.parse_args()

data = open(options.fa,'r')
out = open(options.out,'w')

AATAAA=0
TATAAA=0
CATAAA=0
GATAAA=0
ATTAAA=0
others=0
ACTAAA=0
AGTAAA=0
AAAAAA=0
AACAAA=0
AAGAAA=0
AATTAA=0
AATCAA=0
AATGAA=0
AATATA=0
AATACA=0
AATAGA=0
AATAAT=0
AATAAC=0
AATAAG=0

for i in data:
	i = i.strip()
	if re.search("AATAAA",i):
		AATAAA = AATAAA +1	
	elif re.search("TATAAA",i):
		TATAAA = TATAAA +1
	elif re.search("CATAAA",i):
		CATAAA = CATAAA +1
	elif re.search("GATAAA",i):
		GATAAA = GATAAA +1
	elif re.search("ATTAAA",i):
		ATTAAA = ATTAAA +1
	elif re.search("ACTAAA",i):
                ACTAAA = ACTAAA +1
	elif re.search("AGTAAA",i):
                AGTAAA = AGTAAA +1
	elif re.search("AAAAAA",i):
                AAAAAA = AAAAAA +1
	elif re.search("AACAAA",i):
                AACAAA = AACAAA +1
	elif re.search("AAGAAA",i):
                AAGAAA = AAGAAA +1
	elif re.search("AATTAA",i):
                AATTAA = AATTAA +1
	elif re.search("AATCAA",i):
                AATCAA = AATCAA +1
	elif re.search("AATGAA",i):
                AATGAA = AATGAA +1
	elif re.search("AATATA",i):
                AATATA = AATATA +1
	elif re.search("AATACA",i):
                AATACA = AATACA +1
        elif re.search("AATAGA",i):
                AATAGA = AATAGA +1
        elif re.search("AATAAT",i):
                AATAAT = AATAAT +1
        elif re.search("AATAAC",i):
                AATAAC = AATAAC +1
	elif re.search("AATAAG",i):
                AATAAG = AATAAG +1
	else:
		others = others +1


a = AATAAA + TATAAA + CATAAA + GATAAA + ATTAAA + ACTAAA + AGTAAA + AAAAAA + AACAAA + AAGAAA + AATTAA + AATCAA + AATGAA + AATATA + AATACA + AATAGA + AATAAT + AATAAC + AATAAG + others

out.write("AATAAA" + "\t" + str(AATAAA) + str(AATAAA/a) + "\n")
out.write("TATAAA" + "\t" + str(TATAAA) + str(TATAAA/a) + "\n")
out.write("CATAAA" + "\t" + str(CATAAA) + str(CATAAA/a) + "\n")
out.write("GATAAA" + "\t" + str(GATAAA) + str(GATAAA/a) + "\n")
out.write("ATTAAA" + "\t" + str(ATTAAA) + str(ATTAAA/a) + "\n")
out.write("ACTAAA" + "\t" + str(ACTAAA) + str(ACTAAA/a) + "\n")
out.write("AGTAAA" + "\t" + str(AGTAAA) + str(AGTAAA/a) + "\n")
out.write("AAAAAA" + "\t" + str(AAAAAA) + str(AAAAAA/a) + "\n")
out.write("AACAAA" + "\t" + str(AACAAA) + str(AACAAA/a) + "\n")
out.write("AAGAAA" + "\t" + str(AAGAAA) + str(AAGAAA/a) + "\n")
out.write("AATTAA" + "\t" + str(AATTAA) + str(AATTAA/a) + "\n")
out.write("AATCAA" + "\t" + str(AATCAA) + str(AATCAA/a) + "\n")
out.write("AATGAA" + "\t" + str(AATGAA) + str(AATGAA/a) + "\n")
out.write("AATATA" + "\t" + str(AATATA) + str(AATATA/a) + "\n")
out.write("AATACA" + "\t" + str(AATACA) + str(AATACA/a) + "\n")
out.write("AATAGA" + "\t" + str(AATAGA) + str(AATAGA/a) + "\n")
out.write("AATAAT" + "\t" + str(AATAAT) + str(AATAAT/a) + "\n")
out.write("AATAAC" + "\t" + str(AATAAC) + str(AATAAC/a) + "\n")
out.write("AATAAG" + "\t" + str(AATAAG) + str(AATAAG/a) + "\n")
out.write("others" + "\t" + str(others) + str(others/a) + "\n")

data.close()
out.close()





