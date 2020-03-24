#!/bin/miniconda3/bin/python

"""bed to fa"""

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
parser.add_option('-b','--bed',
                dest='bed',
                help='This should be bed format')

parser.add_option('-g','--genome',
                dest='genome',
                help='genome file')
(options,args) = parser.parse_args()



sense = open("sense.table",'w')
antisense = open("antisense.table",'w')

data = open(options.bed,'r')

for i in data:
	i = i.split('\t')
	strand = i[3]
	if strand == "+":
		a = int(i[4])-1
		sense.write(i[0] + "\t" + str(a-100) + "\t" + str(a+101) + "\n")
	elif strand == "-":
		a = int(i[4])-1
		antisense.write(i[0] + "\t" + str(a-100) + "\t" + str(a+101) + "\n")
sense.close()
antisense.close()

os.system('bedtools getfasta -fi %s -bed sense.table -fo sense.table.fa' % (options.genome))
os.system('bedtools getfasta -fi %s -bed antisense.table -fo antisense.table.fa' % (options.genome))

os.system("sed -i '/^>/d' sense.table.fa")
os.system("sed -i '/^>/d' antisense.table.fa")

with open('sense.table.fa','r') as table:
	with open('sense.table.fa.matrix','w') as f:
		for i in table:
        		i = list(i)
			if "N" in i:
				pass
			else:
        			for j in i:
                			f.write(j + "\t")
table.close()
f.close()


complement_dict = {"A":"T","C":"G","G":"C","T":"A"}

with open('antisense.table.fa','r') as table:
        with open('antisense.table.fa.matrix','w') as f:
                for i in table:
                        i = list(i)
			if "N" in i:
                                pass
                        else:
			#reverse
				i.reverse()
                        	for j in i:
					if j != '\n':
						j = str(j)
						a = complement_dict[j]
                                		f.write(a + '\t')
					else:
						f.write('\n')
table.close()
f.close()

os.system("cat sense.table.fa.matrix antisense.table.fa.matrix > Total.table.fa.matrix")











