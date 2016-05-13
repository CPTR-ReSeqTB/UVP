#! /usr/bin/python

## iterates through coverage files to fine loci with no coverage

import sys
import re
from string import join

dele = False
loci = []
ids = ["Rv0005","Rv0006","Rv0407","Rv0486","Rv0635","Rv0667","Rv0676c","Rv0677c","Rv0678","Rv0682","Rv0701","Rv1137","Rv1258c","Rv1305","rrs","rrl","fabG11/inhA promoter","Rv1483","Rv1484","Rv1694","Rv1704c","Rv1908c","Rv1909c","Rv1988","Rv2043c","pncA promoter","Rv2416c","eis promoter","Rv2428","ahpC promoter","Rv2447c","Rv2671","Rv2763c","Rv2764c","Rv3197A","whiB7 promoter","Rv3261","Rv3262","Rv3547","Rv3792","Rv3793","Rv3794","Rv3795","Rv3806c","Rv3854c","Rv3919c"]


input1 = sys.argv[1]
input2 = sys.argv[2]

fh1 = open(input1 + "_Resistance_Region_Coverage.txt",'r')

for lines in fh1:
    if lines.startswith("Chrom"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    if int(fields[4]) < 10:
       dele = True
       print input2 + "\t" + fields[3] + "\t" + "Complete deletion"
    elif int(fields[5]) < 90:
       dele = True
       print input2 + "\t" + fields[3] + "\t" + "Partial deletion"
    loci.append(fields[3])
for genes in ids:
    if genes not in loci:
       dele = True
       print input2 + "\t" + genes + "\t" + "Complete deletion"
if dele == False:
   print input2 + "\tNo gene deletion inferred"



