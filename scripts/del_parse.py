#!/usr/bin/env python
from __future__ import print_function
import sys
## iterates through coverage files to find loci with no coverage
"""
Author: Matthew Ezewudo

CPTR ReSeqTB Project - Critical Path Institute
"""


dele = False
ids   = []
names = []
loci  = []
start = []
end   = []

input3 = sys.argv[3]
fh3 = open(input3,'r')
for lines in fh3:
    fields = lines.rstrip("\r\n").split("\t")
    start.append(fields[1])
    end.append(fields[2])
    ids.append(fields[4])
    names.append(fields[3])
fh3.close()

input1 = sys.argv[1]
input2 = sys.argv[2]

fh1 = open(input1 + "_genome_region_coverage.txt",'r')
print("Sample ID" + "\t" + "CHROM" + "\t" + "START" + "\t" + "END" + "\t" + "GENENAME" + "\t" + "GENEID" + "\t" + "VARIANT TYPE")

for lines in fh1:
    if lines.startswith("Chrom"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    ind    = ids.index(fields[4])
    if float(fields[5]) < 10.0:
       dele = True
       print(input2 + "\t"  + "NC_000962" + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\t" + fields[4] + "\t" + "Complete deletion")
    elif int(fields[6]) < 96:
       dele = True
       print(input2 + "\t"  + "NC_000962" + "\t" + fields[1] + "\t" + fields[2] + "\t" + fields[3] + "\t" + fields[4] + "\t" + "Partial deletion")
    loci.append(fields[4])
for genes in ids:
    if genes not in loci:
       ind = ids.index(genes)
       dele = True
       print(input2 + "\t"  + "NC_000962" + "\t" + start[ind] + "\t" + end[ind] + "\t" + names[ind] + "\t" + ids[ind] + "\t" + "Complete deletion")
if dele == False:
   print(input2 + "\t"  + "NC_000962" + "\t" + start[0] + "\t" + end[-1] + "\t" + "NA" + "\t" + "NA" + "\t" + "No gene deletion inferred")
