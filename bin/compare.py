#! /usr/bin/python

import sys
from string import join
(pos,gene)   = ([],[])
(pos2,gene2) = ([],[])
(diff,diff2) = ([],[])
count = 0
input1 = sys.argv[1]
input2 = sys.argv[2]

fh1 = open(input1,'r')
for lines in fh1:
    if lines.startswith(";") or lines.startswith("SNP"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    pos.append(fields[3])
    gene.append(fields[14])

fh1.close()

fh2 = open(input2,'r')
for lines in fh2:
    if lines.startswith("CHROM"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    pos2.append(fields[1])
    gene2.append(fields[11])

fh2.close()

for site in range(0,len(pos2)):
    if pos2[site] in pos:
       count += 1
       diff.append(pos2[site])
       diff2.append(gene2[site])
for i in range(0,len(diff)):
    print diff[i] + "\t" + diff2[i] 
print "similarity in variant positions:  " + str(count) 
       

