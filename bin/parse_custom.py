#! /usr/bin/python

import sys
from string import join


genes = []
des = []
input1 = sys.argv[1]
input2 = sys.argv[2]
count = 0

fh1 = open(input1,'r')
for lines in fh1:
    fields = lines.rstrip("\r\n").split("\t")
    genes.append(fields[3])
    if len(fields) > 4:
       des.append(fields[4])
fh1.close()

fh2 = open(input2,'r')
for lines in fh2:
    if lines.startswith(";") or lines.startswith("SNP"):
       continue
    line = lines.rstrip("\r\n")
    folds = line.split("\t")
    if folds[14] in genes or folds[14].find('PE') != -1:
           count +=1
           continue
    print line
  
   
