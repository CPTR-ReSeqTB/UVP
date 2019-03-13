#!/usr/bin/env python
from __future__ import print_function
import sys

## compares raw vcf file to final annotation file to find low confidence variants
"""
Author: Matthew Ezewudo

CPTR ReSeqTB Project - Critical Path Institute
"""

input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]

(id,pos,ref,alt) = ([],[],[],[])

fh1 = open(input1,'r')
fh2 = open(input2,'r')
fh3 = open(input3,'a')


for line in fh1:
  if line.startswith("Sample ID"):
     continue
  fields = line.rstrip("\r\n").split("\t")
  id.append(fields[0])
  pos.append(fields[2])
  ref.append(fields[3])
  alt.append(fields[4])

fh1.close()

for line in fh2:
  if line.startswith("#"):
     continue
  fields = line.rstrip("\r\n").split("\t")
  if fields[1] not in pos and fields[6] != 'PASS':
     print(id[0] + "\t" + fields[1] + "\t" + fields[3] + "\t" + fields[4] + "\t" + fields[6], file=fh3)

fh2.close()
fh3.close()

