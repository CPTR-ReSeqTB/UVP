#! /usr/bin/python

import sys

input1 = sys.argv[1]
fh1 = open(input1,'r')
moree = False
for lines in fh1:
    fields = lines.rstrip("\r\n").split("\t")
    if fields[1] == '50bp_duplicated':
       continue
    if fields[0].startswith("IG"):
       rvs = fields[0].split("_")
       newrvs = rvs[1].split("-")
       fields[0] = newrvs[0]
       nextrv = newrvs[1]
       moree = True
    if moree:
       print fields[0] + "\t" + "NA"
       print nextrv + "\t" + "NA"
    else:
       print fields[0] + "\t" + fields[1]
fh1.close()

