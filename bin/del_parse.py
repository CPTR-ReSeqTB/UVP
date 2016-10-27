#! /usr/bin/python

## iterates through coverage files to find loci with no coverage

import sys
import re
from string import join

dele = False
ids   = []
names = []
loci  = []

input3 = sys.argv[3]
fh3 = open(input3,'r')
for lines in fh3:
    fields = lines.rstrip("\r\n").split("\t")
    ids.append(fields[4])
    names.append(fields[3])
fh3.close()

input1 = sys.argv[1]
input2 = sys.argv[2]

fh1 = open(input1 + "_genome_region_coverage.txt",'r')
print "Sample ID" + "\t" + "CHROM" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "Read Depth" + "\t" + "Quality" + "\t" + "Percent Alt allele" + "\t" +  "Annotation" + "\t" + "Variant Type" + "\t" + "Nucleotide Change" + "\t" + "Position within CDS " + "\t" + "Amino acid change" + "\t" + "REF Amino acid" + "\t" + "ALT Amino Acid" + "\t" + "Codon Position" + "\t" "Gene name" + "\t" + "Gene ID" + "\t" + "Transcript ID" + "\t" + "Annotation details"

for lines in fh1:
    if lines.startswith("Chrom"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    ind    = ids.index(fields[4])
    if float(fields[5]) < 10.0:
       dele = True
       print  input2 + "\t"  + "NC_000962" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "Complete deletion" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + names[ind] + "\t" + fields[4] + "\t" + "NA" + "\t" + "NA"
    elif int(fields[6]) < 96:
       dele = True
       print  input2 + "\t"  + "NC_000962" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "Partial deletion" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + names[ind] + "\t" + fields[4] + "\t" + "NA" + "\t" + "NA"
    loci.append(fields[4])
for genes in ids:
    if genes not in loci:
       ind = ids.index(genes)
       dele = True
       print  input2 + "\t"  + "NC_000962" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "Complete deletion" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + names[ind] + "\t" + genes + "\t" + "NA" + "\t" + "NA"
if dele == False:
   print  input2 + "\t"  + "NC_000962" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "No gene deletion inferred" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\t" + "NA"
