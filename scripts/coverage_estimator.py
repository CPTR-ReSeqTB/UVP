#!/usr/bin/env python
from __future__ import print_function
import sys

""" This reads the output file from samtools depth command """
""" and iterates over all the base position to calculate """
""" the average coverage depth and width  """
"""
Author: Matthew Ezewudo

CPTR ReSeqTB Project - Critical Path Institute
"""



input1 = sys.argv[1]

depth = 0
count = 0

fh1 = open(input1, 'r')
for lines in fh1:
    fields = lines.rstrip("\r\n").split("\t")
    depth += int(fields[2])
    count += 1
fh1.close()
av_depth = depth/count
perc_cov = float((count/4411532.00)*100.00)
perc_cov_str = "{0:.2f}".format(perc_cov)
print("Unified Analysis Pipeline Version: UVPv2.6.0\n")
print("Average Genome Coverage Depth: " + str(av_depth))
print("Percentage of Reference genome covered: " + perc_cov_str)
