#! /usr/bin/python

import sys
from string import join

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
print "Average Genome Coverage Depth: " + str(av_depth)
print "Percentage of Reference genome covered: " + perc_cov_str
