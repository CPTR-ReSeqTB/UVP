#! /usr/bin/python

import sys

input1 = sys.argv[1]
input2 = sys.argv[2]
(start,end,features,length,perc_cov,temp_start) = ([],[],[],[],[],[])
(new_start,new_end,sum_cov) = ("","",0)
SnX = False

fh1 = open(input1,'r')
for lines in fh1:
    fields = lines.rstrip("\r\n").split("\t")
    start.append(fields[1])
    end.append(fields[2])
    features.append(fields[3])
    length.append(fields[6])
    perc_cov.append(fields[7])

fh1.close()

fh2 = open(input2,'r')
for lines in fh2:
    fields = lines.rstrip("\r\n").split("\t")
    if fields[1] in end:
       new_end = fields[1]
       ind = end.index(fields[1])
       avg_cov = sum_cov/int(length[ind])
       avg_cov_str = str(avg_cov)
       covp = int(float(perc_cov[ind]) * 100)
       covp_str = str(covp)
       print "NC_000962" + "\t" + new_start + "\t" + new_end + "\t" + features[ind] + "\t" + avg_cov_str + "\t" + covp_str
       if len(temp_start) > 0:
          new_start = temp_start[0]
          SnX = True
       else:
          SnX = False
       temp_start = []
       sum_cov = 0
    elif fields[1] in start:
       if SnX == False:
          new_start = fields[1]
          SnX = True
       else:
         temp_start.append(fields[1])
       sum_cov += int(fields[2])
    elif SnX == False:
         continue
    else:
         sum_cov += int(fields[2])
fh2.close()


       
