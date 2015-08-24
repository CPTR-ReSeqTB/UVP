#! /usr/bin/python

import sys
from string import join

input1 = sys.argv[1]
input2 = sys.argv[2]

fh1 = open(input1, 'r')
sublinn = ""
(lineage,position,ref,alt) = ([],[],[],[])
prevlin = []
prevsub = []
tribes = ["lineages","Indo Oceanic","East Asian","East African Indian","Euro American","African I","African II","Ethiopian"]
(concord,discord,concord1,discord1,count) = (0,0,0,0,0)
discordance = False
sublinneage = False
linfour = ""
BOV = ""
BOV_AFRI = ""


for lines in fh1:
    if lines.startswith('#'):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    lineage.append(fields[0])
    position.append(fields[1])
    ref.append(fields[2])
    alt.append(fields[3])
fh1.close()

fh2 = open(input2,'r')
for lines in fh2:
    fields = lines.rstrip("\r\n").split("\t")
    if fields[1] == '1759252':
       linfour = fields[1]
    if fields[1] == '2831482':
       BOV = fields[1]
    if fields[1] == '1882180':
       BOV_AFRI = '1882180'
       
    if fields[1] in position:
       ind = position.index(fields[1])
       if alt[ind] == fields[3]:
          if len(lineage[ind]) > 1:
             sublin = lineage[ind]
             prevsub.append(sublin)
             sublinn = prevsub[0]
             print "SNP" + " " + position[ind] + " " + "suggests sub-lineage: " + lineage[ind]
             if prevsub[0] != sublin:
                discord += 1
             else:
                concord +=1
             for i in range(0,len(prevsub)):
                 if len(sublinn) < len(prevsub[i]) : 
                    sublinn = prevsub[i]
                    
          else:
                lin = lineage[ind]
                prevlin.append(lin)
                print "SNP" + " " + position[ind] + " " + "suggests lineage: " + lineage[ind]
                if prevlin[0] != lin:
                   discord1 += 1
                else:
                  concord1 += 1

split_first = ['NA']
if len(prevsub) > 0:
   split_first = sublinn.split(".")
   sublinneage = True
if len(prevlin) == 0:
   if len(BOV) == 0 or len(BOV_AFRI) == 0:
      for i in range(0,len(prevsub)): 
          split_lin = prevsub[i].split(".")
          if split_lin[0] != split_first[0]:
             discordance = True
          if split_lin[1] != split_first[1]:
             discordance = True
      if discordance:
         print "no precise lineage inferred"
      else:
         print "Lineage: " + split_first[0] + " :  " + tribes[int(split_first[0])]
         print "Sub-lineage: " + sublinn
   elif len(prevsub) == 0:
     if len(linfour) > 2:
         print "SNP" + " " + linfour + "  suggests sublineage 4.9"
         print "Lineage: " + "4" + " :  " + "Euro American"
     else:
         print "No Informative SNPs detected"
else:
     for i in range(0,len(prevlin) - 1):
         if prevlin[i - 1] != prevlin[i]:
            discordance = True
     if discordance == True:
        print "no concordance between predicted lineage and sublineage(s)" 
     if len(prevsub) > 0 and prevlin[0] != split_first[0]:
        if len(BOV_AFRI) > 0:
           print "Lineage: "  + prevlin[0] 
           print "Sub-lineage: " + "BOV_AFRI"
        elif len(BOV) > 0:
           print "Lineage: " + "BOV"
        else:
           print "no concordance between predicted lineage and sublineage"
     else:  
        print "Lineage: " + prevlin[0] + " " + tribes[int(prevlin[0])]
        if len(sublinn) > 1:
           print "Sub-lineage: " + sublinn


              
     



         

