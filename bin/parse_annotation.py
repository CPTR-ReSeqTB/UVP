#! /usr/bin/python

import sys
import re
from string import join


""" The script accepts a SnpEff annotated VCF file and the sample ID name (string) as input options """
""" it parses files and creates a final annotation file that is in a ReseqTB mappable format """


input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]


position            = ""
reference           = ""
alternate           = ""
annotation          = ""
variant             = ""
read_depth          = ""
quality             = ""
perc_alt            = ""
nucleotide_change   = ""
nuc_change          = ""
transcript_pos      = ""
amino_acid_change   = ""
orig_aacid          = ""
new_aacid	    = ""
codon_pos	    = ""
gene_name           = ""
gene_id             = ""
transcript          = ""
annotation_details  = ""
position1           = ""
reference1          = ""
alternate1          = ""
annotation1         = ""
variant1            = ""
read_depth1         = ""
quality1            = ""
perc_alt11          = ""
nucleotide_change1  = ""
transcript_pos1     = ""
amino_acid_change1  = ""
orig_aacid1         = ""
new_aacid1	    = ""
codon_pos1	    = ""
gene_name1          = ""
gene_id1            = ""
transcript1         = ""
annotation_details1 = ""
Block               = False
(genez,genezid,start,stop,gene_anot,strand) = ([],[],[],[],[],[])
nuc_change  = ""
dic = {'A':'T','T':'A','C':'G','G':'C'}
ref_comp = ""
alt_comp = ""

fh3 = open(input3,'r')
for lines in fh3:
    lined = lines.rstrip("\r\n").split("\t")
    if lines.startswith("H37Rv"):
       continue
    genez.append(lined[0])
    genezid.append(lined[1])
    start.append(lined[2])
    stop.append(lined[3])
    gene_anot.append(lined[4])
    strand.append(lined[5])
    
fh1 = open(input1,'r')
print "Sample ID" + "\t" + "CHROM" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "Read Depth" + "\t" + "Quality" + "\t" + "Percent Alt allele" + "\t" +  "Annotation" + "\t" + "Variant Type" + "\t" + "Nucleotide Change" + "\t" + "Position within CDS " + "\t" + "Amino acid change" + "\t" + "REF Amino acid" + "\t" + "ALT Amino Acid" + "\t" + "Codon Position" + "\t" "Gene name" + "\t" + "Gene ID" + "\t" + "Transcript ID" + "\t" + "Annotation details"  

for lines in fh1:
    if lines.startswith("#"):
       continue
    fields = lines.rstrip("\r\n").split("\t")
    position = fields[1]
    reference = fields[3]
    alternate = fields[4]
    quality = fields[5]
    rarr = fields[9].split(":")
    num_all_array = rarr[1].split(",")
    num_all_2 = num_all_array[1]
    read_depth = rarr[2]
    num_all = rarr[3]
    perc_alt1 = float(num_all_2)/float(read_depth)*100.0
    if perc_alt1 > 100.0:
       perc_alt1 = 100.0
    perc_alt = "{0:.2f}".format(perc_alt1)
    if float(read_depth) < 10.0:
       continue
    subfields = fields[7].split(";")
    if subfields[-1].startswith("ANN"):
       annot = subfields[-1]
    else:
       annot   = subfields[-2]
    subannot   = annot.split(",")
    smallannot = subannot[0].split("|")
    if smallannot[2] == "MODIFIER":
       for x in range(0,82):
           if (int(start[x]) -1) < int(position) < (int(stop[x]) + 1):
              annotation = gene_anot[x]
              if genez[x] == 'rrs':
                 nuc_change = str((int(position)) - (int(start[x]) - 1))
                 gene_id = 'MTB000019'
                 nucleotide_change = "c." + nuc_change + reference + ">" + alternate
              elif genez[x] == 'rrl':
                 nuc_change = str((int(position)) - (int(start[x]) - 1))
                 gene_id = 'MTB000020'
                 nucleotide_change = "c." + nuc_change + reference + ">" + alternate
              elif genez[x] == 'crfA':
                 nuc_change = str((int(position)) - (int(start[x]) - 1))
                 gene_id = 'crfA'
                 nucleotide_change = "c." + nuc_change + reference + ">" + alternate
              elif strand[x] == 'forward':
                 gene_id = genezid[x]
                 nuc_change = str((int(position)) - (int(stop[x]) + 1))
                 nucleotide_change = "c." + nuc_change + reference + ">" + alternate
              elif strand[x] == 'reverse':
                 ref_comp = ""
                 alt_comp = ""
                 for char in reference:
                     ref_comp += dic[char]
                 for char in alternate:
                     alt_comp += dic[char]
                 gene_id = genezid[x]
                 nuc_change = str((int(start[x]) -1) - int(position))
                 nucleotide_change = "c." + nuc_change + ref_comp + ">" + alt_comp
              gene_name = genez[x]
              amino_acid_change  = 'NA'
              if len(fields[4]) > len(fields[3]):
                 if strand[x] == 'forward':
                    nucleotide_change = "c." + nuc_change + "_" + str(int(nuc_change) + 1) + "ins" + alternate[len(reference):]
                    if genez[x] == 'crfA':
                       transcript_pos = nuc_change + "-" + str(int(nuc_change) + 1)
                 elif strand[x] == 'reverse':
                    nucleotide_change = "c." + str(int(nuc_change) - 1) + "_" + nuc_change + "ins" + alt_comp[len(reference):][::-1]  
                 variant = "Insertion"
              elif len(fields[3]) > len(fields[4]):
                 if strand[x] == 'forward':
                    if len(reference) - len(alternate) == 1:
                       nucleotide_change = "c." + str(int(nuc_change) + len(alternate)) + "del" + reference[len(alternate):]
                       if genez[x] == 'crfA':
                          transcript_pos = str(int(nuc_change) + len(alternate))
                    else:
                       nucleotide_change = "c." + str(int(nuc_change) + len(alternate)) + "_" + str(int(nuc_change) + len(reference) - 1) + "del" + reference[len(alternate):]
                       if genez[x] == 'crfA':
                          transcript_pos = str(int(nuc_change) + len(alternate)) + "-" + str(int(nuc_change) + len(reference) - 1)
                 elif strand[x] == 'reverse':
                    if len(reference) - len(alternate) == 1:
                       nucleotide_change = "c." + str(int(nuc_change) - len(reference) - 1) + "del" + ref_comp[len(alternate):][::-1]
                    else:
                       nucleotide_change = "c." + str(int(nuc_change) - len(reference) - 1) + "_" + str(int(nuc_change) - len(alternate))  + "del" + ref_comp[len(alternate):][::-1]
                 variant = "Deletion"
              else:
                 variant = "SNP"
              transcript         = 'NA'
              if genez[x] == 'crfA':
                 transcript_pos  = nuc_change
              else:
                 transcript_pos  = 'NA'
              orig_aacid         = 'NA'
              new_aacid	         = 'NA'
              codon_pos	         = 'NA'
              annotation_details = ','.join(subannot[0:])
              break               
           else:
                annotation         = 'Non-Coding'
                variant            = 'NA'
                nucleotide_change  = smallannot[9]
                amino_acid_change  = 'NA'
                gene_name          = 'NA'
                gene_id            = 'NA'
                transcript         = 'NA'
                transcript_pos     = 'NA'
                orig_aacid         = 'NA'
                new_aacid	   = 'NA'
                codon_pos	   = 'NA'
                annotation_details = ','.join(subannot[1:])
       if len(position1) != 0:
          if Block == True:
              print  input2 + "\t"  + fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + 'MNV' + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + 'Block_Substitution' + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1 + "\t" + annotation_details1
              Block = False
          else:
              print  input2 + "\t" + fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + variant1 + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + new_aacid1 + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1 + "\t" +  annotation_details1
        
    else:
        
        if smallannot[10][2:5] == smallannot[10][-3:]:
           annotation = 'Synonymous'
        else:
           annotation = 'Non-synonymous'
        nucleotide_change  = smallannot[9]
        if len(smallannot[10]) < 1:
           amino_acid_change  = "NA"
        else:
           amino_acid_change  = smallannot[10]
        gene_name          = smallannot[3]
        gene_id            = smallannot[4]
        transcript         = smallannot[6]
        if gene_name == 'erm_37_':
           gene_name = 'erm(37)'
        annotation_details = ','.join(subannot[1:])
        if 'del' in nucleotide_change or 'del' in amino_acid_change:
           variant = 'Deletion'
        elif 'ins' in nucleotide_change or 'ins' in amino_acid_change:
            variant = 'Insertion'
        elif 'dup' in nucleotide_change or 'dup' in amino_acid_change:
           variant = 'Insertion'
        else:
            variant = 'SNP'
        if variant == 'Insertion' or variant == 'Deletion':
           new_aacid = 'NA'
           if '_' in smallannot[9]:
              array1 = smallannot[9].split("_")
              po1 = array1[0].split(".")
              pos1 = po1[1]
              pos2 = re.findall(r'\d+', array1[1])[0]
              transcript_pos = pos1 + "-" + pos2
           else:
              transcript_pos = re.findall(r'\d+', smallannot[9])[0]
           if '_' in smallannot[10]:
              array2 = smallannot[10].split("_")
              po11 = array2[0].split(".")
              orig_aacid = po11[1][0:3]
              pos11  = po11[1][3:]
              pos12  = re.findall(r'\d+', array2[1])[0]
              codon_pos = pos11 + "-" + pos12
           
           else:
              if len(smallannot[10]) > 0:
                 codon_pos = re.findall(r'\d+', smallannot[10])[0]
                 orig_aacid = smallannot[10][2:5]
              else:
                  codon_pos =  "NA"
                  orig_aacid = "NA"     
        else :
            orig_aacid = smallannot[10][2:5]
            if '*' in smallannot[10] or '?' in smallannot[10] :
               new_aacid  = 'NA'
            else:
               new_aacid  = smallannot[10][-3:]
            transcript_pos = re.findall(r'\d+', smallannot[9])[0]
            codon_pos = re.findall(r'\d+', smallannot[10])[0]

        for x in range(0,82):
           if (int(start[x]) -1) < int(position) < (int(stop[x]) + 1):
              annotation = gene_anot[x]    
              if strand[x] == 'forward':
                 gene_id  =  genezid[x]
                 nuc_change = str((int(position)) - (int(stop[x]) + 1))
                 nucleotide_change = "c." + nuc_change + reference + ">" + alternate
              elif strand[x] == 'reverse':
                 ref_comp = ""
                 alt_comp = ""
                 for char in reference:
                     ref_comp += dic[char]
                 for char in alternate:
                     alt_comp += dic[char]
                 gene_id  =  genezid[x]
                 nuc_change = str((int(start[x]) -1) - int(position))
                 nucleotide_change = "c." + nuc_change + ref_comp + ">" + alt_comp
              gene_name = genez[x]
              amino_acid_change  = 'NA'
              if len(fields[4]) > len(fields[3]):
                 if strand[x] == 'forward':
                    nucleotide_change = "c." + nuc_change + "_" + str(int(nuc_change) + 1) + "ins" + alternate[len(reference):]
                 elif strand[x] == 'reverse':
                    nucleotide_change = "c." + str(int(nuc_change) - 1) + "_" + nuc_change + "ins" + alt_comp[len(reference):][::-1]
                 variant = "Insertion"
              elif len(fields[3]) > len(fields[4]):
                 if strand[x] == 'forward':
                    if len(reference) - len(alternate) == 1:
                       nucleotide_change = "c." + str(int(nuc_change) + len(alternate)) + "del" + reference[len(alternate):]
                    else:
                       nucleotide_change = "c." + str(int(nuc_change) + len(alternate)) + "_" + str(int(nuc_change) + len(reference) - 1) + "del" + reference[len(alternate):]
                 elif strand[x] == 'reverse':
                    if len(reference) - len(alternate) == 1:
                       nucleotide_change = "c." + str(int(nuc_change) - len(reference) - 1) + "del" + ref_comp[len(alternate):][::-1]
                    else:
                       nucleotide_change = "c." + str(int(nuc_change) - len(reference) - 1) + "_" + str(int(nuc_change) - len(alternate))  + "del" + ref_comp[len(alternate):][::-1]
                 variant = "Deletion"
              else:
                 variant = "SNP"
              transcript         = 'NA'
              transcript_pos     = 'NA'
              orig_aacid         = 'NA'
              new_aacid          = 'NA'
              codon_pos          = 'NA'
              annotation_details = ','.join(subannot[0:])
              break

        if len(position1) != 0:
           if codon_pos == codon_pos1 and (int(position) - int(position1)) < 4 and float(perc_alt11) > 98.0 :
              Block = True   
              print  input2 + "\t"  + fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + 'MNV' + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + 'Block_Substitution' + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1 + "\t" +  annotation_details1
           elif Block == True:
              print input2 + "\t" + fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + 'MNV' + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + 'Block_Substitution' + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1 + "\t" +  annotation_details1
              Block = False
           else:   
              print input2 + "\t" + fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + variant1 + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + new_aacid1 + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1 + "\t" +  annotation_details1 
           
    position1           = position
    reference1          = reference
    alternate1          = alternate
    annotation1         = annotation
    variant1            = variant
    read_depth1         = read_depth
    quality1            = quality
    perc_alt11          = perc_alt
    nucleotide_change1  = nucleotide_change
    transcript_pos1     = transcript_pos
    amino_acid_change1  = amino_acid_change
    orig_aacid1         = orig_aacid
    new_aacid1	        = new_aacid
    codon_pos1	        = codon_pos
    gene_name1          = gene_name
    gene_id1            = gene_id
    transcript1         = transcript
    annotation_details1 = annotation_details   
                 
if Block == True:
              print input2 + "\t" +  fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + 'MNV' + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + 'Block_Substitution' + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1 + "\t" +  annotation_details1

else:
    print  input2 + "\t" + fields[0] + "\t" + position1 + "\t" + reference1 + "\t" + alternate1 + "\t" + read_depth1 + "\t" + quality1 + "\t" + perc_alt11 + "\t" + annotation1 + "\t" + variant1 + "\t" + nucleotide_change1 + "\t" + transcript_pos1 + "\t" + amino_acid_change1 + "\t" + orig_aacid1 + "\t" + new_aacid1 + "\t" + codon_pos1 + "\t" + gene_name1 + "\t" + gene_id1 + "\t" + transcript1  + "\t" +  annotation_details1

