#! /usr/bin/python

""" The script accepts a SnpEff annotated VCF file and the sample ID name (string) as input options """
""" it parses files and creates a final annotation file that is in a ReseqTB mappable format """

import sys
import re
from string import join

input1 = sys.argv[1]
input2 = sys.argv[2]

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
    read_depth = rarr[2]
    num_all = rarr[3]
    perc_alt1 = float(num_all)/float(read_depth)*100.0
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
       if "rrs" in annot and 1471845 < int(position) < 1473382 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(int(position) - 1471845)
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'rrs'
           gene_id            = 'MTB000019'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
       elif "rrl" in annot and 1473657 < int(position) < 1476495 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(int(position) - 1473657)
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'rrl'
           gene_id            = 'MTB000020'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
       elif 1673338 < int(position) < 1673440 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(int(position) - 1673440)
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'inhA promoter'
           gene_id            = 'inhA promoter'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
       elif 2289241 < int(position) < 2289282 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(2289241 - int(position))
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'pncA promoter'
           gene_id            = 'pncA promoter'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
       elif 2715332 < int(position) < 2715384 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(2715332 - int(position))
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'eis promoter'
           gene_id            = 'eis promoter'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
       elif 2725991 < int(position) < 2726193 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(int(position) - 2726193)
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'aphC promoter'
           gene_id            = 'aphC promoter'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
       elif 3568679 < int(position) < 3569109 :
           annotation         = 'Intragenic_variant'
           nuc_change         = str(3568679 - int(position))
           if len(fields[4]) > 1: 
              variant = "INDEL"
           else:
              variant = "SNP"
           nucleotide_change  = "c." + nuc_change + reference + ">" + alternate
           amino_acid_change  = 'NA'
           gene_name          = 'whiB7 promoter'
           gene_id            = 'whiB7 promoter'
           transcript         = 'NA'
           transcript_pos     = 'NA'
           orig_aacid         = 'NA'
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
           annotation_details = ','.join(subannot[0:])
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
           new_aacid	  = 'NA'
           codon_pos	  = 'NA'
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
        amino_acid_change  = smallannot[10]
        gene_name          = smallannot[3]
        gene_id            = smallannot[4]
        transcript         = smallannot[6]
        annotation_details = ','.join(subannot[1:])
        if 'ins' in nucleotide_change or 'ins' in amino_acid_change:
           variant = 'Insertion'
        elif 'del' in nucleotide_change or 'del' in amino_acid_change:
            variant = 'Deletion'
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
        
        if len(position1) != 0:
           if codon_pos == codon_pos1 and (int(position) - int(position1)) < 4 :
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

       
