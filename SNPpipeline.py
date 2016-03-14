#! /usr/bin/python
"""

"""
import sys
import subprocess
import os
import types
import gzip
from datetime import datetime

class snp():

    def __init__(self, input, outdir, reference, name, paired,input2, verbose, argString):
        self.name = name
        self.fOut1 = "Output"
        self.fOut = self.fOut1 + "/" + outdir
        self.flog = "Master_log"
        self.input = input
        self.outdir = self.fOut + "/tmp"
        self.tmp = self.outdir + "/tmp"
        self.prinseq = self.fOut + "/prinseq"
        self.qualimap = self.fOut + "/qualimap"
        self.kraken = self.fOut + "/kraken"
        #self.kvarq = self.fOut + "/kvarq"
        self.paired = paired
        self.input2 = input2
        self.verbose = verbose
        self.reference = reference
        self.__finalVCF = ''
        self.__annotation = ''
        self.__final_annotation = ''

        # Make the output directory, and start the log file.
        self.__logged = False
        if not os.path.isfile(self.fOut1):
           self.__CallCommand('mkdir', ['mkdir', self.fOut1])
        self.__CallCommand('mkdir', ['mkdir', self.fOut])
        if not os.path.isfile(self.flog):
           self.__CallCommand('mkdir', ['mkdir', self.flog])
        self.__CallCommand('mkdir', ['mkdir', '-p', self.tmp])
        self.__CallCommand('mkdir', ['mkdir', '-p', self.prinseq])
        self.__CallCommand('mkdir', ['mkdir', '-p', self.qualimap])
        self.__CallCommand('mkdir', ['mkdir', '-p', self.kraken])
        #self.__CallCommand('mkdir', ['mkdir', '-p', self.kvarq])
        self.__log = self.fOut + "/" + self.name + ".log"
        self.__lineage = self.fOut + "/" + self.name + ".lineage_report.txt"
        self.__logFH = open(self.__log, 'w')
        self.__logFH.write(argString + "\n\n")
        #self.__logged = True
        self.__mlog = self.flog + "/" + "master.log"
        self.__logFH2 = open(self.__mlog, 'a')
        #self.__logFH2.write(argString + "\n\n")
        self.__logged = True
				
	# Format Validation
	self.__fastqval        = '/scicomp/home/krt7/binary/fastQValidator'

        #fastq QC
        self.__prinseq         = '/scicomp/home/krt7/binary/prinseq-lite.pl'
        self.__kraken          = '/scicomp/home/krt7/KRAKEN/kraken'
        self.__krakendb        = '/scicomp/home/krt7/KRAKEN/minikraken'
        self.__krakenreport    = '/scicomp/home/krt7/KRAKEN/kraken-report'
        self.__pigz            = '/scicomp/home/krt7/binary/pigz'
        self.__unpigz            = '/scicomp/home/krt7/binary/unpigz'
        
        # Mapping
        self.__bwa      = '/scicomp/home/krt7/binary/bwa'
        self.__samtools = '/scicomp/home/krt7/binary/samtools/samtools'
        self.__qualimap = '/scicomp/home/krt7/binary/qualimap_v2.1.1/qualimap'


        # Picard-Tools
        self.__picard     = '/scicomp/home/krt7/binary/picard/picard.jar'
        self.__sortsam    = '/usr/local/bin/SortSam.jar'
        self.__samformat  = '/usr/local/bin/SamFormatConverter.jar'
        self.__readgroups = '/usr/local/bin/AddOrReplaceReadGroups.jar'
        self.__bamindex   = '/usr/local/bin/BuildBamIndex.jar'
        self.__seqdiction = '/usr/local/bin/CreateSequenceDictionary.jar'
        self.__markdupes  = '/usr/local/bin/MarkDuplicates.jar'

        # SNP / InDel Calling
        self.__gatk = '/scicomp/home/krt7/binary/GenomeAnalysisTK.jar'
        
        # Other
        self.__bcftools       = '/scicomp/home/krt7/binary/bcftools/bcftools'
        self.__vcfannotate    = '/scicomp/home/krt7/binary/vcftools/bin/vcf-annotate'
        self.__vcftools       = '/scicomp/home/krt7/binary/vcftools/bin/vcftools'
        self.__vcfutils       = '/scicomp/home/krt7/binary/bcftools/vcfutils.pl' 
        self.__annotator      = '/scicomp/home/krt7/binary/snpEff/snpEff.jar' 
        self.__parser         = '/scicomp/home/krt7/binary/parse_annotation.py'
        self.__lineage_parser = '/scicomp/home/krt7/binary/lineage_parser.py'
        self.__lineages       = '/scicomp/home/krt7/binary/lineage_markers.txt'
        self.__exclusion      = '/scicomp/home/krt7/binary/loci_filtered.txt'
        self.__excluded       = '/scicomp/home/krt7/binary/excluded_loci.txt'
        self.__coverage_estimator = '/scicomp/home/krt7/binary/coverage_estimator.py'
        self.resloci          = '/scicomp/home/krt7/binary/included_loci.txt'
        self.__bedlist        = '/scicomp/home/krt7/binary/bed_list.txt'     
        self.__resis_parser   = '/scicomp/home/krt7/binary/resis_parser.py'

    """ Shell Execution Functions """
    def __CallCommand(self, program, command):
        """ Allows execution of a simple command. """
        out = ""
        err = ""
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out,err = p.communicate()

        if (type(program) is list):
            o = open(program[1], 'w')
            o.write(out)
            o.close()
            out = ""
            program = program[0]
        
        if (self.__logged):
            self.__logFH.write('---[ '+ program +' ]---\n')
            self.__logFH.write('Command: \n' + ' '.join(command) + '\n\n')
            if out:
                self.__logFH.write('Standard Output: \n' + out + '\n\n')
            if err:
                self.__logFH.write('Standard Error: \n' + err + '\n\n')

        return 1
    
    """ Input Validation """
    def runVali(self):
        self.__ifVerbose("Validating Input file.")
        valiOut = self.fOut + "/validation"
        self.__CallCommand('mkdir', ['mkdir', '-p', valiOut]) 
        """ Validates format and coverage of input fastq files """
    	num_lines = sum(1 for line in gzip.open(self.input))
        lines = int(num_lines)/4.0
        
        if self.paired:
           num_lines2 = sum(1 for line in gzip.open(self.input2))
           lines2 = int(num_lines)/4.0 	 			
	if self.paired:
           self.__CallCommand(['fastQValidator', valiOut + "/result1.out"], [self.__fastqval, '--file', self.input])
           self.__CallCommand(['fastQValidator', valiOut + "/result2.out"], [self.__fastqval, '--file', self.input2])
           output1 = valiOut + "/result1.out"
           output2 = valiOut + "/result2.out"
           self.__CallCommand(['cat', valiOut + "/result.out"], ['cat', output1, output2])
           self.__CallCommand('rm', ['rm', output1, output2 ])
        else:  
	   self.__CallCommand(['fastQValidator', valiOut + '/result.out'], [self.__fastqval, '--file', self.input])
        self.__CallCommand('mv', ['mv', valiOut + '/result.out', valiOut + '/Validation_report.txt'])	
	output = valiOut + "/Validation_report.txt"
        fh2 = open (output, 'r')
        for line in fh2:
            lined=line.rstrip("\r\n")
            if lined.startswith("Returning"):
               comments = lined.split(":")
               if comments[2] != " FASTQ_SUCCESS":
                  self.__CallCommand('mv', ['mv', self.input, self.flog])
                  if self.paired:
                     self.__CallCommand('mv', ['mv', self.input2, self.flog])  
                  self.__logFH.write("Input not in fastq format\n")
                  i = datetime.now()
                  self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:  " + self.input + "\t" + "not in fastq format\n")
                  sys.exit(1)
               if lines < 100000:
                  self.__CallCommand('mv', ['mv', self.input, self.flog])
                  self.__logFH.write("not enough read coverage\n")
                  i = datetime.now()
                  self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.input + "\t" + "not enough read coverage\n")
                  if self.paired:
                     if lines < 100000 or lines2 < 100000:
                        self.__CallCommand('mv', ['mv', self.input2, self.flog])
                        i = datetime.now()
                        self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.input2 + "\t" + "not enough read coverage\n")
                        sys.exit(2)
                  sys.exit(3)
    
    """ Fastq QC """
    def runPrinseq(self):
        self.__ifVerbose("Performing prinseq QC.")
        if self.paired:
           self.__CallCommand(['gunzip', self.prinseq + "/" + self.name + "1.fastq"],['gunzip','-c', self.input])
           self.__CallCommand(['gunzip', self.prinseq + "/" + self.name + "2.fastq"],['gunzip', '-c', self.input2])
           self.__CallCommand('prinseq', [self.__prinseq, '-fastq', self.prinseq + "/" + self.name + '1.fastq', '-fastq2', 
                              self.prinseq + "/" + self.name + '2.fastq', '-min_qual_mean', '20', '-ns_max_n', '0', 
                              '-out_format', '3', '-out_good', self.prinseq + "/" + self.name + "good", 
                              '-out_bad', self.prinseq + "/" + self.name + "bad" ])
        else:
           self.__CallCommand(['gunzip', self.prinseq + "/" + self.name + ".fastq"],['gunzip','-c', self.input])
           self.__CallCommand('prinseq', [self.__prinseq, '-fastq', self.prinseq + "/" + self.name + '.fastq', 
                              '-min_qual_mean', '20', '-ns_max_n', '0', '-out_format', '3',
                              '-out_good', self.prinseq + "/" + self.name + "good", 
                              '-out_bad', 'null', '-log', self.prinseq + "/" + self.name + "log"])
        if self.paired:
           self.__CallCommand('pigz', [self.__pigz, '-p', '10', self.prinseq + "/" + self.name + "good" + "_1.fastq"])
           self.__CallCommand('pigz', [self.__pigz, '-p', '10', self.prinseq + "/" + self.name + "good" + "_2.fastq"])
           self.__CallCommand('rm', ['rm', self.prinseq + "/" + self.name + "bad_1.fastq",
                               self.prinseq + "/" + self.name + "bad_2.fastq",
                               self.prinseq + "/" + self.name + "good_1_singletons.fastq",
                               self.prinseq + "/" + self.name + "good_2_singletons.fastq",
                               self.prinseq + "/" + self.name + '1.fastq', self.prinseq + "/" + self.name + '2.fastq'])
           self.input  = self.prinseq + "/" + self.name + "good" + "_1.fastq.gz"
           self.input2 = self.prinseq + "/" + self.name + "good" + "_2.fastq.gz"
        else:
           self.__CallCommand('gzip', [self.__pigz, '-p', '10', self.prinseq + "/" + self.name + "good" + ".fastq"])
           self.__CallCommand('rm', ['rm', self.prinseq + "/" + self.name + "bad.fastq",
                              self.prinseq + "/" + self.name + "good_singletons.fastq",
                              self.prinseq + "/" + self.name + '.fastq'])
           self.input = self.prinseq + "/" + self.name + "good" + ".fastq.gz"

    """ Species specificity check """
    def runKraken(self):
        self.__ifVerbose("Running Kraken.")
        valiOut = self.fOut + "/validation"
        self.__logFH.write("########## Running Kraken. ##########\n")
        if self.paired:
           self.__CallCommand(['kraken', self.kraken + "/kraken.txt"],[self.__kraken, '--db', 
                               self.__krakendb, '--gzip-compressed', self.input, self.input2,
                               '--paired', '--fastq-input', '--threads', '12', '--classified-out',
                                self.name + "_classified_Reads.fastq"])
           self.__CallCommand(['krakenreport', self.kraken + "/final_report.txt"],[self.__krakenreport, '--db',
                               self.__krakendb, self.kraken + "/kraken.txt"])
        else:
           self.__CallCommand(['kraken', self.kraken + "/kraken.txt"],[self.__kraken, '--db', 
                               self.__krakendb, '--gzip-compressed', self.input, '--fastq-input', 
                               '--threads', '4', '--classified-out', self.name + "_classified_Reads.fastq"])                     
           self.__CallCommand(['krakenreport', self.kraken + "/final_report.txt"],[self.__krakenreport, '--db',
                               self.__krakendb, self.kraken + "/kraken.txt"])
        krakenOut = self.kraken + "/final_report.txt"
        cov = 0
        fh1 = open(krakenOut,'r')
        for lines in fh1:
            fields = lines.rstrip("\r\n").split("\t") 
            if fields[5].find("Mycobacterium tuberculosis") != -1 or fields[5].find("Mycobacterium canettii") != -1:
               cov += float(fields[0])
        if cov < 90:
           self.__CallCommand('mv', ['mv', self.input, self.flog])
           self.__CallCommand('rm', ['rm', self.kraken + "/kraken.txt"])
           self.__logFH.write("not species specific\n")
           i = datetime.now()
           self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.input + "\t" + "not species specific\n")
           if self.paired:
              self.__CallCommand('mv', ['mv', self.input2, self.flog])
              self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.input2 + "\t" + "not species specific\n")
           sys.exit(2) 
    
    """ Aligners """ 
    def runBWA(self, bwa):
        """ Align reads against the reference using bwa."""
        self.__ranBWA = True
        self.__ifVerbose("Running BWA.")
        self.__logFH.write("########## Running BWA. ##########\n")
        bwaOut = self.outdir + "/bwa"
        self.__CallCommand('mkdir', ['mkdir', '-p', bwaOut])
        
        self.__ifVerbose("   Building BWA index.")
        self.__bwaIndex(bwaOut + "/index")
        self.__alnSam = bwaOut + "/bwa.sam"
        self.__bwaLongReads(bwaOut)
        self.__ifVerbose("") 
        self.__processAlignment()
          
    def __bwaIndex(self, out):
        """ Make an index of the given reference genome. """ 
        self.__CallCommand('mkdir', ['mkdir', '-p', out])
        self.__CallCommand('cp', ['cp', self.reference, out + "/ref.fa"])
        self.reference = out + "/ref.fa"
        self.__CallCommand('bwa index', [self.__bwa, 'index', self.reference])
        self.__CallCommand('CreateSequenceDictionary', ['java', '-jar', self.__picard, 
                           'CreateSequenceDictionary', 'R='+self.reference,'O='+ out + "/ref.dict"])
        self.__CallCommand('samtools faidx', ['samtools', 'faidx', self.reference ])

    
    def __bwaLongReads(self, out):
        """ Make use of bwa mem """
        if self.paired:
            self.__ifVerbose("   Running BWA mem on paired end reads.")
            self.__CallCommand(['bwa mem', self.__alnSam], [self.__bwa, 'mem', '-t', '10', '-R', 
                               "@RG\tID:" + self.name + "\tSM:" + self.name + "\tPL:ILLUMINA", 
                                self.reference, self.input, self.input2])
        else:
            self.__ifVerbose("   Running BWA mem on single end reads.")
            self.__CallCommand(['bwa mem', self.__alnSam], [self.__bwa, 'mem', '-t', '10', '-R', 
                               "@RG\tID:" + self.name + "\tSM:" + self.name + "\tPL:ILLUMINA", 
                                self.reference, self.input])       
    
    def __processAlignment(self):
        """ Filter alignment using GATK and Picard-Tools """
        self.__ifVerbose("Filtering alignment with GATK and Picard-Tools.")
        self.__logFH.write("########## Filtering alignment with GATK and Picard-Tools. ##########\n")
        GATKdir = self.outdir + "/GATK"
        self.__CallCommand('mkdir', ['mkdir', '-p', GATKdir])

        """ Convert SAM to BAM"""
        if (self.__ranBWA):
            self.__ifVerbose("   Running SamFormatConverter.")
            self.__CallCommand('SamFormatConverter', ['java', '-Xmx4g', '-jar', self.__picard, 'SamFormatConverter',  
                               'INPUT='+ self.__alnSam, 'VALIDATION_STRINGENCY=LENIENT', 
                               'OUTPUT='+ GATKdir +'/GATK.bam', ])
        else:
            self.__CallCommand('cp', ['cp', self.__alnSam, GATKdir +'/GATK.bam'])


        """ Run mapping Report and Mark duplicates using Picard-Tools"""
        self.__ifVerbose("   Running SortSam.")
        self.__CallCommand('SortSam', ['java', '-Xmx4g', '-Djava.io.tmpdir=' + self.tmp, '-jar', self.__picard, 'SortSam',  
                           'INPUT='+ GATKdir +'/GATK.bam', 'SORT_ORDER=coordinate', 'OUTPUT='+ GATKdir +'/GATK_s.bam', 
                           'VALIDATION_STRINGENCY=LENIENT', 'TMP_DIR=' + self.tmp])
        self.__ifVerbose("   Running Qualimap.")
        self.__CallCommand('qualimap bamqc', [self.__qualimap, 'bamqc', '-bam', GATKdir +'/GATK_s.bam', '-outdir', self.qualimap])
        self.__ifVerbose("   Running MarkDuplicates.")
        self.__CallCommand('MarkDuplicates', ['java', '-Xmx8g', '-jar', self.__picard, 'MarkDuplicates',  
                           'INPUT='+ GATKdir +'/GATK_s.bam', 'OUTPUT='+ GATKdir +'/GATK_sd.bam',
                           'METRICS_FILE='+ GATKdir +'/MarkDupes.metrics', 'ASSUME_SORTED=true', 
                           'REMOVE_DUPLICATES=false', 'VALIDATION_STRINGENCY=LENIENT'])         
        self.__ifVerbose("   Running AddOrReplaceReadGroups.")
        self.__CallCommand('AddOrReplaceReadGroups', ['java', '-Xmx8g', '-jar', self.__picard, 'AddOrReplaceReadGroups', 
                           'INPUT='+ GATKdir +'/GATK_sd.bam', 'OUTPUT='+ GATKdir +'/GATK_sdr.bam',
                           'SORT_ORDER=coordinate', 'RGID=GATK', 'RGLB=GATK', 'RGPL=Illumina', 'RGSM=GATK', 
                           'RGPU=GATK', 'VALIDATION_STRINGENCY=LENIENT'])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand('BuildBamIndex', ['java', '-Xmx8g', '-jar', self.__picard, 'BuildBamIndex',  
                           'INPUT='+ GATKdir +'/GATK_sdr.bam', 'VALIDATION_STRINGENCY=LENIENT'])

        """ Re-alignment around InDels using GATK """
        self.__ifVerbose("   Running RealignerTargetCreator.")
        self.__CallCommand('RealignerTargetCreator', ['java', '-Xmx32g', '-jar', self.__gatk, '-T', 
                           'RealignerTargetCreator', '-I', GATKdir +'/GATK_sdr.bam', '-R', self.reference, 
                           '-o', GATKdir +'/GATK.intervals', '-nt', '8'])
        self.__ifVerbose("   Running IndelRealigner.")
        self.__CallCommand('IndelRealigner', ['java', '-Xmx4g', '-jar', self.__gatk, '-T', 'IndelRealigner', '-l', 
                           'INFO', '-I', GATKdir +'/GATK_sdr.bam', '-R', self.reference, '-targetIntervals', 
                           GATKdir +'/GATK.intervals', '-o', GATKdir +'/GATK_sdrc.bam'])
        self.__ifVerbose("   Running SortSam.")
        self.__CallCommand('SortSam', ['java', '-Xmx4g', '-Djava.io.tmpdir=' + self.tmp, '-jar', self.__picard,'SortSam',  
                           'INPUT='+ GATKdir +'/GATK_sdrc.bam', 'SORT_ORDER=coordinate', 'TMP_DIR=' + self.tmp, 
                           'OUTPUT='+ GATKdir +'/GATK_sdrcs.bam', 'VALIDATION_STRINGENCY=LENIENT'])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand('BuildBamIndex', ['java', '-Xmx4g', '-jar', self.__picard, 'BuildBamIndex', 
                           'INPUT='+ GATKdir +'/GATK_sdrcs.bam', 'VALIDATION_STRINGENCY=LENIENT'])

        """ Filter out unmapped reads """
        self.__finalBam = self.fOut + '/'+ self.name + '_sdrcsm.bam'
        self.__ifVerbose("   Running samtools view.")
        self.__CallCommand('samtools view', ['samtools', 'view', '-bhF', '4', '-o', self.__finalBam, 
                           GATKdir +'/GATK_sdrcs.bam'])
        self.__ifVerbose("   Running BuildBamIndex.")
        self.__CallCommand('BuildBamIndex', ['java', '-Xmx4g', '-jar', self.__picard, 'BuildBamIndex', 'INPUT='+ self.__finalBam, 
                           'VALIDATION_STRINGENCY=LENIENT'])
        self.__ifVerbose("")
        self.__CallCommand('rm', ['rm', '-r', self.tmp])
    
    """ Callers """
    def runSamTools(self):
        """ Call SNPs and InDels using SamTools """
        if os.path.isfile(self.__finalBam):
            self.__ifVerbose("Calling SNPs/InDels with SamTools.")
            self.__logFH.write("########## Calling SNPs/InDels with SamTools. ##########\n")
            samDir = self.outdir + "/SamTools"
            self.__CallCommand('mkdir', ['mkdir', '-p', samDir])

            """ Call SNPs / InDels using mpileup, bcftools, vcfutils. """
            self.__ifVerbose("   Running samtools mpileup.")
            self.__CallCommand(['samtools mpileup', samDir + '/samtools.mpileup'], ['samtools', 'mpileup', '-Q', '20', '-q', '20', '-t', 'DP,DV,DPR', 
                               '-ugf', self.reference, self.__finalBam])
            self.__ifVerbose("   Running bcftools view.")
            self.__CallCommand(['bcftools view', samDir + '/samtools.vcf'], 
                               ['bcftools', 'call', '-vcf', 'GQ', samDir + '/samtools.mpileup'])
            self.__ifVerbose("   Running vcfutils.pl varFilter.")
            self.__CallCommand(['vcfutils.pl varFilter', samDir +'/SamTools.vcf'], 
                               [self.__vcfutils, 'varFilter', '-D1500', samDir + '/samtools.vcf'])
            self.__ifVerbose("   Filtering VCf file using vcftools.")
            self.__CallCommand(['vcf-annotate filter', self.fOut + "/" + self.name +'_SamTools.vcf'], 
                               ['vcf-annotate', '--filter', 'SnpCluster=3,10/Qual=20/MinDP=10/MinMQ=20', samDir +'/SamTools.vcf'])
            self.__CallCommand(['vcftools remove-filtered-all', self.fOut + "/" + self.name +'_SamTools_Resistance_filtered.vcf'], 
                                  [self.__vcftools, '--vcf', self.fOut + "/" + self.name +'_SamTools.vcf',
                                  '--stdout', '--bed', self.resloci, '--remove-filtered-all', '--recode', '--recode-INFO-all'])
            self.__CallCommand(['vcftools remove-filtered-all', self.fOut + "/" + self.name +'_SamTools_filtered.vcf'], 
                                   [self.__vcftools, '--vcf', self.fOut + "/" + self.name +'_SamTools.vcf',
                                   '--stdout', '--exclude-bed', self.__excluded, '--remove-filtered-all', '--recode', '--recode-INFO-all'])
            self.__CallCommand('mv', ['mv', samDir + '/samtools.mpileup', self.fOut + "/" + self.name + '.mpileup'])
            self.__CallCommand(['samtools depth', samDir + '/coverage.txt'],
                                ['samtools','depth', self.__finalBam])
            self.__CallCommand(['bedtools coverage', samDir + '/bed_coverage.txt' ],
                                ['bedtools','coverage', '-abam', self.__finalBam, '-b', self.__bedlist])
            self.__CallCommand(['sort', samDir + '/bed_sorted_coverage.txt' ],
                                ['sort', '-nk', '2', samDir + '/bed_coverage.txt'])                 
            
            """ Set final VCF """
            if not self.__finalVCF: 
                self.__finalVCF = self.fOut + "/" + self.name +'_SamTools_filtered.vcf'
                
        else:
            # print error  
            pass  
    
            
    def annotateVCF(self):
        """ Annotate the final VCF file """
        if self.__finalVCF:
           self.__ifVerbose("Annotating final VCF.")
           self.__CallCommand(['SnpEff', self.fOut + "/" + self.name +'_annotation.txt'],
                                ['java', '-Xmx4g', '-jar', self.__annotator, 'NC_000962', self.__finalVCF])
           self.__annotation = self.fOut + "/" + self.name +'_annotation.txt'
           self.__ifVerbose("parsing final Annotation.")
           self.__CallCommand(['parse annotation', self.fOut + "/" + self.name +'_Final_annotation.txt'],
                              ['python', self.__parser, self.__annotation, self.name])

           self.__CallCommand(['SnpEff', self.fOut + "/" + self.name +'_Resistance_annotation.txt'],
                                ['java', '-Xmx4g', '-jar', self.__annotator, 'NC_000962', self.fOut + "/" + self.name +'_SamTools_Resistance_filtered.vcf'])
           self.__ifVerbose("parsing final Annotation.")
           self.__CallCommand(['parse annotation', self.fOut + "/" + self.name +'_Resistance_Final_annotation.txt'],
                              ['python', self.__parser, self.fOut + "/" + self.name +'_Resistance_annotation.txt', self.name])
            
        else:
            self.__ifVerbose("Use SamTools, GATK, or Freebayes to annotate the final VCF.")
    
    def runLineage(self):
        """ Run lineage Analysis """
        self.__ifVerbose("Running Lineage Analysis")
        self.__final_annotation = self.fOut + "/" + self.name +'_Final_annotation.txt'
        self.__CallCommand(['lineage parsing', self.fOut + "/" + self.name +'_Lineage.txt'],
                              ['python', self.__lineage_parser, self.__lineages, self.__final_annotation, self.__lineage, self.name])
        count1 = 0
        count2 = 0
        fh1 = open(self.__lineage,'r')
        for line in fh1:
            lined = line.rstrip("\r\n")
            i = datetime.now()
            if "No Informative SNPs" in lined:
                self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
            elif "no precise lineage" in lined:
                self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
            elif "no concordance" in lined:
                self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "no clear lineage classification\n")
        fh1.close()
        fh2 = open(self.fOut + "/" + self.name +'_Resistance_Final_annotation.txt','r')
        for lines in fh2:
            lined = lines.rstrip("\r\n").split("\t")
            if lined[16] == "rrs":
               count1 += 1
            elif lined[16] == "rrl":
               count2 += 1
        if count1 > 5 or count2 > 5 :
           self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "mixed species suspected\n")
        fh2.close()
      
    def runCoverage(self):
        """ Run Genome Coverage Statistics """
        cov = ""
        self.__ifVerbose("Running Genome Coverage Statistics")
        samDir = self.outdir + "/SamTools"
        i = datetime.now()
        self.__CallCommand(['coverage estimator', self.fOut + "/" + self.name + '_Coverage.txt'],
                            ['python', self.__coverage_estimator, samDir + '/coverage.txt'])
        self.__CallCommand(['resistance region coverage estimator', self.fOut + "/" + self.name + '_Resistance_Region_Coverage.txt'],
                            ['python', self.__resis_parser, samDir + '/bed_sorted_coverage.txt', samDir + '/coverage.txt'])
        fh2 = open(self.fOut + "/" + self.name + '_Coverage.txt','r')
        for line in fh2:
            if line.startswith("Average"):
               cov_str = line.split(":")
               cov = cov_str[1].strip(" ")
        if int(cov) < 10:
           self.__logFH2.write(i.strftime('%Y/%m/%d %H:%M:%S') + "\t" + "Input:" + "\t" + self.name + "\t" + "low genome coverage depth\n")
        fh2.close()                  
        
    def cleanUp(self):
        """ Clean up the temporary files, and move them to a proper folder. """
        self.__CallCommand('rm', ['rm', '-r', self.outdir])
        self.__CallCommand('rm', ['rm',  self.input])
        self.__CallCommand('rm', ['rm',  self.fOut + "/" + self.name +'_annotation.txt'])
        self.__CallCommand('rm', ['rm',  self.fOut + "/" + self.name +'_Resistance_annotation.txt'])
        self.__CallCommand('rm', ['rm',  self.__finalBam])
        self.__CallCommand('rm', ['rm',  self.fOut + '/'+ self.name + '_sdrcsm.bai'])
        self.__CallCommand('rm', ['rm',  self.kraken + "/kraken.txt"])
       
        
        if self.paired:
           self.__CallCommand('rm', ['rm',  self.input2])
        valiOut = self.fOut + "/validation"
        self.__CallCommand('rm', ['rm', '-r', self.prinseq])

    def __ifVerbose(self, msg):
        """ If verbose print a given message. """
        if self.verbose: print msg
