#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_seq_files``
================================
"""
# import sys
from os.path import splitext,exists,basename,dirname
import pandas as pd
from glob import glob
import subprocess
from Bio import SeqIO
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'

def contains_barcode(barcode_seq, read_seq) : 
    """
    This searches for barcode sequence in given sequence read \
    in windows of +/- 5 nucleotide upstream and downstream.
    
    :param barcode_seq: sequence of barcode (string).
    :param barcode_seq: target sequence (string).
    """
    for i in range(0, 5) : # (12,17) : #start looking around the end of the 15th base (4 char tag + 11 char mid)
        if barcode_seq in str(read_seq[i:]) :
            return True

def Rs2mergedandjoined(fastq_R1_read,fastq_R2_read,fastq_merged_f,fastq_joined_f):
    """
    This merges and joins R1 and R2 reads.
    
    merged files (*_dplxd_merged.fastq) have the R1 and R2 mixed together without joining. 
    eg. R1 R2 R1 R2 ..

    merged files (*_dplxd_joined.fastq) have the R1 and R2 joined together.
    eg.	R1R2 R1R2 ..
            
    :param fastq_R1_read: sequence of R1 read of .fastq (str).
    :param fastq_R2_read: sequence of R2 read of .fastq (str).
    :param fastq_merged_f: file object for output merged .fastq.  
    :param fastq_joined_f: file object for output joined .fastq.
    """
    from dms2dfe.lib.convert_seq import revcom,revers
    #MERGED
    fastq_R1_seq=str(fastq_R1_read.seq)
    fastq_R1_qua=''.join([chr(x + 33) for x in fastq_R1_read.letter_annotations['phred_quality']])
    fastq_writer(fastq_merged_f,fastq_R1_read.id+"R1",fastq_R1_seq,fastq_R1_qua)
    fastq_R2_seq=str(fastq_R2_read.seq)
    fastq_R2_qua=''.join([chr(x + 33) for x in fastq_R2_read.letter_annotations['phred_quality']])
    fastq_writer(fastq_merged_f,fastq_R2_read.id+"R2",fastq_R2_seq,fastq_R2_qua)

    #JOINED
    fastq_seq_joined=fastq_R1_seq+revcom(fastq_R2_seq)
    fastq_qua_joined=fastq_R1_qua+revers(fastq_R2_qua)  
    # fastq_seq_joined,fastq_qua_joined=reorder_ATGtoTAA(fastq_seq_joined,fastq_qua_joined) # # make reads start from ATG
    fastq_writer(fastq_joined_f,fastq_R1_read.id,fastq_seq_joined,fastq_qua_joined)

def fastq2dplx(fastq_R1_fh,fastq_R2_fh,barcode_R1s,barcode_R2s,fastq_fns) :
    """
    This demultiplexes fastq.
    
    :param fastq_R1_fh: path to R1 .fastq file.
    :param fastq_R2_fh: path to R2 .fastq file.
    :param barcode_R1s: barcodes to demultiplex forward reads.
    :param barcode_R1s: barcodes to demultiplex reverse reads.
    :param fastq_fns: file names of demultiplexed .fastq files (saved in same folder as the non-multiplexed input file).
    """
    
    from dms2dfe.lib.convert_seq import cmplmt,revers,revcom

    fastq_R1_reads=SeqIO.parse(fastq_R1_fh, 'fastq')
    fastq_R2_reads=SeqIO.parse(fastq_R2_fh,"fastq")

    fastq_dh=dirname(fastq_R1_fh)
    for fastq_fn in fastq_fns:
        fastq_merged_fh="%s/%s_dplxd_merged.qcd.fastq" % (fastq_dh,fastq_fn)
        fastq_joined_fh="%s/%s_dplxd_joined.qcd.fastq" % (fastq_dh,fastq_fn)
        fastq_merged_f =open(fastq_merged_fh,'w')
        fastq_joined_f =open(fastq_joined_fh,'w')

        unresolved_merged_fh="%s_unresolved_merged.qcd.fastq" % (fastq_R1_fh)
        unresolved_joined_fh="%s_unresolved_joined.qcd.fastq" % (fastq_R1_fh)
        unresolved_merged_f =open(fastq_merged_fh,'w')
        unresolved_joined_f =open(fastq_joined_fh,'w')

        for fastq_R1_read in fastq_R1_reads :
                fastq_R2_read =    fastq_R2_reads.next() #for fastq_R2_read in fastq_R2_reads :
                if fastq_R1_read.id in fastq_R2_read.id : # check contiguiency(?)
                    # print ">>> STATUS    : %s" % fastq_R1_read.id
                    for i in range(len(barcode_R1s)): #zip(barcode_R1s,barcode_R2s,fastq_fns):
                        barcode_R1=barcode_R1s[i]
                        barcode_R2=barcode_R2s[i]
                        fastq_fn    =fastq_fns[i]
                        # print "%s\n%s\n" % (barcode_R1,fastq_R1_read.seq)
                        if contains_barcode(barcode_R1 , str(fastq_R1_read.seq) ) and \
                        contains_barcode(barcode_R2 , str(fastq_R2_read.seq)) :
                            flag=1
                            break
                        elif contains_barcode(barcode_R1 , str(fastq_R2_read.seq) ) and \
                        contains_barcode(barcode_R2 , str(fastq_R1_read.seq)) :
                            flag=1
                            break
                        else:
                            flag=0
                    if flag==1:
                        Rs2mergedandjoined(fastq_R1_read,fastq_R2_read,fastq_merged_f,fastq_joined_f)
                    else:
                        Rs2mergedandjoined(fastq_R1_read,fastq_R2_read,unresolved_merged_f,unresolved_joined_f)
        fastq_merged_f.close()
        fastq_joined_f.close()
        unresolved_merged_f.close()
        unresolved_joined_f.close()
    fastq_R1_reads.close()
    fastq_R2_reads.close()

def fastq_writer(otpt_f,read_id,read_seq,read_qua):
    """
    This appends a sequence reads in fastq format. 
    
    :param otpt_f: file object of output file.
    :param read_id: ID of fastq file (str).
    :param read_seq: sequence (str).
    :param read_qua: quality scores (str).
    """  
    otpt_f.write("@"+read_id+"\n")
    otpt_f.write(read_seq+"\n")
    otpt_f.write("+\n")
    otpt_f.write(read_qua+"\n")
    
def getusablefastqs_list(prj_dh):
    """
    Checks sanity of fastq files mentioned in configuration excel file.
    
    :param prj_dh: path to project directory.
    :returns fastq_fhs: list of paths to .fastq files chich can be used downstream. 
    """
    lbls=pd.read_csv(prj_dh+'/cfg/lbls')
    lbls=lbls.set_index('varname')
    fastq_fhs=[]
    for lbli,lbl in lbls.iterrows() :
        fastq_R1_fh= str(lbl['fhs_1'])
        fastq_R2_fh= str(lbl['fhs_2'])
        ext=str(splitext(fastq_R1_fh)[1])
        if ("fastq" in ext) or ("fq" in ext):
            if pd.isnull(fastq_R2_fh): # single ended or sbam
                if not "bam" in ext:
                    fastq_R2_fh=""
                else:
                    logging.warning("skipping %s " % fastq_R1_fh)
            if not any(list("s.bam" in fh for fh in glob(fastq_R1_fh+"*"))):
                fastq_fhs.append([fastq_R1_fh,fastq_R2_fh])                     
    return fastq_fhs

def fastq2qcd(fastq_fhs_list,trimmomatic_fh): #[fastq_R1_fh,fastq_R2_fh]
    """
    Uses Trimmomatic to filter and trim low quality reads. 
    
    :param fastq_fhs_list: list of input .fastq files.
    :param trimmomatic_fh: path to trimmomatic source file.
    """ 
    fastq_R1_fh=fastq_fhs_list[0]
    fastq_R2_fh=fastq_fhs_list[1]
    if not ".qcd.fastq" in fastq_R1_fh:
        if exists(fastq_R1_fh) and (not exists(fastq_R1_fh+".qcd.fastq")):
            if trimmomatic_fh!='nan':
                log_fh="%s.fastq2qcd.log" % fastq_R1_fh
                log_f = open(log_fh,'a')

                if fastq_R2_fh and exists(fastq_R2_fh): #PE
                    bashCommand = "java -jar %s PE -threads 1 -phred33 -trimlog %s.qcd.fastq.log %s %s %s.qcd.fastq %s.qcd.unpaired %s.qcd.fastq %s.qcd.unpaired LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" \
                    % (trimmomatic_fh,fastq_R1_fh,fastq_R1_fh,fastq_R2_fh,fastq_R1_fh,fastq_R1_fh,fastq_R2_fh,fastq_R2_fh)
                    subprocess.call(bashCommand.split(),shell=True,stdout=log_f, stderr=subprocess.STDOUT)
                else:                
                    bashCommand = "java -jar %s SE -threads 1 -phred33 -trimlog %s.qcd.fastq.log %s %s.qcd.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" \
                    % (trimmomatic_fh,fastq_R1_fh,fastq_R1_fh,fastq_R1_fh)
                    # print bashCommand
                    subprocess.call(bashCommand.split(),shell=True,stdout=log_f, stderr=subprocess.STDOUT)
                log_f.close()
            else:
                sys.exit()
                logging.error("trimmomatic_fh not defined")     
        else:
            logging.info("skipping: %s" % basename(fastq_R1_fh))     
    else:
        logging.info("already processed: %s" % basename(fastq_R1_fh))     
                   
def qcd2sbam(fastq_fhs_list,fsta_fh,alignment_type,bt2_ref_fh,bowtie2_fh,samtools_fh): #[fastq_R1_fh,fastq_R2_fh]
    """
    This aligns the fastq files using bowtie2.
    The resultant .sam files are converted to sorted .bam files (.sbam) for downstream ana1_sam2mutmat module. 
    
    :param fastq_fhs_list: list of input .fastq files.
    :param fsta_fh: path to reference .fasta file.
    :param alignment_type: type of alignment ['loc' : local | 'glob' : global]
    :param bt2_ref_fh: path to bowtie2 reference file.
    :param bowtie2_fh: path to bowtie2 source file.
    :param samtools_fh: path to samtools source file.
    """
    # global alignment_type
    if not ".qcd.fastq" in fastq_fhs_list[0]:
        fastq_R1_fh=fastq_fhs_list[0]+".qcd.fastq"
        fastq_R2_fh=fastq_fhs_list[1]+".qcd.fastq"
    else:
        fastq_R1_fh=fastq_fhs_list[0]
        fastq_R2_fh=fastq_fhs_list[1]

    if "glob" in alignment_type: #global_alignment
        alignment_type=""
    else : #local_alignment
        alignment_type="-local" 
    if exists(fastq_R1_fh) :
        if not exists(fastq_R1_fh+".sam.s.bam"):
            if fastq_R2_fh and exists(fastq_R2_fh): # paired end
                bashCommand = "%s -p 1 --very-sensitive%s --no-discordant --no-mixed %s -1 %s -2 %s -S %s.sam 2> %s.sam.log" % (bowtie2_fh,alignment_type,bt2_ref_fh,fastq_R1_fh,fastq_R2_fh,fastq_R1_fh,fastq_R1_fh)       
            elif exists(fastq_R1_fh): # single end
                bashCommand = "%s -p 1 --very-sensitive%s %s -q %s -S %s.sam 2> %s.samlog" % (bowtie2_fh,alignment_type,bt2_ref_fh,fastq_R1_fh,fastq_R1_fh,fastq_R1_fh)
            else:
                logging.error("can not find : %s or %s" % (basename(fastq_R1_fh),basename(fastq_R2_fh)))
            logging.info("processing: %s " % basename(fastq_R1_fh))

            log_fh="%s.qcd2sbam.log" % fastq_R1_fh
            log_f = open(log_fh,'a')

            if not exists(fastq_R1_fh+".sam"):
                subprocess.call(bashCommand.split(),shell=True,stdout=log_f, stderr=subprocess.STDOUT)
            else :
                logging.info("qcd2sbam : sam already done")
            if not exists(fastq_R1_fh+".sam.bam"): 
                bashCommand= "%s view -bS -T %s %s.sam -o %s.sam.bam" % (samtools_fh,fsta_fh,fastq_R1_fh,fastq_R1_fh)
                subprocess.call(bashCommand.split(),shell=True,stdout=log_f, stderr=subprocess.STDOUT)
            else :
                logging.info("qcd2sbam : bam already done")
            bashCommand= "%s sort %s.sam.bam %s.sam.s" % (samtools_fh,fastq_R1_fh,fastq_R1_fh)
            subprocess.call(bashCommand.split(),shell=True,stdout=log_f, stderr=subprocess.STDOUT)
            bashCommand= "%s index %s.sam.s.bam" % (samtools_fh,fastq_R1_fh)
            subprocess.call(bashCommand.split(),shell=True,stdout=log_f, stderr=subprocess.STDOUT) 
            log_f.close()
        else :  
            logging.info("already processed: %s " % basename(fastq_R1_fh))
    else :
        logging.error("can not find: %s " % basename(fastq_R1_fh))

def fasta_nts2prt(fsta_fh,host='coli',fsta_prt_fh=None):
    from dms2dfe.lib.convert_seq import cds2aas
    from Bio import SeqIO,Seq,SeqRecord
    from Bio.Alphabet import IUPAC
    
    with open(fsta_fh,'r') as fsta_data:
        for fsta_record in SeqIO.parse(fsta_data, "fasta") :
            fsta_id=fsta_record.id
            fsta_seq=str(fsta_record.seq)
            break
    if fsta_prt_fh==None:
        fsta_prt_fh="%s_prt%s" % (splitext(fsta_fh)[0],splitext(fsta_fh)[1])
    fsta_prt_f = open(fsta_prt_fh, "w")
    fsta_seq_prt=cds2aas(fsta_seq,host,stop_codon='*')
    fsta_data_prt = SeqRecord.SeqRecord(Seq.Seq(fsta_seq_prt,IUPAC.protein), id = splitext(basename(fsta_fh))[0]+'_prt', description='')
    SeqIO.write(fsta_data_prt, fsta_prt_f, "fasta")
    fsta_prt_f.close()
    return fsta_seq_prt

def cctmr_fasta2ref_fasta(fsta_fh,cctmr):
    from dms2dfe.lib.convert_seq import cds2aas
    from Bio import SeqIO,Seq,SeqRecord
    from Bio.Alphabet import IUPAC

    fsta_cctmr1_fh="%s_cctmr1.fasta" % (splitext(fsta_fh)[0])
    with open(fsta_fh,'r') as fsta_data:
        for fsta_record in SeqIO.parse(fsta_data, "fasta") :
            fsta_id=fsta_record.id
            fsta_seq=str(fsta_record.seq)
            fsta_cctmr1_seq=fsta_seq[(cctmr[0]-1)*3:(cctmr[1]-1)*3]
            break
    fsta_cctmr1_f = open(fsta_cctmr1_fh, "w")
    fsta_data = SeqRecord.SeqRecord(Seq.Seq(fsta_cctmr1_seq,IUPAC.ExtendedIUPACDNA), id = fsta_id, description='')
    SeqIO.write(fsta_data, fsta_cctmr1_f, "fasta")
    fsta_cctmr1_f.close()
    return fsta_cctmr1_fh

def fasta_writer(otpt_f,read_id,read_seq):
    """
    This appends a sequence reads in fastq format. 
    
    :param otpt_f: file object of output file.
    :param read_id: ID of fastq file (str).
    :param read_seq: sequence (str).
    :param read_qua: quality scores (str).
    """  
    otpt_f.write(">"+read_id+"\n")
    otpt_f.write(read_seq+"\n")
