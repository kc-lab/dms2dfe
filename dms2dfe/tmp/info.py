#!usr/bin/python

# source file for dms2dfe's configuration 

host='coli' #Host name for assigning codon table
Ni_cutoff='8' #Lower cut off for frequency per mutant
Q_cutoff='30' #Lower cut off for Phred score quality
transform_type='log' #Type of transformation of frequecies of mutants
norm_type='none' #Type of normalization of fold changes
alignment_type='loc' #Alignment type
cores='2' #Number of cores to be used
mut_type='single' #Whether input data is of single mutation or double mutations
rescaling='FALSE' #Optional: Position wise rescaling the fold changes by fold changes of synonymous (wild type) mutations
mut_subset='N' #Optional: Subset of mutations to be used for down-stream analysis
ml_input='Fi' #Optional: Whether to use Preferential enrichments or Fitness scores for identification of molecular constraints
clips='nan' #Optional: Clip upstream UPTO and downstream FROM codons (space delimited) rg. 10<SPACE>167
fsta_fh='found in project directory' #Optional: Path to reference fasta file
pdb_fh='found in project directory' #Optional: Path to pdb file
active_sites='nan' #Optional: residue numbers of active sites (space delimited) eg. 68<SPACE>192
cctmr='nan' #Optional: if reference sequence is concatamer (space delimited) eg. 1<SPACE>265<SPACE>268<SPACE>532
trimmomatic_com='nan' #Optional: additional commands to pass to trmmomatic
bowtie2_com='nan' #Optional: additional commands to pass to bowtie2
dssp_fh='nan' #Optional: path to dssp module (dependencies)
trimmomatic_fh='nan' #Optional: path to trimmomatic source (.jar) file (dependencies)
bowtie2_fh='nan' #Optional: path to bowtie2 source file
samtools_fh='nan' #Optional: path to samtools source file
clustalo_fh='nan' #Optional: path to clustal omega source file
msms_fh='nan' #Optional: path to MSMS source file (for calculation of residue depths)
rate4site_fh='nan' #Optional: path to rate4site source file (for calculation of conservation scores)
rscript_fh='~/anaconda/envs/dms2dfe/bin/Rscript' #Optional: path to Rscript (for use of Deseq2) can be identified by executing command 'which Rscript'
