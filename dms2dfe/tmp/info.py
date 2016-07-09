#!usr/bin/python

# source file for dms2dfe's configuration 

host='e_coli' #Host name for assigning codon table [coli | yeast | sapiens]
Ni_cutoff='0' #Cut off for frequency per mutant
Q_cutoff='30' #Cut off for Phred score quality
norm_type='syn' #Type of normalization across samples [wild: wrt wild type | syn : wrt synonymous mutations | none : fold change serves as fitness]
fsta_fh='TEM1_Firnberg_et_al_2014/btem.fasta' #Path to reference fasta file
pdb_fh='TEM1_Firnberg_et_al_2014/btem.pdb' #Path to pdb file
alignment_type='loc' #Alignment type [loc:local | glob:global]
cores='3' #Number of cores to be used
active_sites='68 166' #Optional: residue numbers of active sites (space delimited) eg. 68<SPACE>192
dssp_fh='dms2dfe_dependencies/dssp-2.0.4-linux-amd64' #Optional: path to dssp module
trimmomatic_fh='dms2dfe_dependencies/Trimmomatic-0.33/trimmomatic-0.33.jar' #Optional: path to trimmomatic source (.jar) file
cctmr='nan' #Optional: if reference sequence is concatamer (space delimited) eg. 1<SPACE>265<SPACE>268<SPACE>532" 
bowtie2_fh='dms2dfe_dependencies/bowtie2-2.2.1/bowtie2' #Optional: path to bowtie2 source file
samtools_fh='dms2dfe_dependencies/samtools-0.1.18/samtools' #Optional: path to samtools source file
