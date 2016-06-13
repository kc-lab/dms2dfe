#!usr/bin/python

# source file for dms2dfe's configuration 

host='sapiens' #Host name for assigning codon table [coli | yeast | sapiens]
Ni_cutoff='8' #Cut off for frequency per mutant
Q_cutoff='30' #Cut off for Phred score quality
norm_type='wild' #Type of normalization across samples [wild: wrt wild type | syn : wrt synonymous mutations | none : fold change serves as fitness]
alignment_type='loc' #Alignment type [loc:local | glob:global]
cores='3' #Number of cores to be used
fsta_fh='APH2_Melnikov_et_al_2014/aph_wt_nt_cctmr.fasta' #Optional: Path to reference fasta file
pdb_fh='APH2_Melnikov_et_al_2014/1ND4.pdb' #Optional: Path to pdb file
active_sites='nan' #Optional: residue numbers of active sites (space delimited) eg. 68<SPACE>192
dssp_fh='dms2dfe_dependencies/dssp-2.0.4-linux-amd64' #Optional: path to dssp module (dependencies)
trimmomatic_fh='dms2dfe_dependencies/Trimmomatic-0.33/trimmomatic-0.33.jar' #Optional: path to trimmomatic source (.jar) file (dependencies)
bowtie2_fh='dms2dfe_dependencies/bowtie2-2.2.1/bowtie2' #Optional: path to bowtie2 source file
samtools_fh='dms2dfe_dependencies/samtools-0.1.18/samtools' #Optional: path to samtools source file
cctmr='1 265 268 532' #Optional: if reference sequence is concatamer (space delimited) eg. 1<SPACE>265<SPACE>268<SPACE>532
