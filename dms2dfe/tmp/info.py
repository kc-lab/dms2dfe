#!usr/bin/python

# source file for dms2dfe's configuration 

host='coli' #Host name for assigning codon table [coli | yeast | sapiens]
Ni_cutoff='3' #Cut off for frequency per mutant
Q_cutoff='30' #Cut off for Phred score quality
norm_type='wild' #Type of normalization across samples [wild: wrt wild type | syn : wrt synonymous mutations | none : fold change serves as fitness]
alignment_type='loc' #Alignment type [loc:local | glob:global]
cores='12' #Number of cores to be used
fsta_fh='miseq2/gmr_wt.fasta' #Optional: Path to reference fasta file
clips='10 167' #nan
pdb_fh='miseq2/1BO4.pdb' #Optional: Path to pdb file
active_sites='147' #Optional: residue numbers of active sites (space delimited) eg. 68<SPACE>192
dssp_fh='dms2dfe_dependencies/dssp-2.0.4-linux-amd64' #Optional: path to dssp module (dependencies)
trimmomatic_fh='dms2dfe_dependencies/Trimmomatic-0.33/trimmomatic-0.33.jar' #Optional: path to trimmomatic source (.jar) file (dependencies)
bowtie2_fh='dms2dfe_dependencies/bowtie2-2.2.1/bowtie2' #Optional: path to bowtie2 source file
samtools_fh='dms2dfe_dependencies/samtools-0.1.18/samtools' #Optional: path to samtools source file
cctmr='nan' #Optional: if reference sequence is concatamer (space delimited) eg. 1<SPACE>265<SPACE>268<SPACE>532
clustalo_fh='dms2dfe_dependencies/clustalo-1.2.2-Ubuntu-x86_64' #nan
msms_fh='dms2dfe_dependencies/msms/msms.x86_64Linux2.2.6.1' #nan
rate4site_fh='dms2dfe_dependencies/rate4site/rate4site-3.0.0/src/rate4site/rate4site' #nan
fsta_id='gmr-wt' #nan
fsta_seq='ATGTTACGCAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAGGTGGCTCAAGTATGGGCATCATTCGCACATGTAGGCTCGGCCCTGACCAAGTCAAATCCATGCGGGCTGCTCTTGATCTTTTCGGTCGTGAGTTCGGAGACGTAGCCACCTACTCCCAACATCAGCCGGACTCCGATTACCTCGGGAACTTGCTCCGTAGTAAGACATTCATCGCGCTTGCTGCCTTCGACCAAGAAGCGGTTGTTGGCGCTCTCGCGGCTTACGTTCTGCCCAAGTTTGAGCAGCCGCGTAGTGAGATCTATATCTATGATCTCGCAGTCTCCGGCGAGCACCGGAGGCAGGGCATTGCCACCGCGCTCATCAATCTCCTCAAGCATGAGGCCAACGCGCTTGGTGCTTATGTGATCTACGTGCAAGCAGATTACGGTGACGATCCCGCAGTGGCTCTCTATACAAAGTTGGGCATACGGGAAGAAGTGATGCACTTTGATATCGACCCAAGTACCGCCACCTAA' #nan
fsta_len='534' #nan
