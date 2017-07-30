#!usr/bin/python

# source file for dms2dfe's configuration 

host='coli' #Host name for assigning codon table
Ni_cutoff='8' #Lower cut off for frequency per mutant
Q_cutoff='30' #Lower cut off for Phred score quality
transform_type='log' #Type of transformation of frequecies of mutants
norm_type='GLM' #Type of normalization of fold changes
alignment_type='loc' #Alignment type
cores='8' #Number of cores to be used
mut_type='single' #Whether input data is of single mutation or double mutations
rescaling='FALSE' #Optional: Position wise rescaling the fold changes by fold changes of synonymous (wild type) mutations
mut_subset='N' #Optional: Subset of mutations to be used for down-stream analysis
ml_input='FC' #Optional: Whether to use Preferential enrichments or Fitness scores for identification of molecular constraints
clips='nan' #Optional: Clip upstream UPTO and downstream FROM codons (space delimited) rg. 10<SPACE>167
fsta_fh='gmr_analysis/gmr_wt.fasta' #Optional: Path to reference fasta file
pdb_fh='gmr_analysis/1BO4.pdb' #Optional: Path to pdb file
active_sites='147' #Optional: residue numbers of active sites (space delimited) eg. 68<SPACE>192
cctmr='nan' #Optional: if reference sequence is concatamer (space delimited) eg. 1<SPACE>265<SPACE>268<SPACE>532
trimmomatic_com='nan' #Optional: additional commands to pass to trmmomatic
bowtie2_com='nan' #Optional: additional commands to pass to bowtie2
dssp_fh='dms2dfe_dependencies/dssp-2.0.4-linux-amd64' #Optional: path to dssp module (dependencies)
trimmomatic_fh='dms2dfe_dependencies/Trimmomatic-0.33/trimmomatic-0.33.jar' #Optional: path to trimmomatic source (.jar) file (dependencies)
bowtie2_fh='dms2dfe_dependencies/bowtie2-2.2.1/bowtie2' #Optional: path to bowtie2 source file
samtools_fh='dms2dfe_dependencies/samtools-0.1.20/samtools' #Optional: path to samtools source file
clustalo_fh='dms2dfe_dependencies/clustalo-1.2.2-Ubuntu-x86_64' #Optional: path to clustal omega source file
msms_fh='dms2dfe_dependencies/msms/msms.x86_64Linux2.2.6.1' #Optional: path to MSMS source file (for calculation of residue depths)
rate4site_fh='dms2dfe_dependencies/rate4site/rate4site-3.0.0/src/rate4site/rate4site' #Optional: path to rate4site source file (for calculation of conservation scores)
rscript_fh='/usr/bin/Rscript' #Optional: path to Rscript (for use of Deseq2) can be identified by executing command 'which R'
prj_dh='/home/smfret/test/ms_datasets-0.0.3/analysis/gmr_analysis' #nan
fsta_id='gmr-wt' #nan
fsta_seq='ATGTTACGCAGCAGCAACGATGTTACGCAGCAGGGCAGTCGCCCTAAAACAAAGTTAGGTGGCTCAAGTATGGGCATCATTCGCACATGTAGGCTCGGCCCTGACCAAGTCAAATCCATGCGGGCTGCTCTTGATCTTTTCGGTCGTGAGTTCGGAGACGTAGCCACCTACTCCCAACATCAGCCGGACTCCGATTACCTCGGGAACTTGCTCCGTAGTAAGACATTCATCGCGCTTGCTGCCTTCGACCAAGAAGCGGTTGTTGGCGCTCTCGCGGCTTACGTTCTGCCCAAGTTTGAGCAGCCGCGTAGTGAGATCTATATCTATGATCTCGCAGTCTCCGGCGAGCACCGGAGGCAGGGCATTGCCACCGCGCTCATCAATCTCCTCAAGCATGAGGCCAACGCGCTTGGTGCTTATGTGATCTACGTGCAAGCAGATTACGGTGACGATCCCGCAGTGGCTCTCTATACAAAGTTGGGCATACGGGAAGAAGTGATGCACTTTGATATCGACCCAAGTACCGCCACCTAA' #nan
fsta_len='534' #nan
prt_seq='MLRSSNDVTQQGSRPKTKLGGSSMGIIRTCRLGPDQVKSMRAALDLFGREFGDVATYSQHQPDSDYLGNLLRSKTFIALAAFDQEAVVGALAAYVLPKFEQPRSEIYIYDLAVSGEHRRQGIATALINLLKHEANALGAYVIYVQADYGDDPAVALYTKLGIREEVMHFDIDPSTATX' #nan
