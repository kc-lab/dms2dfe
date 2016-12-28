#!usr/bin/python

# source file for dms2dfe's configuration 

host='e_coli' #Host name for assigning codon table [coli | yeast | sapiens]
Ni_cutoff='8' #Cut off for frequency per mutant
Q_cutoff='30' #Cut off for Phred score quality
norm_type='none' #Type of normalization across samples [wild: wrt wild type | syn : wrt synonymous mutations | none : fold change serves as fitness]
alignment_type='loc' #Alignment type [loc:local | glob:global]
cores='8' #Number of cores to be used
clips='nan' #Optional: Clip upstream UPTO and downstream FROM codons (space delimited) rg. 10<SPACE>167
fsta_fh='b161227_modeller_conf_int/btem.fasta' #Optional: Path to reference fasta file
pdb_fh='b161227_modeller_conf_int/btem.pdb' #Optional: Path to pdb file
active_sites='68 166' #Optional: residue numbers of active sites (space delimited) eg. 68<SPACE>192
cctmr='nan' #Optional: if reference sequence is concatamer (space delimited) eg. 1<SPACE>265<SPACE>268<SPACE>532
trimmomatic_com='nan' #Optional: additional commands to pass to trmmomatic
bowtie2_com='nan' #Optional: additional commands to pass to bowtie2
dssp_fh='dms2dfe_dependencies/dssp-2.0.4-linux-amd64' #Optional: path to dssp module (dependencies)
trimmomatic_fh='dms2dfe_dependencies/Trimmomatic-0.33/trimmomatic-0.33.jar' #Optional: path to trimmomatic source (.jar) file (dependencies)
bowtie2_fh='dms2dfe_dependencies/bowtie2-2.2.1/bowtie2' #Optional: path to bowtie2 source file
samtools_fh='dms2dfe_dependencies/samtools-0.1.18/samtools' #Optional: path to samtools source file
clustalo_fh='dms2dfe_dependencies/clustalo-1.2.2-Ubuntu-x86_64' #Optional: path to clustal omega source file
msms_fh='dms2dfe_dependencies/msms/msms.x86_64Linux2.2.6.1' #Optional: path to MSMS source file (for calculation of residue depths)
rate4site_fh='dms2dfe_dependencies/rate4site/rate4site-3.0.0/src/rate4site/rate4site' #Optional: path to rate4site source file (for calculation of conservation scores)
fsta_id='btem' #nan
fsta_seq='TAATAAATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGCAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA' #nan
fsta_len='867' #nan
prt_seq='**MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW*' #nan
prj_dh='/home/kclabws1/Documents/propro/writ/prjs/1_dms_software/data/datasets/TEM1_Firnberg_et_al_2014/b161227_modeller_conf_int' #nan
