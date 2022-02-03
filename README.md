# GBS data analysis codes
## Python and Bash scripts to run data analysis on Gene data

### Template files can be changed to run gene analysis on new datasets

## template_denovopipe.sh
Will take fastq.gz (gzip compressed) files and create:\n
\t Multi Locus Sequence Type files\n
\t Capsular Serotype\n
\t Antimicrobial resistance gene detection\n
\t Corrected fastq files (gzip compressed)\n
\t de-novo .fasta files\n

## template_treepipe.sh
Will take a selection of isolates fastq files, uses a REFERENCE file and an EXCLUSION file, and create:\n
DEPENDENT on template_bacdate.R and template_hierbaps.R\n
\t VCF Type files against a reference genome to call SNPs\n
\t SNP log file with all positions marked\n
\t Pseudo genomes from SNP and complete REFERENCE genome\n
\t Hierbaps grouping using SNP (overly simplistically: grouping by distance)\n
\t Phylogenetic tree file using RAxML\n
\t Clonal Frame ML from RAxML tree\n
\t Bacdating from CFML output\n
