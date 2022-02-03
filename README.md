# GBS data analysis codes
## Python and Bash scripts to run data analysis on Gene data

### Template files can be changed to run gene analysis on new datasets

## template_denovopipe.sh
Will take fastq.gz (gzip compressed) files and create:
  Multi Locus Sequence Type files
  Capsular Serotype
  Antimicrobial resistance gene detection
  Corrected fastq files (gzip compressed)
  de-novo .fasta files

## template_treepipe.sh
Will take a selection of isolates fastq files, uses a REFERENCE file and an EXCLUSION file, and create:
DEPENDENT on template_bacdate.R and template_hierbaps.R
  VCF Type files against a reference genome to call SNPs
  SNP log file with all positions marked
  Pseudo genomes from SNP and complete REFERENCE genome
  Hierbaps grouping using SNP (overly simplistically: grouping by distance)
  Phylogenetic tree file using RAxML
  Clonal Frame ML from RAxML tree
  Bacdating from CFML output
