# GBS data analysis codes
## Python and Bash scripts to run data analysis on Gene data
<p> These are some of the scripts I have written to perform data analysis on our Group B Streptococcus agalactiae. I am currently adapting these to be usable for any further bacterial gene analysis our lab might require.
</p>
<ol>
  <li>These scripts use several known genetic tools, including:<ul>
      <li>Strainseeker (similar genetic match)</li>
      <li>SPAdes (contig assembly)</li>
      <li>Rapberry (read quality)</li>
      <li>SRST2 (MLST and Capsular Serotype calling)</li>
      <li>RGI (AMR calling)</li>
      <li>RAxML (Phylogenetic tree building)</li>
      <li>Hierbaps (Statistical Grouping of samples by SNP differences)</li>
      <li>CFML (Phylogenetic tree building)</li>
      <li>Bacdating (Historical phylogenetic inference)</li>
    </ul></li>
</ol>


### Template files can be changed to run gene analysis on new datasets

## template_denovopipe.sh
<ol>
  <li>Will take fastq.gz (gzip compressed) files and create:<ul>
      <li>Multi Locus Sequence Type files Capsular Serotype</li>
      <li>Antimicrobial resistance gene detection</li>
      <li>Corrected fastq files (gzip compressed)</li>
      <li>de-novo .fasta files</li>
    </ul></li>
</ol>

## template_treepipe.sh
<ol>
  <li>Will take a selection of isolates fastq files, uses a REFERENCE file and an EXCLUSION file, and create: <ul>
      <li>VCF Type files against a reference genome to call SNPs</li>
      <li>SNP log file with all positions marked</li>
      <li>Pseudo genomes from SNP and complete REFERENCE genome</li>
      <li>Hierbaps grouping using SNP (overly simplistically: grouping by distance)</li>
      <li>Phylogenetic tree file using RAxML</li>
      <li>Clonal Frame ML from RAxML tree</li>
      <li>Bacdating from CFML output</li>
    </ul>*DEPENDENT on template_bacdate.R and template_hierbaps.R</li>
</ol>
