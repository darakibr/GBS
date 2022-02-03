#!/bin/bash

### modify these variables ###
folder="TEMP" #  >>> enter the name you want to give folder, will also dictate output names <<< #
REFERENCE="REFISOLATE" # >>> enter the reference isolate name (should match .fa file name) <<< #

### This will create the vcf files to be used to call SNP ###
mkdir {corr,seeker,denovo,stats}
cd 

for f in  ; do #For XXXX insert strain or accession numbers#
		cp /home/user/Documents/GBS_typeIV/ngbstemp_fastq/"$f"_R1.fastq.gz .
		cp /home/user/Documents/GBS_typeIV/ngbstemp_fastq/"$f"_R2.fastq.gz .
		
		#trim R1 reads to 25k (100k lines) for strainseeker#
		gzip -d "$f"_R1.fastq.gz
		head -n 100000 "$f"_R1.fastq >"$f"_R1-25k.fastq
		gzip "$f"_R1.fastq
		
		#strainseeker#
		perl /home/user/gen-soft/strainseeker/seeker.pl -i "$f"_R1-25k.fastq -d /home/user/gen-soft/strainseeker/ss_db_w32_4324 -o "$f"_ss.txt
		
		### START SPADES de novo assembly; saves corrected reads for use in subsequent steps ###
		spades.py -t 28 -1 "$f"_R1.fastq.gz -2 "$f"_R2.fastq.gz -o "$f"_spades-de-novo 
		cd "$f"_spades-de-novo/corrected/
		mv "$f"_R1.fastq.00.0_0.cor.fastq.gz "$f"_R1-corr.fastq.gz
		mv "$f"_R2.fastq.00.0_0.cor.fastq.gz "$f"_R2-corr.fastq.gz
		mv "$f"_R1-corr.fastq.gz /home/user/Documents/GBS_typeIV/ngbstemp_fastq/corr/
		mv "$f"_R2-corr.fastq.gz /home/user/Documents/GBS_typeIV/ngbstemp_fastq/corr/
		cd ..
		mv contigs.fasta "$f"_contigs.fasta
		mv "$f"_contigs.fasta ..
		mv spades.log "$f"_spades.log
		mv "$f"_spades.log ../stats/
		cd ..
		
		### RENAME contigs ###
		awk '/^>/{print ">""'"$f"'""_contigs-" ++i; next}{print}' < "$f"_contigs.fasta > "$f"_contigs2.fasta
		
		### Remove contigs <4K ###
		perl /home/user/gen-soft/removesmalls.pl 1000 "$f"_contigs2.fasta >"$f"_contigs-1K.fasta
		
		#Raspberry for qc on corrected reads#
		raspberry /home/user/Documents/GBS_typeIV/ngbstemp_fastq/corr/"$f"* >"$f"_readqc.txt
		
		#cleanup#
		mv "$f"_ss.txt ./seeker/
		mv "$f"_contigs-1K.fasta ./denovo/
		mv "$f"_readqc.txt ./stats/
		rm -r "$f"_spades-de-novo
		rm *.fastq
		rm *.fastq.gz
		rm "$f"_contigs2.fasta
		rm "$f"_contigs.fasta
		
done
