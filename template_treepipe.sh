#!/bin/bash

### modify these variables ###
folder="TEMP" #  >>> enter the name you want to give folder, will also dictate output names <<< #
REFERENCE="REFISOLATE" # >>> enter the reference isolate name (should match .fa file name) <<< #

### This will create the vcf files to be used to call SNP ###
mkdir "$folder"
cd "$folder"
mkdir snpswap
mkdir vcf-filtered
mkdir vcf
cp ../"$REFERENCE".fa . &&
cp ../exclusion"$REFERENCE".txt ./"$folder"exclude.txt &&

 # >>> enter strain or accession numbers in place of ISOLATES <<< #
for X in ISOLATES ; do
	if [ ! -f /home/user/Documents/GBS_typeIV/vcf-filtered/"$X"_dp15q10f75-decomp.vcf
	then
		cp /home/user/Documents/GBS_typeIV/corr/"$X"_R1-corr.fastq.gz . &&
		cp /home/user/Documents/GBS_typeIV/corr/"$X"_R2-corr.fastq.gz . &&
		#SMALT code #
			#build ref index - do not need if already done - only needs done once, see above#
			#map reads with help of reference genome#
			smalt map -f sam -n 30 -o "$X"_map.sam ../"$REFERENCE"_ref_k14s8 "$X"_R1-corr.fastq.gz "$X"_R2-corr.fastq.gz &&
		#Freebayes code #
			#convert sam to bam using samtools#
			samtools view -S -b -@ 10 "$X"_map.sam >"$X".bam &&
			#sort bam file#
			samtools sort -@ 10 "$X".bam > "$X"sort.bam &&
			#freebayes call polys against reference genome .fa#
			freebayes -f "$REFERENCE".fa -p 1 "$X"sort.bam >"$X".vcf &&
			#filter VCF#
			bcftools view -i "MIN(INFO/DP)>15 & QUAL>10" "$X".vcf -o "$X"_dp15q10.vcf &&
			bcftools view -i "MIN(INFO/AO/INFO/DP)>0.75" "$X"_dp15q10.vcf -o "$X"_dp15q10f75.vcf &&
			vt decompose_blocksub "$X"_dp15q10f75.vcf -o "$X"_dp15q10f75-decomp.vcf &&
		### Cleanup ###
		rm "$X"_*.fastq.gz &&
		rm "$X"*.bam &&
		rm "$X"*.sam &&
		mv "$X".vcf ./vcf/ &&
		rm "$X"_dp15q10.vcf &&
		rm "$X"_dp15q10f75.vcf &&
		echo "Finished vcf for $X"
	else
		cp /home/user/Documents/GBS_typeIV/vcf-filtered/"$X"_dp15q10f75-decomp.vcf .
	fi
		### Make snpswap files ###
	if [ ! -f /home/user/Documents/GBS_typeIV/snpswap/"$X"-swap.fasta.gz
	then
		prephix.py -batchid "$X" "$X"_dp15q10f75-decomp.vcf &&
		snp_swapper.pl "$REFERENCE".fa "$X".snp &&
		mv "$REFERENCE".fa.swapped "$X"-swap.fa &&
		mv "$X"-swap.fa ./snpswap/ &&
		rm "$X"*.log &&
		rm "$X".indel &&
		rm "$X".ref &&
		rm "$X".snp
	else
		cp /home/user/Documents/GBS_typeIV/snpswap/"$X"-swap.fasta.gz ./"$X"-swap.fa.gz &&
		gunzip "$X"-swap.fa.gz &&
		mv "$X"-swap.fa ./snpswap/
	fi
		echo "Finished SNP-swap for $X"
		### Move vcf-filtered files into a separate folder ###
		mv "$X"_dp15q10f75-decomp.vcf ./vcf-filtered/ &&
		echo "done with $X"
done

### merge all the SNP-swap files into one large file to feed into CFML ###
	cp ./"$REFERENCE".fa merged"$folder"-swap.fa &&
	sed -i "1s/.*/>$REFERENCE/" merged"$folder"-swap.fa &&
	cat ./snpswap/*-swap.fa >> merged"$folder"-swap.fa &&
	sed -i "s/_dp15q10f75-decomp.vcf//g" merged"$folder"-swap.fa &&
	gzip ./snpswap/*-swap.fa &&
	mv ./snpswap/*-swap.fa.gz ../snpswap/ &&
	rm -d snpswap

### Remove Indels, and build SNP only file to use for tree building ###
	cd ./vcf-filtered/ &&
	prephix.py -batchid "$folder"_exclude -exclude ../"$folder"exclude.txt *.vcf &&
	phrecon.py "$folder"_exclude.ref "$folder"_exclude.snp &&
	# Replaces the filename with the Isolate name in the output file
	mv "$folder"_exclude.snp.reconstructed "$folder"_exclude.fa &&
	sed -i "s/_dp15q10f75-decomp.vcf//g" "$folder"_exclude.fa &&
	sed -i "s/REF/$REFERENCE/g" "$folder"_exclude.fa &&

### Run hierpabs R script from template_hierbaps.R file ###
	cp ../../template_hierbaps.R ./"$folder"_hierbaps.R &&
	sed -i "s/FOLDER/$folder/g" "$folder"_hierbaps.R &&
	echo " *** hierbaps is running ..."
	R CMD BATCH "$folder"_hierbaps.R &&
	rm "$folder"_hierbaps.Rout &&

### Build tree using RAxML ###
	raxmlHPC-PTHREADS-AVX -s "$folder"_exclude.fa -n "$folder"_exclude -m GTRCAT -f a -x 123 -N 100 -p 456 -T 32 &&
	mv RAxML* .. &&
	mv ./*.vcf ../../vcf-filtered/ &&
	cd .. &&

### Use RAxML tree for ClonalFrame and get output PDF image ###
	ClonalFrameML ./RAxML_bestTree."$folder"_exclude ./merged"$folder"-swap.fa CFML"$folder" -emsim 100 &&
	Rscript ../../../gen-soft/ClonalFrameML/src/cfml_results.R CFML"$folder" &&
	echo "Done with clonal frame, run bacdate on R" &&

### Run bacdate R scrip from template_bacdate.R file ###
	cp ../template_bacdate.R ./"$folder"_bacdate.R &&
	sed -i "s/FOLDER/$folder/g" "$folder"_bacdate.R &&
	echo " *** Bacdating is running ..." &&
	R CMD BATCH "$folder"_bacdate.R &&
	echo " *** output with Model comparison can be found at the end of $folder_hierbaps.Rout file" &&
	echo "Done with bacdating PDF image results can be found in $folder"
cd ..
