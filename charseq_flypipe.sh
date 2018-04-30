#!/bin/bash

set -o errexit

#usage: bash charseq_flypipe_SE_stranded.sh [/path/to/R01.fastq.gz] [/desired/working/path] [# of cores to use] [dirty or cleanup] [split reads location to avoid reruning bridge splitter eg /07_09_2016_dashtest/Sample1]
#NOTE: you must have super-deduper, bowtie2, and samtools in your path and python 2.7 called by python with the tools os, re, sys, bz2, time, gzip, string, argparse, and bio SeqIO

SE_READS=$1
DATAPATH=$2
CORES=$3


if [ -z "$CORES" ]
then
	echo "USAGE: bash charseq_flypipe.sh [/path/to/R01.fastq.gz] [/complete/desired/output/path/] [# of cores to use] [dirty or cleanup]"
	echo "Where the first argument is the single ended gzipped fastq file"
	echo "Second is where you would like your output files"
	echo "Third is the number of cores you'd like to use [enter 1 if you don't want to parallelize]"
	echo "Fourth is optional: type cleanup to automatically gzip intermediate files after the run to save space"
	echo ""
	echo "Test usage: bash charseq_flypipe.sh /compelete/path-to/data/cl8-test.fastq.gz /outputdirectory 4"
	exit
fi

#sanity things to ensure we can find scripts on others' machines
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
SCRIPTDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

if ! type python
then
	echo "Python is not an executable in your path. Please install Python 2.7 and add it to your path."
	echo "https://www.python.org/downloads/"
	exit
fi
if ! type bowtie2
then
	echo "bowtie2 is not an executable in your path. Please install bowtie2 and add it to your path."
	echo "http://bowtie-bio.sourceforge.net/bowtie2/index.shtml"
	exit
fi
if ! type samtools
then
	echo "samtools is not an executable in your path. Please install samtools and add it to your path."
	echo "http://www.htslib.org/"
	exit
fi
if ! type super_deduper
then
	echo "super_deduper is not an executable in your path. Please install super_deduper and add it to your path."
	echo "https://github.com/dstreett/Super-Deduper"
	exit
fi
if ! $(python -c "import os" &> /dev/null); then
    echo "Lacking python package os"
fi
if ! $(python -c "import re" &> /dev/null); then
    echo "Lacking python package re"
fi
if ! $(python -c "import sys" &> /dev/null); then
    echo "Lacking python package sys"
fi
if ! $(python -c "import bz2" &> /dev/null); then
    echo "Lacking python package bz2"
fi
if ! $(python -c "import time" &> /dev/null); then
    echo "Lacking python package time"
fi
if ! $(python -c "import gzip" &> /dev/null); then
    echo "Lacking python package gzip"
fi
if ! $(python -c "import string" &> /dev/null); then
    echo "Lacking python package string"
fi
if ! $(python -c "import argparse" &> /dev/null); then
    echo "Lacking python package argparse"
fi
if ! $(python -c "from Bio import SeqIO" &> /dev/null); then
    echo "Lacking biopython package SeqIO"
fi
if ! $(python -c "from Bio import AlignIO" &> /dev/null); then
    echo "Lacking biopython package AlignIO"
fi
if ! $(python -c "from itertools import izip" &> /dev/null); then
    echo "Lacking python package itertools"
fi
if ! $(python -c "from optparse import OptionParser" &> /dev/null); then
    echo "Lacking python package optparse"
fi
if ! $(python -c "import os" &> /dev/null); then
    echo "Lacking python package python-levenshtein"
fi

samplefastq=${SE_READS##*/}
SAMPLE=$(echo $samplefastq | cut -f 1 -d '.')

if [ -d $DATAPATH ]
then
    echo "DATAPATH exists"
else
    mkdir $DATAPATH
fi

if [ -d $DATAPATH/$SAMPLE ] 
then
    echo "Sample directory exists."
else
    mkdir $DATAPATH/$SAMPLE
fi
cd $DATAPATH/$SAMPLE

echo "Summary statistics for "$SAMPLE > $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt


	if [ -s $DATAPATH/$SAMPLE/rmpcrdups/$SAMPLE._nodup_PE1.fastq.gz ]
	then
	    echo "PCR Duplicates have already been removed"
	else
	if [ -d $DATAPATH/$SAMPLE/rmpcrdups ] 
	then
   		echo "rmpcrdups directory exists."
	else
    	mkdir $DATAPATH/$SAMPLE/rmpcrdups
	fi

	cd $DATAPATH/$SAMPLE/rmpcrdups

	echo Removing PCR Duplicates $SAMPLE $(date +"%T")

	####CHANGE TO SUPER DEDUPER
	#python /raid1/lab/nteran/code/remove_pcrduplicates.py --SE $SE_READS --NAME $SAMPLE. >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	
	LINELENGTH=$(zcat $SE_READS | head -n 2 | tail -n 1 | wc -c)
	LLM1=$(( ${LINELENGTH} - 1 ))

	super_deduper -U $SE_READS -v -q -l $LLM1 -s 1 -g -q -p $SAMPLE.
	fi

	echo "Simple count of reads" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $SE_READS | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "Simple bridge count" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $SE_READS | grep -c 'ACCGGCGTCCAAG\|CTTGGACGCCGGT' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "Simple count of post-pcr-duplicate-removal-reads" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/rmpcrdups/$SAMPLE._nodup_PE1.fastq.gz | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "Simple bridge post-pcr-duplicate-removal-reads" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/rmpcrdups/$SAMPLE._nodup_PE1.fastq.gz | grep -c 'ACCGGCGTCCAAG\|CTTGGACGCCGGT' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt

	##TRIM
	#takes about 45 min for 8M reads

	if [ -s $DATAPATH/$SAMPLE/trim/$SAMPLE.trim.fastq.gz ]
	then
	    echo "Reads have already been trimmed"
	else
	if [ -d $DATAPATH/$SAMPLE/trim ] 
	then
   		echo "trim directory exists."
	else
    	mkdir $DATAPATH/$SAMPLE/trim
	fi
	cd $DATAPATH/$SAMPLE/trim
	java -jar $SCRIPTDIR/trimmomatic-0.35.jar SE -threads $CORES -phred33 -trimlog logFile.txt $DATAPATH/$SAMPLE/rmpcrdups/$SAMPLE._nodup_PE1.fastq.gz $DATAPATH/$SAMPLE/trim/$SAMPLE.trim.fastq.gz ILLUMINACLIP:$SCRIPTDIR/illuminaClipping_main.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	fi

	echo "Simple count of post-trim reads" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/trim/$SAMPLE.trim.fastq.gz | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "Simple bridge post-trim" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/trim/$SAMPLE.trim.fastq.gz | grep -c 'ACCGGCGTCCAAG\|CTTGGACGCCGGT' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt

	##Run char_bridge on SE reads
	if [ -s $DATAPATH/$SAMPLE/bridge/$SAMPLE.dna.bridgePE.fastq.gz ]
	then
	    echo "Bridge has already been found and reads split"
	else
	mkdir $DATAPATH/$SAMPLE/bridge
	cd $DATAPATH/$SAMPLE/bridge


	python $SCRIPTDIR/char_bridge_trackall.py --FASTQGZ $DATAPATH/$SAMPLE/trim/$SAMPLE.trim.fastq.gz --NAME $DATAPATH/$SAMPLE/bridge/$SAMPLE.
	fi

	echo "Reads read by char_bridge with some kind of bridge" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/bridge/$SAMPLE.bridgePE.bridgeposition.fastq.gz | wc -l >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "Reads that contain passed bridge filters" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/bridge/$SAMPLE.dna.bridgePE.fastq.gz | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "      Failed bridge filters: multiple bridges" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/bridge/$SAMPLE.bridgePE.dupbridge.fastq.gz | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "      Failed bridge filters: no bridge" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/bridge/$SAMPLE.bridgePE.nobridge.fastq.gz | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	echo "      Failed bridge filters: too short RNA or DNA" >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
	zcat $DATAPATH/$SAMPLE/bridge/$SAMPLE.bridgePE.tooshort.fastq.gz | grep -c '^+$' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt

	#set a variable to denote where the split reads are
	SPLIT_DNA=$DATAPATH/$SAMPLE/bridge/$SAMPLE.dna.bridgePE.fastq.gz
	SPLIT_RNA=$DATAPATH/$SAMPLE/bridge/$SAMPLE.rna.bridgePE.fastq.gz

if [ -s $DATAPATH/$SAMPLE/align/bam_rev/$SAMPLE.ChARrna.tRNA.bam ];
then
	    echo "BAM files already exist"
else
	    echo "Aligning bams" $(date +"%T")
	    if [ ! -d $DATAPATH/$SAMPLE/align/ ]; then
	    mkdir $DATAPATH/$SAMPLE/align/
		fi
	    if [ ! -d $DATAPATH/$SAMPLE/align/bam ]; then
	    mkdir $DATAPATH/$SAMPLE/align/bam
	    fi
	    if [ ! -d $DATAPATH/$SAMPLE/align/unal_fwd ]; then
	    mkdir $DATAPATH/$SAMPLE/align/unal_fwd
		fi
		if [ ! -d $DATAPATH/$SAMPLE/align/bam_rev ]; then
	    mkdir $DATAPATH/$SAMPLE/align/bam_rev
		fi
	cd $DATAPATH/$SAMPLE/align
	##Align DNA to genome
	##Bowtie2
	(bowtie2 -p$CORES -D 20 -R 3 -N 1 -L 18 -i S,1,0.50 -x $SCRIPTDIR/genomes/genome -U <(gunzip -c $SPLIT_DNA) | samtools view -bS - -o $DATAPATH/$SAMPLE/align/bam/$SAMPLE.ChARdna.bam) 2>> genome_statistics.txt
	
	(bowtie2 -p$CORES -D 20 -R 3 -N 1 -L 18 -i S,1,0.50 -x $SCRIPTDIR/genomes/genome -U <(gunzip -c $SPLIT_RNA) | samtools view -bS - -o $DATAPATH/$SAMPLE/align/bam/$SAMPLE.ChARrna.genome.bam) 2>> genome_statistics.txt

	cd $DATAPATH/$SAMPLE/align

	for genome in $SCRIPTDIR/genomes/*.fasta
	do
	transcriptome=$(echo $genome | cut -f 1,2 -d '.')
	transcriptome_id=$(echo $transcriptome | cut -f 3 -d '-')
	echo $transcriptome_id $SAMPLE >> transcriptome_statistics.txt
	(bowtie2 -p$CORES -D 20 -R 3 -N 1 -L 18 -i S,1,0.50 --norc --un-gz $DATAPATH/$SAMPLE/align/unal_fwd/$transcriptome_id.fq.gz -x $transcriptome -U <(gunzip -c $SPLIT_RNA) | samtools view -bS - -o $DATAPATH/$SAMPLE/align/bam/$SAMPLE.ChARrna.$transcriptome_id.bam) 2>> transcriptome_statistics.txt
	echo  >> transcriptome_statistics.txt
	done

	#try again with things that didn't align forward
	cd $DATAPATH/$SAMPLE/align
	for genome in $SCRIPTDIR/genomes/*.fasta
	do
	transcriptome=$(echo $genome | cut -f 1,2 -d '.')
	transcriptome_id=$(echo $transcriptome | cut -f 3 -d '-')
	echo $transcriptome_id $SAMPLE >> transcriptome_statistics.txt
	(bowtie2 -p$CORES -D 20 -R 3 -N 1 -L 18 -i S,1,0.50 --nofw -x $transcriptome -U <(gunzip -c $DATAPATH/$SAMPLE/align/unal_fwd/$transcriptome_id.fq.gz) | samtools view -bS - -o $DATAPATH/$SAMPLE/align/bam_rev/$SAMPLE.ChARrna.$transcriptome_id.bam) 2>> transcriptome_statistics.txt
	echo  >> transcriptome_statistics.txt
	done


fi

if [ -s $DATAPATH/$SAMPLE/align/sam/$SAMPLE.ChARdna.sam ];
then
    echo "SAM files already exist"
else
    echo "Producing SAMs" $(date +"%T")
    if [ -d $DATAPATH/$SAMPLE/align/sam/ ]
    then
    cd $DATAPATH/$SAMPLE/align/sam/
    else
    mkdir $DATAPATH/$SAMPLE/align/sam_rev/
    mkdir $DATAPATH/$SAMPLE/align/sam/
    cd $DATAPATH/$SAMPLE/align/sam/
    fi
##Bams back to sams
    for file in $DATAPATH/$SAMPLE/align/bam/*.bam
    do
        fileh=${file##*/}
        samtools view -F4 $file > $DATAPATH/$SAMPLE/align/sam/${fileh/.bam/.sam}
    done
    for file in $DATAPATH/$SAMPLE/align/bam_rev/*.bam
    do 
    fileh=${file##*/}
    samtools view -F4 $file > $DATAPATH/$SAMPLE/align/sam_rev/${fileh/.bam/.sam}
    done
fi 

if [ -d $DATAPATH/$SAMPLE/newo/simple/ ]
	then
	    echo "Simplified folders already exist."
	else
	mkdir $DATAPATH/$SAMPLE/newo/
	mkdir $DATAPATH/$SAMPLE/newo/simple/
	cd $DATAPATH/$SAMPLE/newo/simple/
	##turn the sams into simple .text
	#be sure they're space deliminated
	echo "Making simple text files of read info..."
	for file in $DATAPATH/$SAMPLE/align/sam/*sam;
	do
	fileh=${file##*/}
	name=$(echo $fileh | sed -r 's/(.[^.]+){1}$//g')
	awk '{if($2 == 0) sign="+"; else if ($2 == 16) sign="-"; else sign="."; print $1,$3,$4,$5,$6,length($10),sign}' $file > $DATAPATH/$SAMPLE/newo/simple/$name.simple.txt
	done

	for file in $DATAPATH/$SAMPLE/align/sam_rev/*sam;
	do
	fileh=${file##*/}
	name=$(echo $fileh | sed -r 's/(.[^.]+){1}$//g')
	awk '{if($2 == 0) sign="+"; else if ($2 == 16) sign="-"; else sign="."; print $1,$3,$4,$5,$6,length($10),sign}' $file > $DATAPATH/$SAMPLE/newo/simple/$name.simple.rev.txt
	done

	#from sams aligned to genome, adds DNA adjusted coordinate
	#order of *sizeadjusted.txt: uniqueID, chrom, pos, adj pos, MAPQ, CIGAR, read length, strand on genome
	echo "Running size_adjust to make tab delimited reads file with genome-coordinate DNA info..."
	python $SCRIPTDIR/size_adjust.py $DATAPATH/$SAMPLE/newo/simple/$SAMPLE.ChARdna.simple.txt $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdna.sizeadjusted.txt

	python $SCRIPTDIR/size_adjust.py $DATAPATH/$SAMPLE/newo/simple/$SAMPLE.ChARrna.genome.simple.txt $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARrna.genome.sizeadjusted.txt
fi

if [ -d $DATAPATH/$SAMPLE/newo/matched/ ]
	then
	    echo "Matched folders already exist."
	else
	mkdir $DATAPATH/$SAMPLE/newo/matched
	cd $DATAPATH/$SAMPLE/newo/matched
	##Python script for add_DNA_to_transcriptome.py [rna.simple] [dna.simple] [output.simple]
	#matches DNA coordinate with RNA coordinate
	#order of *ChARdnaRna.txt: uniqueID, DNAchrom, DNApos, DNAadj pos, DNAMAPQ, DNACIGAR, DNA read length, DNA strand on genome, RNA fbgn/tr, RNApos, RNAMAPQ, RNACIGAR, RNA length of read, RNA strand on trascriptome
	echo "Adding DNA to transcriptome to make ChARdnaRna.txt files..."

	for file in $DATAPATH/$SAMPLE/newo/simple/*rna*.simple.txt
	do
	fileh=${file##*/}
	name=$(echo $fileh | sed -r 's/(.[^.]+){2}$//g')
	echo $name >> add_DNA_to_transcriptome_statistics.txt
	echo "forward RNA reads:" >> add_DNA_to_transcriptome_statistics.txt
	wc -l $file >> add_DNA_to_transcriptome_statistics.txt
	python $SCRIPTDIR/add_DNA_to_transcriptome.py $file $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdna.sizeadjusted.txt $DATAPATH/$SAMPLE/newo/matched/$name.ChARdnaRna.txt >> add_DNA_to_transcriptome_statistics.txt
	done

	for file in $DATAPATH/$SAMPLE/newo/simple/*rna*.simple.rev.txt
	do
	fileh=${file##*/}
	name=$(echo $fileh | sed -r 's/(.[^.]+){3}$//g')
	echo $name >> add_DNA_to_transcriptome_statistics.txt
	echo "reverse RNA reads:" >> add_DNA_to_transcriptome_statistics.txt
	wc -l $file >> add_DNA_to_transcriptome_statistics.txt
	python $SCRIPTDIR/add_DNA_to_transcriptome.py $file $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdna.sizeadjusted.txt $DATAPATH/$SAMPLE/newo/matched/$name.ChARdnaRna.rev.txt >> add_DNA_to_transcriptome_statistics.txt
	done

fi

if [ -d $DATAPATH/$SAMPLE/newo/matchedinfosize/ ]
	then
	    echo "Matched folders already exist."
	else
	#######
	## Python script to combine fasta info and simple files
	#usage: python add_gene_info_to_transciptome.py *ChARdnaRna.txt /genomes/dmel-all-$transcriptome-coords.dict ouputfile.info.txt
	#adds transcriptome specific information to ChARdnaRna.txt
	#order of info.txt: uniqueID, DNAchrom, DNApos, DNAadj pos, DNAMAPQ, DNACIGAR, DNA read length, DNA strand on genome, RNA fbgn/tr, RNApos, RNAMAPQ, RNACIGAR, RNA length of read, RNA strand on trascriptome, RNA type, RNA locfield, RNAchrom, RNATSS, RNATTS, RNA strand, RNA splice sites, RNA name, RNA parentgn, RNA parent, RNA length of transcript
	#rna aligned to genome is skipped because no transcriptome info 
	mkdir $DATAPATH/$SAMPLE/newo/matchedinfosize
	mkdir $DATAPATH/$SAMPLE/newo/matchedinfo/
	echo "Making infosize files by looking up transcriptome information in coords.dict file..."
	cd $DATAPATH/$SAMPLE/newo/matchedinfo/
	for file in $DATAPATH/$SAMPLE/newo/matched/*.ChARdnaRna.txt
	do
	fileh=${file##*/}
	transcriptome=$(echo $fileh | cut -f 3 -d '.')
	name=$(echo $fileh | sed -r 's/(.[^.]+){1}$//g')
	echo $name
	echo $transcriptome
	if [ "$transcriptome" = "genome" ];
	then
	        continue
	fi
	python $SCRIPTDIR/add_gene_info_to_transcriptome.py $file $SCRIPTDIR/dmel_coords.dict/dmel-all-$transcriptome-r5.57.coords.dict $DATAPATH/$SAMPLE/newo/matchedinfo/$name.info.txt
	done

	for file in $DATAPATH/$SAMPLE/newo/matched/*.ChARdnaRna.rev.txt
	do
	fileh=${file##*/}
	transcriptome=$(echo $fileh | cut -f 3 -d '.')
	name=$(echo $fileh | sed -r 's/(.[^.]+){2}$//g')
	echo $name
	echo $transcriptome
	python $SCRIPTDIR/add_gene_info_to_transcriptome.py $file $SCRIPTDIR/dmel_coords.dict/dmel-all-$transcriptome-r5.57.coords.dict $DATAPATH/$SAMPLE/newo/matchedinfo/$name.info.rev.txt
	done

	##Pyton script to add RNA adjusted position and create masterfiles
	#usage:python size_adjust_master.py *info.txt outfile.infosize.txt /raid1/lab/nteran/dmel_coords.dict/dmel-all-gene-r5.57.coords.dict
	#order of infosize.txt is same as defined in google doc masterfile
	echo "Adding genome coordinate information to RNA to make infosize file..."
	for file in $DATAPATH/$SAMPLE/newo/matchedinfo/*info.txt
	do 
	fileh=${file##*/}
	name=$(echo $fileh | sed -r 's/(.[^.]+){2}$//g')
	echo $name 
	wc -l $file
	python $SCRIPTDIR/size_adjust_master.py $file $DATAPATH/$SAMPLE/newo/matchedinfosize/$name.infosize.txt $SCRIPTDIR/dmel_coords.dict/dmel-all-gene-r5.57.coords.dict
	done

	for file in $DATAPATH/$SAMPLE/newo/matchedinfo/*info.rev.txt
	do 
	fileh=${file##*/}
	name=$(echo $fileh | sed -r 's/(.[^.]+){3}$//g')
	echo $name 
	wc -l $file
	python $SCRIPTDIR/size_adjust_master.py $file $DATAPATH/$SAMPLE/newo/matchedinfosize/$name.infosize.rev.txt $SCRIPTDIR/dmel_coords.dict/dmel-all-gene-r5.57.coords.dict
	done
fi



##cat them all together
#work in order
#split by RNA strand alignemnt to sense or antisense of transcripts
if [ -s $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.notranscriptome.txt ]
	then
	echo "Merged file already exists"
	else
	echo "Concatenating transcriptome master files to make prioritymerged file..."
	for TRANSCRIPTOME in tRNA miscRNA ncRNA transcript three_prime_UTR five_prime_UTR exon intron miRNA gene gene_extended2000
	do
	awk '{if ($27 == "+") print $0}'  $DATAPATH/$SAMPLE/newo/matchedinfosize/$SAMPLE.ChARrna.$TRANSCRIPTOME.ChARdnaRna.infosize.txt >> $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.sense.txt
	#awk '{if ($27 == "-") print $0}'  $DATAPATH/$SAMPLE/newo/matchedinfosize/$SAMPLE.ChARrna.$TRANSCRIPTOME.ChARdnaRna.infosize.txt >> $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.antisense.txt
	awk '{if ($27 == "-") print $0}'  $DATAPATH/$SAMPLE/newo/matchedinfosize/$SAMPLE.ChARrna.$TRANSCRIPTOME.ChARdnaRna.infosize.rev.txt >> $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.antisense.txt
	done
	
	echo "Removing duplicates to make prioritized and ribominus files..."
	sort -t$'\t' -k4,4 -u $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.sense.txt > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.sense.txt
	# sort -t$'\t' -k4,4 -u $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.antisense.txt | fgrep -x -f <( awk '{print $4}' $SAMPLE.ChARdnaRna.prioritized.sense.txt) -v > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt
	
	sort -k4,4 $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.antisense.txt > $DATAPATH/$SAMPLE/newo/infanti.txt
	sort -k1,1 $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.sense.txt > $DATAPATH/$SAMPLE/newo/infsense.txt
 	python $SCRIPTDIR/sortedfilter.py -i $DATAPATH/$SAMPLE/newo/infanti.txt -v $DATAPATH/$SAMPLE/newo/infsense.txt -c 4 -d 1 > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt

	grep -v -e 'rRNA' -e 'FBgn0085813' $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.sense.txt > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.ribominus.sense.txt
	grep -v -e 'rRNA' -e 'FBgn0085813' $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.ribominus.antisense.txt
	cd $DATAPATH/$SAMPLE/
	#searches for unique IDs found in RNA aligned to genome that were not found aligned to other transcriptomes

	sort -t$'\t' -k4,4 -u $DATAPATH/$SAMPLE/newo/$SAMPLE.ChARdnaRna.prioritymerged.antisense.txt > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt
	grep -v -e 'rRNA' -e 'FBgn0085813' $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.ribominus.antisense.txt
	cd $DATAPATH/$SAMPLE/

	#searches for unique IDs found in RNA aligned to genome that were not found aligned to other transcriptomes
	#fgrep -x -f <( awk '{print $4}' $SAMPLE.ChARdnaRna.prioritized.sense.txt) -v $DATAPATH/$SAMPLE/newo/matched/$SAMPLE.ChARrna.genome.ChARdnaRna.txt | fgrep -x -f <( awk '{print $4}' $SAMPLE.ChARdnaRna.prioritized.antisense.txt) -v > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.notranscriptome.txt
	

	sort -k4,4 $DATAPATH/$SAMPLE/newo/matched/$SAMPLE.ChARrna.genome.ChARdnaRna.txt > $DATAPATH/$SAMPLE/newo/infgenome.txt
	cat $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt >> $DATAPATH/$SAMPLE/newo/infsense.txt
	sort -k1,1 $DATAPATH/$SAMPLE/newo/infsense.txt > $DATAPATH/$SAMPLE/newo/inftransc.txt
 	python $SCRIPTDIR/sortedfilter.py -i $DATAPATH/$SAMPLE/newo/infgenome.txt -v $DATAPATH/$SAMPLE/newo/inftransc.txt -c 4 -d 1 > $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.notranscriptome.txt

 	rm $DATAPATH/$SAMPLE/newo/infanti.txt
 	rm $DATAPATH/$SAMPLE/newo/infsense.txt
 	rm $DATAPATH/$SAMPLE/newo/infgenome.txt
 	rm $DATAPATH/$SAMPLE/newo/inftransc.txt
	#final number of matched reads:
fi

echo "Tabulating final stats..."
echo 'Final number of matched reads' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
echo 'sense reads'>> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
cat $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.sense.txt | wc -l >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
echo 'antisense reads' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
cat $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.prioritized.antisense.txt | wc -l >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
echo 'Final number of ribominus reads' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
echo 'sense reads' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
cat $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.ribominus.sense.txt | wc -l >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
echo 'antisense reads' >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt
cat $DATAPATH/$SAMPLE/$SAMPLE.ChARdnaRna.ribominus.antisense.txt | wc -l >> $DATAPATH/$SAMPLE/$SAMPLE.summary_statistics.txt


####Possible files to clean up:
if [ $4 == "cleanup" ]
then
echo "Cleaning up SAM files and gzipping master files..."
rm $DATAPATH/$SAMPLE/align/sam/* & #these can be generated from the bam files using
gzip $DATAPATH/$SAMPLE/newo/simple/* &
gzip $DATAPATH/$SAMPLE/newo/matched/* &
gzip $DATAPATH/$SAMPLE/newo/matchedinfo/* &
echo "Finished!" $(date +"%T")
else
echo "Finished!" $(date +"%T")
fi

