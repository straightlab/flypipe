#!/usr/bin/python
USAGE="""Nikki Teran 2016-03-21 from program started 2016-1-16
python char_bridge.py --R1 R1.fastq.gz --R2 R1.fastq.gz
python char_bridge.py --ASSEMBLED S0.assembled.fastq.gz

Separates single end reads into RNA and DNA based on an orientation specific bridge.

Input files must be gzipped
"""
import os
import re
import sys
import bz2
import time
import gzip
import string
import argparse
import charbridgetools as cbt
from Bio import SeqIO
from Bio import AlignIO
from itertools import izip
from optparse import OptionParser

#first things first
bridge1='ACCGGCGTCCAAG' #fwd, so RNA, bridge,
#bridge2='CTTGGACGCCGGT' #rev, so DNA, bridge, RNA
b1_len=len(bridge1)
#b2_len=len(bridge2)

#all the arguments we want to parse.
# python char_bridge_trackall.py --FASTQGZ meaningless.fasta.gz --NAME test. --minRNA 2 --minDNA 2
parser = argparse.ArgumentParser(description = USAGE)
parser.add_argument('--FASTQGZ', dest = 'FASTQGZ', default = "R1.bridgePE.fastq.gz", help = 'REQUIRED FOR PAIRED END READ FUNCTIONALITY: Read 1 input file in gzipped fastq format or forward unassembled read from PEAR')

parser.add_argument('--DNA', dest = 'DNA', default = "dna.bridgePE.fastq.gz", help = 'Desired output for the DNA sequences')
parser.add_argument('--RNA', dest = 'RNA', default = "rna.bridgePE.fastq.gz", help = 'Desired output for the RNA sequences')

parser.add_argument('--NB', dest = 'NB', default = "bridgePE.nobridge.fastq.gz", help = 'Desired output for Read 1 lines that do not contain the bridge')

parser.add_argument('--DB', dest = 'DB', default = "bridgePE.dupbridge.fastq.gz", help = 'Desired output for Read 1 lines that contain duplicate bridges')

parser.add_argument('--TS', dest = 'TS', default = "bridgePE.tooshort.fastq.gz", help = 'Desired output for RNA when either the RNA or DNA is too short')

parser.add_argument('--POS', dest = 'POS', default = "bridgePE.bridgeposition.fastq.gz", help = 'Desired output for list of bridge positions, if bridge is found')

parser.add_argument('--minRNA', dest = 'minRNA', default = 18, help = 'The minimum length of RNA, shorter reads will be discarded. Default: 18')
parser.add_argument('--minDNA', dest = 'minDNA', default = 18, help = 'The minimum length of RNA, shorter reads will be discarded. Default: 18')

parser.add_argument('--NAME', dest = 'NAME', default = "", help = 'Sample name to add to start of all output files.')

args = parser.parse_args()
name=args.NAME

R1_fh = args.FASTQGZ

DNA_fh = name+args.DNA
RNA_fh = name+args.RNA

NB_fh = name+args.NB

DB_fh = name+args.DB

TS_fh = name+args.TS

POS_fh = name+args.POS

minDNA = args.minDNA
minRNA = args.minRNA

#open files
rna_out = gzip.open(RNA_fh, 'w')
dna_out = gzip.open(DNA_fh, 'w')

nb_out = gzip.open(NB_fh, 'w')
db_out = gzip.open(DB_fh, 'w')
ts_out = gzip.open(TS_fh, 'w')
pos_out = gzip.open(POS_fh, 'w')

with gzip.open(R1_fh) as read_in:
	linecount=0
	for line in read_in:
		line=line.strip()
		linecount=linecount+1
		if linecount==1:
			#if line =1 save the name
			name=line
		elif linecount==2:
		#if line =2, hunt for the bridge
			ntseq=line
			rnaseq,dnaseq,orientation,position,bridgenum=cbt.bridgehunter(bridge1,line)
			#check length
			#will do it right before outputting
		elif linecount==4:
			#if line =3 its just a stupid plus sign
			#if line =4 split the quality score along where the dna was
			linecount=0
			#remember to flip the line if we found the string backwards
			if orientation=='R':
				line=line[::-1]

			rnaqual=line[0:len(rnaseq)]
			dnaqual=line[-len(dnaseq):]
			#print out the position information for everything. There should be as 1/4 as many lines in pos as in the input file
			pos_out.write(str(position)+'\t'+str(len(line)-position-len(bridge1))+'\n')
			if bridgenum==0:
				nb_out.write(name+'\n'+ntseq+'\n'+'+'+'\n'+line+'\n')
			elif bridgenum>1:
				db_out.write(name+'\n'+ntseq+'\n'+'+'+'\n'+line)
			elif len(rnaseq) < int(minRNA) or len(dnaseq) < int(minDNA):
				ts_out.write(name+'\n'+ntseq+'\n'+'+'+'\n'+line+'\n')
			elif bridgenum==1:
				rna_out.write(name+'\n'+rnaseq+'\n'+'+'+'\n'+rnaqual+'\n')
				dna_out.write(name+'\n'+dnaseq+'\n'+'+'+'\n'+dnaqual+'\n')
			else:
				print 'Error: read not printed to any designated file'

