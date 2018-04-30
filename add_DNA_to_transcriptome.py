#!/usr/bin/python
'''
#python add_DNA_to_transcriptome.py [rna_align.transcriptome.smiple] [dna.align.genome.simple_p1] [rna.dna.merge]
### this program takes RNA reads aligned to certain transcriptome and DNA reads aligned to genome (with adjusted coordiante)
#### and looks for uniqueID and adds DNA info to RNA file
##### 
'''
#The sys module allows us to read from the command line. it should be installed by default
import sys

#initialize empty dictionary
dictionary={}

#optionally keep count of the number of untranslatable IDs
untrans=0

#Retrieve filenames from user open files
rnatranscriptomealignfile=sys.argv[1]
dnagenomealign=sys.argv[2]
outfile=sys.argv[3]

infile = open(rnatranscriptomealignfile, 'r')
dictfile= open(dnagenomealign, 'r')
outfile=open(outfile, 'w')

#populate the dictionary with contents of the translation file
for line1 in dictfile:
	line1=line1.rstrip('\n')
	col1 = line1.split("\t")
	dictionary[col1[0]]='\t'.join(col1[1:])

	#Read through simple file and output  (using dictionary) to outFile
for line2 in infile:
	line2=line2.rstrip('\n')
	col2 = line2.split(" ")
	if col2[0] in dictionary:
		outfile.write(col2[0]+"\t"+dictionary[col2[0]]+"\t"+'\t'.join(col2[1:])+"\n")
		
	else:
		untrans=untrans+1
		

print "There were "+str(untrans)+" untranslated gene ids\n"


