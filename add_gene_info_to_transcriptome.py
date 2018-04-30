#!/usr/bin/python
'''
#python add_gene_info_to_transcriptome.py [rna_align.transcriptome_dna.merge.smiple] [transcriptome_dict] [rna.dna.merge_add_gene_info]
### this program takes RNA reads aligned to certain transcriptome with corresponding DNA reads aligned to genome (with adjusted coordiante)
#### and looks for FBtR and adds gene info based on transcriptome specific dictionary
##### 
'''
#The sys module allows us to read from the command line. it should be installed by default
import sys

#initialize empty dictionary
dictionary={}

#optionally keep count of the number of untranslatable IDs
untrans=0

#Retrieve filenames from user open files
rnatranscriptomednamerge=sys.argv[1]
transcriptomedict=sys.argv[2]
outfile=sys.argv[3]

infile = open(rnatranscriptomednamerge, 'r')
dictfile= open(transcriptomedict, 'r')
outfile=open(outfile, 'w')

#populate the dictionary with contents of the translation file
for line1 in dictfile:
	line1=line1.rstrip('\n')
	col1 = line1.split("\t")
	dictionary[col1[0]]='\t'.join(col1[1:])

	#Read through simple file and output  (using dictionary) to outFile
for line2 in infile:
	line2=line2.rstrip('\n')
	col2 = line2.split("\t")
	if col2[8] in dictionary:
		outfile.write(col2[0]+"\t"+'\t'.join(col2[1:])+"\t"+dictionary[col2[8]]+"\n")
		
	else:
		untrans=untrans+1
		

print "There were "+str(untrans)+" untranslated gene ids\n"


