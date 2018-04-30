#!/usr/bin/python
###this program adjusts the position of the read to account for chromosomes before it
##### as if genome is one long sequence
### order of columns: uniqueID, DNA_chrom, DNA_pos, DNA_adjusted_pos 
###order of chromosomes in dictionary
##python size_adjust.py [inputfile(simple.txt)] [outputfile]

#The sys module allows us to read from the command line. it should be installed by default
import sys

#initialize empty dictionary
dictionary={'chr2L': 0, 'chr2LHet': 23011544, 'chr2R': 23380416, 'chr2RHet':44527124, 'chr3L': 47815885, 'chr3LHet':72359442, 'chr3R': 74914933, 'chr3RHet':102819986, 'chr4': 105337493, 'chrX': 106689350, 'chrXHet': 129112177, 'chrYHet': 129316289, 'chrM':129663327}

#optionally keep count of the number of unadjusted rows
untrans=0

#Retrieve filenames from user open files
simplefile=sys.argv[1]
#transfile=sys.argv[2]
outfile=sys.argv[2]

infile = open(simplefile, 'r')
#dictfile= open(transfile, 'r')
outfile=open(outfile, 'w')

#Read through simple file and output  (using dictionary) to outFile
for line in infile:
	line=line.rstrip('\n')
	col = line.split(" ")
	pos=int(col[2])
	add=int(dictionary[col[1]])
	if col[1] in dictionary:
		outfile.write(str(col[0])+"\t"+str(col[1])+"\t"+str(pos)+"\t"+str(pos+add)+"\t"+str(col[3])+"\t"+str(col[4])+"\t"+str(col[5])+"\t"+str(col[6])+"\n")
	else:
		untrans=untrans+1

print "There were "+str(untrans)+" unadjusted rows"


