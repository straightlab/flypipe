#!/usr/bin/python
###this program adjusts the position of the read to account for chromosomes before it
##### as if genome is one long sequence
###order of chromosomes in dictionary
##python size_adjust_master.py [inputfile(dnaRna.matched.txt)] [outputfile] /raid1/lab/nteran/dmel_coords.dict/dmel-all-gene-r5.57.coords.dict

#The sys module allows us to read from the command line. it should be installed by default
import sys
import time

#initialize size adjust dictionary
dictionary={'chrM':129663327,'chr2L': 0, 'chr2LHet': 23011544, 'chr2R': 23380416, 'chr2RHet':44527124, 'chr3L': 47815885, 'chr3LHet':72359442, 'chr3R': 74914933, 'chr3RHet':102819986, 'chr4': 105337493, 'chrX': 106689350, 'chrXHet': 129112177, 'chrYHet': 129316289, 'chrdmel_mitochondrion_genome':129663327, 'chrU':129682851}

#initialize empty dictionary
genedict={}
genetss={}
genetes={}

#optionally keep count of the number of unadjusted rows
untrans=0

#Retrieve filenames from user open files
matchedfile=sys.argv[1]
transfile=sys.argv[3]
outfile=sys.argv[2]

infile = open(matchedfile, 'r')
dictfile= open(transfile, 'r')
outfile=open(outfile, 'w')

#populate the dictionary with contents of the translation file
for line1 in dictfile:
#	print line
	line1=line1.rstrip('\n')
	col1 = line1.split("\t")
	genedict[col1[9]]=[col1[8],col1[4],col1[5]]

#Read through simple file and output  (using dictionary) to outFile
for line in infile:
	line=line.rstrip('\n')
	col = line.split("\t")
	if col[16] == 'chrM':
		continue
	if col[16] == 'chrdmel_mitochondrion_genome':
		continue
	if col[16] == 'chrU':
		continue 
	DNAchrom=str(col[1])
	uniqID=str(col[0])
	DNAadjpos=int(col[3])
	DNAstrand=str(col[7])
	DNAmapq=int(col[4])
	DNAcigar=str(col[5])
	RNAfb=str(col[8])
	RNAmapq=int(col[10])
	RNAcigar=str(col[11])
	RNAreftype=str(col[14])
	RNArefchrom=str(col[16])
	RNArefname=str(col[21])
	DNApos=int(col[2])
	DNApos1=int(DNApos-1)
	RNATSS=int(col[17])
	RNATTS=int(col[18])
	RNApos=int(col[9])
	RNApos1=int(RNApos-1)
	RNAadd=int(dictionary[col[16]])
	RNAreadstrand=str(col[13])
	RNAtranscriptstrand=str(col[19])
	RNAsplicesites=str(col[20])
	RNAtranscriptlength=int(col[24])
	RNAreadlength=int(col[12])
	RNAparentgn=str(col[22])
	add=int(1)
	RNAgenename=str(genedict[RNAparentgn][0])
	RNAgenetss=str(genedict[RNAparentgn][1])
	RNAgenetes=str(genedict[RNAparentgn][2])
	if RNATSS < RNATTS:
		RNAbedstart = RNATSS
		RNAbedstop = RNATTS
		RNAgenestart = RNAgenetss
		RNAgenestop = RNAgenetes
	else:
		RNAbedstart = RNATTS
		RNAbedstop = RNATSS 
		RNAgenestart = RNAgenetes
		RNAgenestop = RNAgenetss
	RNAgenelength = int(int(RNAgenestop)-int(RNAgenestart))
	if RNATSS < RNATTS:
		if RNAreadstrand == "+":
			RNApos2=int(RNApos1+RNAreadlength)
			RNAchrompos = int(RNATSS+RNApos2)
		else:
			RNApos2=int(RNApos1)
			RNAchrompos = int(RNATSS+RNApos2)
	else:
		if RNAreadstrand == "+":
			RNApos2 = int(RNApos1-RNAreadlength)
			RNAchrompos = int(RNATSS-RNApos2)
		else: 
			RNApos2 = int(RNApos1)
			RNAchrompos = int(RNATSS-RNApos2)
	if RNArefchrom in dictionary:
		outfile.write(DNAchrom+"\t"+str(DNApos1)+"\t"+str(DNApos1+add)+"\t"+uniqID+"\t"+str(DNAadjpos-1)+"\t"+DNAstrand+"\t"+str(DNAmapq)+"\t"+DNAcigar+"\t"+RNAfb+"\t"+str(RNApos2)+"\t"+str(RNAmapq)+"\t"+RNAcigar+"\t"+RNAreftype+"\t"+RNArefchrom+"\t"+str(RNAbedstart)+"\t"+str(RNAbedstop)+"\t"+RNArefname+"\t"+str(RNATSS+RNAadd)+"\t"+RNArefchrom+"\t"+str(RNATSS)+"\t"+str(RNATTS)+"\t"+RNAsplicesites+"\t"+RNArefchrom+"\t"+str(RNAchrompos)+"\t"+str(RNAchrompos+add)+"\t"+str(RNAchrompos+RNAadd)+"\t"+RNAreadstrand+"\t"+RNAgenename+"\t"+RNAtranscriptstrand+"\t"+str(RNAtranscriptlength)+"\t"+str(RNAreadlength)+"\t"+RNAparentgn+"\t"+RNAgenetss+"\t"+RNAgenetes+"\t"+RNArefchrom+"\t"+RNAgenestart+"\t"+RNAgenestop+"\t"+str(RNAgenelength)+"\t"+str(RNApos1)+"\n")
		
	else:
		untrans=untrans+1
		print line

print "There were "+str(untrans)+" unadjusted rows"
