#!/usr/bin/python
import Levenshtein

def bridgehunter(bridge,line):
	position,bridgenum,orientation=bridgepos_multicheck(bridge,line)

	if bridgenum==0:
		#then we don't have a bridge
		position=-99
		dna=''
		rna=''
	if bridgenum==1:
		#jackpot!
		#is it fwd or reverse?
		if orientation=='F':
			#forward
			prna=line[0:position]

			pdna=line[(position+len(bridge)):]

			rna,dna=trimrnadna(prna,pdna)
		elif orientation=='R':
			#reverse
			#so we want to reverse the sequences
			backwardsdna=line[0:position+1]

			backwardsrna=line[(position+len(bridge)):]

			pdna=reverse_complement(backwardsdna)
			prna=reverse_complement(backwardsrna)

			rna,dna=trimrnadna(prna,pdna)
			position=len(line)-position-len(bridge)-1 #so it equals the equivalent if fwd
		else:
			print 'Error: charbridgetools bridgehunter: bridgenum not 1'
	else:
		orientation='F'
		rna=''
		dna=''
		position=-100
	return rna,dna,orientation,position,bridgenum

def trimrnadna(prna,pdna):
	#the only issue with the rna will be the extra 6 a's and possibly more from a polyA thing.
	cleanrna=prna.rstrip("A")

	#dna has more problems
	#set up a fake cleandna to be overwritten if need be
	cleandna=pdna
	dpnii_mismatches=['GATCGATCTAGGATC', 'GATCTAGGATCGATC', 'GATCTAGGATC', 'GATCGATC', 'GATC']
	for i in dpnii_mismatches:
		if fuzz_align_once(i,pdna,0):#could make the mismatch permissions 1 or 0
			position,dist=fuzz_align_once(i,pdna,0)#could make the mismatch permissions 1 or 0
			if position==0:
				#cleandna="GATC"+pdna[(len(i)-1):]
				cleandna=pdna[(len(i)):]
				break
	return cleanrna,cleandna


def fuzz_align_once(s_seq,l_seq,mismatch):
	for i, base in enumerate(l_seq):  # loop through equal size windows
		l_subset = l_seq[i:i+len(s_seq)]
		dist = Levenshtein.distance(l_subset, s_seq)
		if dist <= mismatch:  # find first then break
			return i, dist
			break

def bridgepos_multicheck(bridge,line):
	count=0
	#try finding the bridge
	fuzz=fuzz_align_once(bridge,line,1)
	#if you didnt find it, try the revcomp of the bridge
	if not(fuzz):
		fuzz_rc=fuzz_align_once(reverse_complement(bridge),line,1)
		if not(fuzz_rc):
			#if no bridge either way, just give up.
			return 0,0,'F'
		else:
			#if we found the bridge backwards once, check to make sure we don't see it again
			restofline=line[(fuzz_rc[0]+len(bridge)):]
			fuzz_rc_2=fuzz_align_once(reverse_complement(bridge),restofline,1)
			if not(fuzz_rc_2):
				#so this means we found it once but not twice
				return fuzz_rc[0],1,'R'
			else:
				#this means we found it multiple times
				return fuzz_rc[0],2,'R'
	else:
		#so we found it in the fwd direction
		#check the rest of the read in the fwd direction
		restofline=line[(fuzz[0]+len(bridge)):]
		fuzz_fw=fuzz_align_once(bridge,restofline,1)
		#and the rc
		fuzz_rc=fuzz_align_once(reverse_complement(bridge),restofline,1)
		#if either returns something, it means there's a double bridge
		if fuzz_fw or fuzz_rc:
			return fuzz[0],2,'F'
		else:
			#if there isn't, return what we want
			return fuzz[0],1,'F'

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):	
	for k,v in alt_map.iteritems():
		seq = seq.replace(k,v)
	bases = list(seq) 
	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	for k,v in alt_map.iteritems():
		bases = bases.replace(v,k)
	return bases

