import sys
import os
import argparse
import itertools
import subprocess
import fnmatch
from Bio import SeqIO

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Adjusting reference MELT hit coordinates to new MSA genome coordinates", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#required = parser.add_argument_group('required arguments')
	#Give input coordinate file of MELT hits (CHROM# + POS)
	parser.add_argument('-i', '--input', help='MELT coordinate list to use as input file in the given directory', required=True)
	#Argument of directory containing MSA scaffolds (need full path)
	parser.add_argument('-d', '--directory', type=str, help='Location of directory of the query scaffold files', required=True)
	#Argument of directory containing MSA scaffolds (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file', required=True
	
	args = parser.parse_args()
	COORDS = args.input
	DIR = args.directory
	OUTDIR = args.outdir
	
	return COORDS, DIR
	
COORDS, DIR, OUTDIR = get_args()

#argument sanity checks
#if not args.directory:
#  	sys.exit('You must provide the directory containing your scaffold alignments.')
#if not args.input:
#  	sys.exit('You must provide the list of MELT coordinates.')

print('The input coordinate list is ' + INPUT +'.')
print('The scaffold MSA sequence directory is ' + DIR + '.')

cd OUTDIR

#open an ouput coordinate file
f = open("adj_meltcoords.txt", "w+")

cd DIR

for LINE in open(COORDS):
#store list of chromosomes as dictionary, 
# my_dict = ["CHROM1":[(23,178)(83,290)]}
my_dict = {}

for line in open(COORDS)
	line = line.split()
	CHROM = line[0]
	START = int(line[1])
	STOP = int(line[2])
	if line[0] in my_dict:
		my_dict[line[0]].append((START,STOP))
	else:
		my_dict[CHROM] = [(START, STOP)]
		
#	for every key in the dictionary, open MSA
# need sorted coords list
#edit to list
	for CHROM in my_dict:
		MSA = CHROM + ".fas"
		if not os.path.isfile(MSA):
			#print and continue on list, or print and exit
			print("Cannot find file: {}".format(MSA))
			#sys.exit()
		CHROM_ADJ = 0
		for SEQ in SeqIO.parse(MSA, 'fasta'):
			if SEQ.id == "MyoLuc":
			#argument of SPECIES_NAME
			# if seq.id == SPECIES_NAME:
			#for loop of count up to CHROM[1]
				#Variable = access CHROM list, first insertion on chrom, start pos
				FIRST_START = my_dict[CHROM][0][0]
				#Subseq is up to the start pos of first insert
				FIRST_SUBSEQ = SEQ.seq[:FIRST_START]
				#add +1 for each gap character up to first insertion start pos
				CHROM_ADJ += FIRST_SUBSEQ.count("-")
				#for insertion (Start,stop)
				for INSERTION in my_dict[CHROM]:
					SUBSEQ = SEQ.seq[INSERTION[0]:INSERTION[1]]
					#might need to +/- 1 for counting system
					NUM_DASH = SUBSEQ.count("-")
					CHROM_ADJ += NUM_DASH
					NEW_START = INSERTION[0] + CHROM_ADJ
					NEW_END = INSERTION[1] + CHROM_ADJ
					#tell how to format the output; {} = items
					#f.write("{}\t{}\t{}\n".format(CHROM,NEW_START,NEW_END))
					sys.stdout.write("{}\t{}\t{}\n".format(CHROM,NEW_START,NEW_END))
	else:
		print("Cannot find MSA.")
else:
	print("Unable to read input list. Try again.")
	
#close output file
f.close()