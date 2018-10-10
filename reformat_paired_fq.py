import sys
import os
import argparse
import subprocess
import fnmatch
import itertools
import gzip
from collections import deque

#Created on Tues Sep 4 2018 by Nicole S. Paulat

##Made to work with interleaved paired_end FASTQ files with the assumption that the file is in standard format, with the forward read then reverse read of each pair listed

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Custom formatting of FASTQ headers", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input FASTQ file
	parser.add_argument('-i', '--input', help='FASTQ file basename to use as input within the current directory', required=True)
	#Argument of directory containing the FASTQ file (need full path)
	#parser.add_argument('-d', '--directory', type=str, help='Location of directory of the FASTQ file', required=True)
	#Argument of the output directory
	parser.add_argument('-od','--outdir', type=str, help='Path to the desired output directory', required=True)
	
	args = parser.parse_args()
	FILE = args.input
	#DIR = args.directory
	OUT_DIR = args.outdir
	
	return FILE, OUT_DIR
	#return FILE, DIR

FILE, OUT_DIR = get_args()	
#FILE, DIR = get_args()

#define files, file paths
ZIP = FILE
BASENAME = os.path.basename(FILE).split(".")[0]
#ZIP = os.path.join(DIR,FILE)
FASTQ = BASENAME + ".fastq"
OUTPUT_DIR = OUT_DIR
#OUTPUT_DIR = "/lustre/scratch/npaulat/formatted_fastq/"
OUTPUT = os.path.join(OUTPUT_DIR, FASTQ)

#open zipped input FASTQ file within the Python virtual environment to read as text
d = gzip.open(ZIP, 'rt')

#define function to find length of the zipped input FASTQ file
def find_length():
	size = 0
	with d:
		for line in d:
			size += 1
			pass

	return size

LENGTH = find_length()
d.close()
print(LENGTH)

#f = open(OUTPUT + FASTQ, "a+")
f = open(OUTPUT, "w+")

#define function to generate list of header line index values
def my_range(START, END, STEP):
	i = START
	END = LENGTH - 1
	while i <= END:
		yield i
		i += STEP

#start at first line of file
i = 0
#start at 5th line of file
j = 4
#start numbering for first of pair headers at 1
header1 = 1
#start numbering second of pair headers at 1
header2 = 1

#generate the Header1 index deque/queue
HEADER1_INDEX_GEN = my_range(i, LENGTH, 8)
HEADER1_INDEX = []
for x in HEADER1_INDEX_GEN:
	HEADER1_INDEX.append(x)
HEADER1_QUEUE = deque(HEADER1_INDEX)

#generate the Header2 index deque/queue
HEADER2_INDEX_GEN = my_range(j, LENGTH, 8)
HEADER2_INDEX = []
for x in HEADER2_INDEX_GEN:
	HEADER2_INDEX.append(x)
HEADER2_QUEUE = deque(HEADER2_INDEX)

#Check pairing by comparing Header list lengths
#probably should replace this with more accurate check via BioPython
LIST1_LEN = len(HEADER1_INDEX)
LIST2_LEN = len(HEADER2_INDEX)

if LIST2_LEN != LIST1_LEN:
	print("There appears to be {} unpaired read in your file.".format(str(abs(LIST1_LEN - LIST2_LEN))))
	print(HEADER1_INDEX)
	print(HEADER2_INDEX)

###code check to ensure all pairs are in forward read then reverse read order

#open zipped FASTQ within the python virtual environment for reading as text
#d = gzip.open(ZIP, 'rt')
with gzip.open(ZIP, 'rt') as d:
#enumerate the file lines with corresponding index values
	for index, line in enumerate(d):
		#if Header1 queue is not empty
		if HEADER1_QUEUE:
			#if line index value is in the list of index values of every 8th line
			if index == HEADER1_QUEUE[0]:
				#check that the line starts with the "@" symbol
				if line[0] == "@":
					#replace that line with formatted header for 1/2 of pairs
					f.write("@" + '{:0>12}'.format(header1) + " OP:i:1\n")
					#update Header1 numbering
					header1 += 1
					#remove used index number from Header1 queue
					HEADER1_QUEUE.popleft()
				#if "@" symbol not there, something is wrong, break off
				else:
					print("Did not find header of line {} as expected. Aborted.".format(index))
					sys.exit(1)
			#if line index value is not the first item in Header1 queue
			elif not index == HEADER1_QUEUE[0]:
				#if Header2 queue is not empty
				if HEADER2_QUEUE:
					#if line index value is the first item in Header2 queue
					if index == HEADER2_QUEUE[0]:
						if line[0] == "@":
							#replace that line with formatted header for 1/2 of pairs
							f.write("@" + '{:0>12}'.format(header2) + " OP:i:2\n")
							#update Header2 numbering
							header2 += 1
							#remove used index number from Header2 queue
							HEADER2_QUEUE.popleft()
						#if "@" symbol not there, something is wrong, break off
						else:
							print("Did not find header of line {} as expected. Aborted.".format(index))
							sys.exit(1)
					else:
						f.write(line)
				else:
					f.write(line)
		#if Header1 queue is empty, start at Header2 queue
		elif HEADER2_QUEUE:
			#if line index value is the first item in Header2 queue
			if index == HEADER2_QUEUE[0]:
				#check that the line starts with the "@" symbol
				if line[0] == "@":
					#replace that line with formatted header for 2/2 of pairs
					f.write("@" + '{:0>12}'.format(header2) + " OP:i:2\n")
					header2 += 1
					HEADER2_QUEUE.popleft()
				#if "@" symbol not there, something is wrong, break off
				else:
					print("Did not find header of line {} as expected. Aborted.".format(index))
					sys.exit(1)
			#else Header1 queue empty, Header2 queue first item non-match
			#write out line as is
			else:
				f.write(line)
		#if line index value is not on list, write out as is
		else:
			f.write(line)

f.close()
d.close()