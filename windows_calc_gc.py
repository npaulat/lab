import sys
import os
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
from pylab import savefig
from Bio import SeqIO
from Bio.SeqUtils import GC

##### Creating sliding windows and Calculating GC content #####
###############################################################

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Extracting genome sequence from large scaffolds into windows and calculating the GC content of all windows.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input FASTA file of genome scaffolds
	parser.add_argument('-i', '--input', type=str, help='genome FASTA file to use as input in the given directory', required=True)
	#Give directory containing the input file (need full path if not in current directory)
	parser.add_argument('-d', '--directory', type=str, help='Location of directory of the input file', required=True)
	#Give the output file name [default input file basename]
	parser.add_argument('-p', '--prefix', type=str, help='Option to specify output name, defaults to basename of input file')
	#Give minimum length threshold for scaffolds
	parser.add_argument('-ml', '--minlen', type=int, help='Minimum length of scaffolds to use, defaults to 30kb', default=30000)
	#Give window size
	parser.add_argument('-w', '--window', type=int, help='Size of window, defaults to 10kb', default=10000)
	#Give window size
	parser.add_argument('-s', '--step', type=int, help='Size of step, defaults to 10kb', default=10000)
	#Give the output directory (will make a new subdirectory)[default to current]
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output file(s), defaults to current directory')
	#Argument to produce secondary output file = FASTA of windows used [default = no]
	parser.add_argument('-wo', '--winout', type=str, help='Make window FASTA file, y or n', default = 'n')
	
	args = parser.parse_args()
	FILE = args.input
	DIR = args.directory
	if args.prefix:
		PREFIX = args.prefix
	else:
		BASENAME = os.path.basename(FILE).split(".")[0]
		PREFIX = BASENAME
	MIN_LEN = args.minlen
	WINDOW = args.window
	STEP = args.step
	if args.outdir:
		OUT_DIR = args.outdir
	else:
		OUT_DIR = '.'
	WIN_OUT = args.winout
	
	return FILE, DIR, PREFIX, MIN_LEN, WINDOW, STEP, OUT_DIR, WIN_OUT
	
#FILE, DIR, PREFIX, MIN_LEN, WINDOW, STEP, OUT_DIR, WIN_OUT = get_args()

##### Steps #####
# 1. Only use scaffolds with a minimum length (Ex: 30,000 bp)
# 2. Create sliding windows of given length (Ex: 10,000 bp) using "def" function
# 3. Calculate GC content for each window using "def" function
# 4. Plot frequency histogram 


#for making windows
#for a given scaffold, get sequence length
#calculate the number of complete windows by dividing seq length by step size
#for each window within the range of full-length windows (window # x step size)
#to be exactly non-overlapping, step size must equal window size
#return the sequence of the window

def make_windows(SEQUENCE, WIN_SIZE, STEP_SIZE):
	SEQ_LEN = len(SEQUENCE)
	WIN_NUM = int(SEQ_LEN/STEP_SIZE)
	for WIN in range(0, WIN_NUM * STEP_SIZE, STEP_SIZE):
		yield SEQUENCE[WIN:WIN+WIN_SIZE]

#for calculating GC content of windows
#initialize scaffold counter and window counter
#using the genome file, read each scaffold
#if scaffold sequence length is 30k+, up scaffold count, make windows
#calculate GC content using Bio.SeqUtils GC function, up window count, append to list
#use list of GC % to make histogram

def gc_from_windows():
	FILE, DIR, PREFIX, MIN_LEN, WINDOW, STEP, OUT_DIR, WIN_OUT = get_args()
	print('Genome file is ' + FILE + '.')
	print('Minimum scaffold size = ' + str(MIN_LEN) + '.')
	print('Window size = ' + str(WINDOW) + '.')
	
	GC_LIST = []
	SC_COUNT = 0
	WIN_COUNT = 0
	
	with open(FILE, "r") as GENOME:
		for SCAFFOLD in SeqIO.parse(GENOME, 'fasta'):
			SEQ = SCAFFOLD.seq
			if len(SEQ) >= MIN_LEN:
				SC_COUNT += 1
				for SLICE in make_windows(SEQ, WINDOW, STEP):
					#if WIN_OUT == 'y':
						#OUTPUT2 = PREFIX + "_windows.fa"
						#with open(OUTPUT2,"w+") as WIN_FILE:
							#SeqIO.write(SLICE, WIN_FILE, "fasta")
					WIN_COUNT += 1
					WIN_GC = GC(SLICE)
					GC_LIST.append(WIN_GC)
	
	plt.hist(GC_LIST)
	plt.savefig(PREFIX + "_GC_histogram.png")
	
	#print('Scaffolds used (over {} bases = {}'.format(str(MIN_LEN), str(SC_COUNT))
	#print('Windows used = {}'.format(str(WIN_COUNT)))
	
	
if __name__ =='__main__':gc_from_windows()
