import sys
import os
import argparse
import itertools
import subprocess
from pyfaidx import Faidx
from pyfaidx import Fasta
import time
import re
import stat
import shutil
from shutil import copyfile
import errno
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
#import numpy as np
#import pandas as pd

# Where RepeatMasker is stored
REPEATMASKER = "/lustre/work/daray/software/RepeatMasker"
# Where this script can find liftUp, twoBitInfo and twoBitToFa
BIN_DIR = "/lustre/work/daray/software"

# Define arguments
def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Generate SGE cluster runs for RepeatMasker; built in RepeatMasker parameters are -xsmall [softmasks repetitive regions] -a [.align output file] -gff [generates a GFF format output] -pa [runs in parallel], please see RepeatMasker for details of these run options", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input genome FASTA
	parser.add_argument('-i', '--input', type=str, help='genome file in FASTA format', required=True)
	#Argument of species name
	parser.add_argument('-sp', '--species', type=str, help='Source species of query DNA FASTA', required=False)
	# Desired batch number
	parser.add_argument('-b', '--batch_count', type=int, help='Batch count', default=50)
	# Input genome directory
	parser.add_argument('-dir', '--genome_dir', type=str, help='Path to genome FASTA', required=True)
	# Argument for output directory
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output subdirectory', default='.')
	# Which queue to use
	parser.add_argument('-q', '--queue', type=str, help='Select the queue to run RepeatMasker in [quanah|hrothgar] with the quanah option being the general quanah omni queue, and hrothgar being the communitycluster Chewie queue', choices=['quanah', 'hrothgar'], default='quanah')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-lib', type=str, help='RepeatMasker run parameter custom library "-lib [filename]" option', required=False)
	#Argument of RepeatMasker run parameter
	parser.add_argument('-xsmall', type=str, help='Select a RepeatMasker masking option as lowercase bases [-xsmall], default is to mask as Ns', action='store_true')
#	#Argument of RepeatMasker run parameter
#	parser.add_argument('-engine', type=str, help='RepeatMasker run parameter "-engine <search_engine>" option; select a non-default search engine to use, otherwise RepeatMasker will used the default configured at install time; [crossmatch|abblast|rmblast|hmmer]', choices=['crossmatch', 'abblast', 'rmblast', 'hmmer'], required=False)
#	#Argument of RepeatMasker run parameter
#	parser.add_argument('-inv', type=str, help='RepeatMasker parameter flag "-inv" option; alignments are presented in the orientation of the repeat', action='store_true')
#	#Argument of RepeatMasker run parameter
	parser.add_argument('-nolow', type=str, help='RepeatMasker parameter flag "-nolow" option; does not mask low complexity DNA or simple repeats', action='store_true')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-s', '-speed', type=str, help='RepeatMasker run parameter "-q" or "-s" option; q=quick search; 5-10% less sensitive, 3-4 times faster than default; s=slow search; 0-5% more sensitive, 2.5 times slower than default', choices=['q', 's'], required=False)
#	#Argument of RepeatMasker run parameter
#	parser.add_argument('-div', type=int, help='RepeatMasker run parameter "-div [number]" option; masks only those repeats that are less than [number] percent diverged from the consensus sequence', required=False)
	
	args = parser.parse_args()
	GENOME = args.input
	SPECIES = args.species
	BATCH_COUNT = args.batch_count
	GENOME_DIR = args.genome_dir
	OUTDIR = args.outdir
	QUEUE = args.queue
	LIBRARY = args.lib
	XSMALL = args.xsmall
#	ENGINE = args.engine
#	INV = args.inv
	NOLOW = args.nolow
	SPEED = args.speed
#	DIV = args.div
	
#	return GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV
	return GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, NOLOW, SPEED
	
#GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV = get_args()
GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, NOLOW, SPEED = get_args()

# Sanity checks
print("The query genome is {}.\n".format(GENOME))
print("{} batches will be made.\n".format(str(BATCH_COUNT)))
print("The genome FASTA is located in '{}'.\n".format(GENOME_DIR))
print("The output directory is '{}'.\n".format(OUTDIR))
print("The job queue is {}.\n".format(QUEUE))

if not SPECIES or LIBRARY:
	sys.exit("Must supply value for option 'species' or 'lib'!")
if SPECIES and LIBRARY:
	sys.exit("Only supply a value for one option: 'species' or 'lib'! Not both!")

#FLAGS = [LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV]
FLAGS = [LIBRARY, XSMALL, NOLOW, SPEED]
if not FLAGS:
	print("All default RepeatMasker parameters were used, no custom library.")
else:
	print("Custom parameters used:\n")
	if XSMALL:
		print("-xsmall flag used.\n")
#	if INV:
#		print("-inv flag used.\n")
	if NOLOW:
		print("-nolow flag used.\n")
	if LIBRARY:
		print("-lib flag used. Custom library is '{}'.\n".format(os.path.basename(LIBRARY)))
#	if ENGINE:
#		print("-engine flag used. Changed search engine to {}.\n".format(ENGINE))
	if SPEED:
		print("-{} flag used. Search sensitivity has changed.\n".format(SPEED))
#	if DIV:
#		print("-div flag used. RepeatMasker will mask only repeats that are less than {}% diverged from the consensus sequence.\n".format(str(DIV)))

if not os.path.isdir(GENOME_DIR):
	sys.exit("The given genome directory, '{}', does not exist.".format(GENOME_DIR))

GENOME_FASTA = os.path.join(GENOME_DIR, GENOME)

#if not os.path.isfile(GENOME_FASTA):
#	sys.exit("The given genome file '{}' does not exist.").format(GENOME_FASTA)
#if os.stat(GENOME).st_size==0:
#	sys.exit("The genome file, '{}', is empty.").format(GENOME_FASTA)
#if not os.path.isfile(LIBRARY):
#	sys.exit("The given library file '{}' does not exist.").format(LIBRARY)
#if os.stat(LIBRARY).st_size==0:
#	sys.exit("The library file, '{}', is empty.").format(LIBRARY)

try:
	if not os.path.getsize(GENOME_FASTA) > 0:
		sys.exit("The genome file, '{}', is empty.".format(GENOME_FASTA))
except OSError as e:
	sys.exit("The genome file '{}' does not exist or is inaccessible.".format(GENOME_FASTA))
	
try:
	if not os.path.getsize(LIBRARY) > 0:
		sys.exit("The library file, '{}', is empty.".format(LIBRARY))
except OSError as e:
	sys.exit("The library file '{}' does not exist or is inaccessible.".format(LIBRARY))
	
if not os.path.isdir(OUTDIR):
	sys.exit("The output directory '{}' does not exist.".format(OUTDIR))

PARTITION_DIR = os.path.join(GENOME_DIR, "RMPart")

SLOTS_PER_BATCH = 10
NUM_BATCHES = BATCH_COUNT

check_empty(PARTITION_DIR)
# &checkEmpty($partitionDir);

if LIBRARY:
	copyfile(LIBRARY, PARTITION_DIR)
	LIB_FILE = os.path.basename(LIBRARY)
	LIBRARY = os.path.join(PARTITION_DIR, LIB_FILE)
# my $lib;
# if ( exists $options{'lib'} )
# {
  # system("cp $options{'lib'} $partitionDir");
  # my ( $vol, $dir, $file ) = File::Spec->splitpath( $options{'lib'} );
  # $lib="$partitionDir/$file";
# }

#simple_partition()
#build_DoLift()
# &simplePartition();
# &buildDoLift();

# exit;


############## FUNCTIONS #################

# Subroutine (1)
def check_empty(PARTITION_DIR):
	PARTITION_DIR = os.path.abspath(PARTITION_DIR)
	if not os.path.exists(PARTITION_DIR):
		try:
			os.makedirs(PARTITION_DIR)
		except OSError as e:
			if e.errno != errno.EEXIST:
				raise
		print("Made '{}' directory.".format(PARTITION_DIR))
	else:
		if not os.listdir(PARTITION_DIR):
			print("'{}' is empty. Continuing.".format(PARTITION_DIR))
		else:
			print("'{}' is not empty. Removing contents and continuing.".format(PARTITION_DIR))
		#	os.remove(os.path.join(PARTITION_DIR, '*'))
		## To remove anything in the RMPart folder from a previous run (all symbolic links (not expected here) and files and subdirectories) without deleting the RMPart directory itself
			for FILE in os.listdir(PARTITION_DIR):
				FILE_PATH = os.path.join(PARTITION_DIR, FILE)
				try:
					shutil.rmtree(FILE_PATH)
				except OSError:
					os.remove(FILE_PATH)

# Subroutine (2)
def get_batches(NUM_BATCHES, GENOME_FASTA):
	# Return a 3-level list(ref): partitions -> chunks -> chunk properties (scaffold + coordinates)
	PARTS = []
	GENOME_NAME = os.path.basename(GENOME_FASTA).split(".")[0]
	TOTAL_SIZE = 0
	SEQS = {}
	
	FAIDX = Faidx(GENOME_FASTA)
	FASTA_IDX = GENOME_FASTA + ".fai"
	
	with open(FASTA_IDX) as FILE:
		for LINE in FILE:
			LINE = LINE.rstrip()
			SEQ, SEQ_SIZE, JUNK = LINE.split("\t", 2)
			TOTAL_SIZE += int(SEQ_SIZE)
			SEQS[SEQ] = int(SEQ_SIZE)
	
	if NUM_BATCHES > 0:
		CHUNK_SIZE = int(TOTAL_SIZE / NUM_BATCHES) + 1
	
	BATCHES = []
	CURRENT_BATCH_SIZE = 0
	for SCAFFOLD in SEQS:
		SEQ_SIZE = SEQS[SCAFFOLD]
		SEQ_IDX = 0
		
		while SEQ_SIZE > 0:
			if (CURRENT_BATCH_SIZE + SEQ_SIZE) > CHUNK_SIZE:
				FILL_SIZE = CHUNK_SIZE - CURRENT_BATCH_SIZE
				CHUNK_INFO = str(GENOME_NAME + ":" + SCAFFOLD + ":" + str(SEQ_IDX) + "-" + str(SEQ_SIZE))
				#NOTE: For scaffold size, always refer back to the index dict, not SEQ_SIZE,
				# since SEQ_SIZE changes depending on if the whole scaffold was used in
				# a single batch or not (as in the if statement of this loop)
				PARTS.append([SCAFFOLD, SEQS[SCAFFOLD], SEQ_IDX, FILL_SIZE, CHUNK_INFO])
				BATCHES.append([PARTS])
				PARTS = []
				SEQ_IDX += FILL_SIZE
				SEQ_SIZE -= FILL_SIZE
				CURRENT_BATCH_SIZE = 0
			else:
				CHUNK_INFO = str(GENOME_NAME + ":" + SCAFFOLD + ":" + str(SEQ_IDX) + "-" + str(SEQ_SIZE))
				PARTS.append([SCAFFOLD, SEQS[SCAFFOLD], SEQ_IDX, SEQ_SIZE, CHUNK_INFO])
				CURRENT_BATCH_SIZE += SEQ_SIZE
				SEQ_SIZE = 0
	#unclear if BATCHES will be in the appropriate hierarchy of lists/parts(elements) atm
	# This bit must be outside of the for loop, otherwise each iteration thru the loop 
	# will append the current PARTS list to BATCHES x# of scaffolds in the PARTS list
	if PARTS:
		BATCHES.append([PARTS])
	
	return BATCHES
# Subroutine (3)
def part_path_from_num(PART_NUM, LEVELS):
	MAX_DIR_SIZE = 1000
	# Given a partition ID number and number of levels in directory tree,
	# determine its path and base filename.
	LEAF_ID = PART_NUM % MAX_DIR_SIZE
	PART_MOD = int(PART_NUM / MAX_DIR_SIZE)
	#Use this line if actually include that dumb loop below
	#PATH = str('{:03d}'.format(LEAF_ID))
	PATH = str('{:03d}/'.format(LEAF_ID))
	## So the below loop should never come into play (as far as I can tell),
	## because LEVELS should always = 1, so w/ i=1 i is never less than LEVELS 
	## (unless # BATCHES exceeds MAX_DIR_SIZE)
	## w/n an RMPart dir, PATH folders are always 000 to n, no addnl nums on either side
	#for i in range(1, LEVELS):
	#i = 1
	#while i < LEVELS:
		#PATH = PATH + str('{:03d}/'.format(PART_MOD % MAX_DIR_SIZE))
		#PART_MOD = int(PART_MOD / MAX_DIR_SIZE)
		#i+= 1
	#my partName = $path;
	PART_NAME = PATH
	#$partName =~ s@/@@g ## is telling it to replace "/" with nothing, globally (90% sure)
	#	var	match operator of must contain	s is strict or s is the start of s///g???
	PART_NAME = PART_NAME.replace('/', '')
	
	return(PATH, PART_NAME)

# Subroutine (4)
#partially tested, incomplete
def simple_partition():
	print("Generating list of batches...\n")
	PARTS = get_batches(BATCH_COUNT, GENOME_FASTA)
	NUM_PARTS = len(PARTS)
	# Basically LEVELS always equals 1 until NUM_PARTS >= MAX_DIR_SIZE
	MAX_DIR_SIZE = 1000
	LEVELS = 1 + int(log(NUM_PARTS) / log(MAX_DIR_SIZE))
	
	TRAILING_SIZE = 0
	MIN_SEQ_CNT = -1
	MAX_SEQ_CNT = 0
	for BATCH in PARTS:
		SEQ_SIZE = 0
		for SEQS in BATCH:
			SEQ_CNT = 0
			for SEQ in SEQS:
				SEQ_CNT += 1
				SEQ_SIZE += SEQ[3]
			#print(SEQ_CNT)
			if (MIN_SEQ_CNT < 0) or (MIN_SEQ_CNT > SEQ_CNT):
				MIN_SEQ_CNT = SEQ_CNT
			if (MAX_SEQ_CNT < SEQ_CNT):
				MAX_SEQ_CNT = SEQ_CNT
		if SEQ_SIZE != CHUNK_SIZE:
			TRAILING_SIZE = SEQ_SIZE
		#print(SEQ_SIZE)
	print("Batch Stats:\n")
	print("  - {} batches with between {} and {} sequences.".format(str(len(PARTS)), str(MIN_SEQ_CNT), str(MAX_SEQ_CNT)))
	if TRAILING_SIZE:
		print("  - {} batches with {} bp and one trailing batch with {} bp.".format(str(len(PARTS) - 1), str(CHUNK_SIZE), str(TRAILING_SIZE)))
	else:
		print("  - All batches contain exactly {} bp.".format(str(CHUNK_SIZE)))
	
	### UNTESTED BELOW THIS POINT ###
	
	## Combining my and Dr. Ray's def for making the qsub files
	
	PREFIX == GENOME_NAME
	SLOTS_PER_BATCH == PROC
	
	def generate_rmsubs(BATCH_COUNT, QUEUE, PREFIX, LIBRARY, SPECIES, XSMALL, NOLOW, SPEED, PROC):
		LIB_OR_SPECIES = ""
		if LIBRARY:
			LIB_OR_SPECIES = " -lib " + LIBRARY + " "
		else:
			LIB_OR_SPECIES = " -species " + SPECIES + " "
		
		ADD_PARAMS = str(LIB_OR_SPECIES)
		if NOLOW:
			ADD_PARAMS = ADD_PARAMS + "-nolow "
		if XSMALL:
			ADD_PARAMS = ADD_PARAMS + "-xsmall "
		if SPEED:
			ADD_PARAMS = ADD_PARAMS + "-" + str(SPEED) + " "	
			
		# Build SGE file
		for BATCH in range(BATCH_COUNT):
			BATCHNUMBER = '{:03}'.format(BATCH)
			SGEBATCH_NAME = "batch-" + BATCH_NUMBER + ".sh"
			SGEBATCH_PATH = os.path.join(PARTITION_DIR, PATH)
			SGEBATCH = os.path.join(SGEBATCH_PATH, SGEBATCH_NAME)
			
			if QUEUE = 'hrothgar':
				PROJECT = 'communitycluster'
				KUE = 'Chewie'
			elif QUEUE = 'quanah':
				PROJECT = 'quanah'
				KUE = 'omni'
			else:
				sys.exit('Please choose hrothgar or quanah as the queue to use.')
				
			with open(SGEBATCH) as BATCH_FILE:
			#with open("RMPart/" + BATCH + "/" + QSUBFILENAME, 'w') as THISFILE:
				BATCH_FILE.write("#!/bin/sh\n")
				BATCH_FILE.write('#$ -V\n')
				BATCH_FILE.write('#$ -cwd\n')
				BATCH_FILE.write('#$ -S /bin/bash\n')
				BATCH_FILE.write("#\$ -N {}-batch-{}\n".format(PREFIX, BATCH_NUMBER))
				BATCH_FILE.write("#\$ -o \$JOB_NAME.o.\$JOB_ID\n")
				BATCH_FILE.write("#\$ -e \$JOB_NAME.e.\$JOB_ID\n")
				BATCH_FILE.write('#$ -q ' + KUE + '\n')
				BATCH_FILE.write('#$ -pe sm {}\n'.format(str(PROC)))
				BATCH_FILE.write('#$ -P ' + PROJECT + '\n')
				BATCH_FILE.write('\n')
				BATCH_FILE.write("cd {}\n\n".format(SGEBATCH_PATH))
				BATCH_FILE.write("{}/RepeatMasker{}-inv -a -gff -s -pa {} {}.fa >& run.log\n".format(REPEATMASKER, ADD_PARAMS, str(PROC - 1), BATCH_NUMBER))

			QSUB.write("qsub {}\n".format(SGEBATCH)
