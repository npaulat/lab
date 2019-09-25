import sys
import os
import argparse
import itertools
import subprocess
import fnmatch
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
	#Argument of RepeatMasker run parameter
	parser.add_argument('-engine', type=str, help='RepeatMasker run parameter "-engine <search_engine>" option; select a non-default search engine to use, otherwise RepeatMasker will used the default configured at install time; [crossmatch|abblast|rmblast|hmmer]', choices=['crossmatch', 'abblast', 'rmblast', 'hmmer'], required=False)
	#Argument of RepeatMasker run parameter
	parser.add_argument('-inv', type=str, help='RepeatMasker parameter flag "-inv" option; alignments are presented in the orientation of the repeat', action='store_true')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-nolow', type=str, help='RepeatMasker parameter flag "-nolow" option; does not mask low complexity DNA or simple repeats', action='store_true')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-s', '-speed', type=str, help='RepeatMasker run parameter "-q" or "-s" option; q=quick search; 5-10% less sensitive, 3-4 times faster than default; s=slow search; 0-5% more sensitive, 2.5 times slower than default', choices=['q', 's'], required=False)
	#Argument of RepeatMasker run parameter
	parser.add_argument('-div', type=int, help='RepeatMasker run parameter "-div [number]" option; masks only those repeats that are less than [number] percent diverged from the consensus sequence', required=False)
	
	args = parser.parse_args()
	GENOME = args.input
	SPECIES = args.species
	BATCH_COUNT = args.batch_count
	GENOME_DIR = args.genome_dir
	OUTDIR = args.outdir
	QUEUE = args.queue
	LIBRARY = args.lib
	XSMALL = args.xsmall
	ENGINE = args.engine
	INV = args.inv
	NOLOW = args.nolow
	SPEED = args.speed
	DIV = args.div
	
	return GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV
	
GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV = get_args()

# Sanity checks
print("The species is {}, the query genome is {}.\n").format(SPECIES, GENOME)
print("{} batches will be made.\n").format(str(BATCH_COUNT))
print("The genome FASTA is located in '{}'.\n").format(GENOME_DIR)
print("The output directory is '{}'.\n").format(OUTDIR)
print("The job queue is {}.\n").format(QUEUE)

if not SPECIES or LIBRARY:
	sys.exit("Must supply value for option 'species' or 'lib'!")
if SPECIES and LIBRARY:
	sys.exit("Only supply a value for one option: 'species' or 'lib'! Not both!")

FLAGS = [LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV]
if not FLAGS:
	print("All default RepeatMasker parameters were used, no custom library.")
else:
	print("Custom parameters used:\n")
	if XSMALL:
		print("-xsmall flag used.\n")
	if INV:
		print("-inv flag used.\n")
	if NOLOW:
		print("-nolow flag used.\n")
	if LIBRARY:
		print("-lib flag used. Custom library is '{}'.\n").format(os.path.basename(LIBRARY))
	if ENGINE:
		print("-engine flag used. Changed search engine to {}.\n").format(ENGINE)
	if SPEED:
		print("-{} flag used. Search sensitivity has changed.\n").format(SPEED)
	if DIV:
		print("-div flag used. RepeatMasker will mask only repeats that are less than {}% diverged from the consensus sequence.\n").format(str(DIV))


if not os.path.isdir(GENOME_DIR):
	sys.exit("The given genome directory, '{}', does not exist.").format(GENOME_DIR)

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
		sys.exit("The genome file, '{}', is empty.").format(GENOME_FASTA)
except OSError as e:
	sys.exit("The genome file '{}' does not exist or is inaccessible.").format(GENOME_FASTA)
	
try:
	if not os.path.getsize(LIBRARY) > 0:
		sys.exit("The library file, '{}', is empty.").format(LIBRARY)
except OSError as e:
	sys.exit("The library file '{}' does not exist or is inaccessible.").format(LIBRARY)
	
if not os.path.isdir(OUTDIR):
	sys.exit("The output directory '{}' does not exist.").format(OUTDIR)

PARTITION_DIR = os.path.join(GENOME_DIR, "RMPart")

SLOTS_PER_BATCH = 10
MAX_DIR_SIZE = 1000
NUM_BATCHES = BATCH_COUNT

PARTITION_DIR = os.path.abspath(PARTITION_DIR)
if not os.listdir(PARTITION_DIR):
	print("{} is empty. Continuing.").format(PARTITION_DIR)
else:
	print("{} is not empty. Removing contents and continuing.").format(PARTITION_DIR)
	os.remove(os.path.join(PARTITION_DIR, '*'))

def get_batches():
	# Return a 3-level list(ref): partitions -> chunks -> chunk properties (scaffold + coordinates)
	PARTS = []
	
	## Open up the genome file and tabulate the sequence sizes
	COMMAND = BIN + "/twoBitInfo " + GENOME_FASTA + " stdout |"
	#open(P, $COMMAND)
		#|| die "Couldn't open pipe ($COMMAND): $_\n"
	TOTAL_SIZE = 0
	SEQS = ()
	#while (<P>) {
	#	chomp;
	#	my ($seq, $seqsize) = split("\t");
	# ....
	
	if NUM_BATCHES > 0:
		CHUNK_SIZE = int(TOTAL_SIZE / NUM_BATCHES) + 1
		
	BATCHES = ()
	CURRENT_BATCH_SIZE = 0
	for SEQ in SEQS:
		SEQ_SIZE = SEQ_SIZES{SEQ}
		SEQ_IDX = 0
		while SEQ_SIZE > 0:



