import argparse
import os
import sys
import re
import subprocess
import pandas as pd
import numpy as np
import pyfaidx
from pyfaidx import Fasta
#from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import logging
import time
import datetime
pd.options.mode.chained_assignment = None  # default='warn'

LOGGER = logging.getLogger(__name__)

#############################################################################################
##### extract.align.py - A program to extract TEs from a genome. 						#####
##### Required values: 																	#####
#####		--genome  input genome file (in fasta format) 								#####
#####		--blast  input blast file (tab format [6]) 									#####
#####		--consTEs  file of consensus elements 										#####
#####		--seqbuffer  additional 5'/3' bases in extracted sequence (default 500) 	#####
#####		--seqnum  number of sequences to extract (default 50) 						#####
##### Optional values: 																	#####
#####		--align  aligns the sequences with MUSCLE (default off) 					#####
#####		--index  genome .fai file 													#####
#############################################################################################

##### FUNCTIONS #####

# Set up input arguments
def get_args():
	# What this script does
	parser = argparse.ArgumentParser(description="extract.align.py - a program to extract TEs from a genome. Required values: --genome  input genome file (in fasta format), --blast  input blast file (tab format [6]), --consTEs  file of consensus elements, --seqbuffer  additional 5'/3' bases in extracted sequence (default 500), --seqnum  number of sequences to extract (default 50); Optional values: --align  aligns the sequences with MUSCLE (default off)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	# Argument of the input genome file
	parser.add_argument('-g', '--genome', type=str, help='Name of your input genome file, include path if outside current directory', required=True)
	# Argument for the input genome index file, if not provided, is generated
	parser.add_argument('-in', '--index', type=str, help='Name of your input genome index .fai file, include path if outside current directory', required=False)
	# Argument of the input blast file
	parser.add_argument('-b', '--blast', type=str, help='Name of your input blast file', required=True)
	# Argument for the library of consensus TEs
	parser.add_argument('-c', '--consTEs', type=str, help='Name of your consensus TE library', required=True)
	# Argument for length of sequence buffer to extract, default 500 bp on either side
	#parser.add_argument('-s', '--seqbuffer', type=int, help='Option to specify length of sequence buffer (flanking bp) to extract', default=500)
	# Alternative argument for sequence buffer, upstream only
	parser.add_argument('-lb', '--leftbuffer', type=int, help='Left buffer size. The number of bp of flanking sequence for each hit to be extracted along with the hit. Optional, Default = 1000', default = 1000)
	# Alternative argument for sequence buffer, downstream only
	parser.add_argument('-rb', '--rightbuffer', type=int, help='Right beffer size. The number of bp of flanking sequence for each hit to be extracted along with the hit. Optional, Default = 1000', default = 1000)
	# Argument for number of sequences (top matches) to extract, default 50
	parser.add_argument('-n', '--seqnum', type=int, help='Option to specify number of to sequences to extract', default=50)
	# Argument for MUSCLE alignment, default off
	parser.add_argument('-a', '--align', help='Option to align extracted BLAST hits, default no', action='store_true')
	# Argument for output directory, default to current directory
	parser.add_argument('-od', '--outdir', help='Output directory path, default is current directory', default=".")
	# Set up logger
	parser.add_argument("-log", "--log_level", default="INFO")
	
	args = parser.parse_args()
	GENOME = args.genome
	if args.index:
		INDEX = args.index
	else:
		INDEX = 'No'
	BLAST = args.blast
	CONSTES = args.consTEs
	#BUFFER = args.seqbuffer
	LBUFFER = args.leftbuffer
	RBUFFER = args.rightbuffer
	SEQNUM = args.seqnum
	if args.align:
		ALIGN = 'Yes'
	else:
		ALIGN = 'No'
	OUTDIR = args.outdir
	LOG = args.log_level
	return GENOME, INDEX, BLAST, CONSTES, LBUFFER, RBUFFER, SEQNUM, ALIGN, OUTDIR, LOG
	#return GENOME, INDEX, BLAST, CONSTES, BUFFER, SEQNUM, ALIGN, OUTDIR
#GENOME, INDEX, BLAST, CONSTES, LBUFFER, RBUFFER, SEQNUM, ALIGN, OUTDIR, LOG = get_args()
#GENOME, INDEX, BLAST, CONSTES, BUFFER, SEQNUM, ALIGN, OUTDIR = get_args()

#def REVTRANS(SEQ):
#	SEQUENCE = Seq(SEQ, generic_dna)
#	REVSEQ = SEQUENCE.reverse_complement()
#	return REVSEQ

# Index genome (faidx to generate .fai file)
def INDEX_GENOME(OUTDIR, GENOME_FILE):
	LOGGER.info('Indexing the genome')
	GENOMEIDX = Fasta(GENOME_FILE)
	GENOMEPREFIX = os.path.splitext(GENOME_FILE)[0]
	FAIDX = pd.read_csv(GENOME_FILE + '.fai', sep='\t', names=['SCAFFOLD', 'SCAFF_LENGTH', 'three', 'four', 'five'])
	#FAIDX = FAIDX[['SCAFFOLD', 'SCAFF_LENGTH']]
	FILE = GENOMEPREFIX + '.fai'
	INDEX = os.path.join(OUTDIR, FILE)
	FAIDX.to_csv(INDEX, sep='\t', header=False, index=False)
	return INDEX

# def CREATE_TE_OUT_FILES(OUTDIR, CONS_TES):
	# LOGGER.info('Parsing TE library')
	# for SEQ_RECORD in SeqIO.parse(CONS_TES, "fasta"):
		# SEQ_ID = re.sub('#', '__', SEQ_RECORD.id)
		# SEQ_ID = re.sub('/', '___', SEQ_ID)
		# SEQ = str(SEQ_RECORD.seq)
		# OUT = SEQ_ID + ".fas"
		# OUTFILE = os.path.join(OUTDIR, OUT)
		# with open(OUTFILE, "w+") as OUTPUT:
			# #SeqIO.write(SEQ, OUTPUT, "fasta")
			# OUTPUT.write(">CONSENSUS-" + SEQ_ID + "\n" + SEQ + "\n")

def ORGANIZE_BLAST_HITS(OUTDIR, BLAST_HITS, LBUFFER, RBUFFER, SEQ_NUM, INDEX, GENOME):
	if INDEX == 'No':
		INDEX = INDEX_GENOME(OUTDIR, GENOME)
	LOGGER.info('Organizing BLAST hits')
## Read in blast data
	BLAST_DATA = pd.read_csv(BLAST_HITS, sep='\t', names=['QUERYNAME', 'SCAFFOLD', 'C', 'D', 'E', 'F', 'QUERYSTART', 'QUERYSTOP', 'SCAFSTART', 'SCAFSTOP', 'E-VALUE', 'BITSCORE'])
## Convert to bed format
	BLASTDATA = BLAST_DATA[['SCAFFOLD', 'SCAFSTART', 'SCAFSTOP', 'QUERYNAME', 'E-VALUE', 'BITSCORE']]
	BLASTDATA.insert(4, 'STRAND', '+')
	BLASTDATA.loc[BLASTDATA.SCAFSTOP < BLASTDATA.SCAFSTART, 'STRAND'] = '-'
	BLASTDATA.SCAFSTART, BLASTDATA.SCAFSTOP = np.where(BLASTDATA.SCAFSTART > BLASTDATA.SCAFSTOP, [BLASTDATA.SCAFSTOP, BLASTDATA.SCAFSTART], [BLASTDATA.SCAFSTART, BLASTDATA.SCAFSTOP])
## Generate list of query names
	QUERYLIST = BLASTDATA.QUERYNAME.unique()
## For each query, sort hits
	QUERY_CLUSTER_FILES = []
	BUFFED_QUERY_CLUSTER_FILES = []
	for QUERY in QUERYLIST:
		QUERY_HITS = BLASTDATA[BLASTDATA['QUERYNAME'] == QUERY].sort_values(by=['E-VALUE', 'BITSCORE'], ascending=[True, False])
		QUERY_HITS.reset_index(inplace=True, drop=True)
	## Select top SEQ_NUM hits, output QUERY BED file
		QUERY_HITS = QUERY_HITS.head(SEQ_NUM)
		OUT_FILE1 = QUERY + ".bed"
		OUTPUT1 = os.path.join(OUTDIR, OUT_FILE1)
		QUERY_HITS.to_csv(OUTPUT1, sep='\t', header=False, index=False)
	## Append QUERY output file to list
		QUERY_CLUSTER_FILES.append(OUTPUT1)
	## Add BUFFER to SCAFSTART, SCAFEND, extract sequences, generate BLAST_HITS FASTA file
		#QUERY_BUFFED = QUERY_HITS[['SCAFFOLD', 'SCAFSTART', 'SCAFSTOP', 'QUERYNAME', 'STRAND', 'E-VALUE', 'BITSCORE']]
		#QUERY_BUFFED.insert(2, 'BUFF_START', QUERY_BUFFED.SCAFSTART - LBUFFER)
		#QUERY_BUFFED.loc[QUERY_BUFFED.BUFF_START < 0, 'BUFF_START'] = 0
		#QUERY_BUFFED.insert(4, 'BUFF_STOP', QUERY_BUFFED.SCAFSTOP + RBUFFER)
		QUERY_HITS.SCAFSTART = QUERY_HITS.SCAFSTART - LBUFFER
		QUERY_HITS.SCAFSTOP = QUERY_HITS.SCAFSTOP + RBUFFER
		QUERY_HITS.loc[QUERY_HITS.SCAFSTART < 0, 'SCAFSTART'] = 0
		### Compare BUFF_END to actual END of the scaffold as listed in FAI
		with open(INDEX, "r") as i:
			SCAFFOLD_DICT = {}
			for LINE in i:
				RECORD = LINE.split('\t')
				SCAFFOLD = RECORD[0]
				SCAFF_LEN = int(RECORD[1])
				try:
					SCAFFOLD_DICT[RECORD[0]].append(SCAFF_LEN)
				except KeyError:
					SCAFFOLD_DICT[SCAFFOLD] = [SCAFF_LEN]
		for SCAFF_INFO in SCAFFOLD_DICT:
			#QUERY_BUFFED['BUFF_STOP'] = np.where(((QUERY_BUFFED['SCAFFOLD'] == SCAFF_INFO) & (QUERY_BUFFED['BUFF_STOP'] > SCAFFOLD_DICT[SCAFF_INFO])), SCAFFOLD_DICT[SCAFF_INFO], QUERY_BUFFED['BUFF_STOP'])
			QUERY_HITS['SCAFSTOP'] = np.where(((QUERY_HITS['SCAFFOLD'] == SCAFF_INFO) & (QUERY_HITS['SCAFSTOP'] > SCAFFOLD_DICT[SCAFF_INFO][0])), SCAFFOLD_DICT[SCAFF_INFO][0], QUERY_HITS['SCAFSTOP'])
		#QUERY_BUFFED = QUERY_BUFFED.drop(columns=['SCAFSTART', 'SCAFSTOP'])
		OUT_FILE2 = QUERY + "_buffered.bed"
		OUTPUT2 = os.path.join(OUTDIR, OUT_FILE2)
		#QUERY_BUFFED.to_csv(OUTPUT2, sep='\t', header=True, index=False)
		QUERY_HITS.to_csv(OUTPUT2, sep='\t', header=False, index=False)
		BUFFED_QUERY_CLUSTER_FILES.append(OUTPUT2)
	LOGGER.info('BLAST hits organized, buffers applied')
	return QUERY_CLUSTER_FILES, BUFFED_QUERY_CLUSTER_FILES

def EXTRACT_BLAST_HITS(OUTDIR, BUFFED_QUERY_CLUSTER_FILES, GENOME, CONS_TES):
	LOGGER.info('Extracting BLAST hits')
	TO_ALIGN_LIST = []
	RECORD_DICT = SeqIO.index(GENOME, "fasta")
	for QUERY_FILE in BUFFED_QUERY_CLUSTER_FILES:
		QUERY_DICT = {}
		with open(QUERY_FILE, "r") as q:
			for QUERY in q:
				LINE = QUERY.split('\t')
				SCAFFOLD = LINE[0]
				SCAFSTART = int(LINE[1])
				SCAFSTOP = int(LINE[2])
				QUERY_NAME = LINE[3]
				STRAND = LINE[4]
				try:
					QUERY_DICT[LINE[0]].append((SCAFSTART, SCAFSTOP, QUERY_NAME, STRAND))
				except KeyError:
					QUERY_DICT[SCAFFOLD] = [(SCAFSTART, SCAFSTOP, QUERY_NAME, STRAND)]
		OUT_FILE = os.path.basename(QUERY_FILE).split('.')[0] + ".fas"
		OUTPUT = os.path.join(OUTDIR, OUT_FILE)
		TO_ALIGN_LIST.append(OUTPUT)
		with open(OUTPUT, "w+") as WRITE_FILE:
			for SCAFF in QUERY_DICT:
				RECORD = RECORD_DICT[SCAFF]
				for HIT in range(len(QUERY_DICT[SCAFF])):
					START = QUERY_DICT[SCAFF][HIT][0] - 1
					STOP = QUERY_DICT[SCAFF][HIT][1]
					STRAND = QUERY_DICT[SCAFF][HIT][3]
					if STRAND == '+':
						QUERY_FASTA = RECORD[START:STOP]
						QUERY_FASTA.description = "{}:{}".format(START, STOP)
						WRITE_FILE.write('>' + QUERY_FASTA.id + '\n')
						WRITE_FILE.write(str(QUERY_FASTA.seq) + '\n')
					elif STRAND == '-':
						QUERY_FASTA = RECORD[START:STOP]
						#QUERY_FASTA_RC = QUERY_FASTA.reverse_complement()
						#QUERY_FASTA.description = "{}:{}".format(STOP, START)
						WRITE_FILE.write('>' + QUERY_FASTA.id + '\n')
						WRITE_FILE.write(str(QUERY_FASTA.reverse_complement().seq) + '\n')
					#for RECORD in SeqIO.parse(GENOME, "fasta"):
						#if SCAFF == RECORD.id:
							#if STRAND == '+':
								#QUERY_FASTA = RECORD[START:STOP]
								#QUERY_FASTA.description = "{}:{}".format(START, STOP)
								#WRITE_FILE.write('>' + QUERY_FASTA.id + '\n')
								#WRITE_FILE.write(str(QUERY_FASTA.seq) + '\n')
							#elif STRAND == '-':
							#	QUERY_FASTA = RECORD[START:STOP]
								#QUERY_FASTA_RC = QUERY_FASTA.reverse_complement()
								#QUERY_FASTA.description = "{}:{}".format(STOP, START)
							#	WRITE_FILE.write('>' + QUERY_FASTA.id + '\n')
							#	WRITE_FILE.write(str(QUERY_FASTA.reverse_complement().seq) + '\n')
	#Add TE consensus sequence to the query top hits .fas file
	LOGGER.info('Parsing TE library')
	for SEQ_RECORD in SeqIO.parse(CONS_TES, "fasta"):
		SEQ_ID = re.sub('#', '__', SEQ_RECORD.id)
		SEQ_ID = re.sub('/', '___', SEQ_ID)
		SEQ = str(SEQ_RECORD.seq)
		LOGGER.info('Appending TE consensus sequences to BLAST hit sequence files')
		for UNALIGNED_FILE in TO_ALIGN_LIST:
			TE_NAME = os.path.basename(UNALIGNED_FILE).split('.')[0].split('_buffered')[0]
			if SEQ_ID == TE_NAME:
				with open(UNALIGNED_FILE, "a+") as SEQ_FILE:
					SEQ_FILE.write('>CONSENSUS-' + RECORD.id + '\n')
					SEQ_FILE.write(SEQ + '\n')
	# for UNALIGNED_FILE in TO_ALIGN_LIST:
		# TE_NAME = os.path.basename(UNALIGNED_FILE).split('.')[0]
		# with open(UNALIGNED_FILE, "a+") as SEQ_FILE:
			# for RECORD in SeqIO.parse(CONS_TES, "fasta"):
				# if TE_NAME == RECORD.id:
					# SEQ_FILE.write('>CONSENSUS-' + RECORD.id + '\n')
					# SEQ_FILE.write(str(RECORD.seq) + '\n')
	return TO_ALIGN_LIST

def ALIGN_BLAST_HITS(OUTDIR, TO_ALIGN_LIST):
	SUBDIR = "muscle"
	ALIGN_DIR = os.path.join(OUTDIR, SUBDIR)
	os.mkdir(ALIGN_DIR)
	LOGGER.info('Aligning TE BLAST hits with MUSCLE')
	ALIGNED_FILES = []
	for UNALIGNED_FILE in TO_ALIGN_LIST:
		OUT_FILE = os.path.basename(UNALIGNED_FILE).split('.')[0] + ".muscle.fas"
		ALIGNED_FILE = os.path.join(ALIGN_DIR, OUT_FILE)
		SOFTWARE = "/lustre/work/daray/software/muscle/"
		ALIGN_COMM = "{}muscle -in {} -out {}".format(SOFTWARE, UNALIGNED_FILE, ALIGNED_FILE)
		subprocess.run(ALIGN_COMM, capture_output=True, shell=True)
		ALIGNED_FILES.append(ALIGNED_FILE)
	LOGGER.info('All MUSCLE alignments completed')
	return ALIGNED_FILES

##Consensus generation function
# def CONSENSUSGEN(OUTDIR, ALIGNED_FILES):
	# SUBDIR = "consensusfiles"
	# OUT_DIR = os.path.join(OUTDIR, 'muscle')
	# CONS_DIR = os.path.join(OUT_DIR, SUBDIR)
	# os.mkdir(CONS_DIR)
	# os.chdir(OUT_DIR)
	# for ALIGNED_FILE in os.listdir('.'):
		# FILEPREFIX = os.path.splitext(ALIGNED_FILE)[0]
		# SOFTWARE = '/lustre/work/daray/software/'
		# CONS_FILE = FILEPREFIX + '_cons.fa'
		# OUT_FILE = os.path.join(CONS_DIR, CONS_FILE)
		# subprocess.run(SOFTWARE + 'trimal/source/trimal -in {} -gt 0.6 -cons 60 -fasta -out {}'.format(ALIGNED_FILE, FILEPREFIX + '_trimal.fa'), shell=True)
		# subprocess.run(SOFTWARE + 'EMBOSS-6.6.0/emboss/cons -sequence ' + FILEPREFIX + '_trimal.fa -outseq ' + FILEPREFIX + '_cons.fa -name ' + FILEPREFIX + '_cons -plurality 3 -identity 3', shell=True)
		# subprocess.run('cat {} {} >{}'.format(FILEPREFIX + '_trimal.fa', FILEPREFIX + '_cons.fa', OUT_FILE), shell=True)
		# os.remove(FILEPREFIX + '_cons.fa')
		# os.remove(FILEPREFIX + '_trimal.fa')
		
		
def CONSENSUSGEN(OUTDIR, ALIGNED_FILES):
	SUBDIR = 'consensusfiles'
	OUT_DIR = os.path.join(OUTDIR, 'muscle')
	CONS_DIR = os.path.join(OUT_DIR, SUBDIR)
	os.mkdir(CONS_DIR)
	os.chdir(OUT_DIR)
	#FILE_NUM = len([name for name in os.listdir('.') if os.path.isfile(name)])
	#for ALIGNED_FILE in os.listdir('.'):
	for ALIGNED_FILE in ALIGNED_FILES:
		FILEPREFIX = os.path.splitext(ALIGNED_FILE)[0]
		BASENAME = os.path.basename(ALIGNED_FILE).split(".fas")[0]
		SOFTWARE = '/lustre/work/daray/software/'
		CONS_FILE = BASENAME + '_cons.fa'
		OUT_FILE = os.path.join(CONS_DIR, CONS_FILE)
		retcode = subprocess.call(SOFTWARE + 'trimal/source/trimal -in {} -gt 0.6 -cons 60 -fasta -out {}'.format(ALIGNED_FILE, FILEPREFIX + '_trimal.fa'), shell=True)
		if retcode < 0:
			if retcode == -11:
				LOGGER.info('Segmentation fault error for ' + ALIGNED_FILE)
			else:
				LOGGER.info('Unexpected error in Trimal processing for ' + ALIGNED_FILE)
		else:
			subprocess.run(SOFTWARE + 'EMBOSS-6.6.0/emboss/cons -sequence ' + FILEPREFIX + '_trimal.fa -outseq ' + FILEPREFIX + '_cons.fa -name ' + FILEPREFIX + '_cons -plurality 3 -identity 3', shell=True)
			subprocess.run('cat {} {} > {}'.format(FILEPREFIX + '_trimal.fa', FILEPREFIX + '_cons.fa', OUT_FILE), shell=True)
			os.remove(FILEPREFIX + '_cons.fa')
			os.remove(FILEPREFIX + '_trimal.fa')

def main():
	# Set-up tasks:
	## Get input arguments
	GENOME, INDEX, BLAST, CONSTES, LBUFFER, RBUFFER, SEQNUM, ALIGN, OUTDIR, LOG = get_args()
	
	if OUTDIR == '.':
		OUTDIR = os.getcwd()
	
	# Setup logging and script timing
	handlers = [logging.FileHandler('extract_align.log'), logging.StreamHandler()]
	logging.basicConfig(format='', handlers = handlers)
	logging.getLogger().setLevel(getattr(logging, LOG.upper()))

	start_time = time.time()

	LOGGER.info('#\n# extract_align.py\n#')
	LOGGER.info('Genome file: ' + GENOME)
	LOGGER.info('Blast file: ' + BLAST)
	LOGGER.info('TE library: ' + CONSTES)
	LOGGER.info('Left buffer size: ' + str(LBUFFER))
	LOGGER.info('Right buffer size: ' + str(RBUFFER))
	LOGGER.info('Number of hits evaluated: ' + str(SEQNUM))
	LOGGER.info('Muscle alignment = ' + ALIGN)
	#LOGGER.info('Trimal/Emboss consensus = ' + EMBOSS)
	LOGGER.info('Trimal/Emboss consensus generation')
	LOGGER.info('Log level: ' + LOG)
	
	# Main steps:
	## Step 1: Organizing blast output-------------------------------------
	QUERY_CLUSTER_FILES, BUFFED_QUERY_CLUSTER_FILES = ORGANIZE_BLAST_HITS(OUTDIR, BLAST, LBUFFER, RBUFFER, SEQNUM, INDEX, GENOME)
	LOGGER.info('All buffered BLAST hit coordinate files generated.')

	## Step 2: Creating ouput files with the consensus --------------------
	#CREATE_TE_OUT_FILES(OUTDIR, CONSTES)
	#LOGGER.info('All consensus TE fasta files generated.')

	## Step 3: Extracting top blast hit sequences -------------------------
	TO_ALIGN_LIST = EXTRACT_BLAST_HITS(OUTDIR, BUFFED_QUERY_CLUSTER_FILES, GENOME, CONSTES)
	LOGGER.info('All BLAST hit sequences extracted.')

	## Step 4: Align if desired--------------------------------------------
	if ALIGN == 'Yes':
		ALIGNED_FILES = ALIGN_BLAST_HITS(OUTDIR, TO_ALIGN_LIST)
		LOGGER.info('All BLAST hit sequence files aligned.')
		CONSENSUSGEN(OUTDIR, ALIGNED_FILES)
		LOGGER.info('All new query consensus sequences generated.')
	
	end_time = time.time()
	LOGGER.info('Run time: ' + str(datetime.timedelta(seconds=end_time-start_time)))

#
# Wrap script functionality in main() to avoid automatic execution
# when imported ( e.g. when help is called on file )
#
if __name__ =="__main__":main()
	
