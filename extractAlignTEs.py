import argparse
import os
import sys
import re
import subprocess
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="extract.align.py - a program to extract TEs from a genome. Required values: --genome  input genome file (in fasta format), --blast  input blast file (tab format [6]), --consTEs  file of consensus elements, --seqbuffer  additional 5'/3' bases in extracted sequence (default 500), --seqnum  number of sequences to extract (default 50); Optional values: --align  aligns the sequences with MUSCLE (default off)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Argument of the input genome file
	parser.add_argument('-g', '--genome', help='Name of your input genome file, include path if outside current directory', required=True)
	#Argument for the input genome index file
	parser.add_argument('-in', '--index', help='Name of your input genome index .fai file, include path if outside current directory', required=False)
	#Argument of the input blast file
	parser.add_argument('-b', '--blast', help='Name of your input blast file', required=True)
	#Argument for the library of consensus TEs
	parser.add_argument('-c', '--consTEs', type=str, help='Name of your consensus TE library', required=True)
	#Argument for length of sequence buffer to extract
	parser.add_argument('-s', '--seqbuffer', type=int, help='Option to specify length of sequence buffer (flanking bp) to extract', default=500)
	#Argument for number of sequences (top matches) to extract
	parser.add_argument('-n', '--seqnum', type=str, help='Option to specify number of to sequences to extract', default=50)
	#Argument for x-axis label
	parser.add_argument('-a', '--align', help='Option to specify x-axis label, default uses file header', required=False)
	parser.add_argument('-od', '--outdir', help='Output directory path, default is current directory', default=".")
	
	args = parser.parse_args()
	GENOME = args.genome
	INDEX = args.index
	BLAST = arg.blast
	CONSTES = args.consTEs
	BUFFER = args.seqbuffer
	SEQNUM = args.seqnum
	if args.align:
		ALIGN = 'Yes'
	else:
		ALIGN = 'No'
	OUTDIR = args.outdir
	return GENOME, INDEX, BLAST, CONSTES, BUFFER, SEQNUM, ALIGN
	
GENOME, INDEX, BLAST, CONSTES, BUFFER, SEQNUM, ALIGN = get_args()

if OUTDIR == '.':
	OUTDIR = os.getcwd()

# Step 1: Organizing blast output-------------------------------------
#organize_blast_hits(BLAST, BUFFER, SEQNUM);
ORGANIZE_BLAST_HITS(BLAST, BUFFER, SEQNUM)

# Step 2: Creating ouput files with the consensus --------------------
CREATE_TE_OUT_FILES(CONSENSUS_TES)

# Step 3: Loading genome into memory----------------------------------		
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
	with open(OUTPUT, "w+") as WRITE_FILE:
		for SCAFF in QUERY_DICT:
			for HIT in range(len(QUERY_DICT[SCAFF])):
				START = QUERY_DICT[SCAFF][HIT][0] - 1
				STOP = QUERY_DICT[SCAFF][HIT][1]
				STRAND = QUERY_DICT[SCAFF][HIT][3]
				for RECORD in SeqIO.parse(GENOME, "fasta"):
					if SCAFF == RECORD.id:
						if STRAND == '+':
							QUERY_FASTA = RECORD[START:STOP]
							QUERY_FASTA.description = "{}:{}".format(START, STOP)
							f.write('>' + QUERY_FASTA.id + '\n')
							f.write(str(QUERY_FASTA.seq) + '\n')
						elif STRAND == '-':
							QUERY_FASTA = RECORD[START:STOP]
							#QUERY_FASTA_RC = QUERY_FASTA.reverse_complement()
							QUERY_FASTA_RC.id = QUERY_FASTA.id
							QUERY_FASTA_RC.description = "{}:{}".format(STOP, START)
							f.write('>' + QUERY_FASTA_RC.id + '\n')
							f.write(str(QUERY_FASTA_RC.reverse_complement().seq) + '\n')

# Step 4: Align if desired--------------------------------------------
if ALIGN == 'Yes':
	SUBDIR = "aligned"
	ALIGN_DIR = os.path.join(OUTDIR, SUB_DIR)
	os.mkdir(ALIGN_DIR)
	for UNALIGNED_FILE in TO_ALIGN_LIST:
		OUT_FILE = os.path.basename(UNALIGNED_FILE).split('.')[0] + ".muscle.fas"
		ALIGNED_FILE = os.path.join(ALIGN_DIR, OUT_FILE)
		ALIGN_COMM = "/lustre/work/daray/software/muscle/muscle -in {}.fas -out {}.muscle.fas".format(UNALIGNED_FILE, ALIGNED_FILE)
		subprocess.run(ALIGN_COMM, capture_output=True, shell=True)


##### FUNCTIONS #####

#def REVTRANS(SEQ):
#	SEQUENCE = Seq(SEQ, generic_dna)
#	REVSEQ = SEQUENCE.reverse_complement()
#	return REVSEQ

def CREATE_TE_OUT_FILES(CONS_TES):
	for SEQ_RECORD in SeqIO.parse(CONS_TES, "fasta"):
		SEQ_ID = re.sub('#', '__', SEQ_RECORD.id)
		SEQ_ID = re.sub('/', '___', SEQ_ID)
		SEQ = SEQ_RECORD.Seq
		OUTPUT = SEQ_ID + ".fas"
		with open(OUTFILE, "w+") as OUTPUT:
			#SeqIO.write(SEQ, OUTPUT, "fasta")
			OUTPUT.write(">CONSENSUS-" + SEQ_ID + "\n" + SEQ + "\n")
	
def ORGANIZE_BLAST_HITS(BLAST_HITS, SEQ_BUFFER, SEQ_NUM):
## Read in blast data
	BLAST_DATA = pd.read_csv(BLAST, sep='\t', names=['QUERYNAME', 'SCAFFOLD', 'C', 'D', 'E', 'F', 'QUERYSTART', 'QUERYSTOP', 'SCAFSTART', 'SCAFSTOP', 'E-VALUE', 'BITSCORE'])
## Convert to bed format
	BLASTBED = BLAST_DATA[['SCAFFOLD', 'SCAFSTART', 'SCAFSTOP', 'QUERYNAME', 'E-VALUE', 'BITSCORE']]
	BLASTBED.insert(4, 'STRAND', '+')
	BLASTBED.loc[BLASTBED.SCAFSTOP < BLASTBED.SCAFSTART, 'STRAND'] = '-'
	BLASTBED.SCAFSTART, BLASTBED.SCAFSTOP = np.where(BLASTBED.SCAFSTART > BLASTBED.SCAFSTOP, [BLASTBED.SCAFSTOP, BLASTBED.SCAFSTART], [BLASTBED.SCAFSTART, BLASTBED.SCAFSTOP])
## Generate list of query names
	QUERYLIST = BLASTBED.QUERYNAME.unique()
## For each query, sort hits
	QUERY_CLUSTER_FILES = []
	BUFFED_QUERY_CLUSTER_FILES = []
	for QUERY in QUERYLIST:
		QUERY_HITS = BLASTBED[BLASTBED['QUERYNAME'] == QUERY].sort_values(by=['E-VALUE', 'BITSCORE'], ascending=[True, False])
		QUERY_HITS.reset_index(inplace=True, drop=True)
	## Select top SEQ_NUM hits, output QUERY BED file
		QUERY_HITS = QUERY_HITS.head(SEQ_NUM)
		OUT_FILE1 = QUERY + ".bed"
		OUTPUT1 = os.path.join(OUTDIR, OUT_FILE1)
		QUERY_HITS.to_csv(OUTPUT1, sep='\t', header=True, index=False)
	## Append QUERY output file to list
		QUERY_CLUSTER_FILES.append(OUTPUT1)
	## Add BUFFER to SCAFSTART, SCAFEND, extract sequences, generate BLAST_HITS FASTA file
		QUERY_BUFFED = QUERY_HITS[['SCAFFOLD', 'SCAFSTART', 'SCAFSTOP', 'QUERYNAME', 'STRAND', 'E-VALUE', 'BITSCORE']]
		QUERY_BUFFED.insert(2, 'BUFF_START', QUERY_BUFFED.SCAFSTART - SEQ_BUFFER)
		QUERY_BUFFED.loc[QUERY_BUFFED.BUFF_START < 0, 'BUFF_START'] = 0
		QUERY_BUFFED.insert(4, 'BUFF_STOP', QUERY_BUFFED.SCAFSTOP + SEQ_BUFFER)
		# Compare BUFF_END to actual END of the scaffold as listed in FAI
		with open(INDEX, "r") as i:
			SCAFFOLD_DICT = {}
			for LINE in INDEX:
				SCAFFOLD = LINE[0]
				SCAFF_LEN = LINE[1]
				try:
					SCAFFOLD_DICT[LINE[0]].append((SCAF_LEN))
				except KeyError:
					QUERY_DICT[SCAFFOLD] = [(SCAF_LEN)]
		for SCAFF_INFO in SCAFFOLD_DICT:
			QUERY_BUFFED['BUFF_STOP'] = np.where(((QUERY_BUFFED['SCAFFOLD'] == SCAFF_INFO) & (QUERY_BUFFED['BUFF_STOP'] > SCAFF_DICT[SCAFF_INFO])), SCAFF_DICT[SCAFF_INFO], QUERY_BUFFED['BUFF_STOP'])
		QUERY_BUFFED = QUERY_BUFFED.drop(columns=['SCAFSTART', 'SCAFSTOP'])
		OUT_FILE2 = QUERY + "_buffered.bed"
		OUTPUT2 = os.path.join(OUTDIR, OUT_FILE2)
		QUERY_BUFFED.to_csv(OUTPUT2, sep='\t', header=True, index=False)
		BUFFED_QUERY_CLUSTER_FILES.append(OUTPUT2)
