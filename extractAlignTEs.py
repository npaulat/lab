from Bio import SearchIO
for Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import argparse
import os
import sys
import pandas as pd


def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="extract.align.py - a program to extract TEs from a genome. Required values: --genome  input genome file (in fasta format), --blast  input blast file (tab format [6]), --consTEs  file of consensus elements, --seqbuffer  additional 5'/3' bases in extracted sequence (default 500), --seqnum  number of sequences to extract (default 50); Optional values: --align  aligns the sequences with MUSCLE (default off)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Argument of the input genome file
	parser.add_argument('-g', '--genome', help='Name of your input genome file, include path if outside current directory', required=True)
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
	
	args = parser.parse_args()
	GENOME = args.genome
	BLAST = arg.blast
	CONSTES = args.consTEs
	BUFFER = args.seqbuffer
	SEQNUM = args.seqnum
	if args.align:
		ALIGN = 'Yes'
	else:
		ALIGN = 'No'
	
	return GENOME, BLAST, CONSTES, BUFFER, SEQNUM, ALIGN
	
GENOME, BLAST, CONSTES, BUFFER, SEQNUM, ALIGN = get_args()


# Step 1: Organizing blast output-------------------------------------
#organize_blast_hits(BLAST, BUFFER, SEQNUM);

# Step 2: Creating ouput files with the consensus --------------------
#create_te_out_files(CONSTEs);

# Step 3: Loading genome into memory----------------------------------


# Step 4: Align if desired--------------------------------------------



##### FUNCTIONS #####

def REVTRANS(SEQ):
	SEQUENCE = Seq(SEQ, generic_dna)
	REVSEQ = SEQUENCE.reverse_complement()
	return REVSEQ
	
#def ORGANIZE_BLAST_HITS(BLAST_HITS, SEQ_BUFFER, SEQ_NUM):
	
	
	
	
	
	
#def CREATE_TE_OUT_FILES(CONSENSUS_TES):
	
	
	
	
	
	
	
	
	
	