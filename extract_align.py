import sys
import os
import argparse
import itertools
import subprocess

module_load = 'module load intel rmblast/2.6.0'
subprocess.run(['bash', '-c', module_load])

#1. run blast: blastn –query <TE library> -db <genome> -out <arbitrary name of outfile> -outfmt 6 (is to make it tab delimited)
#2. run extractAlignTEs.pl: extractAlignTEs.pl –genome <genome> –blast <blast output file> –consTEs <TE library> --seqBuffer <# of bps flanking TE> --seqNum <# of closest blast hits> --align

LIBRARY = '/lustre/work/daray/extractAlign_assignment/starting_library.fas'
GENOME = '/lustre/work/daray/extractAlign_assignment/cPer_m.fa.gz'
OUTPUT = '/lustre/scratch/npaulat/extract_align/cPer_m_blast.out'

def run_blastn(LIBARY, GENOME, OUTPUT):
	blasting = 'blastn -q ' + LIBRARY + ' -db ' + GENOME + ' -out ' + OUTPUT + ' -outfmt 6'
	subprocess.run(['bash', '-c', blasting])

BLAST_OUT = OUTPUT
BUFFER = 300
HIT_LIM = 5

def extract_align(GENOME, BLAST_OUT, LIBRARY, BUFFER, HIT_LIM):
	call_script = 'extractAlignTEs.pl -g ' + GENOME + ' -blast ' + BLAST_OUT + '-consTEs ' + LIBRARY + ' --seqBuffer ' + BUFFER + ' --seqNum ' + HIT_LIM + ' --align'
	subprocess.run(['perl', '-c', call_script])

run_blastn(LIBRARY, GENOME, OUTPUT)
extract_align(GENOME, BLAST_OUT, LIBRARY, BUFFER, HIT_LIM)