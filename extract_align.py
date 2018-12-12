import sys
import os
import argparse
import itertools
import subprocess

#1. run blast: blastn –query <TE library> -db <genome> -out <arbitrary name of outfile> -outfmt 6 (is to make it tab delimited)
#2. run extractAlignTEs.pl: extractAlignTEs.pl –genome <genome> –blast <blast output file> –consTEs <TE library> --seqBuffer <# of bps flanking TE> --seqNum <# of closest blast hits> --align

parser=argparse.ArgumentParser()
parser.add_argument('-gen', '--genome', help='data input file - FASTA', required=True)
parser.add_argument('-loc', '--location', type=str, help='full location of working directory, i.e. no period', required=True)
parser.add_argument('-lib', '--library', type=str, help='The TE library to blast your genome with' , required=True)
parser.add_argument('-flank', '--flankbuffer', type=int, help='The number of bps flanking the blasted TEs', required=True)
parser.add_argument('-num', '--seqnum', type=int, help='The number of closest blast hits included in extractAlign output', required=True)
args=vars(parser.parse_args())

GENOME = args["genome"]
WORKING_DIR = args["location"]
LIBRARY = args["library"]
BUFFER = args["flankbuffer"]
HIT_LIM = args["seqnum"]

#WORKING_DIR = '/lustre/scratch/npaulat/extract_align/'
#LIBRARY = '/lustre/work/daray/extractAlign_assignment/starting_library.fas'
#GENOME = '/lustre/work/daray/extractAlign_assignment/cPer_m.fa.gz'
#OUTPUT = '/lustre/scratch/npaulat/extract_align/cPer_m_blast.out'

os.chdir(WORKING_DIR)

module_load = 'module load intel rmblast/2.6.0'
subprocess.run(['bash', '-c', module_load])

#subprocess.check_call('module load intel rmblast/2.6.0', shell=True)

#create BLAST database
subprocess.check_call('makeblastdb -in {} -dbtype nucl'.format(GENOME),shell=True)

#run blastn
subprocess.check_call('blastn -query {} -db {} -out blast_results -outfmt 6'.format(LIBRARY, GENOME), shell=True)

#def run_blastn(LIBRARY, GENOME, OUTPUT):
#	blasting = 'blastn -q ' + LIBRARY + ' -db ' + GENOME + ' -out ' + OUTPUT + ' -outfmt 6'
#	subprocess.run(['bash', '-c', blasting])

#BLAST_OUT = OUTPUT
#BUFFER = 300
#HIT_LIM = 5

#run Extract Align
subprocess.check_call('perl {}/extractAlignTEs.pl -genome {} -blast ExAlign_out -consTEs {} -seqBuffer {} -seqNum {} -align'.format(WORKING_DIR, GENOME, LIBRARY, BUFFER, HIT_LIM), shell=True)

#def extract_align(GENOME, BLAST_OUT, LIBRARY, BUFFER, HIT_LIM):
#	call_script = 'extractAlignTEs.pl -g ' + GENOME + ' -blast ' + BLAST_OUT + '-consTEs ' + LIBRARY + ' --seqBuffer ' + BUFFER + ' --seqNum ' + HIT_LIM + ' --align'
#	subprocess.run(['perl', '-c', call_script])

#run_blastn(LIBRARY, GENOME, OUTPUT)
#extract_align(GENOME, BLAST_OUT, LIBRARY, BUFFER, HIT_LIM)
