from Bio import SeqIO
import statistics as stat
import pandas as pd
import argparse

def get_args():
	parser = argparse.ArgumentParser(description="Python seminar assignment. Will parse transcriptome data.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--fasta', type=str, help='Name of the peptide file', required=True)
	parser.add_argument('-b', '--blastx', type=str, help='Name of the blastx file', required=True)
	parser.add_argument('-p', '--prefix', type=str, help='A prefix to label the output', required=True)

	args = parser.parse_args()
	FASTA = args.fasta
	BLASTX = args.blastx
	PREFIX = args.prefix

	return FASTA, BLASTX, PREFIX
FASTA, BLASTX, PREFIX = get_args()

print("Input fasta file is " + FASTA + ".")
print("Input blastx file is " + BLASTX + ".")
print("Prefix = " + PREFIX + ".")

LENLIST = []
IDLIST = []

with open(PREFIX + '.stats.txt', 'w') as STATS:

	for RECORD in SeqIO.parse(FASTA, 'fasta'):
		LENLIST.append(len(RECORD.seq))
			
	print('Total bases = ' + str(sum(LENLIST)) + '.' )
	STATS.write('Total bases = ' + str(sum(LENLIST)) + '.\n')
	print('Total transcripts = ' + str(len(LENLIST)) + '.')
	STATS.write('Total transcripts = ' + str(len(LENLIST)) + '.\n')
	print('Mean length of transcripts = ' + str(stat.mean(LENLIST)) + '.')
	STATS.write('Mean length of transcripts = ' + str(stat.mean(LENLIST)) + '.\n')
	print('Median length of transcripts = ' + str(stat.median(LENLIST)) + '.')
	STATS.write('Median length of transcripts = ' + str(stat.median(LENLIST)) + '.\n')
	
	LENLIST = sorted(LENLIST)
	MEAN = stat.mean(LENLIST)
	TOTALENTRIES = len(LENLIST)
	LONGEST = max(LENLIST)
	print('The longest transcript = ' + str(LONGEST) + '.') 
	SHORTEST = min(LENLIST)
	print('The shortest transcript = ' + str(SHORTEST) + '.')
	#HALF_ASS = float(len(LENLIST)//2)
	#HALF_ASS = int(sum(LENLIST)//2)
	HALF_ASS = round(sum(LENLIST)//2)
	print('Half of the assembly length = ' + str(HALF_ASS) + '.')
	
	for ENTRY in range(1, TOTALENTRIES):
		RUNNINGTOTAL = sum(LENLIST[0:ENTRY])
#		Why does this step take forever??
		if RUNNINGTOTAL >= HALF_ASS:
			print('Transcript N50 = ' + str(LENLIST[ENTRY-1]) + '.')
			STATS.write('Transcript N50 = ' + str(LENLIST[ENTRY-1]) + '.\n')
			break

	for ENTRY in range(1, TOTALENTRIES):
		RUNNINGTOTAL = sum(LENLIST[0:ENTRY])
#		Why does this seem wrong?
		if RUNNINGTOTAL >= MEAN:
			print('N50 = ' + str(LENLIST[ENTRY-1]) + '.')
			STATS.write('N50 = ' + str(LENLIST[ENTRY-1]) + '.\n')
			break

	for RECORD in SeqIO.parse(FASTA, 'fasta'):
		IDLIST.append(RECORD.id.rsplit('_',1))

	IDS = pd.DataFrame(IDLIST, columns = ['uni', 'iso'])
	COUNTUNI = len(IDS['uni'].unique())
	print('There are ' + str(COUNTUNI) + ' unigenes in the file.')
	STATS.write('There are ' + str(COUNTUNI) + ' unigenes in the file.\n')