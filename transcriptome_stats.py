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
	#HALF_ASS = round(sum(LENLIST)//2)
	HALF_ASS = sum(LENLIST)/2
	print('Half of the assembly length = ' + str(HALF_ASS) + '.')
	
	for ENTRY in range(1, TOTALENTRIES):
		RUNNINGTOTAL = sum(LENLIST[0:ENTRY])
#		Why does this step take forever??
		if RUNNINGTOTAL >= HALF_ASS:
			print('N50 = ' + str(LENLIST[ENTRY-1]) + '.')
			STATS.write('Transcript N50 = ' + str(LENLIST[ENTRY-1]) + '.\n')
			break

	for RECORD in SeqIO.parse(FASTA, 'fasta'):
		IDLIST.append(RECORD.id.rsplit('_',1))

	IDS = pd.DataFrame(IDLIST, columns = ['uni', 'iso'])
	COUNTUNI = len(IDS['uni'].unique())
	print('There are ' + str(COUNTUNI) + ' unigenes in the file.')
	STATS.write('There are ' + str(COUNTUNI) + ' unigenes in the file.\n')
	
	HITLIST = []
	HITSLIST = []
	
	#from Bio.Blast import NCBIXML
	BLAST_FILE = open(BLASTX)
	#BLAST_HITS = NCBIXML.parse(BLAST_FILE)
	
	#if RECORD in FILECONTENTS:
	#	hitcount = hitcount +1
	
	for HIT in BLAST_FILE:
		HITLIST.append(HIT.rsplit('_', 1)[1])
		for LINE in HITLIST:
			field = LINE.split()
			col1 = field[0]
			col2 = field[1]
			col3 = field[2]
				
			HITID = col2, col3
				
			HITSLIST.append(HITID)
		
	HITS = pd.DataFrame(HITSLIST, columns = ['uni', 'iso'])
	
	MATCHLIST = []
	
	for i in IDS:
		if IDS[i]==HITS[i]:
			MATCHLIST.append(i)
	
	print('There are ' + str(len(MATCHLIST)) + ' transcripts with BLAST hits in the reference proteome.\n')
	
	UNIPROT = len(MATCHLIST.unique())
	print('There are ' + str(UNIPROT) + 'unique proteins in the reference proteome with a BLAST hit.')
	
	df = pd.DataFrame(MATCHLIST, columns = ['uni', 'iso'])
	#to get the count of unique lines (ID and isoform) in the match list
	df.groupby(df.domain.str.strip("'"))['ID'].nunique()
	print('Collaspe factor = ' +str(mean((df))) + '\n')
	
	#Number of unigenes with more thatn 1(or 5 or 10) isoforms
	#UNI_ISO = IDS.unique()
	#for IDS['uni'].unique():
	#	if len(IDS['iso'].unique()) > 1:
	#		UNI_ISO.append(IDS)
	#for IDS['uni'].unique() in UNI_ISO:
	#	if len(IDS['iso'].unique()) == 1:
	#		UNIQ_ISO.append(IDS)
	#	elif len(IDS['iso'].unique()) > 1:
	#		MORE_ISO.append(IDS)
	#	elif len(IDS['iso'].unique()) > 5:
	#		PLUS5_ISO.append(IDS)
	#	elif len(IDS['iso'].unique()) > 10:
	#		PLUS10_ISO.append(IDS)
	#print length of each list
	
	#Percentage of cases where each isoform has the same best BLAST hit
	#for ISO in UNI_ISO:
	#	if ISO == BESTHIT:
	#		list.append('YES')
	#	else:
	#		list.append('NO')
	# RECP = # yes lines/len(list)*100
	#print('Cases where each isoform has the same best BLAST hit = ' = str(RECP))
	#Annotation statistics for only the longest isoform of each uni
