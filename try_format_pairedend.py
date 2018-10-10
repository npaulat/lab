###Created on Thurs Sep 6 2018
###@author: jessica_mfs + nicole.paulat@ttu.edu

import os
import argparse
import sys
import subprocess
import itertools

FILE = open(sys.argv[1]).readlines()
FQ_LENGTH = len(FILE)
with open(sys.argv[2],'w') as OUT:
	i = 0
	NUMBER = 1
	LENGTH = 12
	#i,OUT,number,length = 0,[],1,12
	INDEX_NUM = 0
	while i < (FQ_LENGTH):
		if INDEX_NUM %2 == 0:
			OUT.write("@" + '{:0>12}'.format(NUMBER) + "\tOP:i:1\n")
			OUT.write(FILE[i+1])
			OUT.write("+\n")
			OUT.write(FILE[i+3])
			i += 4
			INDEX_NUM += 1
		else:
			OUT.write("@" + '{:0>12}'.format(NUMBER) + "\tOP:i:2\n")
			OUT.write(FILE[i+1])
			OUT.write('+\n')
			OUT.write(FILE[i+3])
			i += 4
			INDEX_NUM += 1
			NUMBER += 1

OUT.close()
	#for line in OUT:
	#	OUT.write(line)