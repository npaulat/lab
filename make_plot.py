import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import seaborn as sns
from pylab import savefig

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="General removal of most MELT duplicate calls and overlapping calls", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#required = parser.add_argument_group('required arguments')
	#Argument of the input data file
	parser.add_argument('-i', '--input', help='f', required=True)
	#Argument of the output file name
	parser.add_argument('-p', '--prefix', help='Option to specify output name, defaults to basename of input file')
	#Argument for the color scheme of the figure
	parser.add_argument('-color', '--color', type=str, help='Color options include seaborn palettes ("deep, muted, bright, pastel, dark, colorblind") or matlibplot palettes', default='colorblind')
	#Argument for figure style
	parser.add_argument('-style', '--style', type=str, help='Option to specify graph style, white, dark, whitegrid, darkgrid, or ticks, default ticks', default='ticks')
	#Argument for figure type, e.g. line graph
	parser.add_argument('-type', '--type', type=str, help='Define type of graph to make, default line graph', default='line')
	#Argument for x-axis label
	parser.add_argument('-x', '--xaxis', help='Option to specify x-axis label, default uses file header', default='default')
	#Argument for y-axis label
	parser.add_argument('-y', '--yaxis', help='Option to specify y-axis label, default uses file header', default='default')
	#dpi
	parser.add_argument('-dpi', '--dpi', type=int, help='Option for specifying dpi of the figure', default=300)
	#legend location
	parser.add_argument('-l', '--legend', type=str, help='Option for figure legend, default False', default='False')
	parser.add_argument('-ll', '--leg_loc', type=str, help='Option for specifying figure legend location in relation to the figure, options: best, (lower/upper/center) left, [lower/upper/center] right, [lower/upper] center', default='out')
	parser.add_argument('-lo', '--leg_out', type=str, help='Option to set legend location outside of the figure, default False', default='False')
	#Argument for figure title
	parser.add_argument('-t', '--title', type=str, help='Option to add a title to the figure, write spaces as underscores, defaults to none', default='default')
	#adjust y-axis
	#save output as other file type
	#Argument of the output directory (need full path)
	parser.add_argument('-od', '--outdir', type=str, help='Give the output directory path, defaults to current directory')
	
	args = parser.parse_args()
	FILE = args.input
	if args.prefix:
		PREFIX = args.prefix
	else:
		BASENAME = os.path.basename(FILE).split(".")[0]
		PREFIX = BASENAME
	COLOR = args.color
	STYLE = args.style
	TYPE = args.type
	XAXIS = args.xaxis
	YAXIS = args.yaxis
	DPI = args.dpi
	LEG = args.legend
	LEG_LOC = args.leg_loc
	TITLE = args.title
	if args.outdir:
		OUT_DIR = args.outdir
	else:
		OUT_DIR = '.'
	
	return FILE, PREFIX, COLOR, STYLE, TYPE, XAXIS, YAXIS, DPI, LEG, LEG_LOC, TITLE, OUT_DIR
	
FILE, PREFIX, COLOR, STYLE, TYPE, XAXIS, YAXIS, DPI, LEG, LEG_LOC, TITLE, OUT_DIR = get_args()

#define files, file paths
BASENAME = os.path.basename(FILE).split(".")[0]
FIGURE = PREFIX + ".png"
OUTPUT = os.path.join(OUT_DIR, FIGURE)

df = pd.read_table(FILE)

column1 = df.columns[0]
column2 = df.columns[1]

sns.set_style(STYLE)
sns.set_palette(COLOR)

#sns.set_context("paper")
#sns.set_context("poster")

if LEG == 'False':
	#fig = sns.relplot(x=column1, y=column2, kind=TYPE, data=df)
	fig = sns.relplot(x=column1, y=column2, kind=TYPE, data=df)
else:
	if LEG_LOC == 'out':
		#fig = sns.relplot(x=column1, y=column2, kind=TYPE, data=df, legend=LEG)
		fig = sns.relplot(x=column1, y=column2, data=df, hue=column2, legend=LEG)
		#plt.legend(loc=LEG_LOC)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		#plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	else:
		fig = sns.relplot(x=column1, y=column2, data=df, hue=column2, legend=LEG)
		plt.legend(loc=LEG_LOC)

sns.despine()

if XAXIS != 'default':
	XAXIS = XAXIS.replace("_"," ")
	fig.set_ylabels(XAXIS)
if YAXIS != 'default':
	YAXIS = YAXIS.replace("_"," ")
	fig.set_xlabels(YAXIS)
	
if TITLE != 'default':
	TITLE = TITLE.replace("_"," ")
	plt.title(TITLE)

fig.savefig(OUTPUT, dpi=DPI)


##sns.relplot(x=None, y=None, hue=None, size=None, style=None, data=None, row=None, col=None, col_wrap=None, row_order=None, col_order=None, palette=None, hue_order=None, hue_norm=None, sizes=None, size_order=None, size_norm=None, markers=None, dashes=None, style_order=None, legend='brief', kind='line', height=5, aspect=1, facet_kws=None, **kwargs)