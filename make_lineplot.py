import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import seaborn as sns
from pylab import savefig

##### DEFINE YOUR ARGUMENTS FOR EVERY VARIABLE YOU WANT TO CHANGE ON YOUR LINE PLOT #####

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Make a line plot from your data file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Argument of the input data file
	parser.add_argument('-i', '--input', help='Name of your input file, include path if outside current directory', required=True)
	#Argument of the output file name [default input file basename]
	parser.add_argument('-p', '--prefix', help='Option to specify output name, defaults to basename of input file')
	#Argument for the color scheme of the figure [default colorblind]
	parser.add_argument('-color', '--color', type=str, help='Color options include seaborn palettes ("deep, muted, bright, pastel, dark, colorblind") or matlibplot palettes', default='colorblind')
	#Argument for figure style [default ticks]
	parser.add_argument('-style', '--style', type=str, help='Option to specify graph style, white, dark, whitegrid, darkgrid, or ticks, default ticks, best options are ticks or white', default='ticks')
	#Argument for figure type [default line]
	parser.add_argument('-type', '--type', type=str, help='Define type of graph to make, default line graph', default='line')
	#Argument for x-axis label [default data column header]
	parser.add_argument('-x', '--xaxis', help='Option to specify x-axis label, default uses file header', default='default')
	#Argument for y-axis label [default data column header]
	parser.add_argument('-y', '--yaxis', help='Option to specify y-axis label, default uses file header', default='default')
	#Argument to change file dpi [default 300]
	parser.add_argument('-dpi', '--dpi', type=int, help='Option for specifying dpi of the figure', default=300)
	#Arguments for including a legend [default no legend]
	parser.add_argument('-l', '--legend', type=str, help='Option to include figure legend (True), default False', default='False')
	#Argument for where to place the legend [default out (outside of figure plot)]
	parser.add_argument('-ll', '--leg_loc', type=str, help='Option for specifying figure legend location in relation to the figure, options: best, (lower/upper/center) left, [lower/upper/center] right, [lower/upper] center', default='out')
	#Argument for including a figure title [default none]
	parser.add_argument('-t', '--title', type=str, help='Option to add a title to the figure, write spaces as underscores, defaults to none', default='default')
	#Argument to manually change x-axis limits from auto-adjusted range
	parser.add_argument('-xl', '--xlim', type=str, help='Option to change x-axis limits from default range, write as pair, e.g. "0, 5", or as single to change only one side "0, "', default='default')
	#Argument to manually change y-axis limits from auto-adjusted range
	parser.add_argument('-yl', '--ylim', type=str, help='Option to change y-axis limits from default range, write as pair, e.g. "0, 5", or as single to change only one side "0,"', default='default')
	#Argument to save image as pdf instead of png [default file saves as png]
	parser.add_argument('-of', '--out_format', type=str, help='Option to save image as pdf instead of png file [pdf]', default='default')
	#Argument of the output directory (need full path) [default to current directory]
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
	XLIM = args.xlim
	YLIM = args.ylim
	OUT_FORMAT = args.out_format
	if args.outdir:
		OUT_DIR = args.outdir
	else:
		OUT_DIR = '.'
	
	return FILE, PREFIX, COLOR, STYLE, TYPE, XAXIS, YAXIS, DPI, LEG, LEG_LOC, TITLE, XLIM, YLIM, OUT_FORMAT, OUT_DIR
	
FILE, PREFIX, COLOR, STYLE, TYPE, XAXIS, YAXIS, DPI, LEG, LEG_LOC, TITLE, XLIM, YLIM, OUT_FORMAT, OUT_DIR = get_args()

INPUT = {}

## Define your output file name and directory ##
# Use your PREFIX argument to name your output file (FIGURE)
FIGURE = PREFIX + ".png"
# Make output directory if needed
if OUT_DIR != '.':
	try:
		os.mkdir(OUT_DIR)
		print("Image file name is " + FIGURE + " located at " + OUT_DIR)
	except FileExistsError:
		print("Image file name is " + FIGURE + " located at " + OUT_DIR)
else:
	print("Image file name is " + FIGURE)

# Use your OUT_DIR argument and FIGURE variable to send your output where you want
OUTPUT = os.path.join(OUT_DIR, FIGURE)

##### BUILDING THE LINE PLOT #####

## Make your file into a pandas dataframe ##
df = pd.read_table(FILE)

## Define your dataframe column headers as variables to use later ##
column1 = df.columns[0]
column2 = df.columns[1]

## Making your figure pretty ##
# Tell seaborn to use your desired graph style and color palette
sns.set_style(STYLE)
sns.set_palette(COLOR)

#sns.set_context("paper")
#sns.set_context("poster")

## Making the figure given various arguments ##
# If you don't ask for a legend
if LEG == 'False':
	#Make the line plot using column1 as your x-axis and column2 as your y-axis
	fig = sns.relplot(x=column1, y=column2, kind=TYPE, data=df)

# If you do ask for a legend
else:
	# If you specified the legend must be outside the figure plot
	if LEG_LOC == 'out':
		# Same as above, except uses the column2 header as the grouping variable in the legend
		#fig = sns.relplot(x=column1, y=column2, kind=TYPE, data=df, legend=LEG)
		fig = sns.relplot(x=column1, y=column2, data=df, hue=column2, legend=LEG)
		# Moves the legend to upper right outside of the figure
		#plt.legend(loc=LEG_LOC)
		plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		#plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

	# If you specified the legend be within the figure plot
	else:
		# Same as above
		fig = sns.relplot(x=column1, y=column2, data=df, hue=column2, legend=LEG)
		# Moves the legend around within the figure plot given your location argument
		plt.legend(loc=LEG_LOC)

# Removes the upper and right border from your figure plot
sns.despine()

## Changing the x-axis label, y-axis label, and figure title from the default data column headers
if XAXIS != 'default':
	XAXIS = XAXIS.replace("_"," ")
	fig.set_ylabels(XAXIS)
if YAXIS != 'default':
	YAXIS = YAXIS.replace("_"," ")
	fig.set_xlabels(YAXIS)
if TITLE != 'default':
	TITLE = TITLE.replace("_"," ")
	plt.title(TITLE)

## Setting the limits of the plot's x and y axes ##
if XLIM != 'default':
	try:
		X1 = float(XLIM.split(", ")[0])
	except ValueError:
		X1 = 'None'
	try:
		X2 = float(XLIM.split(", ")[1])
	except ValueError:
		X2 = 'None'
	if X1 != 'None':
		if X2 != 'None':
			plt.xlim(left=X1, right=X2)
		else:
			plt.xlim(left=X1)
	else:
		plt.xlim(right=X2)

if YLIM != 'default':
	try:
		Y1 = float(YLIM.split(", ")[0])
	except ValueError:
		Y1 = 'None'
	try:
		Y2 = float(YLIM.split(", ")[1])
	except ValueError:
		Y2 = 'None'
	if Y1 != 'None':
		if Y2 != 'None':
			plt.ylim(bottom=Y1, top=Y2)
		else:
			plt.ylim(bottom=Y1)
	else:
		plt.ylim(top=Y2)

## Save to the specified OUTPUT file at the desired dpi level ##
if OUT_FORMAT != 'default':
	FIGURE = PREFIX + ".pdf"
	OUTPUT = os.path.join(OUT_DIR, FIGURE)
	fig.savefig(OUTPUT, dpi=DPI)
else:
	fig.savefig(OUTPUT, dpi=DPI)