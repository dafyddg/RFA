#!/usr/bin/python3

# textgrid2csv.py, D. Gibbon, 2015.02.12
# Convert Praat TextGrid to CSV
import os, re, sys
TAB = '\t'; NL = '\n'; SPACE = " "

# Check input from CLI
if len(sys.argv) < 2:
    print('Usage:',sys.argv[0],'<filename.TextGrid>')
    exit()
textgridfile = sys.argv[1]
if not os.path.isfile(textgridfile):
    print('File', textgridfile, 'does not exist.')
    exit()
csvfile = re.sub('.TextGrid','.csv',textgridfile)

# Read textgrid, preprocess textgrid format
textgrid = open(textgridfile, 'r' ).read().split(NL)
nugrid = []
for line in textgrid:
    line = re.sub(' *$','',line) # Del. fin. spaces
    line = re.sub('^ *','',line) # Del. init. spaces
    line = re.sub('\"','',line)  # Delete quotes
    if line != '':
        nugrid += [line]

# Initialise CSV text string (metadata, column heads)
csvstring = textgridfile + NL
csvstring += TAB.join(
    ['Tier', 'Label', 'xmin', 'xmax', 'xdiff'] ) + NL

# Initialise row identifier flag
rowflag = False

# Loop through lines in preprocessed textgrid
# The Praat indented hierarchy has finite depth
# and can therefore be processed by an iterative
# finite machine - recursion unnecessary:
for line in nugrid:
    if 'name = ' in line:            # Get tier name
        name = line.split(SPACE)[-1]
    if 'intervals [' in line:        # Skip the header
    	rowflag = True
    if rowflag and 'xmin' in line:
        xmin = '%.3f'%float(line.split(SPACE)[-1])
    if rowflag and 'xmax' in line:
        xmax = '%.3f'%float(line.split(SPACE)[-1])
    if rowflag and "text" in line:
        text = line.split(SPACE)[-1]
        xdiff = '%.3f'%(float(xmax)-float(xmin))
        row = TAB.join([name, text, xmin, xmax, xdiff])
        csvstring += row + NL
        rowflag = False
        
# Output CSV file
open(csvfile,'w').write(csvstring)
