#!/usr/bin/python3
# rfa_mult.py
# Dafydd Gibbon
# Created 2021-08-16
# Modified 2022-04-18
# Code first deposited on GitHub 2021-08-16
# http://www.github.com/dafyddg/RFA

# Note: the code is designed for reading
# and is not optimised or particularly pythonic.
# There is no error trapping or input checking.
# See top level README.1st file.

"""
MIT license
Begin license text.
Copyright 2021 Dafydd Gibbon
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
End license text.
"""

"""
Functional specifications and application in this article:

Gibbon, Dafydd (2021). Journal of the International Phonetic Association , First View , pp. 1-33.
DOI: https://doi.org/10.1017/S0025100321000086

This file contains all signal processing and clustering functions.
Sections can be selected for the custom graphics choices in the above article.

Specifications:
Input:
	- mono WAV file, duration at >= 10s
Output:
	CSV file with values:
	- AM LF spectrum frequencies
	- AM LF spectrum magnitudes
	- FM LF spectrum frequencies
	- FM LF spectrum magnitudes
	- AM LF spectrum max magnitude frequencies
	- FM LF spectrum max magnitude frequencies
	- AM LF spectrogram frequencies
	- AM LF spectrogram magnitudes
	- FM LF spectrogram frequencies
	- FM LF spectrogram magnitudes
Methods:
	- AM envelope: full-wave rectification (i.e. np.abs(signal)
	- FM envelope: AMDF (Average Magnitude Difference Function)
	- AM spectrum: FFT
	- FM spectrum: FFT
	- AM spectrogram: moving FFT window
	- FM spectrogram: moving FFT window
Implementation:
	Style: funtional (not OO); readability before pythonicity
	Input: DATA directory
	Output:
		- CSV directory
		- PNG directory
	Main application file:
		- rfa_mult.py
	Configuration file:
		- rfa_mult_conf.file
	Modules:
		- module_F0.py
		- module_spectrogram.py
"""

#===============================================================
#===============================================================
# System and library module import

import sys, re
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wave
from scipy.signal import butter, lfilter
from datetime import datetime
from glob import glob

# RFA custom module import
from module_F0 import *	# FM demodulation (F0 estimation, 'pitch' tracking)
from module_spectrogram import *	# Low frequency spectrogram functions

#===============================================================
#===============================================================
# Date for the CSV output

datetoday = datetime.today().strftime('%Y-%m-%d-%H:%M:%S')

#---------------------------------------------------------------
# Assign configuration parameters to variables (see .conf file)

from rfa_mult_conf import *

#---------------------------------------------------------------
# CSV output

def outputtextlines(text, filename):
	handle = open(filename,'w')
	linelist = handle.write(text)
	handle.close()
	return

def appendtextlines(text, filename):
	handle = open(filename,'a')
	linelist = handle.write(text)
	handle.close()
	return

def csvoutput(i, wavfilebase, vectors, csvfile):
	outputstring = wavfilebase + separator + separator.join(
			[str(mag) for mag in vectors ]) + "\n"
	if i==0: outputtextlines(outputstring, csvfile)
	else: appendtextlines(outputstring, csvfile)
	return

#===============================================================
#===============================================================
# Filenames; WAV file error handling and input

appfilename = re.sub("^.*/","",sys.argv[0])

if len(sys.argv) > 1:
	wavfiledirectory = sys.argv[1]
else:
	wavfiledirectory = "DATA/DATA_RT/"	# demo mode default data
	datasetname = "RFA_N_RT"							# output files prefix

wavfilelist = sorted(glob(wavfiledirectory+"*.wav"))

#===============================================================
#===============================================================
# Cycle through the filenames in the given directory, building lists

wavnamelist = []		# to label nodes in distance map

f0list = []

lfamsmlist = []
lffmsmlist = []

lfammaxflist = []
lfammaxmlist = []
lffmmaxflist = []
lffmmaxmlist = []

lfamsmfile = "CSV/lfamspecmags.csv"
lffmsmfile = "CSV/lffmspecmags.csv"

lfammaxffile = "CSV/lfammaxfreqs.csv"
lfammaxmfile = "CSV/lfammaxmags.csv"

lffmmaxffile = "CSV/lffmmaxfreqs.csv"
lffmmaxmfile = "CSV/lffmmaxmags.csv"

lfamtrajmfile = "CSV/lfamtrajmags.csv"
lfamtrajffile = "CSV/lfamtrajfreqs.csv"
lffmtrajmfile = "CSV/lffmtrajmags.csv"
lffmtrajffile = "CSV/lffmtrajfreqs.csv"

for i, wavfilename in enumerate(wavfilelist):	# Collect spectra for all files
	wavfilebase = re.sub(".*/", "", wavfilename)
	wavefilebase = re.sub(".wav", "", wavfilebase)
	wavnamelist += [ wavfilebase ]
	print(wavfilebase)

	#=======================================================
	# Mono WAV file input and signal time domain properties

	fs, signal = wave.read(wavfilename)	# read sampling frequency and signal
	signallength = len(signal)		# define numerical signal length
	signalseconds = signallength / fs	# define signal length in seconds
	signal = signal / max(abs(signal))	# scale signal: -1 ... 0 ... 1

	#=======================================================
	# Butterworth LP filter (5 Hz: typical upper limit for LF spectrum)

	b, a = butter(5, 5 / (0.5 * fs), btype="low")	# define Butterworth filter

	#===============================================================
	#===============================================================
	# AM demodulation (envelope extraction) by full-wave rectification

	envelope = lfilter(b, a, abs(signal))	# apply filter to create lf envelope
	envelope = envelope-min(envelope) / (max(envelope)-min(envelope))	# scale

	#===============================================================
	# AM low frequency spectral analysis

	# FFT of complete envelope, output magnitude values
	amspecmags = np.abs(np.fft.rfft(envelope))
	amspecmaglen = len(amspecmags)

	# Extraction of low frequency spectrum segment
	lfamspecmaglen = int(round(amspecfreqmax * amspecmaglen / (fs / 2)))
	lfamspecmags = amspecmags[1:lfamspecmaglen]	# DC cutoff
	lfamspmMin = min(lfamspecmags)
	# Scale to 0...1
	lfamspecmags = (lfamspecmags-lfamspmMin) / (np.max(lfamspecmags)-lfamspmMin)

	# Assign LF spectrum frequencies to magnitude values
	lfamspecfreqs = np.linspace(0,fs/2,amspecmaglen)
	lfamspecfreqs = lfamspecfreqs[1:lfamspecmaglen]

	# Identification of highest magnitude spectral frequencies
	amtopmagscount = magscount
	amtopmags = sorted(lfamspecmags)[-amtopmagscount:]
	amtoppos = [ list(lfamspecmags).index(m) for m in amtopmags ]
	amtopfreqs = [ lfamspecfreqs[p] for p in amtoppos ]

	#===============================================================
	# Create AM spectrogram and max magnitude value trajectory

	ammagarray, amfreqarray, amtrajmags, amtrajfreqs = spectrogramarray(
		signal, fs,amspecfreqmin, amspecfreqmax,
		specdownsample, spectrumpower, specwindowsecs, specstrides
		)

	#===============================================================
	#===============================================================
	# FM demodulation (F0 estimation, pitch extraction)
	# by AMDF (Absolute Magnitude Difference Function)

	f0array, framerate, frameduration = f0estimate(signal, fs)
	f0arraylength = len(f0array)

	#===============================================================
	# FFT low frequency spectral analysis of F0 estimation track

	fmspecmags = np.abs(np.fft.rfft(f0array))
	fmspecmaglen = len(fmspecmags)
	fmspecfreqs = np.linspace(0,framerate/2,fmspecmaglen)

	# Extraction of low frequency spectral segment, with magnitude filter
	lffmspecmaglen = int(round(fmspecfreqmax * fmspecmaglen / (framerate / 2)))
	lffmspecmags = fmspecmags[1:lffmspecmaglen]
	lffmspecmags = lffmspecmags / np.max(lffmspecmags)

	# Magnitude filter for LF formant analysis
	lffmformantmags = [ 0 if m <= fmformantlimit else m for m in lffmspecmags ]
	lffmspecfreqs = fmspecfreqs[1:lffmspecmaglen]	# DC cutoff

	# Identification of highest magnitude spectral frequencies
	fmtopmagscount = magscount
	fmtopmags = sorted(lffmspecmags)[-fmtopmagscount:]
	fmtoppos = [ list(lffmspecmags).index(m) for m in fmtopmags ]
	fmtopfreqs = [ lffmspecfreqs[p] for p in fmtoppos ]

	#===============================================================
	# Create FM spectrogram and max value trajectory

	fmmagarray, fmfreqarray, fmtrajmags, fmtrajfreqs = spectrogramarray(
		f0array, framerate, fmspecfreqmin, fmspecfreqmax,
		specdownsample, spectrumpower, specwindowsecs, specstrides
		)

	#===============================================================
	# AM spectrum CSV outputs
	csvoutput(i, wavfilebase, lfamspecmags, lfamsmfile)
	csvoutput(i, wavfilebase, amtopmags, lfammaxmfile)
	csvoutput(i, wavfilebase, amtopfreqs, lfammaxffile)
	# AM spectrogram max trajectory CSV outputs
	csvoutput(i, wavfilebase, amtrajmags, lfamtrajmfile)
	csvoutput(i, wavfilebase, amtrajfreqs, lfamtrajffile)
	# FM spectrum CSV outputs
	csvoutput(i, wavfilebase, lffmspecmags, lffmsmfile)
	csvoutput(i, wavfilebase, fmtopmags, lffmmaxmfile)
	csvoutput(i, wavfilebase, fmtopfreqs, lffmmaxffile)
	# AM spectrogram max trajectory CSV outputs
	csvoutput(i, wavfilebase, fmtrajmags, lffmtrajmfile)
	csvoutput(i, wavfilebase, fmtrajfreqs, lffmtrajffile)

#===============================================================
# EOF
