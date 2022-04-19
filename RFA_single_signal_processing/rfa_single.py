#!/usr/bin/python3
# rfa_single.py
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

This file contains all signal processing graphics functions.
Sections can be selected for the custom graphics choices in the above article.

Specifications:
Input:
	- mono WAV file, duration at >= 10s
Output:
	PNG graph with panels:
	- AM envelope
	- AM LF spectrum
	- FM envelope (F0 estimation, pitch track)
	- FM LF spectrum
	- AM LF spectrogram
	- FM LF spectrogram
	- AM LF spectrogram max magnitude frequency track
	- FM LF spectrogram max magnitude frequency track
	- AM LF formant dendrogram
	- FM LF formant dendrogram
Methods:
	- AM envelope: full-wave rectification (i.e. np.abs(signal)
	- FM envelope: AMDF (Average Magnitude Difference Function)
	- AM spectrum: FFT
	- FM spectrum: FFT
	- AM spectrogram: moving FFT window
	- FM spectrogram: moving FFT window
	- AM dendrogram: nearest neighbour clustering
	- FM dendrogram: nearest neighbour clustering
Implementation:
	Style: funtional (not OO); readability before pythonicity
	Input: DATA directory
	Output:
		- FIGURES directory with PNG files
	Main application file:
		- rfa.py
	Configuration file:
		- rfa_conf.file
	Modules:
		- module_dendrogram.py
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

# RFA custom module import
from module_F0 import *	# FM demodulation (F0 estimation, 'pitch' tracking)
from module_spectrogram import *	# Low frequency spectrogram functions
from module_dendrogram import *	# Spectral dendrogram drawing functions

#===============================================================
#===============================================================
# Assign configuration parameters to variables (see .conf file)

from rfa_single_conf import *

#===============================================================
#===============================================================
# Filenames; WAV file error handling and input

appfilename = re.sub("^.*/","",sys.argv[0])

if len(sys.argv) > 1:
	wavfilename = sys.argv[1]	# get input filename from command line
else:
	print("Usage: %s <yourwavfile.wav>"%appfilename)
	exit()

wavfilebase = re.sub("^.*/","",wavfilename)
figurefilename = "FIGURES/RFA_%s.png"%wavfilebase
csvfilename = "CSV/RFA%s.csv"%wavfilebase

#===============================================================
#===============================================================
# Mono WAV file input and signal time domain properties

fs, signal = wave.read(wavfilename)							

# read sampling frequency and signal
signallength = len(signal)		# define numerical signal length
signalseconds = signallength / fs	# define signal length in seconds
signal = signal / max(abs(signal))	# normalise signal scale: -1 ... 0 ... 1

#===============================================================
# Butterworth low pass filter (5 Hz is a typical upper limit for the LF spectrum)

b, a = butter(5, 5 / (0.5 * fs), btype="low")	# define Butterworth filter

#===============================================================
#===============================================================
# AM demodulation (envelope extraction) by full-wave rectification

envelope = lfilter(b, a, abs(signal))		# apply filter to create lf envelope
envelope = envelope-min(envelope) / (max(envelope)-min(envelope))	# scale 0 ... 1

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

# Redefinition for column chart display
amrhythmbars = lfamspecfreqs
amweightlist = lfamspecmags

#===============================================================
# Create AM spectrogram and max magnitude value trajectory

ammagarray, amfreqarray, ammaxmags, ammaxfreqs = spectrogramarray(
	signal, fs,amspecfreqmin, amspecfreqmax,
	specdownsample, spectrumpower, specwindowsecs, specstrides)

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
lffmformantmags = [ 0 if m <= fmformantlimit else m for m in lffmspecmags ]	# focus on formants
lffmspecfreqs = fmspecfreqs[1:lffmspecmaglen]	# DC cutoff

# Identification of highest magnitude spectral frequencies
fmtopmagscount = magscount
fmtopmags = sorted(lffmspecmags)[-fmtopmagscount:]
fmtoppos = [ list(lffmspecmags).index(m) for m in fmtopmags ]
fmtopfreqs = [ lffmspecfreqs[p] for p in fmtoppos ]

# Redefinition for column bar display
fmrhythmbars = lffmspecfreqs
fmweightlist = lffmspecmags

#===============================================================
# Create FM spectrogram and max value trajectory

fmmagarray, fmfreqarray, fmmaxmags, fmmaxfreqs = spectrogramarray(
	f0array, framerate, fmspecfreqmin, fmspecfreqmax,
	specdownsample, spectrumpower, specwindowsecs, specstrides)

#===============================================================
#===============================================================
# Graphics definition

# Six rows, two columns
fig,((plt01, plt02), (plt03, plt04), (plt05, plt06), (plt07, plt08), (plt09, plt10), (plt11, plt12)) = plt.subplots(nrows=6, ncols=2, figsize=(figwidth, figheight))

plt.suptitle("%s, fs=%d [RFA M]"%(wavfilename, fs), fontweight="bold")

#===============================================================
#===============================================================
# Plot waveform and envelope

xaxistime = np.linspace(0, signalseconds, signallength)
plt01.plot(xaxistime, signal, color="lightgrey")
plt01.plot(xaxistime, envelope, color="red", label="AM envelope")
plt01.set_xlim(0,np.ceil(signalseconds))
plt01.set_xlabel("Time (s)")
plt01.set_ylabel("Amplitude")
plt01.set_title("1. Waveform and AM envelope demodulation", fontsize=fontsize)
plt01.legend(loc="lower right")
plt01.grid()

#===============================================================
# Frequency domain: AM LF spectrum, high magnitude counts, formant pattern chart

plt02.hist(amrhythmbars, weights=amweightlist, bins=bincount, color="lightgreen", density=True, stacked=True)

plt02.axhline(y=amformantlimit, linestyle=":", label="AM LF formant min "+str(amformantlimit))
plt02.plot(lfamspecfreqs, lfamspecmags, color="green", label="AM LF spectrum")
plt02.scatter(amtopfreqs, amtopmags, s=24, color="red")
for f,m in zip(amtopfreqs, amtopmags):
	plt02.text(f, m-0.1, " %.3fHz"%f, fontsize=6)

plt02.set_xlabel("Freq (Hz)", fontsize=fontsize)
plt02.set_ylabel("Magnitude", fontsize=fontsize)
plt02.set_xlim(amspecfreqmin,amspecfreqmax)
plt02.set_title("2. AM LF spectrum", fontsize=fontsize)
plt02.legend(loc="upper right", fontsize=fontsize-2)
plt02.grid()

#===============================================================
#===============================================================
# FM / F0

xaxistime = np.linspace(0, signalseconds, f0arraylength)	# define x axis in seconds
plt03.scatter(xaxistime, f0array, s=1, color="blue", label="FM, F0")		# plot waveform in grey
plt03.set_ylim(f0min, f0max)
plt03.set_xlabel("Time (s)")
plt03.set_ylabel("Freq (Hz)")
plt03.set_xlim(0,np.ceil(signalseconds))
plt03.set_title("3. FM envelope demodulation (F0 estimation)", fontsize=fontsize)
plt03.legend(loc="upper right")

#===============================================================
# Frequency domain: FM LF spectrum, high magnitude counts, formant pattern chart

plt04.hist(fmrhythmbars, weights=fmweightlist, bins=bincount, color="lightblue", density=True, stacked=True)
plt04.axhline(y=fmformantlimit, linestyle=":", label="FM LF formant min "+str(amformantlimit))
plt04.plot(lffmspecfreqs, lffmspecmags, linewidth=1, color="blue", label="FM LF spectrum")
plt04.scatter(fmtopfreqs, fmtopmags, color="red")
for f,m in zip(fmtopfreqs, fmtopmags):
	plt04.text(f, m-0.1, "%.3fHz"%f, fontsize=8)
plt04.set_xlabel("Freq(Hz)", fontsize=fontsize)
plt04.set_ylabel("Magnitude", fontsize=fontsize)
plt04.set_xlim(fmspecfreqmin,fmspecfreqmax)
plt04.set_title("4. FM LF spectrum", fontsize=fontsize)
plt04.legend(loc="upper right", fontsize=fontsize-2)
plt04.grid()

#===============================================================
#===============================================================
# AM LF spectrogram in heatmap format

plotspectrogramheatmap(plt05, amfreqarray, ammagarray, signalseconds, amspecfreqmin, amspecfreqmax, specgramdotsize, specheatmaptype, fontsize)

plt05.set_title("5. AM rhythm spectrogram (heatmap)", fontsize=fontsize)

# Envelope overlay as an alignment aid
if envelopeoverlay:
	x = np.linspace(0, signalseconds, len(signal))
	env = envelope
	env = amspecfreqmax * env - np.min(env)
	plt05.plot(x, env, color="grey", linewidth=0.5, label="AM envelope")
	plt05.legend(loc="upper right")

#===============================================================
# Choice: AM spectrogram waterfall format and formant pattern column chart

if waterfall:
	arraymax = max([ max(row[1:]) for row in ammagarray ])
	arraylen = len(ammagarray)
	incplus = waterfallincplus*arraymax/arraylen
	inc = 0
	x = np.linspace(amspecfreqmin, amspecfreqmax, len(ammagarray[0]))
	dotsx = []
	dotsy = []
	for row in ammagarray:
		row[0] = 0
		row = row + inc
		plt06.plot(x, row, color='lightgreen', linewidth=0.5, zorder=1)
		ymax = max(row)
		position = list(row).index(ymax)
		dotsx += [x[position]]
		dotsy += [ymax]
		inc += incplus
	plt06.scatter(dotsx, dotsy, s=1, color="red", zorder=1000)
	plt06.set_yticks([])
	plt06.set_xlim(fmspecfreqmin, fmspecfreqmax)
	plt06.set_title("6. AM LF spectrogram (waterfall)", fontsize=fontsize)
	plt06.set_xlabel("Freq (Hz)", fontsize=fontsize)
	plt06.set_ylabel("→ Time →", fontsize=fontsize)
	plt06.grid()

else:
	plt06.set_yticks([])
	plt06.spines["top"].set_visible(False)
	plt06.spines["bottom"].set_visible(False)
	plt06.spines["left"].set_visible(False)
	plt06.spines["right"].set_visible(False)

	plt06.set_title("6. AM reduced resolution rhythm formant pattern (%d bins)"%bincount, 	fontsize=fontsize,)

	plt06.hist(amrhythmbars, weights=amweightlist, bins=bincount, color="grey", density=True, 	stacked=True)

	plt06.set_xlim(amspecfreqmin,amspecfreqmax)
	plt06.set_xlabel("Freq (Hz)")

#===============================================================
#===============================================================
# FM LF spectrogram in heatmap format

plotspectrogramheatmap(plt07, fmfreqarray, fmmagarray, signalseconds, fmspecfreqmin, fmspecfreqmax, specgramdotsize, specheatmaptype, fontsize)

plt07.set_title("7. FM LF spectrogram (heatmap)", fontsize=fontsize)

# Envelope overlay as an alignment aid
if envelopeoverlay:
	x = np.linspace(0, signalseconds, len(signal))
	env = envelope
	env = amspecfreqmax * env - np.min(env)
	plt07.plot(x, env, color="grey", linewidth=0.5, label="AM envelope")
	plt07.legend(loc = "upper right")

#===============================================================
# Choice between FM spectrogram in waterfall format and formant pattern column chart

if waterfall:
	arraymax = max([ max(row[1:]) for row in fmmagarray ])
	arraylen = len(fmmagarray)
	incplus = waterfallincplus*arraymax/arraylen
	inc = 0
	x = np.linspace(fmspecfreqmin, fmspecfreqmax, len(fmmagarray[0]))
	dotsx = []
	dotsy = []
	for row in fmmagarray:
		row[0] = 0
		row = row + inc
		plt08.plot(x, row, color='lightblue', linewidth=0.5, zorder=1)
		ymax = max(row)
		position = list(row).index(ymax)
		dotsx += [x[position]]
		dotsy += [ymax]
		inc += incplus
	plt08.scatter(dotsx, dotsy, s=1, color="red", zorder=1000)
	plt08.set_yticks([])
	plt08.set_xlim(fmspecfreqmin, fmspecfreqmax)
	plt08.set_title("8. FM LF spectrogram (waterfall)", fontsize=fontsize)
	plt08.set_xlabel("Freq (Hz)", fontsize=fontsize)
	plt08.set_ylabel("→ Time →", fontsize=fontsize)
	plt08.grid()

else:
	plt08.set_yticks([])
	plt08.spines["top"].set_visible(False)
	plt08.spines["bottom"].set_visible(False)
	plt08.spines["left"].set_visible(False)
	plt08.spines["right"].set_visible(False)

	plt08.set_title("8. FM reduced resolution rhythm formant pattern (%d bins)"%bincount, fontsize=fontsize,)

	plt08.hist(fmrhythmbars, weights=fmweightlist, bins=bincount, color="grey", density=True, stacked=True)

	plt08.set_xlim(fmspecfreqmin,fmspecfreqmax)
	plt08.set_xlabel("Freq (Hz)")

#===============================================================
#===============================================================
# AM LF spectrogram maximum magnitude frequency vector

x = np.linspace(0,signalseconds, len(ammaxfreqs))
plt09.scatter(x, ammaxfreqs, s=8)
plt09.plot(x, ammaxfreqs)

plt09.set_xlim(0,np.ceil(signalseconds))
plt09.set_ylim(amspecfreqmin, amspecfreqmax)
plt09.set_xlabel("Time (s)")
plt09.set_ylabel("Freq (Hz)")
plt09.set_title("9. AM LF spectrogram max magn frequency trajectory", fontsize=fontsize)
plt09.grid()

# Envelope overlay as an alignment aid
if envelopeoverlay:
	x = np.linspace(0, signalseconds, len(signal))
	env = envelope
	env = amspecfreqmax * env - np.min(env)
	plt09.plot(x, env, color="grey", linewidth=0.5, label="AM envelope")
	plt09.legend(loc="upper right")

#===============================================================
# Dendrogram with AM top formant frequencies

plt10.set_xticks([])
plt10.set_yticks([])
plt10.spines["top"].set_visible(False)
plt10.spines["bottom"].set_visible(False)
plt10.spines["left"].set_visible(False)
plt10.spines["right"].set_visible(False)

drawdendrogram(plt10, zip(amtopfreqs, amtopmags), amboxwidth, amboxheight, amboxx, amboxy, fontsize)

plt10.set_title("10. AM LF formant pattern (%d values)"%len(amtopmags), fontsize=fontsize)

#===============================================================
#===============================================================
# FM LF spectrogram maximum magnitude frequency vector

x = np.linspace(0,signalseconds, len(fmmaxfreqs))
plt11.scatter(x, fmmaxfreqs, s=8)
plt11.plot(x, fmmaxfreqs)

plt11.set_xlim(0,np.ceil(signalseconds))
plt11.set_ylim(fmspecfreqmin, fmspecfreqmax)
plt11.set_xlabel("Time (s)")
plt11.set_ylabel("Freq (Hz)")
plt11.set_title("11. FM LF spectrogram max magn frequency trajectory", fontsize=fontsize)
plt11.grid()

# Envelope overlay as an alignment aid
if envelopeoverlay:
	x = np.linspace(0, signalseconds, len(signal))
	env = envelope
	env = amspecfreqmax * env - np.min(env)
	plt11.plot(x, env, color="grey", linewidth=0.5, label="AM envelope")
	plt11.legend(loc="upper right")

#===============================================================
# Dendrogram with FM top formant frequencies

plt12.set_xticks([])
plt12.set_yticks([])
plt12.spines["top"].set_visible(False)
plt12.spines["bottom"].set_visible(False)
plt12.spines["left"].set_visible(False)
plt12.spines["right"].set_visible(False)

drawdendrogram(plt12, zip(fmtopfreqs, fmtopmags), fmboxwidth, fmboxheight, fmboxx, fmboxy, fontsize)
plt12.set_title("12. FM LF formant pattern (%d values)"%len(amtopmags), fontsize=fontsize)

#===============================================================
# Graph output

plt.tight_layout(pad=3, w_pad=1, h_pad=1)
plt.savefig(figurefilename)
if showgraph:
	plt.show()

# EOF
