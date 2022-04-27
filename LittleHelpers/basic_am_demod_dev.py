#!/usr/bin/python3
# basic_am_demod.py
# D. Gibbon, 2021-03-16

figwidth = 14
figheight = 6
maxfreq = 6
peaks = True
peaklevel = 0.9

#================================================

import sys, re
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wave
from scipy.signal import medfilt

#================================================
# Input command line arguments

wavfilename = sys.argv[1]
pngfilestem = re.sub(".wav","",wavfilename)
appfilestem = re.sub(".py","",sys.argv[0])

# Input audio signal

fs, signal = wave.read(wavfilename)
period = 1/fs
siglen = len(signal)
sigsecs = siglen / fs

#================================================
# Extract AM envelope, apply FFT

env = np.abs(signal)	# Minimalist method :)
mags = np.abs(np.fft.rfft(env))
freqs = np.fft.rfftfreq(env.size, period)

#================================================
# Normalise signal and envelope

signal = signal / max(signal)
env = medfilt(env / max(env), 301) # Only odd values!
env = env/max(env)

# Select LF spectrum
xoffset = 4		# Arbitrary offset to exclude 0
foffset = xoffset / fs
maxsample = int(round(maxfreq*sigsecs))
mags = mags[xoffset:maxsample]
mags = mags / max(mags)
freqs = freqs[xoffset:maxsample]

#================================================
# Generate graphics

fig = plt.figure(figsize=(figwidth, figheight))
x = np.linspace(0, sigsecs, siglen)

plt.suptitle("Duration: %.3f, Sample rate: %d, file: %s"%(sigsecs, fs, wavfilename))

plt.subplot(3,1,1)	# waveform
plt.plot(x, signal, color="b")
plt.plot(x, env, color="orange")
plt.xlim(0,sigsecs)

plt.subplot(3,1,2)	# env
plt.plot(x, env, color="orange")
plt.xlim(0,sigsecs)
plt.title("Envelope")

#================================================
# Spectrum

plt.subplot(3,1,3)	# spectrum
plt.xlim(0, maxfreq)
plt.plot(freqs, mags, color="g")
plt.title("LF spectrum")
for freq, peak in zip(freqs, mags):
	if peak >= peaklevel:
		plt.scatter(freq, peak, color="r")
		plt.text(freq, peak, "%.3f"%freq)
plt.axhline(peaklevel, linestyle="--", color="pink")
print
plt.text(foffset, peaklevel, "PL=%s"%peaklevel)

# Save and show graphics
plt.tight_layout(pad=1, w_pad=1, h_pad=1)
plt.savefig("%s_%s.png"%(appfilestem, pngfilestem))
plt.show()
