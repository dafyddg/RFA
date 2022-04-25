#!/usr/bin/python3
# basic_am_demod.py
# D. Gibbon, 2021-03-16

import sys, re
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wave
from scipy.signal import medfilt

wavfilename = sys.argv[1]
pngfilestem = re.sub(".wav","",wavfilename)
appfilestem = re.sub(".py","",sys.argv[0])
samplerate, signal = wave.read(wavfilename)
period = 1/samplerate
siglen = len(signal)
sigsecs = siglen / samplerate

env = np.abs(signal)
mags = np.abs(np.fft.rfft(env))
freqs = np.fft.rfftfreq(env.size, period)

signal = signal / max(signal)
env = medfilt(env / max(env), 301)
mags = mags / max(mags)

x = np.linspace(0, sigsecs, siglen)
plt.subplot(3,1,1)	# waveform
plt.plot(x, signal)
plt.subplot(3,1,2)	# env
plt.plot(x, env)

maxfreq = 6
maxsample = int(round(maxfreq*sigsecs))

plt.subplot(3,1,3)	# env
mags = mags[int(round(sigsecs/2)):maxsample]
freqs = freqs[int(round(sigsecs/2)):maxsample]
plt.plot(freqs, mags)

plt.savefig("%s_%s.png"%(appfilestem, pngfilestem))
plt.show()
