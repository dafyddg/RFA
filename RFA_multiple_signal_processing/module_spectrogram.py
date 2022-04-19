# module_spectrogram.py
# D. Gibbon
# Created 2021-08-16
# Modified 2022-04-18
# Spectrogram array module for rfa.py

#===============================================================

import numpy as np

#===============================================================

def spectrogramarray(signal, fs, specfreqmin, specfreqmax, specdownsample, spectrumpower, specwindowsecs, specstrides):

	# Brute force downsampling, optional
	signal = signal[::specdownsample]
	fs = int(round(fs/specdownsample))
	period = 1/fs
	signallen = len(signal)
	signalsecs = int(round(signallen / fs))

	#============================================
	windowlen = int(round(specwindowsecs * fs))	# window length sec -> sample
	signalleneffective = signallen - windowlen	# first to last stride pos
	stride = int(round(signalleneffective / specstrides))	# time step

	# Moving window
	
	# Stride start and end counters
	counterstart = np.array(range(0,signalleneffective,stride))	# start & end
	counterend = counterstart + windowlen

	magarray = []
	freqarray = []
	for countstart,countend in zip(counterstart,counterend):
		segment = np.abs(signal[countstart:countend])				# window-length segment
		segment = list(segment)
		segment = segment
		segment = np.array(segment)
		mags = abs(np.fft.rfft(segment))						# FFT magnitudes
		freqs = np.abs(np.fft.rfftfreq(segment.size,period))	# FFT frequencies
#		freqs = np.linspace(0, fs/2, len(mags))
		magarray += [ mags ]								# collect FFTs
		freqarray += [ freqs ]

	#============================================
	# Spectrum properties

	spectrummax = int(round(fs/2))
	rowlen = len(freqarray[0])
	elementsperhertz = int(round( rowlen / spectrummax ))

#	print("Elements per hertz:", elementsperhertz)
	xmin = specfreqmin * elementsperhertz
	xmax = specfreqmax * elementsperhertz
	sfmin = np.int(np.floor(xmin))
	sfmax = np.int(np.ceil(xmax))
	magarray = np.array([ x[sfmin:sfmax] for x in magarray ])**spectrumpower
	freqarray = [ x[sfmin:sfmax] for x in freqarray ]

	#Detect maximum vector through spectrogram
	maxmags = np.array([ max(mags[1:]) for mags in magarray ])
#	print(maxmags.shape)

	#============================================
	# Loop to define spectrogram as a spectrum sequence
	maxfreqs = []
	for mags, freqs in zip(magarray, freqarray):
		maxofmags = np.max(mags)
		# An error with a very deep voice (60Hz) threw an error
		if maxofmags == 0.0: maxofmags = 0.0001	# a hack, sorry
		mags = mags/maxofmags
		mags = list(mags[1:])
		maxmag = np.max(mags)
		maxmagpos = mags.index(maxmag)
		maxfreq = freqs[maxmagpos]
		maxmags += [maxmag]
		maxfreqs += [maxfreq]

	return np.array(magarray), np.array(freqarray), maxmags, maxfreqs

#===============================================================

# Rotation of spectrogram array as heatmap

def plotspectrogramheatmap(pltobj, freqarray, magarray, signalsecs, specfreqmin, specfreqmax, specgramdotsize, specheatmaptype, fontsize):

	# y-axis as scale of the range, number of spectra in spectrogram
	y = np.linspace(specfreqmin,specfreqmax,len(magarray[0]))

	# x-axis as signal time range, number of spectra in spectrogram
	x = np.linspace(0,signalsecs,len(magarray))

	# Colormap is derived from magnitudes at each frequency/spectrum
	for i, freqvals, magvals in zip(x, freqarray, magarray):
		freqvals = freqvals[1:]
		magvals = magvals[1:]
		x = [i] * len(freqvals)
		pltobj.scatter(
			x,freqvals, c=magvals, cmap=specheatmaptype,
			marker="s", s=specgramdotsize)

	# Spectrogram properties
	pltobj.set_xlim(0,np.ceil(signalsecs))
	sfmin = np.floor(specfreqmin)
	sfmax = np.ceil(specfreqmax)
	pltobj.set_ylim(sfmin,sfmax)
	pltobj.grid(b=True, which="major", axis="both")
	pltobj.set_xlabel("Time (s)", fontsize=fontsize)
	pltobj.set_ylabel("Freq (Hz)", fontsize=fontsize)

	return

# EOF
