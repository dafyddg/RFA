# module_fm_demodulation.py
# D. Gibbon,
# Created 2021-07-11
# Modified 2022-04-18

"""
FM demodulation (F0 estimation, pitch tracking)

Preprocessing:
	centre and peak clipping
	Butterworth low pass filter

FM demodulation:
	The algorithm is the AMDF (Absolute Magnitude Difference Function).
	Often known as the Average Magnitude Difference Function.
	In this interpretation of AMDF, the sum, not the average is used.
	Four Pythonic and Non-Pythonic implementations are given.
	Algorithm:
	Moving window from frame to frame through the signal
		- frame duration dependent on selected minimum F0
		- skip interval to next frame is selectable
		- moving window through frame and neighbour
			- offset avoids initial zero, depends on selected max F0
			- sum of differences between window and frame
		- min distance: argmin, add offset, convert to period & frequency
- SMDF, Summed Magnitude Difference Function (variant of AMDF, Average ...)
- This implementation simply uses the sum of differences: simpler and faster
- Instead of "abs" "square" could be used, yielding the ASDF
- The difference is minimal, on a quick test ASDF seems a little cleaner.
- The choice is up to the user.

F0 estimation parameter guesstimation (adjust if you know what you are doing):
- Signal pre-forming,  F0 post-forming settings
- The "female" and "male" options are just general settings.
- The f0min and f0max values can be set to whatever is convenient.
- f0min co-determines the AMDF offset
- f0max co-determines the length of the frame
Four algorithm styles are provided, user selectable: A, B, C, D.
"""

#===============================================================

import numpy as np
from scipy.signal import butter, lfilter, medfilt, tukey

#===============================================================
# The voice variable is just a mnemonic convenience.
voice = "male"

if voice == "male":
	# Male default settings
	f0min =100
	f0max = 300
elif voice == "female":
	# Female default settings
	f0min = 110
	f0max = 350
else:
	# General default settings
	f0min = 70
	f0max = 420

f0medfilter = 3

f0framelengthfactor = 0.75							# relative to f0min, > 1
f0frameskipfactor = 0.5								# resolution, default is 1, the frame length

# Frame, offset and window shape definition
f0frameduration =  1 / f0min
f0frameduration = f0framelengthfactor * f0frameduration
framerate = 2 / f0frameduration						# 2 is because of the fs/2 spectrum resolution

f0diffoffsetlengthfactor = 0.1
f0diffoffsetduration = f0diffoffsetlengthfactor / f0max	# duration in seconds, relative to f0max

#===============================================================
# Centre and peak clipping (adjust if you know what you are doing)

centrethresh = 0.1			# Deals with silence and low volume noise
limitthresh = 0.9

#===============================================================
# Butterworth bandpass filter (adjust if you know what you are dong)

fmbutterhigh = f0min * 2
fmbutterhighorder = 5
fmbutterlow = f0max
fmbutterloworder = 2

#===============================================================
# Butterworth filter

def butterworthfilter(signaldata, cutoff, order, fs, type):

	nyqvist = 0.5 * fs
	normal_cutoff = cutoff / nyqvist
	b, a = butter(order, normal_cutoff, btype=type, analog=False)
	filteredsignal = lfilter(b, a, signaldata)

	return filteredsignal

#===============================================================
# Zero clipping and peak clipping

def clipper(sig,thresh,type):
	if type == "centre":
		clipped = (abs(sig) > thresh).astype(int) * sig
	elif type == "limit":
		clipped = (abs(sig) < thresh).astype(int) * sig
	elif type == "lower":
		clipped = (sig > thresh).astype(int) * sig
	elif type == "upper":
		clipped = (sig < thresh).astype(int) * sig
	else:
		print("Unknown type:",type); sys.exit()
	return np.asarray(clipped)

#===============================================================
#===============================================================

def f0amdf(signal, fs, framestart, framelength, f0diffoffsetlength, algo):

	"""
AMDF FM demodulation (F0 estimation, 'pitch' extraction) for one frame

If multiplication is used instead of subtraction, and argmax instead of argmin, the 	implementations are effectively cross-correlation variants.
"""

	# Define the end of the current frame
	# Define the range of the lag window
	# from the frame start plus offset to frame end
	framestop = framestart + framelength
	movingwindowrange = np.arange(framestart+f0diffoffsetlength, framestop)

	# Make a list of differences between the frame and lag frame
	# The comparison starts after an offset
	# to avoid the identity of frame and zero lag frame
	# The sums of the differences are collected
	# in a list using a Python comprehension loop

	#==================================================
	# "There's more than one way to do it!"
	#==================================================

	# Pythonic comprehension, calculations in loop
	if algo == "A":
		diffsums = [
			np.sum(np.abs(signal[framestart:framestop] - signal[winstart:winstop]))
			for winstart, winstop in zip(movingwindowrange, movingwindowrange+framelength) ]
	
#===============================================================
	# Pythonic comprehension, calculations outside loop
# Version B is on average slightly faster than the others.
	elif algo == "B":
		frame = signal[framestart:framestop]
		movingwindows =  zip(movingwindowrange, movingwindowrange+framelength)
		diffsums = [
			np.sum(np.abs(frame - signal[winstart:winstop]))
			for winstart, winstop in movingwindows ]

#===============================================================
	# Classic for-loop with output list construction
	elif algo == "C":
		frame = signal[framestart:framestop]
		movingwindows =  zip(movingwindowrange, movingwindowrange+framelength)
		diffsums = []
		for winstart, winstop in movingwindows:
			diffsums += [
				np.sum(np.abs(frame - signal[winstart:winstop])) ]

#===============================================================
	# Classic for-loop with pre-defined empty output list
	elif algo == "D":
		frame = signal[framestart:framestop]
		movingwindows =  zip(movingwindowrange, movingwindowrange+framelength)
		diffsums = np.zeros(len(movingwindowrange))
		for i, (winstart, winstop) in enumerate(movingwindows):
			diffsums[i] = np.sum(np.abs(frame - signal[winstart:winstop]))

#===============================================================
	"""
- The position of the first smallest difference in the sums of differences is calculated and added to the offset position in order to determine period estimate (in samples).
- The period estimate is divided by the sampling rate to obtain the period estimate in seconds.
- The inverse of the period estimate is calculated in order to find the F0 estimate.
	"""
	f0 = 1 / ( (np.argmin(diffsums) + f0diffoffsetlength) / fs )

	return f0

#===============================================================
#===============================================================
# Move through the signal from frame to frame, calling the AMDF function

def f0estimate(signal,fs):

	"""
AMDF:
- frame duration is defined relative to specified f0min (the frame has to be long enough 	to capture low frequencies with long periods), and a frame length factor for fine tuning
- start offset duration for difference calculation is defined relative to specified f0max; it has to be just far enough from the frame start so as not to be further than the short period of the maximum frequency), and a length factor for fine tuning
	"""

	framelength = int(f0frameduration * fs)
	frameskip = int(framelength * f0frameskipfactor)

	f0diffoffsetlength = int(f0diffoffsetduration * fs)	# samples

	# F0 preprocessing: clip the low amplitude noise between speech units
	signal = clipper(signal,centrethresh,"centre")
	signal = clipper(signal,limitthresh,"limit")
	signal = butterworthfilter(
		signal, fmbutterlow, fmbutterloworder, fs, "low")
	signal = butterworthfilter(
		signal, fmbutterhigh, fmbutterhighorder, fs, "high")
#	windowshape = tukey(framelength, f0tukeyfraction)	# Not used here

	algo = "B"	# Equivalent AMDF implementations A, B, C, D. On average, B is slightly faster.

	#===============================================================
	# Make an array from list of f0 results for all frames.
	# The list is created using a Python comprehension loop

	f0track = np.array([
		f0amdf(signal, fs, framestart, framelength, f0diffoffsetlength, algo)
		for framestart in range(0, len(signal)-3*framelength, frameskip)
		])

	#===============================================================
	# F0 median smoothing and min max cutoff

	f0track = medfilt(f0track, f0medfilter)
	f0track = [ 0 if (f0 < f0min) or (f0 > f0max) else f0 for f0 in f0track ]

	return f0track, framerate, f0frameduration

#----------------------------------------------------------------------

