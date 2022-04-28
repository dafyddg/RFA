#!/usr/bin/python3
# -*- coding: utf-8 -*-'

# V01: PyAMDF.py
# D. Gibbon
# 2018-10-10

# V02: basic_fm_demod.py
# 2022-04-28

# F0 Estimator based on Average Magnitude Difference Function (AMDF)
# Method:
#	move a window (frame) through the signal
#	at each step,
#		copy the frame
#			move the copy through the frame
#			at each step subtract copy from frame
#			average these differences for the frame
#		find the first smallest difference
#		calculate period from start of frame
#		convert to frequency (F0 estimation)
#		collect the frequencies (F0 track estimation)

import sys, re
import numpy as np
from scipy.signal import medfilt, hilbert
import scipy.io.wavfile as wav
import matplotlib.pyplot as plt

#==============================================================================
#==============================================================================
# Command line parameters

try:
	appname = sys.argv[0]
	filename = sys.argv[1]
	voice = sys.argv[2]
	filebase = re.sub('.wav','',filename)
	appbase = re.sub('.py', '', appname)
	figurename = "%s_%s.png"%(appbase, filebase)

except:
	print('Usage: %s wavfilename high|fairlyhigh|mid|fairlylow|low'%appname); exit()

#==============================================================================
# Voice register heuristic settings

try:
	if voice == 'high':
		f0min = 160
		f0max = 450
		frameduration = 0.005	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.005  # AMDF offset, higher: lower frequencies
	
	elif voice == 'fairlyhigh':
		f0min = 140
		f0max = 400
		frameduration = 0.005	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.02  # AMDF offset, higher: lower frequencies

	elif voice == 'mid':
		f0min = 120
		f0max = 300
		frameduration = 0.007	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.02  # AMDF offset, higher: higher frequencies

	elif voice == 'fairlylow':
		f0min = 120
		f0max = 350
		frameduration = 0.02	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.02  # AMDF offset, higher: higher frequencies

	elif voice == 'low':
		f0min = 90
		f0max = 200
		frameduration = 0.01	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.03  # AMDF offset, higher: lower frequencies

#==============================================================================

	elif voice == 'putonghua':
		f0min = 120
		f0max = 300
		frameduration = 0.007	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.01  # AMDF offset, higher: lower frequencies

	elif voice == 'mlk':
		f0min = 120
		f0max = 350
		frameduration = 0.02	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.02  # AMDF offset, higher: higher frequencies

	elif voice == 'onetoseven':
		f0min = 90
		f0max = 200
		frameduration = 0.02	# Higher: lower frequencies
		framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
		frameoffsetfraction = 0.03  # AMDF offset, higher: lower frequencies

	else:
		print('Unknown voice type.'); exit()

#==============================================================================
# Figure settings
	figwidth = 14
	figheight = 6
	envelopeflag = False

# Waveform settings (0, ..., 1)
	centreclip = 0	# Default 0
	peakclip = 1	# Default 1

# Envelope low pass filter (moving median, always)
	peakwin = 20
	envmed = 201

#==============================================================================
# AMDF settings

	"""
	framestepfactor = 0.2		# AMDF frame step divisor (default: 1)
	
# Typical: low voice: 0.02, high voice 0.01
#	frameduration = 0.007	# Higher: lower frequencies

# Typical: low voice: 0.1, high voice 0.01
#	frameoffsetfraction = 0.02  # AMDF offset, higher: higher frequencies
	"""
	
# F0 post-processing with moving median window
	f0medwin = 5	# Must be an odd integer. Default 5

except:
	print('Parameter error.'); exit()

#==============================================================================
#==============================================================================
# WAV file input

try:
	fs, signal = wav.read(filename)

	signallen = len(signal)
	signalduration = int(round(1.0*signallen/fs))
	abssignal = abs(signal)

except:
	print("Error reading signal."); exit()

#==============================================================================
#==============================================================================
# Preprocessing (centre-clipping)

try:
	signal = signal / float(np.max(abs(signal)))  # Normalisation

	signal = np.array([	s if abs(s)>centreclip else 0 for s in signal ])
	signal = np.array([	s if abs(s)<peakclip else 1 for s in signal ])

	abssignal = abs(signal)

	envelope = [ max(abssignal[i-peakwin:i]) for i in range(peakwin,len(abssignal)) ]
	envelope = np.append(medfilt(envelope,envmed),[0]*peakwin)
	envelope = envelope / max(envelope)
	
except:
	print("Error preprocessing."); exit()

#==============================================================================
# F0 estimation

if True:

	sampframeratio = fs * frameduration

	framelen = int(round(fs * frameduration))

	# Offset to ensure no AMDF zero:zero match
	frameoffsetduration = frameduration * frameoffsetfraction # 3ms, 333 Hz
	frameoffsetposition = int(round(fs * frameoffsetduration))

	framestep = int(np.round(framelen * framestepfactor))
	if framestep < 1 or framestep > (signallen/10):	# validity check
		framestep = framelen
	irange = list(range(0, signallen-2*framelen, framestep))	# Note truncation.

	framecount = int(round(signallen / framestep))

#	newsignallen = framecount * framelen
#	signal = signal[:newsignallen]

# Allocate memory for f0list and AMDF list
	f0list = np.zeros(framecount)
	meandiffs = np.zeros(framelen).tolist()

# Move frame window through signal
	for count, framestart in enumerate(irange):
		framestop = framestart + framelen
		frame = signal[framestart:framestop]

# Calculate Average Magnitude Difference Function with moving window
		for lag, winstart in enumerate(range(framestart,framestop)):
			movingwin = signal[winstart:winstart+framelen]
			meandiffs[lag] = np.mean(abs(frame - movingwin))

# Pick smallest absolute difference in frame
		smallestdiff = np.min(meandiffs[frameoffsetposition:])

# Get position of the smallest absolute difference
		smallestdiffposition = meandiffs.index(smallestdiff)

# Divide the sampling rate by the number of samples in the interval
# That is: t = index/fs; f0 = 1/t
		f0 = fs / smallestdiffposition if smallestdiffposition > 0 else 0

# Extend f0 list
		f0list[count] = f0

if False:
	print("F0 estimation error."); exit()

#==============================================================================
# Post-processing

try:
# Remove f0 values outside defined limits
	f0list = (f0list > f0min).astype(int) * f0list
	f0list = (f0list < f0max).astype(int) * f0list

# Smooth F0 contour
	f0list = medfilt(f0list,f0medwin)

except:
	print("Post-processing failed."); exit()

#==============================================================================
#==============================================================================
# Figures

try:
	fig, (sp1,sp2,sp3)  = plt.subplots(3, 1, figsize=(figwidth,figheight))
	fig.suptitle('AM and FM demodulation. DG, 2022-04-28',fontsize=16)

except:
	print("Figure definition error."); exit()

#==============================================================================
# Plot signal

try:
	x = np.linspace(0,len(signal),len(signal))/fs
	signal = medfilt(signal,31)
	signal = signal / float(np.max(abs(signal)))
	sp1.scatter(x,signal,color='g',s=2)

	if envelopeflag:
		sp1.plot(x,envelope,color='r',linewidth=1)

	sp1.set_title('Waveform')
	sp1.set_xlabel('Time (s)')
	signalduration = float(len(signal))/fs
	sp1.set_xlim(0,signalduration)
	sp1.set_ylim(-1,1)
	sp1.grid(which="both")

except:
	print("Signal representation error."); exit()

#==============================================================================
# Plot AM envelope

try:
	signal = medfilt(signal,31)
	signal = signal / max(abs(signal))
	peakwin = 20
	x = np.linspace(0,signalduration,len(envelope))
	sp2.plot(x,envelope,color='r',linewidth=1)
	sp2.set_title('AM demodulation envelope')
	sp2.set_xlabel('Time (s)')
	sp2.set_xlim(0,signalduration)
	sp2.grid(which="both")

except:
	print("AM envelope plotting error."); exit()

#==============================================================================
# AMDF graph

try:
	sp3.set_title('FM demodulation (F0 estimation, \'pitch\' extraction) using AMDF')
	sp3.set_xlabel('Time (s), frame duration: %.1fms, step: %.1fms, frame offset: %.1fms'%(frameduration*1000, 1000*framestep/fs, 1000*frameoffsetfraction*frameduration))

	x = np.linspace(0,signalduration,len(f0list))
	sp3.set_xlim(0,signalduration)
	sp3.set_ylim(f0min, f0max)

	if envelopeflag:
		envf0ratio = int(np.floor(len(envelope)/len(f0list)))
		z = envelope[::envf0ratio]

		# Shorten or pad z to handle rounding errors
		lendiff = len(z) - len(f0list)
		if lendiff > 0: z = z[:-lendiff]
		else: z = list(z) + [0]*abs(lendiff)

		z = z/max(z)
		z = f0min + (z * (f0max-f0min) * 0.9)
		sp3.plot(x, z,linewidth=1,color='r')

	sp3.scatter(x, f0list, s=1, color='b')
	sp3.grid(which="both")

except:
	print("F0 graph error."); exit()

#==============================================================================
# Format and display figure

plt.tight_layout(pad=1, w_pad=1, h_pad=1)
plt.savefig(figurename)
plt.show()

#===========================================================================EOF
