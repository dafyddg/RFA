# rfa_single_conf.py
# Dafydd Gibbon
# Created 2021-08-16
# Modified 2022-04-18
# Code first deposited on GitHub 2021-08-16
# http://www.github.com/dafyddg/RFA
# Variable assignments for rfa_single.py
# See top level README.1st file.

"""
Note that the F0 parameters are defined in the F0 module, not here.
"""

"""
MIT license
Begin license text.
Copyright 2021 Dafydd Gibbon
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
End license text.
"""

# CSV separator
separator = ","

# Figure show and size
showgraph = True
figwidth = 10
figheight = 10

# Dendrogram box dimensions and position (AM)
amboxwidth = 0.45
amboxheight = 0.07
amboxx = 0.53
amboxy = 0.25

# Dendrogram box dimensions and position (FM)

fmboxwidth = 0.45
fmboxheight = 0.07
fmboxx = 0.53
fmboxy = 0.1

# Column chart resolution
bincount = 40

# Envelope superposition flag (spectrograms, formant tracks)
envelopeoverlay = True

# Waterfall or column chart display switch
waterfall = True

# Waterfall parameters
waterfallincplus = 2
specheatmaptype = "YlOrRd"
specgramdotsize = 60
fontsize = 10

# Minimum and maximum spectrum and spectrogram frequencies
amspecfreqmin = 0
amspecfreqmax = 5

fmspecfreqmin = 0
fmspecfreqmax = 5

# Minimum value of spectral magnitude line for formants
amformantlimit = 0.3
fmformantlimit = 0.3

# Number of spectral magnitude peaks marked in spectrum display
magscount = 6

# Optional high magnitude emphasis
spectrumpower = 1

# Optional brute force downsampling factor (spectral analysis)
# This is not true downsampling, but stepping interval selection
specdownsample = 4

# Set spectrogram moving window duration in seconds for FFT analysis
specwindowsecs = 3

# Choose number of equally spaced rows in spectrogram matrix (default 50)
specstrides = 100	# yields same length spectrograms, for comparison

# EOF
