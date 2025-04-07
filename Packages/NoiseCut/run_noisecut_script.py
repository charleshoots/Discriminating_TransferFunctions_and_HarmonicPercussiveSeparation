import math,pd,obspy
import numpy as np
import librosa
import librosa.display
import matplotlib.pyplot as plt
from obspy import Trace,Inventory,Stream,read
from pathlib import Path
import noisecut



# --------

# Author: Charles Hoots, 2025
# GNU Public license v3.0

# The following script will run every SAC file found in a given DataFolder through the Noisecut function and write each output to a SAC output file in OutFolder.
# Options for writing either a 24-hour or 2-hour output are available at the bottom.

# Assuming you ahave the basic packages installed (numpy, pandas, obspy, librosa, matplotlib), 
# You will also need to install one called librosa (pip install librosa works) 
# After that, you can run this script to process all of your data through noisecut.

# --------




DataFolder = Path('PATH_TO_DATA_FOLDER')
OutputFolder = Path('PATH_TO_OUTPUT_FOLDER')

files = list(DataFolder.glob('*.SAC')) #SAC files for 24-hour events. This line will search for all SAC files in the directory
# files = list(DataFolder.rglob('*.SAC')) #This line will search sub-directories too. 

#Two very important parameters to be made aware of:

# 1-The Time-Frequency Resolution Trade-off Curve:
win_length=163.84 #This is the length of the window in seconds that defines how the spectrogram is parameterized. 
# 163.84 second is the original parameterization found to work best with OBS teleseismic data. 
# A shorter length gives greater temporal resolution at the cost of frequency resolution, and vice versa for larger values.
# The function will use this 163.84 value by default if not specified.
# I'd recommend not changing this value unless you have done a fair amount of testing.
# Read Zali et al. (2023) for more details.

# 2-The similarity matrix minimum wait time (in samples):
width=None 
# In the original creators of this implimentation of HPS via NoiseCut (Zali et al. 2023), 
# the authors have this waiting factor set to TWO-HOURS after testing against teleseismic observations at OBS.
# However, the value is given in samples, and must thus scale with the particular observation sample rate, 
# which ultimately scales with resolution of the spectrogram 
# and thus the shape of the spectrogram used must define what this value is set to.
# I have it set to None here as I have built a small snippet of code 
# that will figure out this value for you based on the spectogram shape that best equals the spectrogram's 
# sample width for 2-hours that takes into account that there is a 75% overlap in the STFT.
# -
# Bottom line, the value that is defined here must equal the spectrogram sample width equal to 2-hours 
# which is by no means the same thing as the original time axis since the spectrogram uses overlapping windows.
# Nonetheless, my code takes care of that for you. So just keep that in mind.



# This line does alot of things all at once in order to keep this code very fast and extremely memory efficient.
# This one line will:
# 1-Read the SAC files one at a time
# 2-Give each one to noisecut
# 3-Take the output and save it to the OutputFolder

# 24-HOUR OUTPUTS:
# [noisecut(read(f)[0].copy(),win_length=win_length,width=width).write(OutputFolder/f.name) for f in files]
# Obviously, you can impliment this in a for loop if you wanted more control over the process, but that's up to you.


# 2-HOUR OUTPUTS:
# If you want to trim the output to just the last 2-hours of the trace for your event (event_length parameter).
# use the following lines and run them instead of the above:

event_length = 7200 #seconds
endtimes = [read(f)[0].stats.endtime for f in files]
[noisecut(read(f)[0].copy(),win_length=win_length,width=width).trim(t_end-event_length,t_end).write(OutputFolder/f.name) for t_end,f in zip(endtimes,files)]