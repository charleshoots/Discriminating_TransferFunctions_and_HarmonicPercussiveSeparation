import math
import pandas as pd
import obspy
import numpy as np
import librosa
import librosa.display
import matplotlib.pyplot as plt
from obspy import Trace,Inventory,Stream,read
from pathlib import Path
import noisecut
# -----
# A simple helper code to process traces through NoiseCut
# --------
# Author: Charles Hoots, 2024
# GNU Public license v3.0

# The following script will run every SAC file found in a given DataFolder through the Noisecut function and write each output to a SAC output file in OutFolder.
# Options for writing either a 24-hour or 2-hour output are available at the bottom.

# This script is just a task manager for traces to send through NoiseCut, an HPS de-noising algorithm by Zahra Zali.
# Read Zali et al. (2023) for more details.

# Assuming you ahave the basic packages installed (numpy, pandas, obspy, librosa, matplotlib), 
# You will also need to install one called librosa (pip install librosa works) 
# After that, you can run this script to process all of your data through noisecut.
# --------



DataFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data/rmresp')
# DataFolder = Path('PATH_TO_DATA_FOLDER')
# OutputFolder = Path('PATH_TO_OUTPUT_FOLDER')

files = list(DataFolder.rglob('*.SAC')) #Searches entire directory (and sub-directories) for SAC files.
files = {fi:np.array(files)[np.isin([of.parent.name for of in files],fi)] for fi in np.unique([of.parent.name for of in files])}
for k in files.keys():files[k]={ev:np.array(files[k]) [np.isin(['.'.join(fi.name.split('.')[:4]) for fi in files[k]],ev)] for ev in (['.'.join(fi.name.split('.')[:4]) for fi in files[k]])}

#Three important parameters to be made aware of:

# | 1 |-The Time-Frequency Resolution Trade-off Curve:
win_length=163.84 #This is the length of the window in seconds that defines how the spectrogram is parameterized. 
# 163.84 second is the original parameterization found to work best with OBS teleseismic data. 
# A shorter length gives greater temporal resolution at the cost of frequency resolution, and vice versa for larger values.
# The function will use this 163.84 value by default if not specified.
# I'd recommend not changing this value unless you have done a fair amount of testing.
# Read Zali et al. (2023) for more details on these parameters.

# 2-The similarity matrix minimum wait time (in samples):
width=None 
# When set to None, will force NoiseCut to use defaults (2 hour wait time).
# | 3 | - Trims the output to a 2-hour length.
trim_output=True

# This one line will:
# 1-Read the SAC files one at a time
# 2-Process each one in noisecut
# 3-Take the output ,trim it (if enabled), then save it to the OutputFolder.

if not trim_output:#24 hour outputs
    [noisecut(read(f)[0].copy(),win_length=win_length,width=width).write(OutputFolder/f.name) for f in files]
if trim_output: #2 hour outputs
    event_length = 7200 #seconds
    endtimes = [read(f)[0].stats.endtime for f in files]
    [noisecut(read(f)[0].copy(),win_length=win_length,width=width).trim(t_end-event_length,t_end).write(OutputFolder/f.name) for t_end,f in zip(endtimes,files)]
