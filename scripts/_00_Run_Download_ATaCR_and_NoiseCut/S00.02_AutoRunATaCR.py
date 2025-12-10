### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *
import time
# step time value
step_time = lambda: print(f"{'X'*30}{'\n'}{'X'*30}{'\n'}{'\n'*1}|{ATaCR_Steps[STEP]}|\n|EXECUTION TIME: {(time.time() - start_time)/60 :.2f} minutes{'\n'*2}{'X'*30}{'\n'}{'X'*30}{'\n'}")
# unpacker value
unpacker=locals().update;unpack=lambda Args,keys=None:[unpacker({k:Args[k]}) for k in [keys if keys is not None else list(Args.keys())][0]]
# hps staquery output value
hps_staquery_output = Path(os.getcwd())/'_DataArchive/HPS_Data/sta_query.pkl'
# Args value
Args=AttribDict()
# ATaCR Steps value
ATaCR_Steps = {1:'Metadata ',
2:'Download events',
3:'Downlosd noise',
4:'Day Spectra',
5:'Average Spectra',
6:'Transfer Functions',
7:'Correct Events'}

## =============================================================================== ##
# ============================================ OPTIONS ==============================
## =============================================================================== ##
## Step-1: Station Metadata. Step a0 in ML-ATaCR. Always run this. Note: This step is implicit and will always run whether you ask for it or not.
## Step-2: Download event data. Step a3 in ML-ATaCR.
## Step-3: Download day data. Step a2 in ML-ATaCR.
## Step-4: Daily Spectra. Step b1 in ML-ATaCR.
## Step-5: Clean and Average Daily Spectra. Step b2 in ML-ATaCR.
## Step-6: Calculate transfer functions. Step b3 ML-ATaCR.
## Step-7: Correct events. Step b4 in ML-ATaCR.
## ===============================================================================
icat=catalog.r.copy() #ATaCR is processed on a receiver catalog, not a source-receiver catalog.
# If you are to subset the data catalog to be processed, do it here.
# icat=icat.loc[['7D.J25A']].copy()
Args.Minmag,Args.Maxmag=6.0,8.0
Args.cleanspectra_flags = '--figQC --figAverage --figCoh --figCross --save-fig'
Args.dailyspectra_flags='--figQC --figAverage --figCoh --save-fig'
Args.STEPS = [4,5,6,7] #Steps to run. Choose from 1 to 7. See ATaCR_Steps dict for reference.
Args.fork = False #Whether to use multiprocessing (forking) or not.
Args.event_mode = False #For NoiseCut downloads. Enable this to download event data in 24-hour segments (last 2-hrs being the default trace window) rather than 2-hours (ATaCR default)
Args.event_window = 3600*2 #7200
Args.channels = 'Z,P,12' #Channels to process
Args.ATaCR_Parent = dirs.ATaCR #The parent directory of the ATaCR codebase.
Args.days=10 #Number of days of noise data to download per station.
Args.ovr=True #Whether to overwrite existing outputs or not.



# # Run ATaCR for each station in the catalog at the steps (STEPS) with the options specified.
## =============================================================================== ##
if Args.event_mode:Args.staquery_output = hps_staquery_output
else:Args.staquery_output = './sta_query.pkl'
if 1 in Args.STEPS:Args.STEPS.pop(np.where(np.array(Args.STEPS)==1)[0][0])
# loop over STEPS
for STEP in Args.STEPS:
    Args.STEP = [STEP]
    # start time value
    start_time = time.time()
    for ii,Station in enumerate(icat.iloc):
        Args.log_prefix=Station.StaName
        print('[//////////////////////////]'*2)
        print('----Station: ' + Station.StaName +  ' (' + str(ii+1) + ' of ' + str(len(icat)) + ')')
        Args.catalog = Station.to_frame().T
        print('[//////////////////////////]'*2)
        Args.message = f'Station {ii+1}/{len(icat)}'
        ObsQA.TOOLS.io.Run_ATaCR(Args)
    step_time()
## =============================================================================== ##
