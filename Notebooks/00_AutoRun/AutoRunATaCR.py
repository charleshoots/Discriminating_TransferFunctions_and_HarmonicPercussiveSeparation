import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
from imports import *
unpacker=locals().update;unpack=lambda Args,keys=None:[unpacker({k:Args[k]}) for k in [keys if keys is not None else list(Args.keys())][0]]

hps_staquery_output = Path(os.getcwd())/'_DataArchive/HPS_Data/sta_query.pkl'
Args=AttribDict()

ATaCR_Steps = {1:'Metadata ',
2:'Download events',
3:'Downlosd noise',
4:'Day Spectra',
5:'Average Spectra',
6:'Transfer Functions',
7:'Correct Events'}
step_time = lambda: print(f"{'X'*30}{'\n'}{'X'*30}{'\n'}{'\n'*1}|{ATaCR_Steps[STEP]}|\n|EXECUTION TIME: {(time.time() - start_time)/60 :.2f} minutes{'\n'*2}{'X'*30}{'\n'}{'X'*30}{'\n'}")
# ---------------------------------------------------------------------------------------------------
# ============================================ LOAD DATA ===========================================
# ---------------------------------------------------------------------------------------------------
### ===============================================================================
## -------------------------------------------------------
## Step-1: Station Metadata. Step a0 in ML-ATaCR. Always run this. Note: This step is implicit and will always run whether you ask for it or not.
## Step-2: Download event data. Step a3 in ML-ATaCR.
## Step-3: Download day data. Step a2 in ML-ATaCR.
## Step-4: Daily Spectra. Step b1 in ML-ATaCR.
## Step-5: Clean and Average Daily Spectra. Step b2 in ML-ATaCR.
## Step-6: Calculate transfer functions. Step b3 ML-ATaCR.
## Step-7: Correct events. Step b4 in ML-ATaCR.
## ===============================================================================
## ===============================================================================
# catalog = catalog[catalog.Station.isin(['FN07A','FN14A','J50A','M07A'])]
# # reverse catalog
# catalog = catalog.iloc[list(np.flip([a for a in range(len(catalog))]))]
# catalog = catalog.iloc[np.where(catalog.StaName=='YO.X01')[0][0]:]
# -----------------------------------------------------------------------------------
## =============================================================================== ##
cat = catalog.copy()
# cat = cat.iloc[89:]
# cat = cat.loc[['YL.A10W','X9.BB420','7D.J28C','7D.G34D']]
# cat = cat.loc[['7D.FS08D']]
Args.Minmag,Args.Maxmag=6.0,8.0
Args.cleanspectra_flags = '--figQC --figAverage --figCoh --figCross --save-fig'
Args.dailyspectra_flags='--figQC --figAverage --figCoh --save-fig'
## =============================================================================== ##
# -----------------------------------------------------------------------------------

Args.STEPS = [6]
Args.fork = False
Args.event_mode = False 
Args.event_window = 3600*2 #7200
Args.channels = 'Z,P,12'
# Args.channels = 'P,12'
# Args.channels = 'Z'
Args.ATaCR_Parent = dirs.ATaCR
Args.days=10
Args.ovr=True


# ----=----=----=----
# |
# ----=----=----=----

# -----# -----# -----# -----# -----# -----# -----# -----# -----
import time
if Args.event_mode:Args.staquery_output = hps_staquery_output
else:Args.staquery_output = './sta_query.pkl'
if 1 in Args.STEPS:Args.STEPS.pop(np.where(np.array(Args.STEPS)==1)[0][0])
for STEP in Args.STEPS:
    Args.STEP = [STEP]
    start_time = time.time()
    for ii,Station in enumerate(cat.iloc):
        # if STEP>4:Args.ovr=True
        # if Station.StaName=='7D.G17B':Args.ovr=False
        # if Station.StaName=='7D.G25B':Args.ovr=False
        Args.log_prefix=Station.StaName
        print('[//////////////////////////]'*2)
        print('----Station: ' + Station.StaName +  ' (' + str(ii+1) + ' of ' + str(len(cat)) + ')')
        Args.catalog = Station.to_frame().T
        print('[//////////////////////////]'*2)
        Args.message = f'Station {ii+1}/{len(cat)}'
        ObsQA.TOOLS.io.Run_ATaCR(Args)
    step_time()

## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##

## Step-1: Station Metadata. Step a0 in ML-ATaCR. Always run this. Note: This step is implicit and will always run whether you ask for it or not.
## Step-2: Download event data. Step a3 in ML-ATaCR.
## Step-3: Download day data. Step a2 in ML-ATaCR.
## Step-4: Daily Spectra. Step b1 in ML-ATaCR.
## Step-5: Clean and Average Daily Spectra. Step b2 in ML-ATaCR.
## Step-6: Calculate transfer functions. Step b3 ML-ATaCR.
## Step-7: Correct events. Step b4 in ML-ATaCR.