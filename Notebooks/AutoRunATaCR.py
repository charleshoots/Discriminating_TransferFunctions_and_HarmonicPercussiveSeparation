# # ------------------------------------------------------------------------------------------------------------------------
# from pathlib import Path
# project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp')
# import shutil
# import numpy as np
# import pandas as pd
# import sys
# from obspy import UTCDateTime
# sys.path.append(str(project_path / 'Packages'))
# sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'))
# sys.path.insert(0, str(project_path / 'Packages' / 'CompCode'))
# sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'/ 'OBStools'))
# import ObsQA
# from comp_tools import *
# # ---------------------------------------------------------------------------------------------------
# # ============================================ FOLDERS ==============================================
# # ---------------------------------------------------------------------------------------------------
# ATaCR_DataFolder = str(project_path / '_DataArchive' / 'ATaCR_Data')
# dirs = OBS.TOOLS.io.dir_libraries(ATaCR_DataFolder)
# catalog_full = pd.read_excel(str(project_path / '_DataArchive' / 'utilities' / 'Janiszewski_etal_2023_StationList.xlsx'))
# catalog = pd.read_pickle(dirs.Catalogs / 'sta_catalog_101524.pkl')
from imports import *
hps_staquery_output = Path(os.getcwd())/'_DataArchive/HPS_Data/sta_query.pkl'
hps_staquery_output = dirs.Events_HPS.parent/'sta_query.pkl'
# ---------------------------------------------------------------------------------------------------
# ============================================ LOAD DATA ===========================================
# ---------------------------------------------------------------------------------------------------
### ===============================================================================
## STEPS = [1,2,3,4,5,6,7] ##Absolutely every step - Downloading adds hour(s) or more to the process
## STEPS = [2,3] ##Everything but the download steps - About 4min for six stations.
## STEPS = [4,5,6,7] ##Everything but the download steps - About 4min for six stations
## -------------------------------------------------------
## -------------------------------------------------------
## Step-1: Station Metadata. Step a0 in ML-ATaCR. Always run this.
## Step-2: Download event data. Step a3 in ML-ATaCR.
## Step-3: Download day data. Step a2 in ML-ATaCR.
## Step-4: Daily Spectra. Step b1 in ML-ATaCR.
## Step-5: Clean and Average Daily Spectra. Step b2 in ML-ATaCR.
## Step-6: Calculate transfer functions. Step b3 ML-ATaCR.
## Step-7: Correct events. Step b4 in ML-ATaCR.
## ===============================================================================
## ===============================================================================
Minmag,Maxmag=6.0,8.0
fork = False
event_mode = True
# STEPS = [1,2,3,4,5,6,7]
# STEPS = [1,3,4,5,6,7]
# STEPS = [1,4,5,6,7]
# STEPS = [7]
# STEPS = [1,4,5,6,7]
STEPS = [5]
# STEPS,event_mode = [3],True
# days = 10
# ...For testing...
# days = ['2012.061','2012.062','2012.063','2012.064']
# catalog = catalog[catalog.Station.isin(['FN07A','FN14A','J50A','M07A'])]
# catalog = catalog[catalog.Network.isin(['ZA','XF','XO','YL'])]
# reverse catalog
# catalog = catalog.iloc[list(np.flip([a for a in range(len(catalog))]))]
catalog = catalog.iloc[np.where(catalog.StaName=='ZA.B04')[0][0]:]

## =============================================================================== ## =============================================================================== ##
cat = catalog.copy() 
# =============================================================================== #
## =============================================================================== ## =============================================================================== ##
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## =============================================================================== ## =============================================================================== ##
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## =============================================================================== ## =============================================================================== ##
if event_mode:
    staquery_output = hps_staquery_output
else:
    staquery_output = './sta_query.pkl'
if 1 in STEPS:
    STEPS.pop(np.where(np.array(STEPS)==1)[0][0])

# event_window = 7200
event_window = 3600*4
channels = 'P,12'
ATaCR_Parent = dirs.ATaCR
for STEP in STEPS:
    for ii,Station in enumerate(cat.iloc):
        ## StaFolder = Path(dirs['Py_RawDayData']) / Station.StaName
        ## Files = list(StaFolder.glob('*.SAC'))
        staname = Station.StaName
        subfolder = staname + '/'
        print('[//////////////////////////]'*2)
        print('----Station: ' + staname +  ' (' + str(ii+1) + ' of ' + str(len(cat)) + ')')
        icatalog = Station.to_frame().T
        print('[//////////////////////////]'*2)
        message = f'Station {ii+1}/{len(cat)}'
        ObsQA.TOOLS.io.Run_ATaCR(icatalog,fork=fork,message=message,staquery_output=staquery_output,event_mode=event_mode, ATaCR_Parent = ATaCR_Parent,STEPS=[STEP],log_prefix=Station.StaName,Minmag=Minmag,Maxmag=Maxmag,event_window=event_window,channels=channels)

## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
