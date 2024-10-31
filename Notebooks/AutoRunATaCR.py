# ------------------------------------------------------------------------------------------------------------------------
from pathlib import Path
project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp')
import shutil
import numpy as np
import pandas as pd
import sys
from obspy import UTCDateTime
sys.path.append(str(project_path / 'Packages'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'))
sys.path.insert(0, str(project_path / 'Packages' / 'CompCode'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'/ 'OBStools'))
import ObsQA
from comp_tools import *
# ---------------------------------------------------------------------------------------------------
# ============================================ FOLDERS ==============================================
# ---------------------------------------------------------------------------------------------------
project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp')
ATaCR_DataFolder = str(project_path / '_DataArchive' / 'ATaCR_Data')
dirs = OBS.TOOLS.io.dir_libraries(ATaCR_DataFolder)[1]
catfolder =  Path(dirs['Py_DataParentFolder']) / 'Catalogs'
datafolder = dirs['Py_DataParentFolder']
eventsfolder = dirs['Py_CorrectedTraces']
eventsfolder = dirs['Py_CorrectedTraces']
ATaCR_Parent = dirs['Py_DataParentFolder']
catalog_full = pd.read_excel(str(project_path / '_DataArchive' / 'utilities' / 'Janiszewski_etal_2023_StationList.xlsx'))
catalog = pd.read_pickle(catfolder / 'sta_catalog_101524.pkl')

evaudit = pd.read_pickle(Path(ATaCR_Parent) / 'Catalogs' / 'event_record_audit.pkl')
hps_staquery_output = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data/sta_query.pkl')

# ---------------------------------------------------------------------------------------------------
# ============================================ LOAD DATA ===========================================
# ---------------------------------------------------------------------------------------------------
# ________________________________________________________________________________________________________________________________________________________________________
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# catalog = catalog[catalog.Station=='M07A']
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ________________________________________________________________________________________________________________________________________________________________________
## ===============================================================================
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
event_mode = False
# STEPS = [1,2,3,4,5,6,7]
# STEPS = [1,3,4,5,6,7]
# STEPS = [1,4,5,6,7]
# STEPS = [7]
# STEPS = [1,4,5,6,7]
STEPS = [7]
# STEPS,event_mode = [3],True
# days = 10
# ...For testing...
# days = ['2012.061','2012.062','2012.063','2012.064']
# days = [UTCDateTime.strptime('2012.061',format='%Y.%j') + i*3600*24 for i in range(30)]
# catalog = catalog[catalog.Station.isin(['M07A','M08A',])]
# catalog = catalog.reset_index()
# catalog = update_event_catalog(catalog,eventsfolder,['2012.069.07.09','2012.181.21.07'])
# catalog = catalog.iloc[np.where(catalog.Station=='J44C')[0][0]].to_frame().T
catalog = catalog[catalog.Station.isin(['FN07A','FN14A','J50A','M07A'])]

## =============================================================================== ## =============================================================================== ##
cat = catalog.copy() # =============================================================================== #
## =============================================================================== ## =============================================================================== ##
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## =============================================================================== ## =============================================================================== ##
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## =============================================================================== ## =============================================================================== ##
dlopy_data = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/DLOPY_Data/sta_query.pkl')
if event_mode:
    staquery_output = hps_staquery_output
else:
    staquery_output = './sta_query.pkl'
if 1 in STEPS:
    STEPS.pop(np.where(np.array(STEPS)==1)[0][0])

# event_window = 7200
event_window = 3600*4
channels = 'Z,P,12'
staquery_output = dlopy_data
ATaCR_Parent = dlopy_data.parent
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
        ObsQA.TOOLS.io.Run_ATaCR(icatalog,fork=fork,staquery_output=staquery_output,event_mode=event_mode, ATaCR_Parent = ATaCR_Parent,STEPS=[STEP],log_prefix=Station.StaName,Minmag=Minmag,Maxmag=Maxmag,event_window=event_window,channels=channels)

## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
