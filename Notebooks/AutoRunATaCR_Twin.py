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
datafolder = dirs['Py_DataParentFolder']
eventsfolder = dirs['Py_CorrectedTraces']
eventsfolder = dirs['Py_CorrectedTraces']
ATaCR_Parent = dirs['Py_DataParentFolder']
catalog_full = pd.read_excel(str(project_path / '_DataArchive' / 'utilities' / 'Janiszewski_etal_2023_StationList.xlsx'))
catalog = pd.read_pickle(Path(ATaCR_Parent) / 'Catalogs' / 'sta_catalog_proxima_test.pkl')
# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(Path(ATaCR_Parent) / 'Catalogs' / 'event_record_audit.pkl')
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
# catalog = catalog.iloc[np.intersect1d(catalog.Station,['M07A'],return_indices=True)[1]]
cat = catalog.copy()
# cat = cat.reset_index()
# cat.loc[0,'Events'] = ['2012.069.07.09']
# display(cat)
event_mode = False
Minmag,Maxmag=6.0,8.0
fork = False
# STEPS = [1,2,3,4,5,6,7]
# STEPS = [1,5,6,7]
STEPS = [1,5]
# STEPS = [1,4]
# STEPS = [1,4,5,6,7]
days = 10

# ...For testing...
# days = ['2012.061','2012.062','2012.063','2012.064']
# days = [UTCDateTime.strptime('2012.061',format='%Y.%j') + i*3600*24 for i in range(30)]
# catalog = catalog[catalog.Station.isin(['M07A','M08A',])]
# catalog = catalog.reset_index()
# catalog = update_event_catalog(catalog,eventsfolder,['2012.069.07.09','2012.181.21.07'])

# ---
evf = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
tff = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/ATaCR_Data/ATaCR_Python/TF_STA')
# catalog = catalog.iloc[np.array([len(list((evf/ s.StaName).glob('*Z.SAC'))) for s in catalog.iloc])<10]

catalog = catalog.iloc[:np.where(catalog.Station=='S01W')[0][0]+1]
catalog = catalog.iloc[np.array([len(list((tff/ s.StaName).glob('*.pkl'))) for s in catalog.iloc])==0]
catalog = catalog.iloc[np.where(catalog.Station=='LD41')[0][0]:]
# catalog = catalog.iloc[np.where(catalog.Station=='J50A')[0][0]+1:]
# catalog = catalog.iloc[:np.where(catalog.Station=='CC05')[0][0]+1]
# catalog = catalog.iloc[np.where(catalog.Station=='J46C')[0][0]+1:]
cat = catalog.copy()
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
## =============================================================================== ## =============================================================================== ##|
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
## =============================================================================== ## =============================================================================== ##|
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
## =============================================================================== ## =============================================================================== ##|
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
STEPS.pop(np.where(np.array(STEPS)==1)[0][0])
for STEP in STEPS:
    for ii,Station in enumerate(cat.iloc):
        if Station.Station=='LD41':
            continue
        ## StaFolder = Path(dirs['Py_RawDayData']) / Station.StaName
        ## Files = list(StaFolder.glob('*.SAC'))
        staname = Station.StaName
        subfolder = staname + '/'
        print('[//////////////////////////]'*2)
        print('----Station: ' + staname +  ' (' + str(ii+1) + ' of ' + str(len(cat)) + ')')
        icatalog = Station.to_frame().T
        print('[//////////////////////////]'*2)
        ObsQA.TOOLS.io.Run_ATaCR(icatalog,days=days,fork=fork,event_mode=event_mode, ATaCR_Parent = ATaCR_Parent,STEPS=[STEP],log_prefix=Station.StaName,Minmag=Minmag,Maxmag=Maxmag)
    if event_mode & (STEP==3):
        Origins = Station.Origin
        Starts = Origins
        if isinstance(Origins[0],obspy.core.event.origin.Origin):
                Starts = [e.time for e in Origins]
        dateformat = '%Y.%j.%H.%M'
        hps_data_folder = (Path(dirs['Py_RawDayData']) / Station.StaName / 'HPS_Data')
        hps_data_folder.mkdir(exist_ok=True)
        days = list(np.unique([s.strftime('%Y.%j') for s in Starts]))
        [[shutil.move(fi,hps_data_folder / fi.name) for fi in list((Path(dirs['Py_RawDayData']) / Station.StaName).glob(d + '*.SAC'))] for d in days]
        # shutil.
        print('....done')
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##