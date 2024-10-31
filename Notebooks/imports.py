# ===================================================================================================
# ============================================ IMPORTS ==============================================
# ===================================================================================================
from pathlib import Path;import shutil,sys,os;import numpy as np;import pandas as pd
parent = Path('/Users/charlesh/Documents/Codes');sys.path.insert(1,str(parent))
sys.path.insert(1,str(Path(parent/'OBS_Methods/NOISE')))
sys.path.append(str(Path(parent/'OBS_Methods/NOISE/METHODS')))
project_path = parent /'OBS_Methods'/'NOISE'/'ATACR_HPS_Comp'
sys.path.insert(1,str(project_path))
sys.path.append(str(project_path / 'Packages'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'))
sys.path.insert(0, str(project_path / 'Packages' / 'CompCode'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'/ 'OBStools'))
from OrientPy import *
import math,scipy
import scipy
import numpy as np
from scipy.signal import stft, detrend
os.environ['PYDEVD_WARN_SLOW_RESOLVE_TIMEOUT'] = '2'
from obspy import Trace
import librosa,librosa.display
import matplotlib.pyplot as plt
import matplotlib.image as img
import matplotlib.gridspec as gridspec
import matplotlib
import obspy
import pickle as pkl
import glob as g
from obspy.clients.fdsn import Client
import datetime,re
from numpy import linalg as eigen
import matplotlib.colors as mcolors
import matplotlib.cm as cm2
from scipy.stats import norm
import scipy.stats as stats
from scipy import fft
from cmcrameri import cm
from scipy.interpolate import RBFInterpolator, InterpolatedUnivariateSpline #<----Experimental
import ObsQA;import ObsQA as ob
from ObsQA.TOOLS import io
from ObsQA.TOOLS.io import get_Noise
from ObsQA.OBSM.classes import OBSMetrics as OBSM
from ObsQA.TOOLS.io import *
from NoiseCut.Source.src import *
# from ObsQA.plots import qtp
from ObsQA import *
import obstools as obs
import cmath
from comp_tools import *
from pathlib import Path
from scipy.signal import csd as _csd
from scipy.signal import stft, detrend
from helper_functions import *
DataFolder = project_path / '_DataArchive'/ 'ATaCR_Data'
archive =  project_path / '_DataArchive'
plotfolder = project_path / '_FigureArchive' / '_GEN5'
# ===================================================================================================
# ============================================  LOAD DATA ===========================================
# ===================================================================================================
dirs = io.dir_libraries(str(DataFolder))
# catfolder = Path(dirs['Py_DataParentFolder']) / 'Catalogs'
# eventsfolder = dirs['Py_CorrectedTraces']
# depreciated: catalog = pd.read_pickle(catfolder / 'event_catalog_updated.pkl')
# depreciated: catalog = pd.read_pickle(catfolder / 'sta_catalog_proxima_test.pkl')
catalog = pd.read_pickle(dirs.Catalogs / 'sta_catalog_101524.pkl')

Folder = Path(plotfolder) / 'MeetingFigs'
Folder.mkdir(exist_ok=True)

# catalog = catalog[catalog.Station.isin(['M08A','M07A'])]
# catalog = catalog.reset_index()
# catalog = update_event_catalog(catalog,eventsfolder,['2012.069.07.09','2012.181.21.07'])
k=1


def write_pickle(file,var):
    import pickle
    with open(str(file), 'wb') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Saved to :' + str(file))
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        b = pickle.load(handle)
    return b