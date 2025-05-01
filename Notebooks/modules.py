# ====================================================================================
# ===================================== IMPORTS ======================================
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| PATH STUFF |||||||||||||||||||||||||||||||||||||
from pathlib import Path;import shutil,sys,os;import numpy as np;import pandas as pd
import warnings,fnmatch,operator,itertools
from IPython.display import clear_output
parent = Path('/Users/charlesh/Documents/Codes');sys.path.insert(1,str(parent))
sys.path.insert(1,str(Path(parent/'OBS_Methods/NOISE')))
sys.path.append(str(Path(parent/'OBS_Methods/NOISE/METHODS')))
project_path = parent /'OBS_Methods'/'NOISE'/'Research'
sys.path.insert(1,str(project_path))
sys.path.append(str(project_path / 'Packages'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'))
sys.path.insert(0, str(project_path / 'Packages' / 'CompCode'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'/ 'OBStools'))
sys.path.insert(0, str(project_path / 'Notebooks'))
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| BASIC PLOT STUFF |||||||||||||||||||||||||||||||
from obspy.core import AttribDict
import matplotlib.pyplot as plt;import matplotlib.image as img
import matplotlib.gridspec as gridspec;import matplotlib
from cmcrameri import cm
cm.Categorical=AttribDict({k:cm.__dict__['batlow'] for k in cm._cmap_names_categorical})
cm.Sequential=AttribDict({k:cm.__dict__['batlow'] for k in cm._cmap_names_sequential})
import matplotlib.colors as mcolors;import matplotlib.cm as cm2
import matplotlib as mpl
# from matplotlib.collections import PatchCollection
# from matplotlib.patches import Rectangle
# from branca.element import Template, MacroElement
# import cartopy
# import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||| BASIC MATH AND DSP STUFF |||||||||||||||||||||||||||||
import math;import numpy as np
os.environ['PYDEVD_WARN_SLOW_RESOLVE_TIMEOUT'] = '2'
import librosa,librosa.display
from numpy import linalg as eigen
from scipy.signal import stft, detrend;import scipy
from scipy.stats import norm;import scipy.stats as stats;from scipy import fft
from scipy.interpolate import RBFInterpolator, InterpolatedUnivariateSpline
from scipy.signal import csd as _csd
import cmath
# ____________________________________________________________________________________
# |||||||||||||||||||||||||||||||||| OBSPY STUFF |||||||||||||||||||||||||||||||||||||
from obspy import Trace,Inventory
import obspy;import pickle as pkl
from obspy.clients.fdsn import Client;import datetime,re
from obspy.geodetics import locations2degrees
from obspy.geodetics import degrees2kilometers
from obspy.core.inventory.inventory import read_inventory
from obspy import read
from obspy import Catalog
from obspy import Stream

# ____________________________________________________________________________________
# ||||||||||||||||||||||| BRANCHED COMMUNITY GITs/PACKAGES |||||||||||||||||||||||||||
from OrientPy import * #DLOPy
from NoiseCut.src import * #Noisecut
import obstools as obs #ATaCR
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| MY CODES |||||||||||||||||||||||||||||||||||||||
import ObsQA;import ObsQA as ob
from ObsQA.TOOLS import io
from ObsQA.TOOLS import TaskMaster
from ObsQA.TOOLS.io import get_noise_metrics
from ObsQA.OBSM.classes import OBSMetrics as OBSM
from ObsQA.OBSM.classes import NoiseWrapper
from ObsQA.TOOLS.io import *
from ObsQA import *;from comp_tools import *
import local_tools as lt
from quick_class import *
from get_reports import *
from helper_functions import *
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||| USEFUL FOLDERS TO KNOW |||||||||||||||||||||||||||||||