# ====================================================================================
# ===================================== IMPORTS ======================================
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| PATH STUFF |||||||||||||||||||||||||||||||||||||
from pathlib import Path;import shutil,sys,os;import numpy as np;import pandas as pd
import warnings,fnmatch,operator,itertools
import numbers
from IPython.display import clear_output
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| BASIC PLOT STUFF |||||||||||||||||||||||||||||||
from obspy.core import AttribDict
import matplotlib.pyplot as plt;import matplotlib.image as img
import matplotlib.gridspec as gridspec;import matplotlib
import matplotlib as mpl;from matplotlib.colors import ListedColormap
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
def np2(n):return int(round(2**np.ceil(np.log2(n))))
# ____________________________________________________________________________________
# ||||||||||||||||||||||| BRANCHED COMMUNITY GITs/PACKAGES |||||||||||||||||||||||||||
# from OrientPy import * #DLOPy
from NoiseCut.src import * #Noisecut
import obstools as obs #ATaCR
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| MY CODES |||||||||||||||||||||||||||||||||||||||
from local_tools import ObsQA as ob
from local_tools.ObsQA.TOOLS import io
# from local_tools.ObsQA.TOOLS import TaskMaster
# from local_tools.ObsQA.TOOLS.io import get_noise_metrics
from local_tools.ObsQA.OBSM.classes import OBSMetrics as OBSM
from local_tools.ObsQA.OBSM.classes import NoiseWrapper
from local_tools.ObsQA.TOOLS.io import *
from local_tools.ObsQA import *
from local_tools.comp_tools import *
import local_tools
import local_tools as lt

from local_tools import math
from local_tools import cat
from local_tools import io
from local_tools import plots
from local_tools import LabeledMatrix
from local_tools import get_reports
from local_tools import helper_functions
from local_tools import quick_class



from local_tools import io as io
from local_tools.math import *
from local_tools.quick_class import *
from local_tools.get_reports import *
from local_tools.helper_functions import *
from local_tools.plots import *
from local_tools.LabeledMatrix import AggregateMeasurements as LM
# from octave_average_ansi import octave_average_ansi
# ____________________________________________________________________________________