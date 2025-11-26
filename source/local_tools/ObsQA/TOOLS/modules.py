# ====================================================================================
# ===================================== IMPORTS ======================================
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| PATH STUFF |||||||||||||||||||||||||||||||||||||
from pathlib import Path;import shutil,sys,os;import numpy as np;import pandas as pd
import warnings,fnmatch,operator,itertools
from IPython.display import clear_output
#This function yields the standard path network called on for everything in this project.
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| BASIC PLOT STUFF |||||||||||||||||||||||||||||||
import matplotlib.pyplot as plt;import matplotlib.image as img
import matplotlib.gridspec as gridspec;import matplotlib
from cmcrameri import cm
import matplotlib.colors as mcolors;import matplotlib.cm as cm2
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
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
from obspy.core.inventory.inventory import read_inventory
from obspy import read
# ____________________________________________________________________________________
# ||||||||||||||||||||||| BRANCHED COMMUNITY GITs/PACKAGES |||||||||||||||||||||||||||
# from modules import *
from NoiseCut.src import * #Noisecut
import obstools as obs #ATaCR
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||| MY CODES |||||||||||||||||||||||||||||||||||||||
import local_tools.ObsQA as ObsQA;import local_tools.ObsQA as ob
from local_tools.ObsQA.TOOLS import io
from local_tools.ObsQA.TOOLS import TaskMaster
from local_tools.ObsQA.TOOLS.io import get_noise_metrics
from local_tools.ObsQA.OBSM.classes import OBSMetrics as OBSM
from local_tools.ObsQA.TOOLS.io import *
from local_tools.ObsQA import *;from local_tools.comp_tools import *
from local_tools.helper_functions import *
# from local_plotter import *  # local_plotter module not present in this project
from local_tools import *
from local_tools.quick_class import *
# ____________________________________________________________________________________
# ||||||||||||||||||||||||||||| USEFUL FOLDERS TO KNOW |||||||||||||||||||||||||||||||