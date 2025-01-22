from modules import *
os.system('cls' if os.name == 'nt' else 'clear')
# kk
# # ====================================================================================
# # ===================================== IMPORTS ======================================
# # ____________________________________________________________________________________
# # ||||||||||||||||||||||||||||||||||| PATH STUFF |||||||||||||||||||||||||||||||||||||
# from pathlib import Path;import shutil,sys,os;import numpy as np;import pandas as pd
# import warnings,fnmatch,operator,itertools
# from IPython.display import clear_output
# parent = Path('/Users/charlesh/Documents/Codes');sys.path.insert(1,str(parent))
# sys.path.insert(1,str(Path(parent/'OBS_Methods/NOISE')))
# sys.path.append(str(Path(parent/'OBS_Methods/NOISE/METHODS')))
# project_path = parent /'OBS_Methods'/'NOISE'/'Research'
# sys.path.insert(1,str(project_path))
# sys.path.append(str(project_path / 'Packages'))
# sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'))
# sys.path.insert(0, str(project_path / 'Packages' / 'CompCode'))
# sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'/ 'OBStools'))
# sys.path.insert(0, str(project_path / 'Notebooks'))
# # ____________________________________________________________________________________
# # ||||||||||||||||||||||||||||||||||| BASIC PLOT STUFF |||||||||||||||||||||||||||||||
# import matplotlib.pyplot as plt;import matplotlib.image as img
# import matplotlib.gridspec as gridspec;import matplotlib
# from cmcrameri import cm
# import matplotlib.colors as mcolors;import matplotlib.cm as cm2
# import matplotlib as mpl
# from matplotlib.collections import PatchCollection
# from matplotlib.patches import Rectangle
# # ____________________________________________________________________________________
# # ||||||||||||||||||||||||||||| BASIC MATH AND DSP STUFF |||||||||||||||||||||||||||||
# import math;import numpy as np
# # os.environ['PYDEVD_WARN_SLOW_RESOLVE_TIMEOUT'] = '2'
# import librosa,librosa.display
# from numpy import linalg as eigen
# from scipy.signal import stft, detrend;import scipy
# from scipy.stats import norm;import scipy.stats as stats;from scipy import fft
# from scipy.interpolate import RBFInterpolator, InterpolatedUnivariateSpline
# from scipy.signal import csd as _csd
# import cmath
# # ____________________________________________________________________________________
# # |||||||||||||||||||||||||||||||||| OBSPY STUFF |||||||||||||||||||||||||||||||||||||
# from obspy import Trace,Inventory
# import obspy;import pickle as pkl
# from obspy.clients.fdsn import Client;import datetime,re
# from obspy.geodetics import locations2degrees
# from obspy.core.inventory.inventory import read_inventory
# from obspy import read
# # ____________________________________________________________________________________
# # ||||||||||||||||||||||| BRANCHED COMMUNITY GITs/PACKAGES |||||||||||||||||||||||||||
# from OrientPy import * #DLOPy
# from NoiseCut.Source.src import * #Noisecut
# import obstools as obs #ATaCR
# # ____________________________________________________________________________________
# # ||||||||||||||||||||||||||||||||||| MY CODES |||||||||||||||||||||||||||||||||||||||
# import ObsQA;import ObsQA as ob
# from ObsQA.TOOLS import io
# from ObsQA.TOOLS import TaskMaster
# from ObsQA.TOOLS.io import get_noise_metrics
# from ObsQA.OBSM.classes import OBSMetrics as OBSM
# from ObsQA.TOOLS.io import *
# from ObsQA import *;from comp_tools import *
# from helper_functions import *
# from local_plotter import *
# from local_tools import *
# from quick_class import *
# # ____________________________________________________________________________________
# # ||||||||||||||||||||||||||||| USEFUL FOLDERS TO KNOW |||||||||||||||||||||||||||||||
project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research')
# project_path = Path('/Volumes/ECHO/Research')
os.chdir(project_path)
dirs = io.dir_libraries(str(project_path))

archive =  project_path / '_DataArchive'
plotfolder = project_path / '_FigureArchive' / '_GEN6'


dirs.Plots = plotfolder;dirs.Archive=archive
dirs.Events_HPS = dirs.Archive / 'HPS_Data' / 'Data'
dirs.Analysis=dirs.Archive/'Analysis'
# ____________________________________________________________________________________
# |||||||||||||||||||||||||||||||||||| Catalogs ||||||||||||||||||||||||||||||||||||||

# HJan23=pd.read_excel(dirs.Catalogs/'Janiszewski_etal_2023_StationList.xlsx')
HJan23=pd.read_pickle(dirs.Catalogs/'Janiszewski_etal_2023_StationList.pkl')
HJan23['Good_Channels']=HJan23.T[-4:].sum().T==4


# catalog = pd.read_pickle(dirs.Catalogs / 'sta_catalog_111524c.pkl')
catalog = pd.read_pickle(dirs.Catalogs /'Catalog_010325.pkl')

# catalog_inventory=Inventory()
# for c in catalog.Inventory:catalog_inventory+=c
# disp=display
# Tweaking the all-caps for experiment names
catalog.Experiment.replace('CASCADIA KECK','KECK',inplace=True)
catalog.Experiment.replace('CASCADIA INITIATIVE','Cascadia',inplace=True)
catalog.Experiment.replace('ALBACORE','Albacore',inplace=True)
catalog.Experiment.replace('MARIANA','Mariana',inplace=True)
catalog.Experiment.replace('PAPUA','Papua',inplace=True)
catalog.Experiment.replace('LAU','Lau',inplace=True)
catalog.Experiment.replace('NOMELT','No Melt',inplace=True)
catalog.Experiment.replace('BLANCO','Blanco',inplace=True)
catalog.set_index('StaName', inplace=True,drop=False)
clear_output(wait=False)

ColorStandard=AttribDict()
ColorStandard.instrument = {'B2': '#e31a1c','KE': '#b2df8a','AB': '#a6cee3','BA': '#cab2d6','AR': '#ff7f00','TRM': '#1f78b4','BG': '#33a02c','BD': '#6a3d9a'}
ColorStandard.seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
ColorStandard.components = {'ZP':'#0c51a6','ZZ':'#2d8297','Z1':'#70cbc0','Z2':'#4f86c5',
'ZZ':'#0c51a6','11':'#2a7e93','22':'#7370cb','PP':'#4f86c5'}
ColorStandard.network = {'2D': '#d2ad90','7A': '#94530c','7D': '#64783c','X9': '#c04797','XF': '#893949',
'XO': '#b25fdf','YL': '#4553b1','YO': '#9797ff','Z6': '#a28463','ZA': '#abc098','ZN': '#4f7a8e'}
PyGMT_PLT_Scatter_Translator = {'o':'c','x':'x','^':'t','s':'s'}

# n_nets=catalog.Network.unique().shape[0]
# # [display(c) for c in [cm.__dict__[e].resampled(70).resampled(n_nets) for e in ['oslo_categorical','nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]]
# sets=['lajolla_categorical','tokyo_categorical','bamako_categorical','buda_categorical','glasgow_categorical','acton_categorical','lapaz_categorical','devon_categorical','hawaii_categorical','nuuk_categorical','imola_categorical']
# _=[display(c) for c in [mcolors.ListedColormap(cm.__dict__[e].resampled(4).colors.mean(axis=0),name=n) for e,n in zip(sets,catalog.Network.unique())]]
# [{n:mcolors.to_hex(mcolors.ListedColormap(cm.__dict__[e].resampled(4).colors.mean(axis=0),name=n).colors)} for e,n in zip(sets,catalog.Network.unique())]
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# def write_pickle(file,var):
#     import pickle
#     with open(str(file), 'wb') as handle:
#         pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     print('Saved to :' + str(file))
# def load_pickle(file):
#     import pickle
#     with open(file, 'rb') as handle:
#         b = pickle.load(handle)
#     return b

lnm=np.array([[0.10,-162.36,5.64]
,[0.17,-166.70,0.00]
,[0.40,-170.00,-8.30]
,[0.80,-166.40,28.90]
,[1.24,-168.60,52.48]
,[2.40,-159.98,29.81]
,[4.30,-141.10,0.00]
,[5.00,-71.36,-99.77]
,[6.00,-97.26,-66.49]
,[10.00,-132.18,-31.57]
,[12.00,-205.27,36.16]
,[15.60,-37.65,-104.33]
,[21.90,-114.37,-47.10]
,[31.60,-160.58,-16.28]
,[45.00,-187.50,0.00]
,[70.00,-216.47,15.70]
,[101.00,-185.00,0.00]
,[154.00,-168.34,-7.61]
,[328.00,-217.43,11.90]
,[600.00,-258.28,26.60]
,[10000.00,-346.88,48.75]])

hnm=np.array([[0.10,-108.73,-17.23]
,[0.22,-150.34,-80.50]
,[0.32,-122.31,-23.87]
,[0.80,-116.85,32.51]
,[3.80,-108.48,18.08]
,[4.60,-74.66,-32.95]
,[6.30,0.66,-127.18]
,[7.90,-93.37,-22.42]
,[15.40,73.54,-162.98]
,[20.00,-151.52,10.01]
,[354.80,-206.66,31.63]])
