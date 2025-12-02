from pathlib import Path;import shutil,sys,os;
import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent))
# --------------------[IMPORTANT] Paths--------------------
# project_path = Path(os.getcwd()) #Generic path (current working directory)
project_path = Path(__file__).parent
# project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research') #Local machine path
# project_path = Path('/Volumes/ECHO/Research') #External drive path

#This function yields the standard path network called on for everything in this project.
os.chdir(project_path.parent);sys.path.append(str(project_path.parent)) #Do not touch this line. All functions/scripts assume the current working directory is the project path defined above.
sys.path.insert(1,str(project_path.parent))
sys.path.append(str(project_path / 'packages'))
sys.path.insert(0, str(project_path / 'packages' / 'scripts'))
sys.path.insert(0, str(project_path / 'packages' / 'NoiseCut'))
sys.path.insert(0, str(project_path / 'packages' / 'ATaCR'))
sys.path.insert(0, str(project_path / 'packages' / 'ATaCR'/ 'OBStools'))


# --------------------Imports--------------------
from paths import *
from modules import *
io.dir_libraries=dir_libraries
import itertools
from collections.abc import Iterable
import local_tools.dataspace as ds

dirs=io.dir_libraries(mkdir=True)
octavg = lt.math.octave_average


primary_models = [dirs.Data/'SNR_Models'/'SNR_acausul.filter_V04_5s_bandwidth_100_bands.pkl', #SNR organized into a dataframe
dirs.Analysis/'BulkLoad.SR.Coherences_092625.pkl', #Coherence organized into a dataframe
dirs.Data/'SNR_Models'/'BulkHold.pkl', #Coherence and SNR in a single dictionary
]
import zipfile
for fi in primary_models:
    if not (fi).exists():
        file=fi.name;fold = fi.parent
        zip_path = fold/f'{file}.zip';out_dir=zip_path.parent
        if zip_path.exists():
            with zipfile.ZipFile(zip_path, "r") as zf:zf.extractall(out_dir)


# -------------------- Registers a more versatile and intuitive way to --------------------
# -------------------- parse a pandas DataFrame under a generic search ability called .loc.--------------------
# Type whatever you want into .log (e.g. df.loc('TRM') and it will (extremely fast) 
# find the relevant key for you and parse the dataset down to rows with that attribute.
os.system('cls' if os.name == 'nt' else 'clear')
def unravel(nested_list):
    flat_list = []
    for item in nested_list:
        test=(isinstance(item, Iterable)) and (not isinstance(item, str))
        if test:flat_list.extend(unravel(item))
        else:flat_list.append(item)
    return flat_list
def loc_any(df: pd.DataFrame, key):
    """
    Return rows whose MultiIndex contains `key` in any level, without level=.
    Assumes `key` occurs in exactly one index level.
    If key not found, returns an empty DataFrame with same columns and index names.
    """
    idx = df.index
    if not isinstance(idx, pd.MultiIndex):
        # Normal Index
        return df.loc[[key]] if key in idx else df.iloc[0:0]
    # Fast path: already in the first level
    if key in idx.levels[0]:
        return df.sort_index(level=0).loc[key]
    # Search other levels; if found, swap it to front
    for i in range(1, idx.nlevels):
        if key in idx.levels[i]:
            df2 = df.swaplevel(0, i).sort_index(level=0)
            try:return df2.loc[key]
            except:return df.iloc[0:0]
    # Not found anywhere → empty DataFrame
    return df.iloc[0:0]
@pd.api.extensions.register_dataframe_accessor("aloc")
class AnyLoc:
    def __init__(self, pandas_obj):
        self._obj = pandas_obj

    def __getitem__(self, key):
        return loc_any(self._obj, key)


# --------------------Data--------------------
# How to use dataspace

# -- Publication setup
# catalog = ds.dataspace() 
#Without any arguments, the dataspace will assume all data for the entire default catalog is present and ready for us.
#If this is not true, this will return an error.


# -- Demo setup
# Constrains catalog to only one station (demo='7D.M07A' for demonstation purposes) 
# and does not pipe in any data to the catalog (aggregate=False ,ie noise data, coherence, snr, or traces)
catalog = ds.dataspace(aggregate=True,demo='7D.M08A')
# Note: The demo argument can be set to a single station name or a list of station names.




# --------------------Some.useful.constants--------------------
ColorStandard=AttribDict()
ColorStandard.instrument = {'B2': '#e31a1c','KE': '#b2df8a','AB': '#a6cee3','BA': '#cab2d6','AR': '#ff7f00','TRM': '#1f78b4','BG': '#33a02c','BD': '#6a3d9a'}
ColorStandard.seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
ColorStandard.components = {'ZP':'#0c51a6','ZZ':'#2d8297','Z1':'#70cbc0','Z2':'#4f86c5',
'ZZ':'#0c51a6','11':'#2a7e93','22':'#7370cb','PP':'#4f86c5'}
ColorStandard.network = {'2D': '#d2ad90','7A': '#94530c','7D': '#64783c','X9': '#c04797','XF': '#893949',
'XO': '#b25fdf','YL': '#4553b1','YO': '#9797ff','Z6': '#a28463','ZA': '#abc098','ZN': '#4f7a8e'}
PyGMT_PLT_Scatter_Translator = {'o':'c','x':'x','^':'t','s':'s'}