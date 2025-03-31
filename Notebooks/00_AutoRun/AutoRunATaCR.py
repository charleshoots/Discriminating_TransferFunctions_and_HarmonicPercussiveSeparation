import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
from imports import *
unpacker=locals().update;unpack=lambda Args,keys=None:[unpacker({k:Args[k]}) for k in [keys if keys is not None else list(Args.keys())][0]]

hps_staquery_output = Path(os.getcwd())/'_DataArchive/HPS_Data/sta_query.pkl'
Args=AttribDict()




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
# cat = pd.read_pickle(dirs.Catalogs / 'Catalog_Test_DensityIncreased.ShalllowIncreased.pkl')


# cat = pd.read_pickle(dirs.Catalogs / 'Catalog_022325.pkl')
# events = lt.cat.unravel_cat(cat[cat.Network=='7D'])
# cat = cat[cat.StaName=='7D.J25A'];stanm=cat.iloc[0].StaName
# cat.Events[0] = Catalog([e for e in events if len(list((dirs.Events/'corrected'/stanm).glob(f'*{e.Name}*.SAC')))==4])
# cat=cat[cat.StaName.isin(['2D.OBS03','2D.OBS06','7D.G17B','7D.G25B','YL.B09W','YO.X10','Z6.16'])]
# cat=cat[cat.StaName.isin(['2D.OBS23'])]
# cat=catalog.copy()

Args.Minmag,Args.Maxmag=6.0,8.0
Args.cleanspectra_flags = '--figQC --figAverage --figCoh --figCross --save-fig'
Args.dailyspectra_flags='--figQC --figAverage --figCoh --save-fig'
## =============================================================================== ##
# -----------------------------------------------------------------------------------

Args.STEPS = [4,5,6,7]
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
if Args.event_mode:Args.staquery_output = hps_staquery_output
else:Args.staquery_output = './sta_query.pkl'
if 1 in Args.STEPS:Args.STEPS.pop(np.where(np.array(Args.STEPS)==1)[0][0])
for STEP in Args.STEPS:
    Args.STEP = [STEP]
    for ii,Station in enumerate(cat.iloc):
        Args.ovr=False
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

## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##


