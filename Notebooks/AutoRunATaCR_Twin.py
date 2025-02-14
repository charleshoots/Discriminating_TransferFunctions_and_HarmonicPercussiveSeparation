from imports import *
hps_staquery_output = Path(os.getcwd())/'_DataArchive/HPS_Data/sta_query.pkl'

# ---------------------------------------------------------------------------------------------------
# ============================================ LOAD DATA ===========================================
# ---------------------------------------------------------------------------------------------------
### ===============================================================================
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
# catalog = catalog[catalog.Station.isin(['FN07A','FN14A','J50A','M07A'])]
# # reverse catalog
# catalog = catalog.iloc[list(np.flip([a for a in range(len(catalog))]))]
# catalog = catalog.iloc[np.where(catalog.StaName=='YO.X01')[0][0]:]
# -----------------------------------------------------------------------------------
## =============================================================================== ##
STEPS = [5,6,7]
fork = False;event_mode = False
# cat = catalog.copy()
cat = pd.read_pickle(dirs.Catalogs / 'Catalog_test.pkl')
cat=cat[cat.StaName.isin(['2D.OBS03','2D.OBS06','7D.G17B','7D.G25B','YL.B09W','YO.X10','Z6.16'])]
Minmag,Maxmag=6.0,8.0
cleanspectra_flags = '--figQC --figAverage --figCoh --figCross --save-fig'
dailyspectra_flags='--figQC --figAverage --figCoh --save-fig'
## =============================================================================== ##
# -----------------------------------------------------------------------------------
if event_mode:
    staquery_output = hps_staquery_output
else:
    staquery_output = './sta_query.pkl'
if 1 in STEPS:
    STEPS.pop(np.where(np.array(STEPS)==1)[0][0])

# event_window = 7200
event_window = 3600*4
channels = 'Z,P,12'
ATaCR_Parent = dirs.ATaCR
days=11
ovr=True
for STEP in STEPS:
    for ii,Station in enumerate(cat.iloc):
        print('[//////////////////////////]'*2)
        print('----Station: ' + Station.StaName +  ' (' + str(ii+1) + ' of ' + str(len(cat)) + ')')
        icatalog = Station.to_frame().T
        print('[//////////////////////////]'*2)
        message = f'Station {ii+1}/{len(cat)}'
        ObsQA.TOOLS.io.Run_ATaCR(icatalog,
        fork=fork,
        message=message,
        staquery_output=staquery_output,
        event_mode=event_mode,
        ATaCR_Parent = ATaCR_Parent,
        STEPS=[STEP],log_prefix=Station.StaName,
        Minmag=Minmag,Maxmag=Maxmag,
        event_window=event_window,
        channels=channels,
        cleanspectra_flags=cleanspectra_flags,
        dailyspectra_flags=dailyspectra_flags,
        days=days,
        ovr=ovr)

## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
