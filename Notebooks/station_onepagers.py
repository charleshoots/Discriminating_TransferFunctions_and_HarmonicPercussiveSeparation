from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from obspy import read
from obspy.geodetics import locations2degrees
# ========================================================================================================================================================
from IPython.display import clear_output

instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}

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
def mirror_events(reports):
    nkeys = [n for n in list(reports[0].__dict__.keys()) if not n=='f']
    mirror = dict()
    for ni,n in enumerate(nkeys):
        skeys = list(reports[0][n].__dict__.keys())
        for si,s in enumerate(skeys):
            stanm = '.'.join([n,s]).replace('n','')
            ev0 = [k.replace('.','') for k in reports[0][n][s].events]
            ev1 = [k.replace('.','') for k in reports[1][n][s].events]
            mirror[stanm] = np.intersect1d(ev0,ev1)
    return mirror
# ================================================================================================================================
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
reportfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences')
bands = ['1-10','10-30','30-100']
# ================================================================================================================================
# ================================================================CODE SNIPPETS===========================================================================
# for a in cm._cmap_names_categorical:
#     display(cm.__dict__[a].resampled(4))
# for (Event,Station,Metrics,Comp) in OBS_Generator(catalog,dirs['Py_DataParentFolder']):
#     print(Event)
# for i,(Event,Station,Metrics,Comp) in zip(range(1),OBS_Generator(catalog,dirs['Py_DataParentFolder'])):
#     print(Event)
def smooth(d,k=10):return np.convolve(d, np.ones(k) / k, mode='same')
NoiseColors = [mcolors.to_hex(m) for m in [cm.__dict__[e].resampled(30).resampled(6).colors for e in ['devon_categorical']][0]]
# [display(c) for c in [cm.__dict__[e].resampled(70).resampled(5) for e in ['nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]]
# np.array([[mcolors.to_hex(c) for c in r] for r in [cm.__dict__[e].resampled(70).resampled(5).colors for e in ['nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]])
print('Stations: ' + str(len(catalog)))
# from local_tools import *
import local_tools as lt
from local_tools.cat import shared_events
from local_tools.io import get_station_events_hps,get_station_events
from local_tools.plots import station_event_page_averages,station_event_page
# # ======================================================================================================================================================



# type='metrics'# type='stream'
plotfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_FigureArchive/_GEN6')
hpsfold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
atacrfold = dirs.Events
averages = True
ovr = True
columns =['Coherence']
minsta=10
cat = catalog.copy()
event_catalog = lt.cat.unravel_cat(cat)
for type in ['metrics','stream']:
# for type in ['stream','metrics']:
    # csd_pairs=[('ZP','#0c51a6')]
    # csd_pairs=[('ZZ','#2a7e93')]
    # csd_pairs=[('Z1','#7370cb')]
    # csd_pairs=[('Z2','#4f86c5')]
    pairs = [[('ZP','#0c51a6')],
            [('ZZ','#2a7e93')],
            [('Z1','#7370cb')],
            [('Z2','#4f86c5')]]
    
    if type=='stream':pairs = [pairs[0]]
    for stai,sta in enumerate(cat.iloc):
        for csd_pairs in pairs:
            for run in [0,1]:
                if run==0:method='ATaCR';evdir = atacrfold;tf='sta.ZP-21.HZ.SAC';mirror_fold = hpsfold
                elif run==1:method = 'Noisecut';tf='HZ.SAC';evdir = [hpsfold,dirs.Events];mirror_fold = atacrfold
                
                if (method=='Noisecut') and (type=='metrics'):continue

                # clear_output(wait=False);os.system('clear')
                print('[*]'*30)
                print('|'.join([str([i for i in range(len(cat))][stai]+1),'/',str(len(catalog)),'|'+sta.StaName]))
                print('[*]'*30)
                if method.lower()=='noisecut':
                    # evmeta = shared_events(event_catalog,dirs,minsta=minsta)
                    # evmeta = sta.Events
                    evmeta = shared_events(event_catalog,dirs,minsta=minsta,stanm=sta.StaName)
                    st_hold,evmeta = get_station_events_hps(sta.StaName,evdir,tf=tf,mirror_fold=mirror_fold,type=type,evmeta=evmeta)
                else:
                    # evmeta = shared_events(event_catalog,dirs,minsta=minsta)
                    # evmeta = sta.Events
                    evmeta = shared_events(event_catalog,dirs,minsta=minsta,stanm=sta.StaName)
                    st_hold,evmeta = get_station_events(sta.StaName,evdir,tf=tf,mirror_fold=mirror_fold,type=type,evmeta=evmeta)
                    raw_reference = st_hold.select(location='*Raw*').copy()

                # clear_output(wait=False);os.system('clear')
                print('[*]'*30)
                print('|'.join([str([i for i in range(len(cat))][stai]+1),'/',str(len(catalog)),'|'+sta.StaName]))
                print('[*]'*30)
                print('Load complete')
                fold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_FigureArchive/_GEN6/StationEventPages/low_rez')
                # np.array([len(list((fold/c.StaName).glob('*.png')))/10 for c in catalog.iloc])
                # np.array([len(list((fold/c.StaName).glob('*events*.png')))/2 for c in catalog.iloc]) #EVENTS=DONE
                bands = [(1,10),(10,30),(30,100)]
                note=''
                folder = plotfolder / 'StationEventPages' / 'low_rez' / sta.StaName;folder.mkdir(parents=True,exist_ok=True)
                if (type.lower()=='metrics') & (averages==True):
                    status = [float(len(list((fold/c.StaName).glob('*AVG*.png')))/2) for c in catalog.iloc]
                    print('Status |'+str(np.round(100*sum(status)/len(status),2))+'%|'+'--'*37)
                    print(np.array(status));print('--'*40);print('Begin plots')
                    file = sta.StaName+'.AVG.'+columns[0]+'.'+method.lower()+'.png'
                    if ((folder / file).exists()) & (not ovr):print('File exists. Skipping.');continue
                    fig = station_event_page_averages(st_hold,sta,evmeta,method,type=type,csd_pairs=csd_pairs,columns=columns,raw_reference=raw_reference)
                else:
                    status = [float(len(list((fold/c.StaName).glob('*coph*.png')))/8) for c in catalog.iloc]
                    print('Status |'+str(np.round(100*sum(status)/len(status),2))+'%|'+'--'*37)
                    print(np.array(status));print('--'*40);print('Begin plots')
                    file = sta.StaName+'.'+note+'.'+method.lower()+'.png'
                    if ((folder / file).exists()) & (not ovr):print('File exists. Skipping.');continue
                    fig = station_event_page(st_hold,sta,evmeta,method,type=type,csd_pairs=csd_pairs)
                if type.lower()=='stream':subfold='Events'
                else:subfold=csd_pairs[0][0]
                note = type.replace('metrics','coph').replace('stream','events')
                if type.lower()=='metrics':note=csd_pairs[0][0]+note
                # folder = plotfolder / 'StationEventPages' / 'low_rez' / sta.StaName;folder.mkdir(parents=True,exist_ok=True)
                save_tight(folder / file,fig,dpi=200)
                plt.close('all')
                k = 1
                k = 1
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX