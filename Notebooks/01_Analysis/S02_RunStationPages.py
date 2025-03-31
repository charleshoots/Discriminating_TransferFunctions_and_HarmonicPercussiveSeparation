import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
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


def sta_metrics(report,sta,
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]],**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = {'Z1':'#0c51a6','ZP':'#2a7e93','ZZ':'#7370cb','Z2':'#4f86c5'}
    note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    # else:note = '\nVariance (shaded)'
    stanm = sta.StaName
    tf=''
    stastr = ' | '.join([stanm+tf,sta.Experiment,'Depth: '+str(int(abs(sta.StaDepth)))+'m, Notch: '+str(int(1/fnotch(1000*abs(sta.StaDepth))))+'s',note])
    nrows = 1;ncols=2
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    fig,axes = plt.subplots(nrows=ncols,ncols=nrows,layout='constrained',sharex='all',squeeze=True,figsize=(8,7))
    axes = axes.reshape(-1)
    fig.suptitle(stastr)
    fn = 1/fnotch(sta.StaDepth)
    methods=['ATaCR','NoiseCut']
    for bi,b in enumerate(columns):
        method = methods[bi]
        metric = b[1][0]
        pair = b[1][1]
        ax = axes[bi]
        for tri in range(len([1])):ax.set_title(f'{method} {pair} {metric}')
        x = report.f;ind=x<=1;x=x[ind]
        y=report[method][stanm].coh[:,ind]
        [ax.scatter(x,yy,label=':'.join([pair]),s=5,facecolor=args.csd_pairs[pair],alpha=0.05) for yy in y]
        [ax.plot(x,yy,label=':'.join([pair]),linestyle=':',linewidth=0.05,alpha=0.05,color=args.csd_pairs[pair]) for yy in y]
        ax.set_xlim(1/500,1);ax.set_xscale('log')
        if metric.lower()=='phase':ax.set_ylim(0,180)
        if metric.lower()=='coherence':ax.set_ylim(0,1.01)
        y_avg = np.mean(y,axis=0)
        y_var = np.std(y,axis=0)**2
        ax.scatter(x,y_avg,label=':'.join([pair]),s=0.5,color=args.csd_pairs[pair])
        ax.plot(x,y_avg,label=':'.join([pair]),linewidth=0.8,alpha=0.8,color=args.csd_pairs[pair])
        ax.fill_between(x,y_avg-y_var,y_avg+y_var, alpha=0.2,color=args.csd_pairs[pair])
        ax.set_xlabel('frequency (Hz)')
        # else:ax.set_xlabel('seconds')
        ax.axvline(1/fn,alpha=0.4,linewidth=1,color='k',linestyle='-.')
        ax.text(1/fn,1,str(int(np.round(fn)))+'s',verticalalignment='top',horizontalalignment='right')
    return fig



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
NoiseColors = [mcolors.to_hex(m) for m in [cm.Categorical.__dict__[e].resampled(30).resampled(6).colors for e in ['devonS']][0]]
# [display(c) for c in [cm.__dict__[e].resampled(70).resampled(5) for e in ['nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]]
# np.array([[mcolors.to_hex(c) for c in r] for r in [cm.__dict__[e].resampled(70).resampled(5).colors for e in ['nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]])
print('Stations: ' + str(len(catalog)))
# from local_tools import *
import local_tools as lt
from local_tools.cat import shared_events
from local_tools.io import get_station_events_hps,get_station_events
from local_tools.plots import station_event_page_averages,station_event_page
# # ======================================================================================================================================================


get_reports

# type='metrics'# type='stream'
plotfolder = dirs.P01.S02
hpsfold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
atacrfold = dirs.Events
averages = True
ovr = True
columns =['Coherence']
minsta=10
cat = catalog.copy()
cat.sort_values(by='StaDepth',inplace=True)
event_catalog = lt.cat.unravel_cat(cat)
# types = ['metrics','stream']
for type in ['metrics']:
# for type in ['stream','metrics']:
    # csd_pairs=[('ZP','#0c51a6')]
    # csd_pairs=[('ZZ','#2a7e93')]
    # csd_pairs=[('Z1','#7370cb')]
    # csd_pairs=[('Z2','#4f86c5')]
    # pairs = [[('ZP','#0c51a6')],
    #         [('ZZ','#2a7e93')],
    #         [('Z1','#7370cb')],
    #         [('Z2','#4f86c5')]]
    pairs = [['ZZ','#2a7e93']]    
    if type=='stream':pairs = [pairs[0]];runs=[0,1]
    else: runs=[99]
    for stai,sta in enumerate(cat.iloc):
        for csd_pairs in pairs:
            report=get_reports(csd_pairs[0],cat,dirs.Archive,dirs,AVG=False,methods=['ATaCR','NoiseCut'])
            for run in runs:
                if run==0:method='ATaCR';evdir = atacrfold;tf='sta.ZP-21.HZ.SAC';mirror_fold = hpsfold
                elif run==1:method = 'Noisecut';tf='HZ.SAC';evdir = [hpsfold,dirs.Events];mirror_fold = atacrfold
                else:method='Both'

                if (type.lower()=='metrics') & (averages==True):
                    file = sta.StaName+'.'+method.lower()+'.'+columns[0]+'.png'
                    folder = plotfolder / ('Coherences' if type.lower()=='metrics' else 'Events')
                else:
                    folder = plotfolder / 'Events'
                    file = sta.StaName+'.'+note+'.'+method.lower()+'.png'
                if ((folder / file).exists()) & (not ovr):print('File exists. Skipping.');continue


                # clear_output(wait=False);os.system('clear')
                print('[*]'*30)
                print('|'.join([f'{method}  | Overall:{np.round((1/len(runs))*100*((stai+1)/len(cat)),1)}% |{str(stai+1)}/{str(len(cat))}',sta.StaName]))
                print('-Load data')
                # if type=='stream':
                #     if method.lower()=='noisecut':
                #         evmeta = shared_events(event_catalog,dirs,minsta=minsta,stanm=sta.StaName)
                #         st_hold,evmeta = get_station_events_hps(sta.StaName,evdir,tf=tf,type=type,evmeta=evmeta)
                #     else:
                #         evmeta = shared_events(event_catalog,dirs,minsta=minsta,stanm=sta.StaName)
                #         st_hold,evmeta = get_station_events(sta.StaName,evdir,tf=tf,type=type,evmeta=evmeta)
                #         raw_reference = st_hold.select(location='*Raw*').copy()
                # else:
                #         evmeta = shared_events(event_catalog,dirs,minsta=minsta,stanm=sta.StaName)
                #         tf='HZ.SAC';evdir = [hpsfold,dirs.Events];mirror_fold = atacrfold
                #         st_hold_noisecut,evmeta = get_station_events_hps(sta.StaName,evdir,tf=tf,type=type,evmeta=evmeta)
                #         evdir = atacrfold;tf='sta.ZP-21.HZ.SAC';mirror_fold = hpsfold
                #         st_hold_atacr,evmeta = get_station_events(sta.StaName,evdir,tf=tf,type=type,evmeta=evmeta)
                #         st_hold = st_hold_atacr.select(location='*Raw*').copy() + st_hold_atacr.select(location='*ATaCR*').copy() + st_hold_noisecut.select(location='*NoiseCut*').copy()
                #         raw_reference = st_hold_atacr.select(location='*Raw*').copy()
                #         del st_hold_noisecut,st_hold_atacr

                # clear_output(wait=False);os.system('clear')
                print('-Plot data')
                # print('|'.join([f'{method} ({run+1}/2) | Overall:{np.round(100*((stai+1)/len(cat)),1)}% |{str(stai+1)}/{str(len(cat))}',sta.StaName]))
                print('[*]'*30)
                print('Load complete')
                # fold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_FigureArchive/_GEN6/StationEventPages/low_rez')
                # np.array([len(list((fold/c.StaName).glob('*.png')))/10 for c in catalog.iloc])
                # np.array([len(list((fold/c.StaName).glob('*events*.png')))/2 for c in catalog.iloc]) #EVENTS=DONE
                bands = [(1,10),(10,30),(30,100)]
                note=''
                if (type.lower()=='metrics') & (averages==True):
                    fig = sta_metrics(report,sta,csd_pairs=csd_pairs,columns=columns)
                else:
                    fig = station_event_page(st_hold,sta,evmeta,method,type=type,csd_pairs=csd_pairs)
                if type.lower()=='stream':subfold='Events'
                else:subfold=csd_pairs[0]
                note = type.replace('metrics','coph').replace('stream','events')
                if type.lower()=='metrics':note=csd_pairs[0]+note
                save_tight(folder / file,fig,dpi=200)
                plt.close('all')
                # del st_hold

lt.cat.banner('S02 COMPLETE!')


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX