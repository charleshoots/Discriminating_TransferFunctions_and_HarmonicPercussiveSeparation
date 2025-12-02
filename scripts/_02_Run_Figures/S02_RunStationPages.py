### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

from local_tools.quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from obspy import read
from obspy.geodetics import locations2degrees
# ========================================================================================================================================================
from IPython.display import clear_output

# instrument colors value
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
# variable
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
# seismometer marker value
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}

# dirs value
dirs = io.dir_libraries()

# write pickle
def write_pickle(file,var):
    import pickle
    with open(str(file), 'wb') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Saved to :' + str(file))
# load pickle
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        # b value
        b = pickle.load(handle)
    return b
# function mirror events
def mirror_events(reports):
    # nkeys value
    nkeys = [n for n in list(reports[0].__dict__.keys()) if not n=='f']
    # mirror value
    mirror = dict()
    for ni,n in enumerate(nkeys):
        # skeys value
        skeys = list(reports[0][n].__dict__.keys())
        for si,s in enumerate(skeys):
            # stanm value
            stanm = '.'.join([n,s]).replace('n','')
            # ev0 value
            ev0 = [k.replace('.','') for k in reports[0][n][s].events]
            # ev1 value
            ev1 = [k.replace('.','') for k in reports[1][n][s].events]
            mirror[stanm] = np.intersect1d(ev0,ev1)
    return mirror
# ================================================================================================================================


# function sta metrics
def sta_metrics(srcat,sta,
    # columns value
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]],**args):
    # defargs value
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    # args value
    args = defargs
    args.csd_pairs = {'Z1':'#0c51a6','ZP':'#2a7e93','ZZ':'#7370cb','Z2':'#4f86c5'}
    # note value
    note=''
    # note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+') \nVariance (shaded)'
    # else:
    note = '\nVariance (shaded)'
    # stanm value
    stanm = sta.StaName
    # tf value
    tf=''
    # stastr value
    stastr = ' | '.join([stanm+tf,sta.Experiment,'Depth: '+str(int(abs(sta.StaDepth)))+'m, Notch: '+str(int(1/fnotch(1000*abs(sta.StaDepth))))+'s',note])
    # nrows value
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

        if method=='ATaCR':y=np.array([sr.Data.Coherence()[method].zp_21.coh for sr in srcat.iloc]).squeeze()[:,ind]
        else:y=np.array([sr.Data.Coherence()['HPS'][b[1][1].lower()].coh for sr in srcat.iloc]).squeeze()[:,ind]
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
reportfolder = dirs.Data/'Analysis'/'NetworkCoherences'
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
print('Stations: ' + str(len(catalog.r)))
# from local_tools import *
import local_tools as lt
from local_tools.cat import shared_events
from local_tools.io import get_station_events_hps,get_station_events
from local_tools.plots import station_event_page_averages,station_event_page
# # ======================================================================================================================================================



# type='metrics'# type='stream'
plotfolder = dirs.P01.S02
hpsfold = dirs.Data/'HPS_Data'/'Data'
atacrfold = dirs.Events
averages = True
ovr = True
columns =['Coherence']

cat = catalog.r.copy()
cat.sort_values(by='StaDepth',inplace=True)

# cat = cat.loc[['ZA.B02','Z6.16','Z6.11','YO.X10','7D.G34D','7D.G26D','7D.G25B','7D.G17B']]

event_catalog = lt.cat.unravel_cat(cat)
# types = ['metrics','stream']
for type in ['metrics']:
    pairs = [['ZZ','#2a7e93']]    
    if type=='stream':pairs = [pairs[0]];runs=[0,1]
    else: runs=[99]
    for stai,sta in enumerate(cat.iloc):
        srcat = catalog.sr[catalog.sr.StaName==sta.StaName].copy()
        for csd_pairs in pairs:
            report=get_reports(methods=['ATaCR','NoiseCut'])
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
                print('[*]'*30)
                print('|'.join([f'{method}  | Overall:{np.round((1/len(runs))*100*((stai+1)/len(cat)),1)}% |{str(stai+1)}/{str(len(cat))}',sta.StaName]))
                print('-Load data')
                print('-Plot data')
                print('[*]'*30)
                print('Load complete')
                bands = [(1,10),(10,30),(30,100)]
                note=''
                if (type.lower()=='metrics') & (averages==True):
                    fig = sta_metrics(srcat,sta,csd_pairs=csd_pairs,columns=columns)
                else:
                    fig = station_event_page(st_hold,sta,evmeta,method,type=type,csd_pairs=csd_pairs)
                if type.lower()=='stream':subfold='Events'
                else:subfold=csd_pairs[0]
                note = type.replace('metrics','coph').replace('stream','events')
                if type.lower()=='metrics':note=csd_pairs[0]+note
                folder.mkdir(exist_ok=True,parents=True)
                save_tight(folder / file,fig,dpi=200)
                plt.close('all')

lt.cat.banner('S02 COMPLETE!')
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX