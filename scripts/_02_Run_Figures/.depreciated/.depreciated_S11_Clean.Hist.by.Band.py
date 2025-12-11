# Author: Charles Hoots
# This code was developed as part of my PhD research in the
# Department of Earth Sciences, University of Hawai‘i at Mānoa.
# Unless otherwise noted, the code is my own original work.
# External libraries and standard research software packages are used as cited.

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

from scipy.stats import iqr
from local_tools.quick_class import *
from local_tools.math import spectra
# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
import statsmodels.api as sm
import local_tools.dataspace as ds
# Noise Spectra
f=cat.r.iloc[0].Data.Noise.Averaged().f
# faxis value
faxis=(f>0)&(f<=1)
# f value
f=f[faxis]
# noise f value
noise_f=f
cat.r['Noise']=[AttribDict({'f':f,
'Z':PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__['cZZ'][faxis]),
'P':PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__['cPP'][faxis]),
'H':np.mean([PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__[c][faxis]) for c in ['c11','c22']],axis=0)
}) for s in cat.r.iloc]
# function custom cmap
def custom_cmap(ind=0,nbins=5):
    if ind==0:cmap = cm.cmaps['glasgow'].reversed().resampled(nbins)
    if ind==1:cmap = cm.cmaps['batlow'].reversed().resampled(nbins)
    return cmap
# figs value
figs = lambda r=3,c=1,f=(5,6),x='all',y='all':plt.subplots(r,c,figsize=f,sharex='all',sharey='all',layout='constrained')




# opts value
opts=AttribDict({})

# --------------------------------------------------------------------------------
# ------------------------------------CODE----------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
width=4;height=4 #Defaults
# # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# # --------------------------------------------------------------------------------
req = None #Standard setup
# figsize value
figsize=(width,height) #width,height
# sets value
sets = [None]
# transpose value
transpose=False
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
######### # --------------------------------------------------------------------------------
######### # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


# req value
req = 'Instrument_Design'
# figsize value
figsize=(6,1);transpose=True #width,height
# sets value
sets=['AB','AR','B2','BA','BG','KE','TRM']
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
# # ####### --------------------------------------------------------------------------------
req = 'Seismometer'
# figsize value
figsize=(6,1);transpose=True #width,height
# sets value
sets=['Guralp CMG3T 120','Trillium 240','Trillium Compact']
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
# ####### --------------------------------------------------------------------------------
req = 'Magnitude'
# figsize value
figsize=(6,1);transpose=True #width,height
# sets value
sets=[[6.0,7.0],[7.0,8.01]]
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
# ####### --------------------------------------------------------------------------------
# ####### --------------------------------------------------------------------------------
req = 'Pressure_Gauge'
# figsize value
figsize=(6,1);transpose=True #width,height
# sets value
sets=['DPG','APG']
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
density=False;notched=False

# justpairs value
justpairs=True

# octav value
octav=False #Large compute

# weighted value
weighted = False

# cumulative value
cumulative=False;stacked=False
# bins value
bins=np.arange(0,1.1,.1)
# bins value
bins = [0,1]
# ls=(5,(3,1))
# ls=(0,(1,.5))
ls='-'
# lw value
lw=1
# sigma value
sigma = False

# reduced value
reduced = True;run_reduce=True

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
safekeys=['Name', 'StaName', 'Station', 'Event', 'Network', 'LaLo', 'Distance',  'Magnitude', 'Stations',
'Latitude', 'Longitude', 'Experiment', 'Environment', 'Pressure_Gauge','StaDepth', 'Start', 'End', 'NoiseAverage', 'Seismometer',
'Sediment_Thickness_m', 'Instrument_Design', 'Distance_from_Land_km','Distance_to_Plate_Boundary_km', 'Surface_Current_ms',
'Crustal_Age_Myr', 'Deployment_Length_days', 'Inventory', 'Coherence','PearsonCC']
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# SR_std=load_pickle(dirs.Analysis/'N_SR.pkl')
# N_SR=load_pickle(dirs.Analysis/'N_SR.pkl')
# NN_SR_std=load_pickle(dirs.Analysis/'N_SR.pkl')
# NN_SR=load_pickle(dirs.Analysis/'N_SR.pkl')

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Data aggregation--------------------------------------------------------

# f value
f=cat.sr.iloc[0].Data.Coherence().f
# notch lambda value
notch_lambda = lambda z,b,f=f:((f<(fnotch(z) if notched else 1)))&(((1/f)<max(b)) & ((1/f)>min(b)))
bands = [[1,100],[1,10],[10,30],[30,100]];methods = ['TF','HPS_Z','HPS_H'];band_keys=AttribDict()
def reduce_data(N_SR,SR_std,band_keys,f,notched=True,justpairs=True,octav=False):
    flim = lambda b,f: ((1/f)<=max(b))&((1/f)>min(b)) #Assumes b is period
    notch_lambda = lambda z,f:f<fnotch(z) if notched else True
    lims = lambda b,f,z:flim(b,f)&notch_lambda(z,f)
    for si,(s,s_std) in enumerate(zip(N_SR.iloc,SR_std.iloc)):
        print(f'Collecting data: {si+1}/{len(N_SR)}')
        s.Coherence.update({k :s.Coherence[k.split('.')[0]][:,lims(band_keys[k],f,s.StaDepth)] for k in band_keys.keys()})
        s_std.Coherence.update({k :s_std.Coherence[k.split('.')[0]][:,lims(band_keys[k],f,s_std.StaDepth)] for k in band_keys.keys()})
        if octav:
            s.Coherence.update({k:octavg(s.Coherence[k],f[lims(band_keys[k],f,s.StaDepth)])[1] if sum(lims(band_keys[k],f,s.StaDepth))>0 else None for k in band_keys.keys()})
            s_std.Coherence.update({k:octavg(s_std.Coherence[k],f[lims(band_keys[k],f,s_std.StaDepth)])[1] if sum(lims(band_keys[k],f,s_std.StaDepth))>0 else None for k in band_keys.keys()})
        if justpairs:
            s.Coherence.update({k :s.Coherence[k].mean() if not np.any(s.Coherence[k]==None) else None for k in band_keys.keys()})
            # for k in band_keys.keys():s.Coherence[k] = s.Coherence[k].mean() if not np.any(s.Coherence[k]==None) else None
            s_std.Coherence.update({k :s_std.Coherence[k].std() if not np.any(s_std.Coherence[k]==None) else None for k in band_keys.keys()})
    return N_SR,SR_std
for m in methods:
    for b in bands:k=f'{m}.{min(b)}to{max(b)}';band_keys[k] = b
if run_reduce:
    file='SR.pkl'
    N_SR,SR_std,NN_SR,NN_SR_std = load_pickle(dirs.Analysis/file),load_pickle(dirs.Analysis/file),load_pickle(dirs.Analysis/file),load_pickle(dirs.Analysis/file)
    N_SR,SR_std = reduce_data(N_SR,SR_std,band_keys,f,notched=True,justpairs=justpairs,octav=octav)
    NN_SR,NN_SR_std = reduce_data(NN_SR,NN_SR_std,band_keys,f,notched=False,justpairs=justpairs,octav=octav)
else:
    if sigma:file=file.replace('.pkl','_std.pkl')
    if reduced:file=file.replace('.pkl','_Reduced.pkl')
    NN_SR = load_pickle(dirs.Analysis/file)
    N_SR = load_pickle(dirs.Analysis/('Notched_'+file))
# # bulk_cat_version='50125'
# # if not notched:
# #     N_SR=load_pickle(dirs.P01.S10/bulk_cat_version/'SR_NoNotch.pkl')
# #     SR_std=load_pickle(dirs.P01.S10/bulk_cat_version/'SR_std_NoNotch.pkl')
# # else:
# #     N_SR=load_pickle(dirs.P01.S10/bulk_cat_version/'N_SR.pkl')
# #     SR_std=load_pickle(dirs.P01.S10/bulk_cat_version/'SR_std.pkl')

# --------------------------------------------------------

status = lambda:print(f'JustPairs:{justpairs}|Octav:{octav}|Density:{density}|Notched:{notched}|{'Summary' if req==None else req} | {'Plotting' if req==None else set}')

# combinations = list(itertools.product([True, False], repeat=2))
# for combi,(density,notched) in enumerate(combinations):


for req in opts.keys():
    figsize,sets,transpose=opts[req]['figsize'],opts[req]['sets'],opts[req]['transpose']
    for set in sets:
        status()
        if req==None:title='All data'
        else:
            interval=f'{min(set)} to {max(set)}' if isinstance(set,list) else f'{set}'
            title=f'{req.replace('_',' ')}, {interval}'
        # DATA SUBSET-----------
        if req==None:
            # All data ----------------
            iSR=N_SR.copy();NN_iSR=NN_SR.copy()
            weighting=N_SR.copy()
        else:
            # Data subsets ----------------
            if isinstance(set,str):partition=N_SR[req]==set
            else:partition=(N_SR[req]>=min(set))&(N_SR[req]<=max(set))
            iSR=N_SR[partition].copy()
            NN_iSR=NN_SR[partition].copy()
            weighting=N_SR[partition].copy()

        # DATA WEIGHTS-----------
        si=0;weights = {}
        weights_nevents_per_sta = np.array([sum(weighting.StaName==s) for s in weighting.StaName])/len(weighting) #[0-1]#Number of events per station relative to the entire bin size.
        weights_sr_per_bin = np.ones(len(weighting))*len(weighting)/len(iSR) #[0-1]#Bin size fraction of the entire data set. weights_sr_per_bin is effectively a constant in each bin and has no effect.
        weights_nfrequencies_per_sta = np.array([sum(f<fnotch(s.StaDepth)) for s in weighting.iloc])/sum(f<fnotch(weighting.StaDepth.min())) #[0-1]#Fraction of frequencies used in this bin contributed to the average from each pair.
        iweights = weights_nevents_per_sta + weights_sr_per_bin + weights_nfrequencies_per_sta
        weights.update({si:iweights})
        if weighted:wgt = [weights[k] for k in [si]]
        else:wgt = [weights[k]*0+1 for k in [si]]
        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        fig,axes=figs(r=4,f=(4,6))

        bands = [[1,100],[1,10],[10,30],[30,100]]
        special_colors = {'tf':cm.cmaps['bilbao'](1/6),'area':cm.cmaps['oslo'](5/6),'hpsh':cm.cmaps['bamako'](3/6),'hpsz':cm.cmaps['imola'](1/6)}
        if not justpairs:weighted=False
        for bi,b in enumerate(bands):
            ax = axes[bi];b=bands[bi]
            k=f'{'TF'}.{min(b)}to{max(b)}'
            for di,D in enumerate([NN_iSR,iSR]):
                if di==0:alpha=0.4
                else:alpha=1.0
                if not notched:alpha=1.0
                keeps = np.atleast_2d(np.array([len(np.atleast_2d(s.Coherence[k]).reshape(-1)) for s in D.iloc])>0)
                if di==0:
                    if justpairs:
                        m='TF';k=f'{m}.{min(b)}to{max(b)}';tf = np.atleast_2d(np.array([s.Coherence[k] for s in D.iloc]).squeeze())
                        m='HPS_Z';k=f'{m}.{min(b)}to{max(b)}';hpsz = np.atleast_2d(np.array([s.Coherence[k] for s in D.iloc]).squeeze())[keeps]
                        m='HPS_H';k=f'{m}.{min(b)}to{max(b)}';hpsh = np.atleast_2d(np.array([s.Coherence[k] for s in D.iloc]).squeeze())[keeps]
                    else:
                        m='TF';k=f'{m}.{min(b)}to{max(b)}';tf = np.atleast_2d(np.array([s.Coherence[k] for s in D.iloc]).squeeze())[keeps[0,:],:]
                        m='HPS_Z';k=f'{m}.{min(b)}to{max(b)}';hpsz = np.atleast_2d(np.array([s.Coherence[k] for s in D.iloc]).squeeze())[keeps[0,:],:]
                        m='HPS_H';k=f'{m}.{min(b)}to{max(b)}';hpsh = np.atleast_2d(np.array([s.Coherence[k] for s in D.iloc]).squeeze())[keeps[0,:],:]
                else:
                    m='TF';k=f'{m}.{min(b)}to{max(b)}';tf = np.hstack([np.atleast_2d(s.Coherence[k]).reshape(-1) for s in D.iloc])
                    m='HPS_Z';k=f'{m}.{min(b)}to{max(b)}';hpsz = np.hstack([np.atleast_2d(s.Coherence[k]).reshape(-1) for s in D.iloc])
                    m='HPS_H';k=f'{m}.{min(b)}to{max(b)}';hpsh = np.hstack([np.atleast_2d(s.Coherence[k]).reshape(-1) for s in D.iloc])

                if weighted:
                    if justpairs:wgt=weights[si][keeps] #*0+1
                    else:wgt = np.atleast_2d(tf) *0 + np.atleast_2d(weights[si][keeps]).T
                y=[tf,hpsz,hpsh];color=[special_colors['tf'],special_colors['hpsz'],special_colors['hpsh']]
                if ((di>0)&notched) or ((di==0)&(not notched)):[ax.axvline(i.mean(),c=c,linewidth=5,alpha=0.7,zorder=-100,ls=ls) for i,c in zip(y,color)]
                for yy,cc in zip(y,color):
                    if not weighted:wgt=yy*0+1
                    yy,wgt = yy.reshape(-1),wgt.reshape(-1)
                    if (di>0)&(not notched):_=ax.hist([yy],bins=bins,weights=[wgt],cumulative=cumulative,density=density,stacked=stacked,histtype='bar',alpha=alpha,linewidth=1.5,zorder=1000,facecolor='w',edgecolor='k')


            ax2=ax.twinx();ax2.set_ylabel(f'{min(b)} to {max(b)}s', color='k',rotation=-90,labelpad=10,fontweight='bold');ax2.set_yticks([])
            if density:ax.set_yticklabels(np.round(ax.get_yticks()*np.diff(bins)[0],2))
            else:ax.set_yscale('log')
        fig.suptitle(title,y=1.1)
        if req==None:file=f'Summary.Bands.{'by_SR' if justpairs else 'by_Fq'}.png'
        else:file=f'{req}.{set if isinstance(set,str) else f'{min(set)}to{max(set)}'}.Bands.{'by_SR' if justpairs else 'by_Fq'}.png'
        fold=dirs.P01.S10/'Sets'
        if octav:fold=fold/'Octav'
        else:fold=fold/'NotOctav'
        if notched:fold=fold/'Notched'
        else:fold=fold/'NoNotch'
        # file=file.replace('.png',('.Weighted' if weighted else '.Unweighted')+'.png')
        file =('PDF.' if density else 'Counts.')+file
        if not req==None:fold=fold/req
        else:fold=fold/'_Summary'
        fold.mkdir(parents=True, exist_ok=True)
        save_tight(fold/file,fig,dpi=800)
        plt.close('all')
        ki=1

ki=1

# Probability Mass Histogram
# Binned Probability Histogram
# Normalized Histogram with Unit Area
# Derivative of the CDF