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
from local_tools.math import cohstats
from scipy.ndimage import gaussian_filter
from local_tools.plots import ax_sta_metrics
import matplotlib.colors as mcolors
from local_tools.quick_class import *
from mne_connectivity import spectral_connectivity_time
# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
from matplotlib.ticker import PercentFormatter
import time;start=time.time()
# runtime value
runtime=lambda:int(time.time()-start)










# -------------------------------------------------------------------------------------------------------------------------------------
SR = cat.sr.copy()
# -------------------------------------------------------------------------------------------------------------------------------------
meta_wins=AttribDict()
# meta_wins.Sediment_Thickness_m=[[i,i+500] for i in np.arange(0,7320+500,500)]
# meta_wins.Magnitude=[[i,i+.5] for i in np.arange(6,8,.5)]
# meta_wins.Network=['2D','7A','7D','X9','XF','YL','YO','ZA','ZN']
# meta_wins.Seismometer=['Guralp CMG3T 120', 'Trillium 240', 'Trillium Compact']
# meta_wins.Instrument_Design=['AB', 'AR', 'B2', 'BA', 'BG', 'KE', 'TRM']
# meta_wins.Environment=['North Atlantic', 'North Pacific', 'Solomon Sea', 'South Pacific']
# meta_wins.Pressure_Gauge=['APG', 'DPG']
# meta_wins.StaDepth=np.array([[i,i+500] for i in np.arange(0,6000,500)])
meta_wins.Band=np.array([[1,10],[10,30],[30,300]])
# -------------------------------------------------------------------------------------------------------------------------------------
# ---
cumulative=True;density=True;outtype='CDF' #Makes a CDF
# cumulative=False;density=True;outtype='PDF' #Makes a PDF
# cumulative=False;density=False;outtype='Hist' #Makes a Histogram
# ---
notched=True #Band limit to only periods sensitive to infragravity surface waves.
# ---
octav=True #Octave averaging. Adds 20s to run time for each category/figure.
# ---
# If true, justpairs reduces coherence values given to the histogram from single frequency measurements for each source-receiver 
# to a single average for each source-receiver. Recommended to keep this at True.
# ---
justpairs=True #Recommended to always keep this at True.
# ---
sigma=True
# ---
relative=False
# ---
stacked=False
# ---
orientation='horizontal' #colorbar
# orientation='vertical' #colorbar
weighted=True #
# ---
norm_pdf=False
# figsize value
figsize=[6,3.15]
if sigma:justpairs=False
# -------------------------------------------------------------------------------------------------------------------------------------


# f value
f = cat.r.Data[0].Coherence().f
# foct value
foct=octavg(cat.sr.Data[0].Coherence().ATaCR.zp_21.coh,f)[0]
# SR value
SR = cat.sr.__deepcopy__()
# SR std value
SR_std=ds.dataspace().sr.copy() #independent load of data for std analysis. Adds 300mb and 20s to the compute, but it's worth the headache.

# SR std value
SR_std=SR_std[SR_std.Magnitude<7.0].copy()
# SR value
SR=SR[SR.Magnitude<7.0].copy();note='__m6_m7'

# SR_std = SR_std[SR_std.Magnitude>=7.0].copy()
# SR = SR[SR.Magnitude>=7.0].copy();note='__m7_m8'

for si,(s,s_std) in enumerate(zip(SR.iloc,SR_std.iloc)):
    print(f'Collecting data: {si+1}/{len(SR)}')
    s.Coherence.update({k :s.Coherence[k][:,f<fnotch(s.StaDepth) if notched else True] for k in ['TF','HPS_Z','HPS_H']})
    s_std.Coherence.update({k :s_std.Coherence[k][:,f<fnotch(s_std.StaDepth) if notched else True] for k in ['TF','HPS_Z','HPS_H']})
    if octav:
        s.Coherence.update({k:octavg(s.Coherence[k],f[f<fnotch(s.StaDepth) if notched else True])[1] for k in ['TF','HPS_Z','HPS_H']})
        s_std.Coherence.update({k:octavg(s_std.Coherence[k],f[f<fnotch(s_std.StaDepth) if notched else True])[1] for k in ['TF','HPS_Z','HPS_H']})
    if justpairs:
        s_std.Coherence.update({k :s_std.Coherence[k].std() for k in ['TF','HPS_Z','HPS_H']})
        s.Coherence.update({k :s.Coherence[k].mean() for k in ['TF','HPS_Z','HPS_H']})

# state value
state=lambda:f'C{int(cumulative)}.D{int(density)}.S{int(sigma)}.T{int(stacked)}.R{int(relative)}.N{int(norm_pdf)}'

# combinations = list(itertools.product([True, False], repeat=5)) #4**2=16 combinations
# for combi,(cumulative,density,sigma,relative,norm_pdf) in enumerate(combinations):
# print(f'Settings----{combi+1}/{len(combinations)}----')
print(state())
print(f"Elapsed time: {(runtime())/60:.2f} minutes")
print('---'*30)
# cumulative=True;density=True;outtype='CDF' #Makes a CDF
# cumulative=False;density=True;outtype='PDF' #Makes a PDF
# cumulative=False;density=False;outtype='Hist' #Makes a Histogram




if (cumulative) & (density):outtype='CDF'
if (not cumulative) & (density):outtype='PDF'
if (not cumulative) & (not density):outtype='Hist'

# xf value
xf = foct if octav else f
for key in meta_wins.keys():
    # sets value
    sets=meta_wins[key]
    # Y TF value
    Y_TF={};Y_HPSZ={};Y_HPSH={};weights={}
    for si,set in enumerate(sets):
        print(f'Slicing: [{key}] {si+1}/{len(sets)}')
        if sigma:iSR=SR_std.copy()
        else:iSR=SR.copy()
        if not key=='Band':
            if isinstance(set,list)or(type(set)==type(np.array([]))):ind=(iSR[key]>=min(set))&(iSR[key]<max(set))
            else:ind=iSR[key]==set
            # icat value
            icat=iSR[ind].copy()
            # weights nevents per sta value
            weights_nevents_per_sta = np.array([sum(icat.StaName==s) for s in icat.StaName])/len(icat) #[0-1]#Number of events per station relative to the entire bin size.
            # weights sr per bin value
            weights_sr_per_bin = np.ones(len(icat))*len(icat)/len(iSR) #[0-1]#Bin size fraction of the entire data set. weights_sr_per_bin is effectively a constant in each bin and has no effect.
            # weights nfrequencies per sta value
            weights_nfrequencies_per_sta = np.array([sum(f<fnotch(s.StaDepth)) for s in icat.iloc])/sum(f<fnotch(icat.StaDepth.min())) #[0-1]#Fraction of frequencies used in this bin contributed to the average from each pair.
            # iweights value
            iweights = weights_nevents_per_sta + weights_sr_per_bin + weights_nfrequencies_per_sta
            weights.update({si:iweights})
            if sigma:
                Y_HPSH.update({si:np.array([s.Coherence.HPS_H for s in icat.iloc])})
                Y_HPSZ.update({si:np.array([s.Coherence.HPS_Z for s in icat.iloc])})
                Y_TF.update({si:np.array([s.Coherence.TF for s in icat.iloc])})
            else:
                Y_HPSH.update({si:np.array([s.Coherence.HPS_H for s in icat.iloc])})
                Y_HPSZ.update({si:np.array([s.Coherence.HPS_Z for s in icat.iloc])})
                Y_TF.update({si:np.array([s.Coherence.TF for s in icat.iloc])})
        else:
            # icat value
            icat=iSR.copy()



            [s.Coherence().TF[:,xf<fnotch(s.StaDepth)] for s in icat.iloc]

            # band value
            band=lambda:(xf[xf<fnotch(s.StaDepth)]<=max(1/set))&(xf[xf<fnotch(s.StaDepth)]>=min(1/set))

            # weights nevents per sta value
            weights_nevents_per_sta = np.array([sum(icat.StaName==s) for s in icat.StaName])/len(icat) #[0-1]#Number of events per station relative to the entire bin size.
            # weights sr per bin value
            weights_sr_per_bin = np.ones(len(icat))*len(icat)/len(iSR) #[0-1]#Bin size fraction of the entire data set. weights_sr_per_bin is effectively a constant in each bin and has no effect.
            # weights nfrequencies per sta value
            weights_nfrequencies_per_sta = np.array([sum(f<fnotch(s.StaDepth)) for s in icat.iloc])/sum(f<fnotch(icat.StaDepth.min())) #[0-1]#Fraction of frequencies used in this bin contributed to the average from each pair.
            # iweights value
            iweights = weights_nevents_per_sta + weights_sr_per_bin + weights_nfrequencies_per_sta
            weights.update({si:iweights})
            if sigma:
                Y_HPSH.update({si:np.array([s.Coherence.HPS_H for s in icat.iloc])})
                Y_HPSZ.update({si:np.array([s.Coherence.HPS_Z for s in icat.iloc])})
                Y_TF.update({si:np.array([s.Coherence.TF for s in icat.iloc])})
            else:
                Y_HPSH.update({si:np.array([s.Coherence.HPS_H for s in icat.iloc])})
                Y_HPSZ.update({si:np.array([s.Coherence.HPS_Z for s in icat.iloc])})
                Y_TF.update({si:np.array([s.Coherence.TF for s in icat.iloc])})

    # nrows = 1 if key=='StaDepth' else len(sets)
    layout='tight' if orientation=='horizontal' else 'none'
    # ncol value
    ncol=3;nrows=1
    if relative:
        fig,axes = plt.subplots(nrows=nrows,ncols=1,
        # figsize value
        figsize=figsize,
        # sharex value
        sharex='all',sharey='all',
        # layout value
        layout=layout) #if orientation=='horizontal' else 'none') #if orientation=='horizontal' else 'none'
    else:
        fig,axes = plt.subplots(nrows=nrows,ncols=ncol,
        # figsize value
        figsize=figsize,
        # sharex value
        sharex='all',sharey='all',
        # layout value
        layout=layout) #if orientation=='horizontal' else 'none') #if orientation=='horizontal' else 'none'
    # cmap value
    cmap=cm.cmaps['glasgow'].reversed().resampled(len(sets))
    if key=='Instrument_Design':cmap=ListedColormap([ColorStandard.instrument[s] for s in sets], name='custom_cmap')
    if key=='Network':cmap=ListedColormap([ColorStandard.network[s] for s in sets], name='custom_cmap')
    # axes value
    axes=np.atleast_2d(axes)
    # yttl value
    yttl = lambda lbl:fr"$\underset{{{lbl}}}{{\gamma\;\;\;\;\;\;\;}}$"
    sigma_yttl = lambda c:fr"$\sigma(\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}})$"
    if sigma:ttl=sigma_yttl
    else:ttl=yttl
    ylabel='' #Makes a (non-cumulative) histogram
    if cumulative:ylabel=ylabel+ 'Cumulative'
    if stacked:ylabel='Relative '+ylabel.lower()
    if density:ylabel=ylabel+' density' if len(ylabel)>0 else 'Density'
    else:ylabel=ylabel + ' count' if len(ylabel)>0 else 'Count'
    ylabel=ylabel+'\nof source-receiver pairs' #Makes a (non-cumulative) histogram
    ylabel=ylabel[0].upper()+ylabel[1:]
    dcoh=0.1;bins=np.arange(0,1+dcoh,dcoh) if not sigma else np.arange(0,0.5+(dcoh/2),(dcoh/2))

    keys=Y_TF.keys();keys=[k for k in keys if len(weights[k])] #remove empty bins.
    if weighted:wgt = [weights[k] for k in keys]
    else:wgt = [weights[k]*0+1 for k in keys]
    tf=[Y_TF[si] for si in keys]
    hpsz=[Y_HPSZ[si] for si in keys]
    hpsh=[Y_HPSH[si] for si in keys]
    color=[cmap(si/len(sets)) for si in keys]

    if outtype=='PDF':
        hi=1
    raxes=axes[0,:]
    lw=3
    nn=[]
    ax=raxes[0]
    if relative:
        ax.set_title(ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z'),y=1.1)
        dcoh=0.1;bins=np.arange(-1,1+dcoh,dcoh)
        y=[(a-b)/b for a,b in zip(tf,hpsz)]
        n,b,p=ax.hist(y,weights=wgt,bins=bins,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
        nn.append(n.max())
    else:
        ax.set_title(ttl('TF.Z'),y=1.1)
        n,b,p=ax.hist(tf,bins=bins,weights=wgt,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
        nn.append(n.max())
        ax=raxes[1]
        ax.set_title(ttl('HPS.Z'),y=1.1)
        n,b,p=ax.hist(hpsz,bins=bins,weights=wgt,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
        nn.append(n.max())
        ax=raxes[2]
        ax.set_title(ttl('HPS.H'),y=1.1)
        n,b,p=ax.hist(hpsh,bins=bins,weights=wgt,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
        nn.append(n.max())
    if norm_pdf:
        if np.max(nn)>1:
            for ax in raxes:ax.set_ylim(0,1.05*np.max(nn));yticks=np.arange(0,1.25,.25)*np.max(nn);ax.set_yticks(yticks);ax.set_yticklabels([f"{y/np.max(nn):.2f}" for y in yticks])
            # raxes[0].set_ylabel('Normalized probability density')
        ylabel=ylabel.replace('Density','probability density').replace('density','probability density')

    raxes[0].set_ylabel(ylabel.replace('  ',' '))
    if orientation=='vertical':fig.subplots_adjust(left=0, right=0.75) 
    if orientation=='horizontal':ax=raxes[1 if not relative else -1];cbar_ax = fig.add_axes([.15, 0.000, 0.8, 0.03])
    else:ax=raxes[-1];cbar_ax = fig.add_axes([.8, 0.15, 0.05, 0.7])
    if cumulative:
        for ax in raxes:ax.set_yticks(np.arange(0,1.25,.25));ax.grid(alpha=0.3)
    if key in ['Instrument_Design', 'Network', 'Pressure_Gauge', 'Environment', 'Seismometer']:
        ncat = len(sets)
        boundaries=np.arange(ncat + 1)
        norm=mpl.colors.BoundaryNorm(boundaries, ncat)
        sm=mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar_ticks=boundaries[:-1] + 0.5  # centers
        cbar_ticklabels = [str(s) for s in sets]
        label=key.replace('_', ' ')
    else:
        # For continuous or binned variables (e.g., StaDepth, Magnitude)
        boundaries = np.sort(np.unique(sets))
        norm = mpl.colors.Normalize(vmin=boundaries.min(), vmax=boundaries.max())
        sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar_ticks = boundaries
        cbar_ticklabels = [str(s) for s in boundaries]
        label = 'Water depth, m' if key == 'StaDepth' else 'Magnitude, Mw'
        if key == 'Sediment_Thickness_m':label = 'Sediment thickness, m'
    cbar=fig.colorbar(sm, cax=cbar_ax, boundaries=boundaries, orientation=orientation, label=label,shrink=0.7, aspect=30)
    cbar.set_ticks(cbar_ticks if not (key in ['StaDepth','Sediment_Thickness_m']) else cbar_ticks[::2])
    cbar.set_ticklabels(cbar_ticklabels if not (key in ['StaDepth','Sediment_Thickness_m']) else cbar_ticklabels[::2])
    fold=dirs.P01.S10
    file=f'{'Octav' if octav else 'NoOctav'}.{'Notched' if notched else 'NoNotch'}_{key}.png'
    # if cumulative:file = 'Cumulative.' + file
    # if density:file = 'Density.' + file
    if stacked:file='Stacked.' + file
    if sigma:file='Sigma.'+file
    else:file='Mu.'+file
    if weighted:file='Weighted.'+file
    else:file='Unweighted.'+file
    if norm_pdf:file='Normed.'+file
    file=f'{outtype}.'+file
    if relative:fold=fold/'Relative'
    if stacked:fold=fold/'Stacked'
    else:fold=fold/'NotStacked'
    fold=fold/outtype
    fold.mkdir(parents=True, exist_ok=True)

    file=file.replace('.png',note+'.png')

    save_tight(fold/file,fig,dpi=700)
    print(f'{key} - Saved')
    plt.close()

print(f"Elapsed time: {(runtime())/60:.2f} minutes")