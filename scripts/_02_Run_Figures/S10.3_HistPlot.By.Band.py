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
# meta_wins.StaDepth=np.array([[i,i+500] for i in np.arange(0,6000,500)])
bands=np.array([[1,10],[10,30],[30,300]])


# meta_wins.Magnitude=[[6,7.0],[7,8.0]] # meta_wins.Magnitude=[[i,i+.5] for i in np.arange(6,8,.5)]
# meta_wins.Sediment_Thickness_m=[[i,i+500] for i in np.arange(0,7320+500,500)]
# meta_wins.Seismometer=['Guralp CMG3T 120', 'Trillium 240', 'Trillium Compact']
# meta_wins.Instrument_Design=['AB', 'AR', 'B2', 'BA', 'BG', 'KE', 'TRM']
# meta_wins.Pressure_Gauge=['APG', 'DPG']
# meta_wins.Environment=['North Atlantic', 'North Pacific', 'Solomon Sea', 'South Pacific']
# meta_wins.Network=['2D','7A','7D','X9','XF','YL','YO','ZA','ZN']

# M value
M=AttribDict()
# M.Magnitude=[[6,7.0],[7,8.0]] # meta_wins.Magnitude=[[i,i+.5] for i in np.arange(6,8,.5)]
# M.Sediment_Thickness_m=[[i,i+500] for i in np.arange(0,7320+500,500)]
# M.Seismometer=['Guralp CMG3T 120', 'Trillium 240', 'Trillium Compact']
# M.Instrument_Design=['AB', 'AR', 'B2', 'BA', 'BG', 'KE', 'TRM']
# M.Pressure_Gauge=['APG', 'DPG']
# M.Environment=['North Atlantic', 'North Pacific', 'Solomon Sea', 'South Pacific']
# M.Network=['2D','7A','7D','X9','XF','YL','YO','ZA','ZN']
M.StaDepth=np.array([[i,i+500] for i in np.arange(0,6000,500)])


# nsets value
nsets = sum([len(i) for i in list(meta_wins.values())])
# -------------------------------------------------------------------------------------------------------------------------------------
# ---
# cumulative=1;density=True;outtype='CDF' #Makes a CDF
cumulative=-1;density=True;outtype='EDF' #Makes a CDF
# cumulative=False;density=True;outtype='PDF' #Makes a PDF
# cumulative=False;density=False;outtype='Hist' #Makes a Histogram
run_reduce = False #
# Re-compiles all data from the raw-files. 
# Takes about 1-min. If disabled will load required data from backups.
notched=True #Band limit to only periods sensitive to infragravity surface waves.
# octav value
octav=True #Octave averaging. Adds 20s to run time for each category/figure.
# If true, justpairs reduces coherence values given to the histogram from single frequency measurements for each source-receiver 
# to a single average for each source-receiver. Recommended to keep this at True.
justpairs=True #Recommended to always keep this at True.
# sigma value
sigma=True
# relative value
relative=False
# stacked value
stacked=False
# unwrapped value
unwrapped = False
# transpose value
transpose = False
# orientation value
orientation='horizontal' #colorbar
# orientation='vertical' #colorbar
weighted=True #
# compact value
compact = False
# norm pdf value
norm_pdf=False
# figsize value
figsize=[6,3.15]
# nkeys value
nkeys=len(meta_wins)>1
if sigma:justpairs=False
if nkeys:unwrapped=False
if not compact:unwrapped=False
# band keys value
band_keys=unravel([[f'{k}.{min(b)}to{max(b)}'  for k in ['TF','HPS_Z','HPS_H']] for b in bands])
# band keys value
band_keys={k:k.split('.')[-1].split('to') for k in band_keys}
# band keys value
band_keys={k:[int(band_keys[k][0]),int(band_keys[k][1])] for k in band_keys.keys()}
# -------------------------------------------------------------------------------------------------------------------------------------
f=cat.r.Data[0].Coherence().f
# foct value
foct=octavg(cat.sr.Data[0].Coherence().ATaCR.zp_21.coh,f)[0]
# safekeys value
safekeys=['Name', 'StaName', 'Station', 'Event', 'Network', 'LaLo', 'Distance',  'Magnitude', 'Stations',
'Latitude', 'Longitude', 'Experiment', 'Environment', 'Pressure_Gauge','StaDepth', 'Start', 'End', 'NoiseAverage', 'Seismometer',
'Sediment_Thickness_m', 'Instrument_Design', 'Distance_from_Land_km','Distance_to_Plate_Boundary_km', 'Surface_Current_ms',
'Crustal_Age_Myr', 'Deployment_Length_days', 'Inventory', 'Coherence','PearsonCC']
# function reduce data
def reduce_data(SR,SR_std,band_keys,f):
    # flim value
    flim = lambda b,f: ((1/f)<=max(b))&((1/f)>min(b)) #Assumes b is period
    # notch lambda value
    notch_lambda = lambda z,f:f<fnotch(z) if notched else True
    # lims value
    lims = lambda b,f,k,z:flim(band_keys[k],f)&notch_lambda(s.StaDepth,f)
    for si,(s,s_std) in enumerate(zip(SR.iloc,SR_std.iloc)):
        print(f'Collecting data: {si+1}/{len(SR)}')
        s.Coherence.update({k :s.Coherence[k.split('.')[0]][:,lims(band_keys[k],f,k,s.StaDepth)] for k in band_keys.keys()})
        s_std.Coherence.update({k :s_std.Coherence[k.split('.')[0]][:,lims(band_keys[k],f,k,s_std.StaDepth)] for k in band_keys.keys()})
        if octav:
            s.Coherence.update({k:octavg(s.Coherence[k],f[lims(band_keys[k],f,k,s.StaDepth)])[1] if sum(lims(band_keys[k],f,k,s.StaDepth))>0 else None for k in band_keys.keys()})
            s_std.Coherence.update({k:octavg(s_std.Coherence[k],f[lims(band_keys[k],f,k,s_std.StaDepth)])[1] if sum(lims(band_keys[k],f,k,s_std.StaDepth))>0 else None for k in band_keys.keys()})
        if justpairs:
            s.Coherence.update({k :s.Coherence[k].mean() if not np.any(s.Coherence[k]==None) else None for k in band_keys.keys()})
            s_std.Coherence.update({k :s_std.Coherence[k].std() if not np.any(s_std.Coherence[k]==None) else None for k in band_keys.keys()})
    return SR,SR_std
# -------------------------------------------------------------------------------------------------------------------------------------



if run_reduce: #Takes ~2min
    # SR value
    SR=cat.sr.__deepcopy__()
    # SR std value
    SR_std=ds.dataspace().sr.copy() #independent load of data for std analysis. Adds 300mb and 20s to the compute, but it's worth the headache.
    SR,SR_std=reduce_data(SR,SR_std,band_keys,f)
    print(f"Reduction elapsed time: {(runtime())/60:.2f} minutes")
else: #Takes ~2seconds
    # bulk_cat_version='71725'
    bulk_cat_version='50125'
    if not notched:
        # SR value
        SR=load_pickle(dirs.P01.S10/bulk_cat_version/'SR_NoNotch.pkl')
        # SR std value
        SR_std=load_pickle(dirs.P01.S10/bulk_cat_version/'SR_std_NoNotch.pkl')
    else:
        # SR value
        SR=load_pickle(dirs.P01.S10/bulk_cat_version/'SR.pkl')
        # SR std value
        SR_std=load_pickle(dirs.P01.S10/bulk_cat_version/'SR_std.pkl')

# state value
state=lambda:f'C{int(cumulative)}.D{int(density)}.S{int(sigma)}.T{int(stacked)}.R{int(relative)}.N{int(norm_pdf)}'

# combinations = list(itertools.product([True, False], repeat=5)) #4**2=16 combinations
# for combi,(cumulative,density,sigma,relative,norm_pdf) in enumerate(combinations):
# print(f'Settings----{combi+1}/{len(combinations)}----')

# xf value
xf = foct if octav else f


for mi in M.keys():
    # meta wins value
    meta_wins = AttribDict()
    meta_wins[mi]=M[mi]


    for ki,key in enumerate(meta_wins.keys()):
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
                iweights = (1/weights_nevents_per_sta) * weights_sr_per_bin * weights_nfrequencies_per_sta # (sta/sr) * (sr/bin) * (samples/sta) = samples/bin
                weights.update({si:iweights})
                Y_HPSH.update({f'{k}__{si}':np.array([s.Coherence[k] for s in icat.iloc]) for k in band_keys.keys() if k.split('.')[0]=='HPS_H'})
                Y_HPSZ.update({f'{k}__{si}':np.array([s.Coherence[k] for s in icat.iloc]) for k in band_keys.keys() if k.split('.')[0]=='HPS_Z'})
                Y_TF.update({f'{k}__{si}':np.array([s.Coherence[k] for s in icat.iloc]) for k in band_keys.keys() if k.split('.')[0]=='TF'})
            else:
                # icat value
                icat=iSR.copy()


                # weights nevents per sta value
                weights_nevents_per_sta = np.array([sum(icat.StaName==s) for s in icat.StaName])/len(icat) #[0-1]#Number of events per station relative to the entire bin size.
                # weights sr per bin value
                weights_sr_per_bin = np.ones(len(icat))*len(icat)/len(iSR) #[0-1]#Bin size fraction of the entire data set. weights_sr_per_bin is effectively a constant in each bin and has no effect.
                # weights_nfrequencies_per_sta = np.array([sum(f<fnotch(s.StaDepth)) for s in icat.iloc])/sum(f<fnotch(icat.StaDepth.min())) #[0-1]#Fraction of frequencies used in this bin contributed to the average from each pair.
                iweights = weights_nevents_per_sta + weights_sr_per_bin #+ weights_nfrequencies_per_sta
                weights.update({si:iweights})
                Y_HPSH.update({f'{k}.{si}':np.array([s.Coherence[k] for s in icat.iloc]) for k in band_keys.keys() if k.split('.')[0]=='HPS_H'})
                Y_HPSZ.update({f'{k}.{si}':np.array([s.Coherence[k] for s in icat.iloc]) for k in band_keys.keys() if k.split('.')[0]=='HPS_Z'})
                Y_TF.update({f'{k}.{si}':np.array([s.Coherence[k] for s in icat.iloc]) for k in band_keys.keys() if k.split('.')[0]=='TF'})

        if ki==0:  #Initial figure setup
            # layout value
            layout='tight' if orientation=='horizontal' else 'none'
            # ncol value
            ncol=3;nrows=3
            # if compact:nrows,ncols=1
            if compact:
                fig,axes = plt.subplots(nrows=nrows,ncols=1,
                # figsize value
                figsize=figsize,
                # sharex value
                sharex='all',sharey='all',
                # layout value
                layout=layout) #if orientation=='horizontal' else 'none') #if orientation=='horizontal' else 'none'            
                axes=np.atleast_2d(axes).T
            elif relative:
                fig,axes = plt.subplots(nrows=nrows,ncols=1,
                figsize=figsize,
                sharex='all',sharey='all',
                layout=layout) #if orientation=='horizontal' else 'none') #if orientation=='horizontal' else 'none'
            else:
                fig,axes = plt.subplots(nrows=nrows,ncols=ncol,
                figsize=figsize,
                sharex='all',sharey='all',
                layout=layout) #if orientation=='horizontal' else 'none') #if orientation=='horizontal' else 'none'
            cmap=cm.cmaps['glasgow'].reversed().resampled((len(meta_wins.keys())+1) if nkeys else (len(sets)+1))
            # cmap=cm.cmaps['batlow'].reversed().resampled(len(meta_wins.keys()))
            if key=='Instrument_Design':cmap=ListedColormap([ColorStandard.instrument[s] for s in sets], name='custom_cmap')
            if key=='Network':cmap=ListedColormap([ColorStandard.network[s] for s in sets], name='custom_cmap')
            axes=np.atleast_2d(axes)
            yttl = lambda lbl:fr"$\underset{{{lbl}}}{{\gamma\;\;\;\;\;\;\;}}$"
            sigma_yttl = lambda c:fr"$\sigma(\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}})$"
            if sigma:ttl=sigma_yttl
            elif (cumulative & density):
                # ttl = lambda lbl:fr"$P\left(\underset{{{lbl}}}{{\gamma\;\;\;\;\;\;\;}}\right)$"
                # ttl = lambda lbl:fr"$P\left(\underset{{{lbl}}}{{P\left(\gamma\;\;\;\;\;\;\;}}\right)$"
                ttl = lambda lbl:fr"$\underset{{{lbl}}}{{\gamma\;\;\;\;\;\;\;}}$"
                # ttl = lambda lbl:fr"$P\left(\underset{{{lbl}}}{{\gamma\phantom{{xxxxx}}}}\right)$"
            else:ttl=yttl
            ylabel='' #Makes a (non-cumulative) histogram
            if cumulative:ylabel=ylabel+ 'Cumulative'
            if stacked:ylabel='Relative '+ylabel.lower()
            if density:ylabel=ylabel+' density' if len(ylabel)>0 else 'Density'
            else:ylabel=ylabel + ' count' if len(ylabel)>0 else 'Count'
            ylabel=ylabel+'\nof source-receiver pairs' #Makes a (non-cumulative) histogram
            ylabel=ylabel[0].upper()+ylabel[1:]
            dcoh=0.1;bins=np.arange(0,1+dcoh,dcoh) if not sigma else np.arange(0,0.5+(dcoh/2),(dcoh/2))
        special_colors = {'tf':cm.cmaps['bilbao'](1/6),'area':cm.cmaps['oslo'](5/6),'hpsh':cm.cmaps['bamako'](3/6),'hpsz':cm.cmaps['imola'](1/6)}
        for k in special_colors.keys():(special_colors[k][0],special_colors[k][1],special_colors[k][2],0.3)

        # -----

        for bi,b in enumerate(bands):
            if not transpose:raxes=axes[bi,:]
            else:raxes=axes[:,bi]
            band=lambda b:f'{min(b)}to{max(b)}'
            keys=list(Y_TF.keys());keys=[k.split('.')[-1] for k in keys if len(weights[int(k.split('__')[-1])])] #remove empty bins.
            if weighted:wgt = [weights[i] for i,_ in enumerate(sets)]
            else:wgt = [weights[i]*0+1 for i,_ in enumerate(sets)]
            tf=[Y_TF[f'TF.{band(bands[bi])}__{si}'] for si,_ in enumerate(sets)]
            hpsz=[Y_HPSZ[f'HPS_Z.{band(bands[bi])}__{si}'] for si,_ in enumerate(sets)]
            hpsh=[Y_HPSH[f'HPS_H.{band(bands[bi])}__{si}'] for si,_ in enumerate(sets)]
            # color=[cmap((si+ki)/nsets) for si,_ in enumerate(sets)]
            if nkeys:color=[cmap((ki)/len(list(meta_wins.keys()))) for si,_ in enumerate(sets)] #Colors by category type (e.g. Pressure), not sub-category (e.g. APG)
            else:color=[cmap(si/len(sets)) for si,_ in enumerate(sets)]
            nonetest = lambda x:x==None

            lw=0.6
            nn=[]
            ax=raxes[0];ax0=[]
            for si,_ in enumerate(sets):
                hpsz[si]=hpsz[si][~(tf[si]==None)];hpsh[si]=hpsh[si][~(tf[si]==None)]
                wgt[si]=wgt[si][~(tf[si]==None)];tf[si]=tf[si][~(tf[si]==None)]
            if unwrapped:
                hpsz=[unravel(hpsz)];hpsh=[unravel(hpsh)]
                tf=[unravel(tf)];wgt=[unravel(wgt)]
                color=[color[-1]]
            if relative:
                if (bi==0)&transpose:ax.y_label(ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z'),y=1.1)
                elif (bi==0)&(not transpose):ax.set_title(ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z'),y=1.1)
                dcoh=0.1;bins=np.arange(-1,1+dcoh,dcoh)
                y=[(a-b)/b for a,b in zip(tf,hpsz)]
                n,b,p=ax.hist(y,weights=wgt,bins=bins,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
                ax0.append([n,b,p])
                nn.append(n.max())
            else:
                if not compact:
                    n,b,p=ax.hist(tf,bins=bins,weights=wgt,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
                    ax0.append([n,b,p]);nn.append(n.max())
                    # n,b,p=ax.hist(tf,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1.0,linewidth=lw*1.5,zorder=-1000,color=None,edgecolor=color)
                    ax=raxes[1]
                    n,b,p=ax.hist(hpsz,bins=bins,weights=wgt,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
                    ax0.append([n,b,p]);nn.append(n.max())
                    # n,b,p=ax.hist(hpsz,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1.0,linewidth=lw*1.5,zorder=-1000,color=None,edgecolor=color)
                    ax=raxes[2]
                    n,b,p=ax.hist(hpsh,bins=bins,weights=wgt,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=color) #facecolor=None,edgecolor='w'
                    ax0.append([n,b,p]);nn.append(n.max())
                    # n,b,p=ax.hist(hpsh,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1.0,linewidth=lw*1.5,zorder=-1000,color=None,edgecolor=color)
                else:
                    # edgeC=[special_colors['hpsh']]
                    edgeC=[special_colors['hpsh']]
                    # n,b,p=ax.hist(hpsh,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,color=color)
                    if len(meta_wins)==1:n,b,p=ax.hist([unravel(hpsh)],bins=bins,weights=[unravel(wgt)],facecolor=['w'],density=density,cumulative=cumulative,stacked=stacked,histtype='barstacked',alpha=1,linewidth=.5,zorder=-1000,color=edgeC,edgecolor='k') #facecolor=None,edgecolor='w'
                    if len(meta_wins)==1:n,b,p=ax.hist([unravel(hpsh)],bins=bins,weights=[unravel(wgt)],facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=1.2,zorder=-500,color=edgeC,edgecolor=edgeC) #facecolor=None,edgecolor='w'
                    ax0.append([n,b,p]);nn.append(n.max())

                    edgeC=[special_colors['hpsz']]
                    # n,b,p=ax.hist(hpsz,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,color=color)
                    if len(meta_wins)==1:n,b,p=ax.hist([unravel(hpsz)],bins=bins,weights=[unravel(wgt)],facecolor=[special_colors['area']],density=density,cumulative=cumulative,stacked=stacked,histtype='barstacked',alpha=1,linewidth=.5,zorder=-1000,color=edgeC,edgecolor='k') #facecolor=None,edgecolor='w'
                    if len(meta_wins)==1:n,b,p=ax.hist([unravel(hpsz)],bins=bins,weights=[unravel(wgt)],facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=1.2,zorder=-500,color=edgeC,edgecolor=edgeC) #facecolor=None,edgecolor='w'
                    ax0.append([n,b,p]);nn.append(n.max())
                    # if ki==(len(meta_wins)-1):n,b,p=ax.hist(hpsz,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1.0,linewidth=1.2,color=None,edgecolor=color)

                    edgeC=[special_colors['tf']]
                    # n,b,p=ax.hist(tf,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,color=color)
                    if len(meta_wins)==1:n,b,p=ax.hist([unravel(tf)],bins=bins,weights=[unravel(wgt)],facecolor=['w'],density=density,cumulative=cumulative,stacked=stacked,histtype='barstacked',alpha=1,linewidth=.5,zorder=-1000,color=edgeC,edgecolor='k') #facecolor=None,edgecolor='w'
                    if len(meta_wins)==1:n,b,p=ax.hist([unravel(tf)],bins=bins,weights=[unravel(wgt)],facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=1.2,zorder=-500,color=edgeC,edgecolor=edgeC) #facecolor=None,edgecolor='w'
                    ax0.append([n,b,p]);nn.append(n.max())
                    # if ki==(len(meta_wins)-1):n,b,p=ax.hist(tf,bins=bins,weights=wgt,facecolor=None,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1.0,linewidth=1.2,color=None,edgecolor=color)
            if norm_pdf:
                if np.max(nn)>1:
                    for ax in raxes:ax.set_ylim(0,1.05*np.max(nn));yticks=np.arange(0,1.25,.25)*np.max(nn);ax.set_yticks(yticks);ax.set_yticklabels([f"{y/np.max(nn):.2f}" for y in yticks])
                    # raxes[0].set_ylabel('Normalized probability density')
                # ylabel=ylabel.replace('Density','probability density').replace('density','probability density')
            # if (ki==0)&(bi==0):raxes[1].set_ylabel(ylabel.replace('  ',' '))
            if ki==0:
                if transpose:raxes[0].set_title(f'{min(bands[bi])} to {max(bands[bi])}s')
                else:raxes[0].set_ylabel(f'{min(bands[bi])} to {max(bands[bi])}s')


            # special_colors = {'tf':cm.cmaps['devon'](1/6),'hpsz':cm.cmaps['bilbao'](1/6),'hpsh':cm.cmaps['oslo'](5/6)}
            # for axi,(ax,aa) in enumerate(zip(raxes,ax0)):
            #     n, bins, patches = aa
            #     n = np.atleast_2d(n)
            #     if unwrapped:n=np.atleast_2d(n).T
            #     for bin_idx in range(n.shape[1]):
            #         group_idx = np.argmin(n[:, bin_idx])
            #         for stack_idx, patch in enumerate([p[bin_idx] for p in patches]):
            #             if stack_idx == group_idx:patch.set_facecolor(special_colors[axi])  # your desired color
            #             else:patch.set_facecolor('white')

            # LaTEXT is breaking the ()
            if not transpose:raxes=axes[0,:]
            else:raxes=axes[:,0]
            # for ax,tl in zip(raxes,[ttl('TF.Z'),ttl('HPS.Z'),ttl('HPS.H')]):
            #     if relative:
            #         if (bi==0)&transpose:ax.y_label(f'P({ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z')})',y=1.1)
            #         elif (bi==0)&(not transpose):ax.set_title(ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z'),y=1.1)
            #     else:
            #         if (bi==0)&transpose:ax.y_label(tl,y=1.1)
            #         elif (bi==0)&(not transpose):ax.set_title(f'$P({tl})$',y=1.1)                
            if bi==0:
                for ax,tl in zip(raxes,[ttl('TF.Z'),ttl('HPS.Z'),ttl('HPS.H')]):
                    if relative:
                        if transpose:ax.y_label(ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z'),y=1.1)
                        elif (not transpose):ax.set_title(ttl('TF.Z')+'-'+ttl('HPS.Z')+'/'+ttl('HPS.Z'),y=1.1)
                    else:
                        if transpose:ax.y_label(tl,y=1.1)
                        elif (not transpose):ax.set_title(tl,y=1.1)

            
    _=[ax.set_yticks([0.0,0.5,1.0]) for ax in axes.reshape(-1)]
    _=[ax.grid(alpha=0.3) for ax in axes.reshape(-1)]

    if orientation=='vertical':fig.subplots_adjust(left=0, right=0.75) 
    # if cumulative:
    #     for ax in raxes:ax.set_yticks(np.arange(0,1.25,.25));ax.grid(alpha=0.3)
    if not unwrapped:
        # ---------------------------------------------------------------
        # COLORBAR-------
        if (orientation=='horizontal'):cbar_ax=fig.add_axes([.05, 0.000, 0.9, 0.02])
        else:cbar_ax=fig.add_axes([.8, 0.15, 0.05, 0.7])
        # if key in ['Instrument_Design', 'Network', 'Pressure_Gauge', 'Environment', 'Seismometer']:
        ncat,nsets=len(meta_wins.keys()),len(sets)
        n=ncat if nkeys else nsets
        labels = meta_wins.keys() if nkeys else sets
        boundaries=np.arange(n+1)
        norm=mpl.colors.BoundaryNorm(boundaries, n)
        sm=mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar_ticks=boundaries[:-1] + 0.5  # centers
        # cbar_ticklabels = [str(s) for s in sets]
        if key=='Magnitude':cbar_ticklabels=[f'{min(i)} to {max(i)}' for i in list(labels)]
        elif key in ['StaDepth','Sediment_Thickness_m']:cbar_ticklabels=[f'{np.mean(i)}' for i in list(labels)]
        else:cbar_ticklabels=[f'{i}' for i in list(labels)]
        label=key.replace('_m',', m').replace('_', ' ').replace('StaDepth','Water Depth')
        # else:
        #     # For continuous or binned variables (e.g., StaDepth, Magnitude)
        #     boundaries = np.sort(np.unique(sets))
        #     norm = mpl.colors.Normalize(vmin=boundaries.min(), vmax=boundaries.max())
        #     sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
        #     sm.set_array([])
        #     cbar_ticks = boundaries
        #     cbar_ticklabels = [str(s) for s in boundaries]
        #     label = 'Water depth, m' if key == 'StaDepth' else 'Magnitude, Mw'
        #     if key == 'Sediment_Thickness_m':label = 'Sediment thickness, m'
        cbar=fig.colorbar(sm, cax=cbar_ax, boundaries=boundaries, orientation=orientation, label=label,shrink=0.7, aspect=30)
        if not (key in ['StaDepth','Sediment_Thickness_m']):cbar.set_ticks(cbar_ticks);cbar.set_ticklabels(cbar_ticklabels)
        else:cbar.set_ticks(cbar_ticks[::2]);cbar.set_ticklabels(cbar_ticklabels[::2])
        # if nkeys:cbar.set_ticklabels(cbar_ticklabels,rotation=-45)
        # else:cbar.set_ticklabels(cbar_ticklabels)
        # if nkeys:cbar.set_ticklabels([i.replace('_','\n').replace('StaDepth','Water Depth').replace('\nm','').replace(' ','\n') for i in list(meta_wins.keys())],rotation=0)
        # else:cbar.set_ticklabels([i.replace('_','\n').replace('StaDepth','Water Depth').replace('\nm','').replace(' ','\n') for i in sets],rotation=0)
        # ------------------------------------------
        # ---------------------------------------------------------------

    fold=dirs.P01.S10
    file=f'{'Octav' if octav else 'NoOctav'}.png'
    # if cumulative:file = 'Cumulative.' + file
    # if density:file = 'Density.' + file
    if stacked:file='Stacked.' + file

    if unwrapped:file='Unwrapped.'+file
    if compact:file='Compact.'+file
    if sigma:file='Sigma.'+file
    else:file='Mu.'+file
    if weighted:file='Weighted.'+file
    else:file='Unweighted.'+file
    if norm_pdf:file='Normed.'+file
    file=f'{outtype}.{key}.{'Notched' if notched else 'NoNotch'}'+file
    if relative:fold=fold/'Relative'
    if stacked:fold=fold/'Stacked'
    else:fold=fold/'NotStacked'
    fold=fold/outtype
    fold = fold/('Mu' if not sigma else 'Sigma')
    fold.mkdir(parents=True, exist_ok=True)

    if nkeys:note = 'CategoricalOverview';file=file.replace(f'{key}',note)

    save_tight(fold/file,fig,dpi=700)
    print(f'{key} - Saved')
    plt.close()
    print(f"Elapsed time: {(runtime())/60:.2f} minutes")