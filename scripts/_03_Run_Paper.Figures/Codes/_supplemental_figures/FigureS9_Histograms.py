### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
from matplotlib.ticker import PercentFormatter
import time;start=time.time()
# runtime value
runtime=lambda:int(time.time()-start)


# plotfolder value
plotfolder=dirs.Ch1/'_supplemental_figures'/'FigureS9_Histograms';plotfolder.mkdir(parents=True,exist_ok=True)



# icat value
icat=catalog.sr.copy()
# usnr value
usnr=unpack_metrics(icat)
# yttl value
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# yttl eta value
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"
# -------------------------------------------------------------------------------------------------------------------------------------
meta_wins=AttribDict()
meta_wins.StaDepth=np.array([[i,i+500] for i in np.arange(0,6000,500)])
meta_wins.Sediment_Thickness_m=[[i,i+500] for i in np.arange(0,7320+500,500)]
meta_wins.Magnitude=[[i,i+.5] for i in np.arange(6,8,.5)]
meta_wins.Network=['2D','7A','7D','X9','XF','YL','YO','ZA','ZN']
meta_wins.Seismometer=['Guralp CMG3T 120', 'Trillium 240', 'Trillium Compact']
meta_wins.Instrument_Design=['AB', 'AR', 'B2', 'BA', 'BG', 'KE', 'TRM']
# meta_wins.Environment=['North Atlantic', 'North Pacific', 'Solomon Sea', 'South Pacific']
meta_wins.Pressure_Gauge=['APG', 'DPG']
# -------------------------------------------------------------------------------------------------------------------------------------
# cumulative=1;density=True;outtype='CDF' #Makes a CDF
# cumulative=-1;density=True;outtype='EDF' #Makes a CDF
cumulative=False;density=False;outtype='Hist' #Makes a simple Histogram
# stacked value
stacked=True
# orientation value
orientation='horizontal' #colorbar
# orientation='vertical' #colorbar
figsize=[6,3.15]
# density value
density = False
# norm pdf value
norm_pdf = False
# -------------------------------------------------------------------------------------------------------------------------------------
bands = ['1_10','10_30','30_100',]
# phases value
phases = ['P','S','Rg']
# layout value
layout='tight' if orientation=='horizontal' else 'none'
# ncols value
ncols=4;nrows=3
# methods value
methods = ['TF Z','HPS Z','HPS 1','HPS 2']
# mtr value
mtr = 'coh'
figlist = phases if mtr=='snr' else ['Coherence']
greek = yttl_eta if mtr=='snr' else yttl
for key in meta_wins.keys():
    for pi,p in enumerate(figlist):
        fig,axes = plt.subplots(nrows=len(bands),ncols=ncols,figsize=figsize,sharex='col',sharey='row',layout=layout)
        axes = np.atleast_2d(axes)
        for mi,method in enumerate(methods):
            if method=='HPS H':
                mset=usnr.__dict__[mtr].__dict__['HPS_1']
                mset.D=np.array([usnr.__dict__[mtr].__dict__['HPS_1'].D,usnr.__dict__[mtr].__dict__['HPS_2'].D]).mean(axis=0)
            else:mset=usnr.__dict__[mtr].__dict__[method.replace(' ','_')]
            if mtr=='snr':mset=mset.R()[p]

            sets=meta_wins[key]
            cmap=cm.cmaps['glasgow'].reversed().resampled(len(sets))
            if key=='Instrument_Design':cmap=ListedColormap([ColorStandard.instrument[s] for s in sets], name='custom_cmap')
            if key=='Network':cmap=ListedColormap([ColorStandard.network[s] for s in sets], name='custom_cmap')
            Y_DAT={}
            for bi,b in enumerate(bands):
                raxes=axes[bi,:]
                ax=raxes[mi]
                print(f'{b} | {p} [{key}]')
                axtitle=p
                ylabel='source-receiver pairs'
                for si,set in enumerate(sets):
                    iSNR=icat.copy()
                    if isinstance(set,list)or(type(set)==type(np.array([]))):ind=(iSNR[key]>=min(set))&(iSNR[key]<max(set))
                    else:ind=iSNR[key]==set
                    cat=icat[ind].copy()
                    Y_DAT.update({si:mset.Average(tuple(list(map(int,b.split('_')))))[ind]})
                keys=Y_DAT.keys()
                Y=[Y_DAT[si].squeeze() for si in keys]
                color=[cmap(si/len(sets)) for si in keys]

                for si,(y,clr) in enumerate(zip(Y,color)):
                    if len(y)==0:continue
                    if bi==0:ax.set_title(greek(method),y=1.2)
                    lw=1
                    bins=np.linspace(y.min(),y.max(),20)
                    n,bns,h=ax.hist(y,bins=bins,density=density,cumulative=cumulative,stacked=stacked,histtype='step',alpha=1,linewidth=lw,zorder=-1000,color=clr) #facecolor=None,edgecolor='w'
                    if (mi==0)&(bi==1):ax.set_ylabel(ylabel.replace('  ',' '),labelpad=7)

                    if (pi+1)==len(phases):
                        rightlabel=f'{b.replace('_',' to ')}s'
                        ax_right=ax.twinx()
                        ax_right.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=True)
                        ax_right.set_ylabel(f'{rightlabel}',fontweight='bold',fontsize=9,rotation=-90,labelpad=15)
                        ax_right.set_yticks([])

            if orientation=='vertical':fig.subplots_adjust(left=0, right=0.75) 
        if orientation=='horizontal':cbar_ax = fig.add_axes([.15, 0.000, 0.8, 0.03])
        else:ax=raxes[-1];cbar_ax = fig.add_axes([.8, 0.15, 0.05, 0.7])
        if key in ['Instrument_Design', 'Network', 'Pressure_Gauge', 'Environment', 'Seismometer']:
            ncat=len(sets)
            boundaries=np.arange(ncat + 1);norm=mpl.colors.BoundaryNorm(boundaries, ncat)
            sm=mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar_ticks=boundaries[:-1]+0.5 # centers
            cbar_ticklabels=[str(s) for s in sets]
            label=key.replace('_', ' ')
        else:
            # For continuous or binned variables (e.g., StaDepth, Magnitude)
            boundaries=np.sort(np.unique(sets));norm=mpl.colors.Normalize(vmin=boundaries.min(), vmax=boundaries.max())
            sm=mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar_ticks = boundaries
            cbar_ticklabels = [str(s) for s in boundaries]
            label = 'Water depth, m' if key == 'StaDepth' else 'Magnitude, Mw'
            if key == 'Sediment_Thickness_m':label = 'Sediment thickness, m'
        cbar_ticks=cbar_ticks if not (key in ['StaDepth','Sediment_Thickness_m']) else cbar_ticks[::2]
        cbar_ticklabels=cbar_ticklabels if not (key in ['StaDepth','Sediment_Thickness_m']) else cbar_ticklabels[::2]
        cbar=fig.colorbar(sm, cax=cbar_ax, boundaries=boundaries, orientation=orientation, label=label,shrink=0.7, aspect=30)
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels(cbar_ticklabels)
        
        file=f'{outtype}.{key}.{mtr.upper()}.png'
        if stacked:file=file+'.Stacked'
        else:file=file+'.Mu'
        if norm_pdf:file=file+'.Normed'
        file=file+'.png'
        save_tight(plotfolder/file,fig,dpi=700)
        print(f'{key} - Saved')
        plt.close()

print(f"Elapsed time: {(runtime())/60:.2f} minutes")