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
cat.sort_values(by='StaDepth',ascending=True,inplace=True)
# ------------------------------
def coherence_distribution_byfreq_plot(iDAT,xf,BAND_WINS,cat,nbins=10,bynetwork=False,notched=True):
    # bins value
    bins=np.linspace(0,1.0,nbins+1)
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=[9,5],sharex='col',sharey='all')
    # f ind value
    f_ind = ((xf>0) & (xf<1))
    # icat value
    icat = cat.copy()
    for i,win in enumerate(BAND_WINS):
        # nets value
        nets = icat.Network.unique()
        # d TF Z value
        d_TF_Z = iDAT[i].TF
        # d HPS Z value
        d_HPS_Z = iDAT[i].HPS_Z
        # d HPS H value
        d_HPS_H = iDAT[i].HPS_H
        # -----
        sets = [d_TF_Z,d_HPS_Z,d_HPS_H]
        # sets ttl value
        sets_ttl = ['TF.Z','HPS.Z','HPS.H']
        for j,d in enumerate(sets):
            ax = axes[i,j]
            if j==0:ax.set_ylabel(f'{int(1/max(win))} to {int(1/min(win))}s')
            if (i==2)&(j==1):ax.set_xlabel('Coherence')
            # ttl value
            ttl = fr"$\underset{{{sets_ttl[j]}}}{{\gamma\;\;\;\;\;\;\;}}$"
            if i==0:ax.set_title(ttl,y=1.15,fontsize=12,fontweight='bold')
            # cumulative value
            cumulative=True;density = True;ylabel='source-receiver pairs, \ncumulative fraction' #Makes a CDF
            # cumulative=True;density = False;ylabel='source-receiver pairs, \ncumulative counts' #Makes a cumulative histogram
            # cumulative=False;density = False;ylabel='source-receiver pairs, \ncounts' #Makes a (non-cumulative) histogram
            if bynetwork:
                # d value
                d = ([np.hstack([d[si] for si in np.where(icat.Network==n)[0]]) for n in nets]);ylabel = ylabel.replace('pairs', 'pairs by network')
                # pl value
                pl = [ax.hist(di,density=density,cumulative=cumulative,stacked=True,bins=bins,label=n,alpha=1,histtype='step',color=mycmaps.categorical.Networks[n]) for di,n in zip(d,nets)]
            else:
                # d value
                d = np.hstack(d);ax.hist(d,density=density,cumulative=cumulative,histtype='bar',bins=bins,alpha=1,facecolor='k',edgecolor='w',linewidth=2,zorder=-1000)
            ax.grid(True,which='both',alpha=.3)
            if j==2:ax.tick_params(axis='y', which='both', labelright=True, labelleft=False)
            else:ax.tick_params(axis='y', which='both', labelright=False, labelleft=False)
            if i==2:
                if len(bins)<5:ax.set_xticks(bins)#ax.set_xticklabels(['%1s'% b for b in bins])
                else:ax.set_xticks([0,.25,.5,.75,1])
    axes[1,2].tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
    ax_right = axes[1,2].twinx();ax_right.set_ylim(axes[1,2].get_ylim());ax_right.set_ylabel(ylabel, fontweight='bold')
    ax_right.tick_params(axis='y', which='both', left=False, right=True, labelleft=False, labelright=True)
    if bynetwork:
        h=[axes[-1,1].scatter(np.nan,np.nan,marker='s',c=mycmaps.categorical.Networks[n],label=n) for n in cat.Network.unique()]
        axes[-1,1].legend(ncol=3,handles=h,loc='lower center',bbox_to_anchor=(0.5, -.9),frameon=False,ncols=10,markerscale=1.5)
    return fig




s = cat.iloc[0]
f=s.Data.Coherence().f
foct=octavg(s.Data.Coherence().ATaCR.zp_21.coh,f)[0]
BAND_WINS = [[1/10,1/30],[1/30,1/100],[1/100,1/300]]
notched = True;octave_av = True
DAT = [0,0,0]
nrows = len(BAND_WINS)
cols = ['TF','HPS','Horizontals']
ncols = len(cols)
for notched in [True]:
    for octave_av in [True]:
        DAT = [0,0,0]
        for i,win in enumerate(BAND_WINS):
            DAT[i] = AttribDict()
            if octave_av:
                xf=foct;f_ind = (xf>=min(win))&(xf<max(win))
                DAT[i].TF = [octavg(s.Data.Coherence().ATaCR.zp_21.coh,f)[1][:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                DAT[i].TF = [y.mean(axis=1) if y.size>1 else np.array(np.nan) for y in DAT[i].TF]
                DAT[i].HPS_Z = [octavg(s.Data.Coherence().HPS.zz.coh,f)[1][:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                DAT[i].HPS_Z = [y.mean(axis=1) if y.size>1 else np.array(np.nan) for y in DAT[i].HPS_Z]
                h1 = [octavg(s.Data.Coherence().HPS['11'].coh,f)[1][:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                h1 = [y.mean(axis=1) if y.size>1 else np.array(np.nan) for y in h1]
                h2 = [octavg(s.Data.Coherence().HPS['22'].coh,f)[1][:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                h2 = [y.mean(axis=1) if y.size>1 else np.array(np.nan) for y in h2]
                DAT[i].HPS_H = [np.nanmean(np.array([ha,hb]),axis=0) for ha,hb in zip(h1,h2)]
            else:
                xf=f;f_ind = (xf>=min(win))&(xf<max(win))
                DAT[i] = AttribDict()
                DAT[i].TF = [s.Data.Coherence().ATaCR.zp_21.coh[:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                DAT[i].TF = np.array([y.mean(axis=0).mean() if y.size>1 else np.array(np.nan) for y in DAT[i].TF])
                DAT[i].HPS_Z = [s.Data.Coherence().HPS.zz.coh[:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                DAT[i].HPS_Z = np.array([y.mean(axis=0).mean() if y.size>1 else np.array(np.nan) for y in DAT[i].HPS_Z])
                h1 = [s.Data.Coherence().HPS['11'].coh[:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                h1 = np.array([y.mean(axis=0).mean() if y.size>1 else np.array(np.nan) for y in h1])
                h2 = [s.Data.Coherence().HPS['22'].coh[:,f_ind & (xf<fnotch(s.StaDepth) if notched else True)] for s in cat.iloc]
                h2 = np.array([y.mean(axis=0).mean() if y.size>1 else np.array(np.nan) for y in h2])
                DAT[i].HPS_H = np.nanmean(np.vstack([h1,h2]),axis=0)
        # ------------------------------
        iDAT = DAT.copy()
        for bynetwork in [True]:
            if octave_av:xf=foct
            else: xf=f
            fig = coherence_distribution_byfreq_plot(iDAT,xf,BAND_WINS,cat,bynetwork=bynetwork,notched=notched)
            #----Save plot
            options = ''.join(
            ['.byNetwork' if bynetwork else '',
            '.notched' if notched else '',
            '.octav'if octave_av else ''])
            file = f'Coherence_Distribution_by_Freq_{'_'.join(cols)}'+options+'.png'
            fold = dirs.P01.S10/('notched' if notched else 'not.notched')/('octav' if octave_av else 'not.octav')
            fold.mkdir(parents=True, exist_ok=True)
            save_tight(fold/file,dpi=700,fig=fig)
            print(f'----\nBy network: {bynetwork}\nNotched: {notched}\nOctave averaged: {octave_av}')