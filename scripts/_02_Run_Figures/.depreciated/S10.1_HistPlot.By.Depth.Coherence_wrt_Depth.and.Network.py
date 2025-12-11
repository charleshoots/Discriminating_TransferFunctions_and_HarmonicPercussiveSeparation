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

# DEPTH WINS value
DEPTH_WINS = [[cat.StaDepth.min()-1,1500],[1500,3000],[3000,cat.StaDepth.max()+1]]
# nrows value
nrows = len(DEPTH_WINS)
# cols value
cols = ['TF','HPS','Horizontals']
# ncols value
ncols = len(cols)

# f value
f = cat.Data[0].Coherence().f;f_ind = ((f>0) & (f<1))
# zslice value
zslice = lambda cat,win:cat[cat.StaDepth.between(win[0],win[1])].copy()

# --------
# ---Comment out after running once. Takes ~2-min to run.
DAT = AttribDict()
for i,win in enumerate(DEPTH_WINS):
    # icat value
    icat = zslice(cat,win)
    DAT[i] = AttribDict()
    DAT[i].TF = {si:s.Data.Coherence().ATaCR.zp_21.coh for si,s in enumerate(icat.iloc)}
    DAT[i].HPS_Z = {si:s.Data.Coherence().HPS['zz'].coh for si,s in enumerate(icat.iloc)}
    DAT[i].HPS_H = {si:np.dstack([s.Data.Coherence().HPS['11'].coh,s.Data.Coherence().HPS['22'].coh]).mean(axis=2) for si,s in enumerate(icat.iloc)}
# # --------


# xf value
xf = f;iDAT = DAT.copy();octave_av = False
# octave av value
octave_av = True
if octave_av:
    # iDAT value
    iDAT = DAT.copy()
    # foct value
    foct = octavg(iDAT[0].TF[0][0,:],f)[0]
    for i in range(len(iDAT)):
        for si in range(len(iDAT[i].TF)):
            iDAT[i].TF[si] = np.array([octavg(iDAT[i].TF[si][ei,:],f)[1].reshape(-1) for ei in range(iDAT[i].TF[si].shape[0])])
            iDAT[i].HPS_Z[si] = np.array([octavg(iDAT[i].HPS_Z[si][ei,:],f)[1].reshape(-1) for ei in range(iDAT[i].HPS_Z[si].shape[0])])
            iDAT[i].HPS_H[si] = np.array([octavg(iDAT[i].HPS_H[si][ei,:],f)[1].reshape(-1) for ei in range(iDAT[i].HPS_H[si].shape[0])])
    # xf value
    xf = foct
    # iDAT oct value
    iDAT_oct = iDAT.copy()
    # iDAT value
    iDAT = DAT.copy()
xf = foct


def coherence_distribution_bydepth_plot(iDAT,xf,DEPTH_WINS,cat,nbins=10,bynetwork=False,notched=True):
    bins=np.linspace(0,1.0,nbins+1)
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,figsize=[9,5],sharex='col',sharey='all')
    f_ind = ((xf>0) & (xf<1))
    for i,win in enumerate(DEPTH_WINS):
        icat = zslice(cat,win)
        nets = icat.Network.unique()
        d_TF_Z = [iDAT[i].TF[si][:,f_ind].mean(axis=1) for si,s in enumerate(icat.iloc)]
        d_HPS_Z = [iDAT[i].HPS_Z[si][:,f_ind].mean(axis=1) for si,s in enumerate(icat.iloc)]
        d_HPS_H = [iDAT[i].HPS_H[si][:,f_ind].mean(axis=1) for si,s in enumerate(icat.iloc)]
        # -----
        # Left of notch
        if notched:
            d_TF_Z = [iDAT[i].TF[si][:,((xf>0) & (xf<fnotch(s.StaDepth)))].mean(axis=1) for si,s in enumerate(icat.iloc)]
            d_HPS_Z = [iDAT[i].HPS_Z[si][:,((xf>0) & (xf<fnotch(s.StaDepth)))].mean(axis=1) for si,s in enumerate(icat.iloc)]
            d_HPS_H = [iDAT[i].HPS_H[si][:,((xf>0) & (xf<fnotch(s.StaDepth)))].mean(axis=1) for si,s in enumerate(icat.iloc)]

        sets = [d_TF_Z,d_HPS_Z,d_HPS_H]
        sets_ttl = ['TF.Z','HPS.Z','HPS.H']
        for j,d in enumerate(sets):
            ax = axes[i,j]
            if j==0:ax.set_ylabel(f'{int(win[0])}-{int(win[1])}m')
            if (i==2)&(j==1):ax.set_xlabel('Coherence')
            ttl = fr"$\underset{{{sets_ttl[j]}}}{{\gamma\;\;\;\;\;\;\;}}$"
            if i==0:ax.set_title(ttl,y=1.15,fontsize=12,fontweight='bold')

            cumulative=True;density = True;ylabel='source-receiver pairs, \ncumulative fraction' #Makes a CDF
            # cumulative=True;density = False;ylabel='source-receiver pairs, \ncumulative counts' #Makes a cumulative histogram
            # cumulative=False;density = False;ylabel='source-receiver pairs, \ncounts' #Makes a (non-cumulative) histogram
            if bynetwork:
                d = ([np.hstack([d[si] for si in np.where(icat.Network==n)[0]]) for n in nets]);ylabel = ylabel.replace('pairs', 'pairs by network')
                pl = [ax.hist(di,density=density,cumulative=cumulative,stacked=True,bins=bins,label=n,alpha=1,histtype='step',color=mycmaps.categorical.Networks[n]) for di,n in zip(d,nets)]
            else:
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


bynetwork=False;notched=False
for bynetwork in [True,False]:
    for notched in [True,False]:
        for octave_av in [False,True]:
            if octave_av:
                datplotable = iDAT_oct.copy();x=foct
            else:
                datplotable = DAT.copy();x=f
            fig = coherence_distribution_bydepth_plot(datplotable,x,DEPTH_WINS,cat,bynetwork=bynetwork,notched=notched)
            #----Save plot
            options = ''.join(
            ['.byNetwork' if bynetwork else '',
            '.notched' if notched else '',
            '.octav'if octave_av else ''])
            file = f'Coherence_Distribution_by_Depth_{'_'.join(cols)}'+options+'.png'
            fold = dirs.P01.S10/('notched' if notched else 'not.notched')/('octav' if octave_av else 'not.octav')
            fold.mkdir(parents=True, exist_ok=True)
            save_tight(fold/file,dpi=700,fig=fig)
            print(f'----\nBy network: {bynetwork}\nNotched: {notched}\nOctave averaged: {octave_av}')