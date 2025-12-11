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

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent.parent))
import os,sys
# --------------------------------------------------------------------------------------------------------------
from imports import * #Standard imports for doing anything in this project. Approx. ~19 seconds.
# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
# plotfolder value
plotfolder=dirs.Ch1/'_supplemental_figures'/'FigureS9_CoherenceHistograms';plotfolder.mkdir(parents=True,exist_ok=True)
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# colors value
colors = AttribDict({'tf':cm.cmaps['bilbao'](1/6),'area':cm.cmaps['oslo'](5/6),'hpsh':cm.cmaps['bamako'](3/6),'hpsz':cm.cmaps['imola'](1/6)})
# figs value
figs = lambda r=3,c=1,f=(5,6),x='all',y='all':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout='constrained')
# width value
width=4;height=4 #Defaults
# opts value
opts=AttribDict({})
# # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# # --------------------------------------------------------------------------------
req = None #Standard setup
# figsize value
figsize=(6,2) #width,height
# sets value
sets = [None];transpose=False
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose,'title':'All data'})})
######### # --------------------------------------------------------------------------------
######### # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
req = 'Magnitude'
# figsize value
figsize=(6,2);transpose=True #width,height
# sets value
sets=[[6.0,7.0],[7.0,8.0]]
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose,'title':req.replace('_',' ') + ' (mW) '})})

# req value
req = 'Instrument_Design'
# figsize value
figsize=(6,2);transpose=True #width,height
# sets value
sets=['AB','AR','B2','BA','BG','KE','TRM']
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose,'title':req.replace('_',' ')})})
# # # ####### --------------------------------------------------------------------------------
req = 'Seismometer'
# figsize value
figsize=(6,2);transpose=True #width,height
# sets value
sets=['Guralp CMG3T 120','Trillium 240','Trillium Compact']
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose,'title':req.replace('_',' ')})})
# # ####### --------------------------------------------------------------------------------
req = 'Pressure_Gauge'
# figsize value
figsize=(6,2);transpose=True #width,height
# sets value
sets=['DPG','APG']
opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose,'title':req.replace('_',' ')})})


# sigma value
sigma=False
f=cat.sr.iloc[0].Data.Coherence().f
bands=[[1,10],[10,30],[30,100]]
justpairs=True
nbins = 10
density=False
octav=False
notched=False
plot_mean_vline=False
icat=cat.sr.copy()
usnr=unpack_metrics(icat)

yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"


mtr = 'coh'
# mtr = 'snr'

greek = yttl if mtr=='coh' else yttl_eta
combinations = list(itertools.product([True, False], repeat=2)) #4**2=16 combinations
for combi,(justpairs,octav) in enumerate(combinations):
    for req in opts.keys():
        print(f'{combi+1}/{len(combinations)} | {req}')
        figsize,sets,transpose,title=opts[req]['figsize'],opts[req]['sets'],opts[req]['transpose'],opts[req]['title']
        for seti,set in enumerate(sets):
            if req==None:partition=icat.StaDepth>0;figtitle=f'{title}' 
            else:
                if isinstance(set,str):partition=icat[req]==set;figtitle=f'{title}, {set}' 
                else:partition=(icat[req]>=min(set))&(icat[req]<=min(set));figtitle=f'{title}, {min(set)} to {max(set)}'
            partition=partition.to_numpy()


            # -------------------------
            DAT = AttribDict()
            if mtr=='coh':
                f=1/usnr[mtr].bands
                TF = usnr.__dict__[mtr].TF_Z.D
                HPSZ = usnr.__dict__[mtr].HPS_Z.D
                HPSH = np.array([usnr.__dict__[mtr].HPS_1.D,usnr.__dict__[mtr].HPS_2.D]).mean(axis=0) #Mean H1/H2 together
            else:
                f=1/usnr.__dict__[mtr].bands
                TF = usnr.__dict__[mtr].TF_Z.R().D
                HPSZ = usnr.__dict__[mtr].HPS_Z.R().D
                HPSH = np.array([usnr.__dict__[mtr].HPS_1.R().D,usnr.__dict__[mtr].R().HPS_2.D]).mean(axis=0) #Mean H1/H2 together
            if octav:
                xf,TF = octavg(TF,f)
                xf,HPSZ = octavg(HPSZ,f)
                xf,HPSH = octavg(HPSH,f)
            else:xf=f
            # -------------------------
            TF=TF[partition,:]
            HPSZ=HPSZ[partition,:]
            HPSH=HPSH[partition,:]


            fig,axes=figs(1,3,f=figsize,x=False,y='all')
            for bi,b in enumerate(bands):
                ax=axes[bi]
                f_ind=(xf>=(1/max(b)))&(xf<=(1/min(b)))
                bTF=TF[:,f_ind]
                bHPSZ=HPSZ[:,f_ind]
                bHPSH=HPSH[:,f_ind]
                if justpairs:bTF,bHPSZ,bHPSH=bTF.mean(axis=1),bHPSZ.mean(axis=1),bHPSH.mean(axis=1)
                bTF,bHPSZ,bHPSH=bTF.reshape(-1),bHPSZ.reshape(-1),bHPSH.reshape(-1)
                binmin,binmax=np.min([bTF,bHPSZ,bHPSH]),np.max([bTF,bHPSZ,bHPSH])
                # binmax=np.round(binmax,1)
                binmin=np.round(binmin,2)
                if binmin<.1:binmin=0
                if binmax>.9:binmax=1.0
                bins=np.linspace(binmin,binmax,nbins)
                # Plotting
                yvals=[bTF,bHPSZ,bHPSH]
                labels = [greek('TF Z'),greek('HPS Z'),greek('HPS H'),]
                colorvals=[colors.tf,colors.hpsz,colors.hpsh]
                for y,c,lbl in zip(yvals,colorvals,labels):
                    ax.hist([y],bins=bins,weights=None,
                    cumulative=False,
                    density=density,
                    stacked=False,
                    histtype='step',
                    alpha=1.0,
                    linewidth=1,
                    zorder=1e3,
                    facecolor='w',
                    edgecolor=[c],
                    label=lbl)
                    if plot_mean_vline:ax.axvline(y.mean(),c=c,linewidth=3,alpha=0.5,zorder=-1e3,ls='-')
                # ------
                xlim=[binmin,binmax]
                ax.set_xlim(xlim)
                if not density:ax.set_yscale('log')
                ax.set_title(f'{min(b)} to {max(b)}s')
                ax.set_xticks(np.linspace(xlim[0],xlim[-1],3))
                ax.set_ylabel('Sample counts' if not justpairs else 'Source-receiver pairs')
                if (seti==0) & (bi==0):ax.legend()
            fig.suptitle(figtitle,y=1.1)
            for ax in axes:ax.set_xticklabels([f'{i:.2f}' for i in ax.get_xticks()])
            if req==None:file=f'Summary.Bands.png'
            else:file=f'{req}.{set if isinstance(set,str) else f'{min(set)}to{max(set)}'}.Bands.png'
            file =('PDF.' if density else 'Counts.')+file
            if justpairs:file='SRPairs.' + file
            else:file='Sample.' + file
            if req==None:file='_'+file
            file = f'{mtr.upper()}.'+file
            save_tight(plotfolder/file,fig,dpi=800)
            plt.close()