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

from local_tools.math import octave_average


# function custom cmap
def custom_cmap(ind=0,nbins=5):
    if ind==0:cmap = cm.cmaps['glasgow'].reversed().resampled(nbins)
    if ind==1:cmap = cm.cmaps['batlow'].reversed().resampled(nbins)
    return cmap

# ----------------------------------------------------------------------------------------------------
bins = np.array([[i,i+500] for i in np.arange(0,6000,500)])
# ----------------------------------------------------------------------------------------------------
boundaries = np.hstack([0,bins[:,1].reshape(-1)])
# binedges value
binedges = np.array([max(b) for b in bins])
# binedges = np.array([min(b) for b in bins]) #min will put a 0-tick on the colorbar
nbins = len(bins)
# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
cmap=custom_cmap(1,nbins)
# cmap value
cmap=custom_cmap(0,nbins)

# instrument colors value
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
# variable
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
# seismometer marker value
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
# colors value
colors=np.array([instrument_colors[sta.Deployment.Instrument_Design] for sta in catalog.r.iloc])
# markers value
markers=np.array([seismometer_marker[sta.Deployment.Seismometer] for sta in catalog.r.iloc])
# Colors InstrumentDesigns value
Colors_InstrumentDesigns=np.unique([sta.Deployment.Instrument_Design for sta in catalog.r.iloc])
# Markers Seismometers value
Markers_Seismometers=np.unique([sta.Deployment.Seismometer for sta in catalog.r.iloc])
# legmarkersize value
legmarkersize=400
# markersize value
markersize=40
# ncols value
ncols=2
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
runs = {}
############################################################
req='Magnitude'
# sets value
sets=[[6,6.5],[6.5,7.0],[7.0,7.5],[7.5,8.0]]
# set ttl value
set_ttl = lambda set:f'Mw {float(min(set))} - {float(max(set))}'
# figsize value
figsize=(6,6)
# ncols value
ncols=2
runs.update({req:[ncols,sets,set_ttl,figsize]})
############################################################
req='StaDepth'
# sets value
sets=np.array([[0,500],[500,2500],[2500,4000],[4000,6000]])
# figsize value
figsize=(6,6)
# ncols value
ncols=2
# set ttl value
set_ttl=lambda set:f'{int(min(set))} - {int(max(set))}m'
runs.update({req:[ncols,sets,set_ttl,figsize]})
############################################################
req='Instrument_Design'
# sets value
sets=['AB','AR','BA','TRM','B2','BG','KE']
figsize=(6,6)
ncols=3
set_ttl=lambda set:f'{set}'
runs.update({req:[ncols,sets,set_ttl,figsize]})
###########################################################
req='Seismometer'
sets=['Trillium Compact','Trillium 240','Guralp CMG3T 120']
figsize=(6,6)
ncols=2
set_ttl=lambda set:f'{set}'
runs.update({req:[ncols,sets,set_ttl,figsize]})
###########################################################
req='Pressure_Gauge'
sets=['DPG','APG']
figsize=(5,5)
ncols=2
set_ttl=lambda set:f'{set}'
runs.update({req:[ncols,sets,set_ttl,figsize]})



# --------------------------------------------------------
# --------------------------------------------------------
sigma=False ####Sets the code to plot in standard deviation mode rather than mean mode.
# --------------------------------------------------------
# --------------------------------------------------------
ymethod='TF_Z';xmethod='HPS_Z'
# ymethod='HPS_H';xmethod='HPS_Z'


cat=catalog.sr.copy()
usnr=unpack_metrics(cat)
f=1/usnr.coh.bands

plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_supplemental_figures'/'FigureS10_Coh.vs.Coh.Scatter';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'


ttl=lambda tl:fr"$\underset{{{tl}}}{{\gamma{'\;\;'*len(tl)}}}$"
f_ind=(f>0)&(f<=1);xf=f[f_ind]
y=0.95;fontweight='bold'
oct_av=True
plot_individual_pairs=False


if sigma:stat=np.std;stat_kwargs={'ddof':1} #Sample standard deviation
else:stat=np.mean;stat_kwargs=None

for req in runs.keys():
    ncols,sets,set_ttl,figsize=runs[req]
    nrows=int(np.ceil(len(sets)/2))
    print(f'Running {req}...')
    for fraci,fraction in enumerate([8]):
        fig,axes=plt.subplots(nrows,ncols,figsize=figsize,sharex='all',sharey='all')
        inds = np.atleast_2d(np.arange(axes.size).reshape(axes.shape)).reshape(-1)
        for i in inds.reshape(-1)[len(sets):]:inds[i]=-1
        inds = np.atleast_2d(inds.reshape(axes.shape))
        leftinds=inds[:,0];leftinds=leftinds[leftinds>=0]
        bottominds=inds.max(axis=0);bottominds=bottominds[bottominds>=0]
        axes=axes.reshape(-1)
        for ax in axes:ax.set_box_aspect(1)
        for axi,(set,ax) in enumerate(zip(sets,axes)):
            print(f'{axi+1}/{len(sets)} | {set}')
            if isinstance(set,str):set_ind=(cat[req]==set)
            else:set_ind=(cat[req]>=min(set))&(cat[req]<=max(set))
            icat=cat[set_ind]
            X=usnr.coh.__dict__[xmethod].D.copy()
            Y=usnr.coh.__dict__[ymethod].D.copy()
            X=X[:,f_ind];Y=Y[:,f_ind]
            X=X[set_ind,:];Y=Y[set_ind,:]
            if oct_av:foct,Y=octavg(Y,xf);foct,X=octavg(X,xf)
            X=np.array([stat(xx[foct<fn]) for xx,fn in zip(X,fnotch(icat.StaDepth))])
            Y=np.array([stat(yy[foct<fn]) for yy,fn in zip(Y,fnotch(icat.StaDepth))])
            markers=[seismometer_marker[k] for k in icat.Seismometer]
            colors=[instrument_colors[k] for k in icat.Instrument_Design]

            # -----Plot individual source-receiver pairs
            if plot_individual_pairs:
                [ax.scatter(xx,yy,c='k',marker=mkr,s=markersize*.7,edgecolor='k',alpha=0.02) for xx,yy,clr,mkr in zip(X,Y,colors,markers) if mkr=='x']
                [ax.scatter(xx,yy,c=clr,marker=mkr,s=markersize*.7,edgecolor=None if mkr=='x' else 'k',alpha=0.02) for xx,yy,clr,mkr in zip(X,Y,colors,markers)]
            # -----Plot the station average
            X=np.array([np.mean(X[icat.StaName==s]) for s in icat.StaName.unique()])
            Y=np.array([np.mean(Y[icat.StaName==s]) for s in icat.StaName.unique()])
            markers=[seismometer_marker[icat[icat.StaName==s].iloc[0].Seismometer] for s in icat.StaName.unique()]
            colors=[instrument_colors[icat[icat.StaName==s].iloc[0].Instrument_Design] for s in icat.StaName.unique()]
            zclr = np.array([cmap(catalog.r.loc[n].StaDepth[0]/max(binedges)) for n in icat.StaName.unique()])

            X=np.array(X);Y=np.array(Y)
            print(f' |   {req.upper()}   |  {set}')
            ax.plot([0,1],[0,1],c='k',lw=.7,alpha=0.7,ls='-.')
            for xx,yy,clr,mkr,zc in zip(X,Y,colors,markers,zclr):
                if mkr=='x':ax.scatter(xx,yy,c=zc,marker=mkr,s=markersize,edgecolor=zc,linewidth=.5)
                else:ax.scatter(xx,yy,c=clr,marker=mkr,s=markersize,edgecolor=zc,linewidth=1)

            if isinstance(set,str):ax.set_title(set_ttl(set),fontweight=fontweight,y=1.01 if req=='Seismometer' else 1.0)
            else:ax.set_title(set_ttl(set),fontweight=fontweight,y=0.99)


            if np.isin(axi,leftinds):ax.set_ylabel(ttl(ymethod.replace('_','.')) if not sigma else (f'{ttl(ymethod.replace('_','.'))} variance'),fontweight=fontweight)
            if np.isin(axi,bottominds):ax.set_xlabel(ttl(xmethod.replace('_','.')) if not sigma else (f'{ttl(xmethod.replace('_','.'))} variance'),fontweight=fontweight)

            if not sigma:ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_xticks([0,0.5,1.0]);ax.set_yticks([0,0.5,1.0])
            else:ax.set_xlim(0,.5);ax.set_ylim(0,.5);ax.set_xticks([0,0.25,.5]);ax.set_yticks([0,0.25,.5])

        for ax in axes[len(sets):]:ax.set_visible(False);ax.set_box_aspect(1)
        # ------------------------------------------------------------------------------------------------------------------------
        plt.subplots_adjust(wspace=0.2, hspace=0.2)
        plt.tight_layout()
        # -----Save
        fname=f'S10.cohVcoh_{xmethod.replace('_','')}_{ymethod.replace('_','')}{'_OctaveAv' if oct_av else None}_by{req}.{save_format}'
        if sigma:
            fname=fname.replace('S10','S10.sigma')
            (plotfolder/'sigma').mkdir(parents=True, exist_ok=True)
            save_tight(plotfolder/'sigma'/('S04.FigureS10_'+fname),fig,dpi=200)
        else:save_tight(plotfolder/('S04.FigureS10_'+fname),fig,dpi=200)
        plt.close('all');del fig,axes
        print('Done')