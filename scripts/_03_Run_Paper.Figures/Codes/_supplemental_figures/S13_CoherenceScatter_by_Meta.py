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

from imports import * #Standard imports for doing anything in this project. Approx. ~19 seconds.
# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
# Noise Spectra
f=cat.r.iloc[0].Data.Noise.Averaged().f
# faxis value
faxis=(f>0)&(f<=1)
# f value
f=f[faxis]
# noise f value
noise_f=f
# yttl value
yttl = lambda lbl:fr"$\underset{{{lbl}}}{{\gamma\;\;\;\;\;\;\;}}$"
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
# notch lambda value
notch_lambda = lambda z,b,f=f:((f<(fnotch(z) if notched else 1)))&(((1/f)<max(b)) & ((1/f)>min(b)))
# SR value
SR = cat.sr.copy();note='' #Everything
# SR=SR[SR.Magnitude<7.0];note='Mw6_7' #Only Mw6.0-7.0
# SR=SR[SR.Magnitude>=7.0];note='Mw7_8' # Only >=Mw7.0
tfzax,hpszax,hpshax = 0,1,2
if True: ## OPTIONS------------------------------------------------------------------------
    # plotfolder = dirs.P01.Parent
    plotfolder=dirs.Ch1/'_supplemental_figures'/'FigureS13_CoherenceScatter_by_Meta';plotfolder.mkdir(parents=True,exist_ok=True)


    # Subsets to only periods longer than the IGP
    notched=False
    # Applies a 1/8 octave averaging to coherence values. 
    # Adds atleast a minute to compute time (normally <20s)
    octav=True #Prepends file name
    # Add experiment names to dots on the scatter plots
    named = False
    # sta avg value
    sta_avg = True #Huge compute if disabled. Don't run unless you have time to spare.
    # mean value
    mean=True
    # deviation value
    deviation=False
    # save value
    save=True
    # stats value
    stats=[]
    # stat = AttribDict({'func':np.std,'kw':{'ddof':1},'YL':[0,0.5],'title':'Sigma'});stats.append(stat)
    # stat = AttribDict({'func':np.mean,'kw':{},'YL':[0,1.02],'title':'Mu'});stats.append(stat)
    stat = AttribDict({'func':np.mean,'kw':{},'YL':None,'title':'Mu'});stats.append(stat)

# -------------------------------####------------------vvvvv CODE vvvvv------------------###-------------------------------
def scatter(ax,x,y,snm,ColorStandard): #Base plot function for everything below
    for s in np.unique(snm):
        # ii value
        ii=snm==s
        # color value
        color=ColorStandard.instrument[cat.r.loc[s].Instrument_Design[0]]
        # marker value
        marker=ColorStandard.seismometer_marker[cat.r.loc[s].Seismometer[0]]
        # isx value
        isx=marker=='x'
        ax.scatter(x[ii],y[ii],c=color,marker=marker,
        # edgecolors value
        edgecolors='k',linewidths=1 if isx else 0.2,
        alpha=1.0 if sta_avg else 0.1)


icat=cat.sr.copy()
usnr=unpack_metrics(icat)
f=1/usnr.coh.bands

z = cat.sr.StaDepth
bands = [[1,10],[10,30],[30,100]]
# bands = [[1,100]]

keys=['StaDepth','Magnitude']
titles=['Water depth,m','Magnitude, Mw']
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
dbin = lambda x:np.array([x[:-1],x[1:]]).T


for fn in [None,'IG']:
    note=''
    if octav:note=note+'octav.'
    if not fn:note=note+'NoNotch.'
    else:note=note+'Notched.'

    for key,title in zip(keys,titles):

        if len(bands)==1:fig,axes=figs(f=(3.5,5),y='cols');axes=np.atleast_2d(axes).T
        else:fig,axes=figs(3,len(bands),f=(3*len(bands),5),y='cols');axes=np.atleast_2d(axes)

        print(key)

        for bi,b in enumerate(bands):
            for stat in stats:
                kw=stat.kw;func=stat.func

                vals=SR[key] #The indexable parameter to sort and aggregate the data by (ie water depth, magnitude)
                snm=SR.StaName.to_numpy()

                f=1/usnr.coh.bands
                H=np.array([usnr.coh.HPS_1.Average(tuple(b),fn=fn),usnr.coh.HPS_2.Average(tuple(b),fn=fn)]).mean(axis=0)
                ZHPS=usnr.coh.HPS_Z.Average(tuple(b),octave=octav,fn=fn)
                ZTF=usnr.coh.TF_Z.Average(tuple(b),octave=octav,fn=fn)
                is_a_receiver_param = key not in ['Distance','Magnitude']
                if sta_avg:
                    if is_a_receiver_param:
                        H=np.array([np.nanmean(H[z==zi],axis=0) for zi in np.unique(vals)])
                        ZHPS=np.array([np.nanmean(ZHPS[z==zi],axis=0) for zi in np.unique(vals)])
                        ZTF=np.array([np.nanmean(ZTF[z==zi],axis=0) for zi in np.unique(vals)])
                        vals=np.unique(vals)
                        snm=np.unique(snm)
                    else:
                        ibins=np.histogram(vals)[1];bins=np.vstack([dbin(ibins),dbin(ibins[::2])]) #Bins
                        ZTF=np.array(unravel([[np.nanmean(ZTF[(snm==s)&(vals>=b[0])&(vals<=b[1])]) for s in np.unique(snm[(vals>=b[0])&(vals<=b[1])])] for b in bins]))
                        ZHPS=np.array(unravel([[np.nanmean(ZHPS[(snm==s)&(vals>=b[0])&(vals<=b[1])]) for s in np.unique(snm[(vals>=b[0])&(vals<=b[1])])] for b in bins]))
                        H=np.array(unravel([[np.nanmean(H[(snm==s)&(vals>=b[0])&(vals<=b[1])]) for s in np.unique(snm[(vals>=b[0])&(vals<=b[1])])] for b in bins]))
                        _vals=np.array(unravel([[vals[(snm==s)&(vals>=b[0])&(vals<=b[1])].mean() for s in np.unique(snm[(vals>=b[0])&(vals<=b[1])])] for b in bins]))
                        _snm=np.array(unravel([[s for s in np.unique(snm[(vals>=b[0])&(vals<=b[1])])] for b in bins]))
                        vals=_vals;snm=_snm

                ind=np.argsort(vals)
                H=H[ind]
                ZHPS=ZHPS[ind]
                ZTF=ZTF[ind]
                vals=vals[ind]
                snm=snm[ind]
                # ------------------------
                ax=axes[hpshax,bi]
                x=vals;y=H
                scatter(ax,x,y,snm,ColorStandard)
                if named:[ax.text(xx,yy,cat.r.loc[n].Experiment,fontsize=3) for xx,yy,n in zip(x,y,snm)]
                if bi==0:ax.set_ylabel(yttl('HPS.H'))
                # ------------------------
                ax=axes[hpszax,bi]
                x=vals;y=ZHPS
                scatter(ax,x,y,snm,ColorStandard)
                if named:[ax.text(xx,yy,cat.r.loc[n].Experiment,fontsize=3) for xx,yy,n in zip(x,y,snm)]
                if bi==0:ax.set_ylabel(yttl('HPS.Z'))
                # ------------------------
                ax=axes[tfzax,bi]
                x=vals;y=ZTF
                scatter(ax,x,y,snm,ColorStandard)
                if bi==0:ax.set_ylabel(yttl('TF.Z'))
                if stat.YL is None:ax.set_ylim(auto=True)
                else:ax.set_ylim(stat.YL)
                # if named:[ax.text(xx,yy,cat.r.loc[n].Experiment,fontsize=3) for xx,yy,n in zip(x,y,snm)]
                # fig.suptitle('Coherence Average',y=1.05)
                # ------------------------------------------------------------------------
                ax=axes[0,bi];ax.set_title(f'{min(b)} to {max(b)}s')
                if stat.title=='Sigma':ax.set_ylim(0,.4)
                if stat.title=='Sigma':
                    for ax in axes.reshape(-1):ax.yaxis.tick_right();ax.yaxis.set_label_position("right")
                ax=axes.reshape(-1)[-1];ax.set_xlabel(title)
            # if key=='StaDepth':
            #     partitions=np.arange(0,6500,500);alpha=0.25;nbins=100
            #     for ax in axes.reshape(-1):[ax.axvspan(i, i+500, color=custom_cmap(0,nbins=nbins)(i/6000), alpha=alpha, zorder=-10) for i in partitions]
        if save: 
            ##SAVING------------------------------------------------------------------------
            file=f'{key}.{'IG_Sensitive.' if fn=='IG' else 'Regardless.of.IG.'}.{stat.title}.{note}png'
            if named:file='Named.'+file
            save_tight(plotfolder/file,fig,dpi=800)
            plt.close()
