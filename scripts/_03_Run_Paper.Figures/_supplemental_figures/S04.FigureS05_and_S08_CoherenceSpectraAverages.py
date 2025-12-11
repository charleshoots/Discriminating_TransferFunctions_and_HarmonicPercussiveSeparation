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

import matplotlib as mpl;from matplotlib.colors import ListedColormap
import time;start=time.time()

# octavg value
octavg=lt.math.octave_average;transpose=False
# opts value
opts=AttribDict({})

# --------------------------------------------------------------------------------
# ------------------------------------CODE----------------------------------------
# --------------------------------------------------------------------------------
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_supplemental_figures'/'FigureS5_CoherenceSpectraAverages';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'

# width value
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


# req = 'Instrument_Design'
# figsize=(6,1);transpose=True #width,height
# sets=['AB','AR','B2','BA','BG','KE','TRM']
# opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
# # ####### --------------------------------------------------------------------------------
# req = 'Seismometer'
# figsize=(6,1);transpose=True #width,height
# sets=['Guralp CMG3T 120','Trillium 240','Trillium Compact']
# opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
# ####### --------------------------------------------------------------------------------
# req = 'Magnitude'
# figsize=(6,1);transpose=True #width,height
# sets=[[6.0,7.0],[7.0,8.01]]
# opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})
# ####### --------------------------------------------------------------------------------
# ####### --------------------------------------------------------------------------------
# req = 'Pressure_Gauge'
# figsize=(6,1);transpose=True #width,height
# sets=['DPG','APG']
# opts.update({req:AttribDict({'req':req,'figsize':figsize,'sets':sets,'transpose':transpose})})

# stats value
stats=[]
# stat value
stat=AttribDict();stat.func=np.mean;stat.title='Mu';stats.append(stat)
# stat value
stat=AttribDict();stat.func=np.std;stat.title='Sigma';stats.append(stat)

# octave av value
octave_av=True
# notched value
notched=False
# plot notch value
plot_notch=False
# ________________________________________________________________________________________________________________________________________________________________
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
cat=catalog.sr.copy()
# usnr value
usnr=unpack_metrics(cat)

# si value
si=0;s = cat.iloc[si];f=s.Data.Coherence().f;foct=octavg(s.Data.Coherence().ATaCR.zp_21.coh,f)[0]
# zmin value
zmin=0;zmax=6000;binspacing=500
# bins value
bins = np.array([[i,i+binspacing] for i in np.arange(zmin,zmax,binspacing)])
# bins = np.flip(bins,axis=0)
boundaries=np.unique(bins);binedges=np.array([max(b) for b in bins]);nbins=len(bins)
# function custom cmap
def custom_cmap(ind=0,nbins=5):
    if ind==0:cmap = cm.cmaps['glasgow'].resampled(nbins+3)
    if ind==1:cmap = cm.cmaps['batlow'].resampled(nbins+3)
    return cmap
# cmap value
cmap=custom_cmap(1,nbins).reversed()
# cmap value
cmap=custom_cmap(0,nbins).reversed()
# darker value
darker=1.25;clrs=cmap(bins.mean(axis=1)/(bins.max()*darker))
# clrs=np.flip(clrs,axis=0) #colors
cmap = ListedColormap(clrs)
# ----------------------------------------------------------------------------------------------------



# reqs value
reqs=list(opts.keys())
# status value
status = lambda:f'R:{ri+1}/{len(reqs)}  |-|  S:{set_i+1}/{len(sets)}  |:|  {req}:{set}'
for ri,req in enumerate(reqs):
    # figsize value
    figsize=opts[req].figsize
    # sets value
    sets=opts[req].sets
    # transpose value
    transpose=opts[req].transpose
    # loop over stats
    for stat in stats:
        for set_i,set in enumerate(sets):
            print(status()+' | Start')
            icat = cat.copy()
            if req==None:pass
            else:
                if isinstance(set,list):icat=icat[(icat[req]>=min(set))&(icat[req]<(max(set)))].copy()
                else:icat=icat[icat[req]==set].copy()
            # -----
            # DATA ------------------------------------------------------------------------------------------

            ind=0;DAT=[0]
            DAT = AttribDict()
            f=1/usnr.coh.bands
            DAT.TF = usnr.coh.TF_Z.D
            DAT.HPS_Z = usnr.coh.HPS_Z.D
            DAT.HPS_H = np.array([usnr.coh.HPS_1.D,usnr.coh.HPS_2.D]).mean(axis=0) #Mean H1/H2 together
            if octave_av:
                DAT.TF = octavg(DAT.TF,f)[1]
                DAT.HPS_Z = octavg(DAT.HPS_Z,f)[1]
                DAT.HPS_H = octavg(DAT.HPS_H,f)[1]
                xf=foct
            else:xf=f
            faxis=(xf>0)&(xf<1)

            iDAT_std = DAT.copy();Z=icat.StaDepth
            iDAT_std.TF = np.array([np.std(iDAT_std.TF[(Z>b0)&(Z<=b1),:],axis=0) for b0,b1 in bins])
            iDAT_std.HPS_Z = np.array([np.std(iDAT_std.HPS_Z[(Z>b0)&(Z<=b1),:],axis=0) for b0,b1 in bins])
            iDAT_std.HPS_H = np.array([np.std(iDAT_std.HPS_H[(Z>b0)&(Z<=b1),:],axis=0) for b0,b1 in bins])

            iDAT = DAT.copy();Z=icat.StaDepth
            iDAT.TF = np.array([stat.func(iDAT.TF[(Z>b0)&(Z<=b1),:],axis=0) for b0,b1 in bins])
            iDAT.HPS_Z = np.array([stat.func(iDAT.HPS_Z[(Z>b0)&(Z<=b1),:],axis=0) for b0,b1 in bins])
            iDAT.HPS_H = np.array([stat.func(iDAT.HPS_H[(Z>b0)&(Z<=b1),:],axis=0) for b0,b1 in bins])

            print(status()+' | Data ')
            # ------------------------------------------------------------------------------------------
            # PLOT ------------------------------------------------------------------------------------------
            flim=[1/300,.97] #hz
            sets_ttl=['TF.Z','HPS.Z','HPS.H']
            linewidth=1.5
            # plot_sigma=False #Variance, unlike sigma, does not have the same units as the random variable (coherence).
            print(status()+' | Plot ')
            if transpose:fig,axes=plt.subplots(nrows=1,ncols=3,figsize=figsize,sharex='all',sharey='all')
            else:fig,axes=plt.subplots(nrows=3,ncols=1,figsize=figsize,sharex='all',sharey='all')
            for axi,(ax,ttl,d,sig) in enumerate(zip(axes,sets_ttl,[iDAT.TF,iDAT.HPS_Z,iDAT.HPS_H],[iDAT_std.TF,iDAT_std.HPS_Z,iDAT_std.HPS_H])):
                # if plot_sigma: #------#------#------#------#------#------#------#------
                #     C=lambda b,xf,k=2:(k*b/(1/xf))
                #     power=1;sig=sig**power;k=.5
                #     _=[ax.plot(C(b.mean(),xf)[xf<fnotch(b.mean(),k)],
                #     y[xf<fnotch(b.mean(),k)] ,c=c,lw=linewidth) 
                #     for c,y,b in zip(clrs,d,bins)]
                #     [ax.fill_between(C(b.mean(),xf)[xf<fnotch(b.mean(),k)] ,
                #     (y-isig)[xf<fnotch(b.mean(),k)] ,(y+isig)[xf<fnotch(b.mean(),k)],
                #     color=c, alpha=0.15, lw=0,zorder=-1e4) 
                #     for c,y,isig,b in zip(clrs,d,sig,bins)]

                if True: #------#------#------#------#------#------#------#------
                    for di,(c,y,b) in enumerate(zip(clrs,d,bins)):
                        kw = {'s':10,'facecolor':c,'c':c,'marker':'s','linewidth':0.1}
                        if plot_notch: 
                            zfn=cat.StaDepth[(cat.StaDepth>=min(b))&(cat.StaDepth<=max(b))].mean();fn=fnotch(zfn)
                            # ax.scatter( xf[xf<=fn][-1] ,y[xf<=fn][-1],color=c,linewidth=10,marker='_',s=25,)
                            ax.axvline(xf[xf<=fn][-1],color=c,linewidth=2,zorder=-1000,alpha=1)
                        ax.scatter(xf,y,**kw,)
                    ax.set_xlim([1/200,1]);ax.set_xscale('log')
                ttl = fr"$\underset{{{ttl}}}{{\gamma{'\;\;'*len(ttl)}}}$"
                if transpose:ax.set_title(ttl,y=1.2)
                else:ax.set_ylabel(ttl)
                ax.grid('major',alpha=0.3,zorder=-1e5)
            ax=axes[-1]
            # ax.set_xlim([1/300,1])
            if transpose:axes[1].set_xlabel('period, s')
            if not transpose:ax.set_xlabel('period, s')
            periods=[100,30,10,1]
            _=ax.set_xticks([1/t for t in periods]);ax.set_xticklabels([t for t in periods])
            # norm = mpl.colors.Normalize(vmin=np.min(bins), vmax=np.max(bins))
            norm = mpl.colors.BoundaryNorm(boundaries,cmap.N+1)
            sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm);sm.set_array([])
            fig.subplots_adjust(right=0.85,top=0.85)  # Make space on the right
            cbar_ax = fig.add_axes([
            0.88,0.1 if transpose else 0.13, 
            0.03,0.8 if transpose else 0.74])  # [left, bottom, width, height]

            # cb = fig.colorbar(sm, cax=cbar_ax)
            cb = fig.colorbar(sm, cax=cbar_ax, boundaries=boundaries, ticks=boundaries[::2])
            cb.set_label('water depth, m', fontweight='bold')
            cb.ax.tick_params(labelsize=10)
            cb.ax.invert_yaxis()
            _=cb.ax.set_yticks(boundaries[::2])
            ax.set_ylim(0,1 if stat.title=='Mu' else 0.5)
            ax.set_yticks(np.unique(np.array([0,.5,1.0])/(1 if stat.title=='Mu' else 2)))
            if transpose:
                if isinstance(set,list):axes[0].set_ylabel(f'Mw {(float(int(min(set))))} - {(float(int(max(set))))}')
                else:axes[0].set_ylabel(set)
            else:
                if isinstance(set,list):axes[0].set_title(f'Mw {(float(int(min(set))))} - {(float(int(max(set))))}')
                else:axes[0].set_title(set)
            # ------------------------------------------------------------------------------------------
            # options = '.'.join([('notched' if notched else 'not.notched'),('octav' if octave_av else 'not.octav')])
            options='.'.join([('octav' if octave_av else 'not.octav'),''])
            if isinstance(set,list):file=f'{req}.{min(set)}.to.{max(set)}.CoherenceSpectrabyDepth_{options}.{save_format}'.replace('..','.')
            else:
                if req==None:file=f'AllData.{stat.title}CoherenceSpectrabyDepth_{options}.{save_format}'.replace('..','.')
                else:file=f'{req}.{stat.title}.{set.replace(' ','.')}.CoherenceSpectrabyDepth_{options}.{save_format}'.replace('..','.')
            fold=plotfolder
            ax.set_xlim([ax.get_xticks().min(),ax.get_xticks().max()])
            if req is not None:fold=fold/req
            fold.mkdir(parents=True, exist_ok=True)
            save_tight(fold/('S04.FigureS05S08_'+file),fig=fig,dpi=700) #dpi=700,
            print(status()+' | ----Completed----')
            plt.close()

end=time.time()
print(f"Elapsed time: {(end - start)/60:.2f} minutes")