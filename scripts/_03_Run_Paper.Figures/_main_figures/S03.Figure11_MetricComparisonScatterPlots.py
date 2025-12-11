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
# figs value
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',w=None,layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout,width_ratios=np.ones(c) if w is None else w)

# cat value
cat = catalog.copy()
# dirs value
dirs=dir_libraries()



# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure11_MetricComparisons';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'


# tilt bin edges value
tilt_bin_edges = [0.08267542, 0.22943596, 0.50754827, 0.88716358, 0.9808266 ]
# yttl value
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# yttl eta value
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"
# mtrs value
mtrs=['coh','snr']
magwins,markers=[[6,7],[7,8]],['v','^']
# scl value
scl=1.4
# figheight value
figheight=5*scl
# figwidth value
figwidth=3*scl


# phband value
phband={'P':[1,10],'S':[10,30],'Rg':[30,100]};phband_og=phband.copy()
# phcolor value
phcolor={'P':'firebrick','S':'skyblue','Rg':'forestgreen'}
# phases value
phases=['Rg','S','P']
# axes = axes.T
msize = 20
# mthd value
mthd = 'HPS_Z'
# mthd value
mthd = 'TF_Z'
# nonavgalpha value
nonavgalpha=0.05
# mthds value
mthds = ['HPS_Z','TF_Z']
# mthds = ['HPS_Z']
icat=cat.sr.copy()
# usnr value
usnr=unpack_metrics(icat)
# loop over mthds
for mthd in mthds:
    fig,axes=figs(3,1,f=(figwidth,figheight),x=False,y='col')
    # phband value
    phband={'P':[1,10],'S':[1,10],'Rg':[1,10]}
    # phband value
    phband={'P':[10,30],'S':[10,30],'Rg':[10,30]}
    # phband value
    phband={'P':[30,100],'S':[30,100],'Rg':[30,100]}

    # mtrs value
    mtrs = [['LT','ST']]
    # names value
    names={'LT':'Noise'.lower(),'ST':'Signal'.lower(),'snr':'snr'}
    # bands value
    bands = [[1,10],[10,30],[30,100]]
    # resx value
    resx=[];resy=[];rcoh=[]
    # fn='IG'
    fn=None
    for axi,(ax,b) in enumerate(zip(axes,bands)):
        # prezx value
        prezx=[];prezy=[];pcoh=[]
        # phband value
        phband={'P':b,'S':b,'Rg':b}
        # phband og value
        phband_og=phband.copy()
        # loop over mtrs
        for mtr in mtrs:
            # mrezx value
            mrezx=[];mrezy=[];mcoh=[]
            for mg,marker in zip(magwins,markers):
                # idx value
                idx=icat.Magnitude.between(mg[0],mg[1],inclusive='both')
                # sn value
                sn = icat.StaName[idx]
                # xmtr value
                xmtr='LT'
                # ymtr value
                ymtr='ST'
                xphaseinvariant,yphaseinvariant=xmtr in ['LT','coh'],ymtr in ['LT','coh']
                xratiod,yratiod=not (xmtr=='coh'),not (ymtr=='coh')
                if not xphaseinvariant:x=[usnr.__dict__[xmtr].__dict__[mthd].R().__dict__[ph].Average(phband[ph],fn=fn) for ph in phases]
                else:
                    if xratiod:x=[usnr.__dict__[xmtr].__dict__[mthd].R().Average(phband[ph],fn=fn) for ph in phases]
                    else:x=[usnr.__dict__[xmtr].__dict__[mthd].Average(phband[ph],fn=fn) for ph in phases]
                if not yphaseinvariant:y=[usnr.__dict__[ymtr].__dict__[mthd].R().__dict__[ph].Average(phband[ph],fn=fn) for ph in phases]
                else:
                    if yratiod:y=[usnr.__dict__[ymtr].__dict__[mthd].R().Average(phband[ph],fn=fn) for ph in phases]
                    else:y=[usnr.__dict__[ymtr].__dict__[mthd].Average(phband[ph],fn=fn) for ph in phases]
                # coh value
                coh=[usnr.__dict__['coh'].__dict__[mthd].Average(phband_og[ph],fn=fn) for ph in phases]
                # coh value
                coh=[c[idx] for c in coh]
                x=[h[idx] for h in x]
                y=[t[idx] for t in y]
                for avg,alpha in zip([False,True],[nonavgalpha,.7]):
                    if not avg:continue
                    if avg:
                        # pl coh value
                        pl_coh = [[np.nanmean(c[sn==s]) for s in np.unique(sn)] for c in coh]
                        # pl y value
                        pl_y=[[np.nanmean(t[sn==s]) for s in np.unique(sn)] for t in y]
                        pl_x=[[np.nanmean(h[sn==s]) for s in np.unique(sn)] for h in x]
                    else:pl_x=x.copy();pl_y=y.copy();pl_coh=coh.copy()
                    pl_x=[10**np.array(xx) for xx in pl_x]
                    pl_y=[10**np.array(yy) for yy in pl_y] #R() computes the log ratio, take the log to return to linear ratio
                    pl_x=[xx.clip(0,1) for xx in pl_x]
                    pl_y=[yy.clip(0,1) for yy in pl_y]
                    [ax.scatter(xx,yy,marker=marker,s=msize,ec='k',lw=0.2,c=phcolor[ph],alpha=alpha) for xx,yy,ph in zip(pl_x,pl_y,phases)]
                    xlabel=rf'$\Delta$ {names[xmtr]}';ylabel=rf'$\Delta$ {names[ymtr]}'
                    if axi==2:ax.set_xlabel(f'{xlabel} ({mthd.split('_')[0]})')
                    if axi==1:ax.set_ylabel(f'{ylabel} ({mthd.split('_')[0]})')
                    # ax.axhline(0,alpha=0.4,ls=':',c='k',zorder=-1e3,lw=0.9);ax.axvline(0,alpha=0.4,ls=':',c='k',zorder=-1e3,lw=0.9)
                mrezx.append(pl_x);mrezy.append(pl_y);mcoh.append(pl_coh)
            prezx.append(mrezx);prezy.append(mrezy);pcoh.append(mcoh)
            if mtr=='coh':lims=[0,1.03]
        resx.append(prezx);resy.append(prezy);rcoh.append(pcoh)
        ax.tick_params(size=7,labelsize=7)

    resx=np.array(resx).squeeze()
    resy=np.array(resy).squeeze()
    rcoh=np.array(rcoh).squeeze()

    # resx=np.array(resx).squeeze()
    # resy=np.array(resy).squeeze()
    # rcoh=np.array(rcoh).squeeze()

    for ax,b in zip(axes,bands):
        lims=[np.min([ax.get_xlim(),ax.get_ylim()]),np.max([ax.get_xlim(),ax.get_ylim()])]
        dx=0.0000
        lims=[lims[0]-dx,lims[1]+dx]
        ax.set_xlim(lims);ax.set_ylim(lims)
        ax.plot(lims,lims,color='k',ls='-',lw=0.5,alpha=0.5,zorder=-1e4)
        # ax.tick_params(size=7,labelsize=7)
        hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c=phcolor[ph],label=f'{ph.replace('Rg','Rayleigh')} ({b[0]}-{b[1]}s)',lw=0.2) for ph in ['P','S','Rg']]
        ax.legend(handles=hdls,loc='lower right',borderaxespad=0,frameon=False)



    fold = plotfolder
    
    file = f'_01_{mthd}.Signal.Nose.Reduction.{save_format}'
    _=save_tight(fold/file,fig,dpi=700)


    # # ====================================================================================================
    # # ====================================================================================================
    # # ====================================================================================================
    # # ====================================================================================================
    # # ====================================================================================================
    # # ====================================================================================================


    magwins=[[6,7],[7,8]]
    phases=['Rg','S','P']
    phband={'P':[1,10],'S':[10,30],'Rg':[30,100]}
    bands = [[1,10],[10,30],[30,100]]
    pearson = lambda mgx,mgy:np.corrcoef(mgy,mgx)[0,1]
    fig,axes=figs(3,1,f=(5,figheight),x=False,y=False)
    ax=axes
    coh=rcoh
    for bi,(b,rx,ry,by,ax) in enumerate(zip(bands,resx,resy,coh,axes)):
        for mg,mmx,mmy,mgy in zip(magwins,rx,ry,by):
            for ph,ppx,ppy,py in zip(phases,mmx,mmy,mgy):
                # x=np.nanmean(b)
                # x=bi+1
                # y=pearson(py,px)
                x=np.array(ppy)/np.array(ppx);y=np.array(py)
                ax.scatter(x,y,c=phcolor[ph],marker='^' if max(mg)==8.0 else 'v',ec='k',lw=0.5)
        ax.axvline(1.0,ls=':',c='k',lw=0.9)
    axes[1].set_ylabel(yttl(mthd.replace('_',' ')))
    axes[-1].set_xlabel(f'{ylabel} / {xlabel} ({mthd.split('_')[0]})')
    for ax,b in zip(axes,bands):
        ax.grid(alpha=0.3,zorder=-1e3)
        ax.tick_params(size=7,labelsize=7)
        hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c=phcolor[ph],label=f'{ph.replace('Rg','Rayleigh')} ({b[0]}-{b[1]}s)',lw=0.2) for ph in ['P','S','Rg']]
        ax.legend(handles=hdls,loc='best',borderaxespad=0,frameon=False)


    fold = plotfolder
    
    file = f'_02_{mthd}.Coherence.with.Delta.Ratio.{save_format}'
    _=save_tight(fold/file,fig,dpi=700)

    waterdepth=np.array([cat.r.aloc[s].iloc[0].StaDepth for s in np.unique(sn)])
    for bi,(b,rx,ry,by,ax) in enumerate(zip(bands,resx,resy,coh,axes)):
        for mg,mmx,mmy,mgy in zip(magwins,rx,ry,by):
            for ph,ppx,ppy,py in zip(phases,mmx,mmy,mgy):
                # x=np.nanmean(b)
                # x=bi+1
                # y=pearson(py,px)
                x=np.array(ppy)/np.array(ppx)
                # y=np.array(py)
                y=waterdepth
                ax.scatter(x,y,c=phcolor[ph],marker='^' if max(mg)==8.0 else 'v',ec='k',lw=0.5)
        ax.axvline(1.0,ls=':',c='k',lw=0.9)
    axes[1].set_ylabel('Water depth, m')
    axes[-1].set_xlabel(f'{ylabel} / {xlabel} ({mthd.split('_')[0]})')
    for ax,b in zip(axes,bands):
        ax.grid(alpha=0.3,zorder=-1e3)
        ax.tick_params(size=7,labelsize=7)
        hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c=phcolor[ph],label=f'{ph.replace('Rg','Rayleigh')} ({b[0]}-{b[1]}s)',lw=0.2) for ph in ['P','S','Rg']]
        ax.legend(handles=hdls,loc='best',borderaxespad=0,frameon=False)


    fold = plotfolder
    
    file = f'_06_{mthd}.WaterDepth.with.Delta.Ratio.{save_format}'
    _=save_tight(fold/file,fig,dpi=700)



    for wi,w in enumerate(['Noise'.lower(),'Signal'.lower()]):
        magwins=[[6,7],[7,8]]
        phases=['Rg','S','P']
        phband={'P':[1,10],'S':[10,30],'Rg':[30,100]}
        bands = [[1,10],[10,30],[30,100]]
        pearson = lambda mgx,mgy:np.corrcoef(mgy,mgx)[0,1]
        fig,axes=figs(3,1,f=(5,figheight),x=False,y=False)
        ax=axes
        coh=rcoh
        for bi,(b,rx,ry,by,ax) in enumerate(zip(bands,resx,resy,coh,axes)):
            for mg,mmx,mmy,mgy in zip(magwins,rx,ry,by):
                for ph,ppx,ppy,py in zip(phases,mmx,mmy,mgy):
                    # x=np.nanmean(b)
                    # x=bi+1
                    # y=pearson(py,px)
                    if w=='Signal'.lower():x=np.array(ppy)
                    else:x=np.array(ppx)
                    y=np.array(py)
                    ax.scatter(x,y,c=phcolor[ph] if w=='Signal'.lower() else 'gainsboro',marker='^' if max(mg)==8.0 else 'v',ec='k',lw=0.5)
            ax.axvline(1.0,ls=':',c='k',lw=0.9)
            ax.grid(alpha=0.3,zorder=-1e3);ax.tick_params(size=7,labelsize=7)
        axes[1].set_ylabel(yttl(mthd.replace('_',' ')))
        axes[-1].set_xlabel(f'{ylabel if w=='Signal'.lower() else xlabel} ({mthd.split('_')[0]})')
        ax.tick_params(size=7,labelsize=7)
        for ax,b in zip(axes,bands):
            if w=='Signal'.lower():hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c=phcolor[ph],label=f'{ph.replace('Rg','Rayleigh')} ({b[0]}-{b[1]}s)',lw=0.2) for ph in ['P','S','Rg']]
            else:hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c='gainsboro',label=f'{b[0]}-{b[1]}s',lw=0.2)]
            ax.legend(handles=hdls,loc='best',borderaxespad=0,frameon=False)


        fold = plotfolder
        
        file = f'_{str(3+wi).zfill(2)}_{mthd}.Coherence.with.Delta.{w}.{save_format}'
        _=save_tight(fold/file,fig,dpi=700)


    for wi,w in enumerate(['Noise'.lower(),'Signal'.lower()]):
        magwins=[[6,7],[7,8]]
        phases=['Rg','S','P']
        phband={'P':[1,10],'S':[10,30],'Rg':[30,100]}
        bands = [[1,10],[10,30],[30,100]]
        pearson = lambda mgx,mgy:np.corrcoef(mgy,mgx)[0,1]
        fig,axes=figs(3,1,f=(5,figheight),x=False,y=False)
        ax=axes
        coh=rcoh
        for bi,(b,rx,ry,by,ax) in enumerate(zip(bands,resx,resy,coh,axes)):
            for mg,mmx,mmy,mgy in zip(magwins,rx,ry,by):
                for ph,ppx,ppy,py in zip(phases,mmx,mmy,mgy):
                    # x=np.nanmean(b)
                    # x=bi+1
                    # y=pearson(py,px)
                    if w=='Signal'.lower():x=np.array(ppy)
                    else:x=np.array(ppx)
                    # y=np.array(py)
                    y=waterdepth
                    ax.scatter(x,y,c=phcolor[ph] if w=='Signal'.lower() else 'gainsboro',marker='^' if max(mg)==8.0 else 'v',ec='k',lw=0.5)
            ax.axvline(1.0,ls=':',c='k',lw=0.9)
            ax.grid(alpha=0.3,zorder=-1e3);ax.tick_params(size=7,labelsize=7)
        axes[1].set_ylabel('Water depth, m')
        axes[-1].set_xlabel(f'{ylabel if w=='Signal'.lower() else xlabel} ({mthd.split('_')[0]})')
        ax.tick_params(size=7,labelsize=7)
        for ax,b in zip(axes,bands):
            if w=='Signal'.lower():hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c=phcolor[ph],label=f'{ph.replace('Rg','Rayleigh')} ({b[0]}-{b[1]}s)',lw=0.2) for ph in ['P','S','Rg']]
            else:hdls=[ax.scatter(np.nan,np.nan,marker='s',ec='k',c='gainsboro',label=f'{b[0]}-{b[1]}s',lw=0.2)]
            ax.legend(handles=hdls,loc='best',borderaxespad=0,frameon=False)


        fold = plotfolder
        
        file = f'_{str(7+wi).zfill(2)}_{mthd}.WaterDepth.with.Delta.{w}.{save_format}'
        _=save_tight(fold/file,fig,dpi=700)




    fig,axes=figs(1,2,f=(3,3))
    ax=axes
    coh=rcoh
    for axi,ax in enumerate(axes):
        for bi,(b,rx,ry,by) in enumerate(zip(bands,resx,resy,coh)):
            for mgi,(mg,mmx,mmy,mgy) in enumerate(zip(magwins,rx,ry,by)):
                for pi,(ph,ppx,ppy,py) in enumerate(zip(phases,mmx,mmy,mgy)):
                    # x=np.array(ppy)/np.array(ppx)
                    x=np.array(py)
                    if axi==0:y=np.array(ppx)
                    else:y=np.array(ppy)
                    y=pearson(x,y)
                    x=bi+1
                    ax.scatter(x,y,s=100/(pi+1),alpha=0.3,c=phcolor[ph],marker='^' if max(mg)==8.0 else 'v',ec='k',lw=0.5,zorder=1/(100*(pi+1)))

        if axi==1:
            ax.tick_params(axis='y', which='both', labelright=True, labelleft=False,
            right=True, left=False)
            ax.set_ylabel(yttl(mthd.replace('_',' ')))
            ax.tick_params(size=7,labelsize=7)
            ax.set_xlim(0.5,3.5)
            ax.set_xticks([1.5,2.5])
            ax.set_xticklabels([10,30])
            ax.axhline(0,ls=':',c='k')
            ax.yaxis.set_label_position('right')   # then ax.set_ylabel('Your label')
            ax.set_ylabel(rf'{yttl(mthd.replace('_',' '))} correlation with $\Delta$ Signal')
            # ax.set_ylabel('Correlation coefficient')
            # ax.set_xlabel('Period band')
        else:
            # ax.tick_params(axis='y', which='both', labelright=True, labelleft=False,
            # right=True, left=False)
            ax.set_ylabel(yttl(mthd.replace('_',' ')))
            ax.tick_params(size=7,labelsize=7)
            ax.set_xlim(0.5,3.5)
            ax.set_xticks([1.5,2.5])
            ax.set_xticklabels([10,30])
            ax.axhline(0,ls=':',c='k')
            # ax.yaxis.set_label_position('right')   # then ax.set_ylabel('Your label')
            ax.set_ylabel(rf'{yttl(mthd.replace('_',' '))} correlation with $\Delta$ Noise')
            # ax.set_ylabel('Correlation coefficient')
            spaces='\;'*7
            ax.set_xlabel(rf'${spaces}$Period band, s',ha='left')

    fold = plotfolder
    
    file = f'S03.Figure11_{mthd}.Correlation.{save_format}'
    _=save_tight(fold/file,fig,dpi=700)
