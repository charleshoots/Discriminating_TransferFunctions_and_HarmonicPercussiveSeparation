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

import statsmodels.api as sm

# dirs value
dirs = io.dir_libraries()
# icat value
icat=catalog.sr.copy()
# usnr value
usnr=unpack_metrics(icat)
# cat value
cat = catalog.copy()


# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_supplemental_figures'/'FigureS14_NarrowSymmetryAnalysisPlots';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'
# Noise Spectra
f=cat.r.iloc[0].Data.Noise.Averaged().f
# faxis value
faxis=(f>0)&(f<=1)
# f value
f=f[faxis]
# noise f value
noise_f=f
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
figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
from obspy.signal.trigger import classic_sta_lta,carl_sta_trig,recursive_sta_lta
# stalta methods value
stalta_methods={'classic':classic_sta_lta,'carl':carl_sta_trig,'recurssive':recursive_sta_lta}
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
# function color master
def color_master(key,sets=np.nan):
    if isinstance(sets,int):sets=np.arange(0,sets,1)
    elif (np.shape(np.atleast_1d(np.array(sets)))[0]-1)==0:sets=np.arange(0,100,1)
    #reverses the colormap of any key in this variable, regardless of other conditions.
    color_overrides=['Pressure_Gauge','Magnitude','Sediment_Thickness_m']
    #Luminnance model, on a scale from 0 (black) to 1 (white)
    darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)],name=cmap.name).resampled(100)
    # luminance value
    luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
    # cmap value
    cmap = darken(cm.cmaps['bamako'].copy(),frac=.8)
    if key=='Instrument_Design':cmap=ListedColormap([ColorStandard.instrument[s] for s in sets], name='custom_cmap')
    elif key=='Network':cmap=ListedColormap([ColorStandard.network[s] for s in sets], name='custom_cmap')
    elif key=='Magnitude':cmap=darken(cm.cmaps['lipari'],.8)
    elif key=='StaDepth':cmap=cm.cmaps['glasgow'].reversed()
    elif key=='Distance':cmap=darken(cm.cmaps['tokyo'].copy())
    elif key=='Sediment_Thickness_m':cmap=darken(cm.cmaps['lapaz'].copy(),.9).reversed()
    elif key=='Coherence':cmap=darken(cm.cmaps['lajolla'].copy(),.9)
    else:cmap=darken(darken(cm.cmaps['glasgow'].resampled(100),0.7).reversed(),0.9)
    # cmap value
    cmap=cmap.resampled(len(sets))
    if (luminance(cmap(0))<luminance(cmap(1e3)))&(not isinstance(sets[0],str)):cmap=cmap.reversed()
    if key in color_overrides:cmap=cmap.reversed()
    # color value
    color=[cmap(si/len(sets)) for si,_ in enumerate(sets)]
    #suggested zorder based on luminance, will put the darker colors in the back
    zorder = np.array([luminance(c) for c in color])
    return color,zorder,cmap
# function gbounds
def gbounds(ix,iy,xscl=0.04,yscl=0.02):
    # ix value
    ix=np.array(ix);iy=np.array(iy)
    # yb value
    yb=iy;yl=np.array([np.min(yb),np.max(yb)])
    # xb value
    xb=(ix);yl=np.array([np.min(xb),np.max(xb)])
    xl,yl=np.array([np.min(xb),np.max(xb)]),np.array([np.min(yb),np.max(yb)])
    # xl value
    xl=xl+np.array([-abs(np.diff(xl))*xscl,abs(np.diff(xl))*xscl]).reshape(-1)
    # yl value
    yl=yl+np.array([-abs(np.diff(yl))*yscl,abs(np.diff(yl))*yscl]).reshape(-1)
    return xl,yl
# function fig setup
def fig_setup(figsize=(6, 4),width_ratios=[.4, .5,.5, .4],height_ratios=[.2, .5],debug=False):
    fig=plt.figure(figsize=figsize)
    # gs value
    gs = gridspec.GridSpec(nrows=2, ncols=4, width_ratios=width_ratios, height_ratios=height_ratios)
    # ax top center left value
    ax_top_center_left = fig.add_subplot(gs[0, 1])
    # ax top center right value
    ax_top_center_right = fig.add_subplot(gs[0, 2],sharey=ax_top_center_left)
    # ax bottom center left value
    ax_bottom_center_left = fig.add_subplot(gs[1, 1],sharex=ax_top_center_left)
    # ax bottom center right value
    ax_bottom_center_right = fig.add_subplot(gs[1, 2],sharex=ax_top_center_right)#,sharey=ax_bottom_center_left)
    # ax bottom left value
    ax_bottom_left = fig.add_subplot(gs[1, 0])#,sharey=ax_bottom_center_left)
    # ax bottom right value
    ax_bottom_right = fig.add_subplot(gs[1, 3])#,sharey=ax_bottom_left)
    if debug:
        for ax, label in zip(
            [ax_top_center_left, ax_top_center_right,
            ax_bottom_left, ax_bottom_center_left, ax_bottom_center_right, ax_bottom_right],
            ['Top Center Left', 'Top Center Right', 'Bottom Left', 'Bottom Center Left', 'Bottom Center Right', 'Bottom Right']):
            ax.set_title(label)
    for ax in [ax_top_center_left, ax_top_center_right]: #, ax_bottom_left, ax_bottom_right]:
        ax.tick_params(labelbottom=False)
    for ax in [ax_bottom_center_left,ax_top_center_right, ax_bottom_right, ax_bottom_center_right]:
        ax.tick_params(labelleft=False)
    # left value
    left=AttribDict({});left.top=ax_top_center_left;left.side=ax_bottom_left;left.center=ax_bottom_center_left
    # right value
    right=AttribDict({});right.top=ax_top_center_right;right.side=ax_bottom_right;right.center=ax_bottom_center_right
    left.name='left';right.name='right'
    fig.subplots_adjust(wspace=.35,hspace=0.1)
    return fig,left,right
# centers value
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
# dbin value
dbin=lambda x:np.array([x[:-1],x[1:]]).T
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Constants and helper functions, don't touch these.
bands=['1_10','10_30','30_100'];phases=['Rg','P','S']
methods=['ATaCR','NoiseCut'];mnames={'ATaCR':'TF_Z','NoiseCut':'HPS_Z'}
f=catalog.sr.iloc[0].Data.Coherence().f
mu_coh = lambda b,method,s: s.Coherence[method].reshape(-1)[((1/f)>=min(b))&((1/f)<=max(b))].mean()
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ------------------------------| OPTIONS |-----------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# Code takes about 90 seconds to save 
# 81 different plots covering the following categories. 
# ---------------------------------------------------------------------------------------
meta_wins=AttribDict() #don't touch this one 
# Comment out ones you don't want to use


# meta_wins.Sediment_Thickness_m=[[i,i+1500] for i in np.arange(0,7320+500,1500)][:-1]
# meta_wins.Network=['2D','7A','7D','X9','XF','YL','YO','ZA','ZN']
# meta_wins.Environment=['North Atlantic', 'North Pacific', 'Solomon Sea', 'South Pacific']
# -----
meta_wins.Magnitude=[[i,i+.5] for i in np.arange(6,8+.5,.5)]
meta_wins.Seismometer=['Guralp CMG3T 120', 'Trillium 240', 'Trillium Compact']
meta_wins.Instrument_Design=['AB', 'AR', 'B2', 'BA', 'BG', 'KE', 'TRM']
meta_wins.Pressure_Gauge=['APG', 'DPG']
meta_wins.StaDepth=np.array([[i,i+500] for i in np.arange(0,6500,500)])
meta_wins.Distance=dbin(np.arange(30,180+30,30))


color_overrides=['Pressure_Gauge','Magnitude','Sediment_Thickness_m'] #reverses the colormap of any key in this variable, regardless of other conditions.
lw=1.4 #histogram linewidth
sz=3 #scatter plot size
NSTD=-1 #Masks high variance data from plot, set to -1 (default) to disable. Will not plot data that exceeds NSTD standard deviations from the mean.



RATIO = ratio = True


x_nbins=20;y_nbins=20 #number of bins for histograms
binrange={'1_10':[-1,3],'10_30':[-1,5],'30_100':[-1,6]}
if ratio:binrange={'1_10':[-0.5,0.5],'10_30':[-1,1],'30_100':[-2,5]}

symmetry=True #forces the left side of the plot to reverse x-axes. visually creating a mirror effect.



# -------------------------
filtertype='acausul';note='V04';fold=dirs.Data/'SNR_Models';file =f'SNR_{filtertype}.filter_{note}.pkl'
SNR=load_pickle(fold/file)
SR=catalog.sr.copy();SR.sort_values(by=['Name','StaName'],inplace=True)
SNR.sort_values(by=['Name','StaName'],inplace=True)
SNR = SNR[np.isin(np.array(SNR.Name + '-' + SNR.StaName),np.array(SR.Name + '-' + SR.StaName))]
assert sum((SR.Name==SNR.Name)&(SR.StaName==SNR.StaName))==len(SNR), 'failed start sets'
# snr=unpacksnr(SNR.copy(),methods=['ATaCR','NoiseCut','Original'])
usnr=unpack_metrics(icat)
# -------------------------
# icoh=meancoh(icat)
# isnr=unpacksnr(icat,ratio=ratio)
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
status=lambda:print(f'[{b}] [ {ki+1}/{len(meta_wins.keys())} ] {key} ')
cumulative=False;density=False;stacked=True #Makes a Histogram

scen=None
for ki,key in enumerate(meta_wins.keys()):
    sets=meta_wins[key]
    color,zorder,cmap = color_master(key,sets)
    for bi,b in enumerate(bands):
        status()
        brange=binrange[b]
        for pi,p in enumerate(phases):
            fig,left,right=fig_setup()
            ysides=[right.side,left.side,left.center,right.center] # axes that must always share the same y-lims
            xsides=[left.center,left.top,right.center,right.top] # axes that must always share the same x-lims
            left_xbins=[];left_ybins=[];right_xbins=[];right_ybins=[]
            for si,set in enumerate(sets):
                if isinstance(set,str):ind=icat[key]==set;zc=0 
                else:ind=(icat[key]>=min(set))&(icat[key]<=max(set));zc=-10*min(set)
                if sum(ind)==0:continue
                for mi,(m,caxes,xb,yb) in enumerate(zip(methods,[left,right],[left_xbins,right_xbins],[left_ybins,right_ybins])):
                    filt = list(map(int,b.split('_')))
                    x=usnr.coh.__dict__[mnames[m]].Average(filt)
                    y=usnr.snr.__dict__[mnames[m]].R()[p].Average(filt)
                    x,y=x[ind].copy(),y[ind].copy()
                    nind=y!=None;x,y=x[nind].copy(),y[nind].copy()
                    if ((len(x)==0) or (len(y)==0)):continue

                    if NSTD>0:
                        cutoff=np.array([(y.mean()-NSTD*y.std()),(y.mean()+NSTD*y.std())])
                        cutaxis=(y>=cutoff.min())&(y<=cutoff.max())
                        x,y=x[cutaxis],y[cutaxis]
                    else:cutaxis=np.ones_like(y,dtype=bool)

                    if scen:pind=Scenario[scen]()[ind][nind][cutaxis]
                    else:pind=[True for _ in range(len(x))]
                    ax=caxes.top;d=x
                    bins=np.linspace(d.min(),d.max(),x_nbins);xb.append(bins)
                    # print(f'{m}: ');print(f'{bins.min()} - {bins.max()}')
                    ax.hist(d[pind],bins=bins,density=density,cumulative=cumulative,
                    stacked=stacked,histtype='step',alpha=1,linewidth=lw,color=color[si],rwidth=1,zorder=zorder[si])
                    ax=caxes.side;d=y
                    bins=np.linspace(d.min(),d.max(),y_nbins);yb.append(bins)
                    ax.hist(d[pind],bins=bins,density=density,cumulative=cumulative,orientation='horizontal',
                    stacked=stacked,histtype='step',alpha=1,linewidth=lw,color=color[si],rwidth=1,zorder=zorder[si])# edgecolor='k'
                    markerupscale=(1-zorder[si])*sz
                    # markerupscale=0
                    ax=caxes.center
                    ax.scatter(x,y, s=sz+markerupscale, c=color[si], cmap=cmap, lw=0.1,edgecolors='k',zorder=zorder[si],alpha=1)
                    ax.set_xlabel(yttl(f'{mnames[m].split('_')[0]}.Z'))
                    if scen:ax.scatter(x[pind],y[pind], s=sz+markerupscale, c=color[si], cmap=cmap, lw=0.1,edgecolors='r',zorder=1e10,alpha=1)

            for ax in [right.side,right.center,left.center]:ax.set_yticks([]);ax.set_yticklabels([]);ax.yaxis.set_minor_locator(plt.NullLocator())
            ax=right.side
            ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False, labelright=False)
            ax_right=ax.twinx()
            ax_right.tick_params(axis='y', which='both', left=False, right=True, labelleft=False, labelright=True)
            ysides.append(ax_right)
            for ax in [right.side,left.side]:ax.set_xlabel('counts')
            for ax in [left.top]:ax.set_ylabel('counts')         
            fig.suptitle(f'{p}, {b.replace('_',' to ')}s',y=.94)
            if scen:fig.suptitle(f'Scenario  {scen.replace('_','.').upper()} | {fig.get_suptitle()}',y=.94)
            cbar_ax=fig.add_axes([.125, -.1, 0.8-.025, 0.03])
            # cbar_ax=fig.add_axes([1.2, 0, 0.03, 0.8])

            if key in ['Instrument_Design', 'Network', 'Pressure_Gauge', 'Environment', 'Seismometer']:
                ncat=len(sets)
                boundaries=np.arange(ncat + 1);norm=mpl.colors.BoundaryNorm(boundaries, ncat)
                cbar_ticks=boundaries[:-1]+0.5 # centers
                cbar_ticklabels=[str(s) for s in sets]
                label=key.replace('_', ' ')
            else:
                # For continuous or binned variables (e.g., StaDepth, Magnitude)
                ncat=len(sets)
                boundaries=sorted(np.unique(sets))
                norm=mpl.colors.BoundaryNorm(boundaries, ncat)
                # boundaries=np.sort(np.unique(sets));norm=mpl.colors.Normalize(vmin=boundaries.min(), vmax=boundaries.max())
                cbar_ticks = boundaries
                cbar_ticklabels = [str(s) for s in boundaries]
                label = {'StaDepth':'Water depth, m','Magnitude':'Magnitude, Mw',
                'Distance':'Distance, °',
                'Sediment_Thickness_m':'Sediment thickness, m'}[key]
            sm=mpl.cm.ScalarMappable(cmap=cmap, norm=norm);sm.set_array([])
            cbar_ticks=cbar_ticks if not (key in ['StaDepth','Sediment_Thickness_m','Distance']) else cbar_ticks[::2]
            cbar_ticklabels=cbar_ticklabels if not (key in ['StaDepth','Sediment_Thickness_m','Distance']) else cbar_ticklabels[::2]
            cbar=fig.colorbar(sm, cax=cbar_ax, boundaries=boundaries, orientation='horizontal', label=label.lower().replace('mw','Mw'),shrink=0.7, aspect=30)
            cbar.set_ticks(cbar_ticks)
            cbar.set_ticklabels(cbar_ticklabels)
            if True:
                xscl=0.04;yscl=0.02
                xl,yl=gbounds([left_xbins,right_xbins],[left_ybins,right_ybins])
                if not ratio:
                    for ax in [left,right]:ax.side.set_ylim(yl);ax.center.set_ylim(yl)
                xl,yl=gbounds([left_xbins],[left_ybins]);left.center.set_xlim(xl);left.top.set_xlim(xl)
                xl,yl=gbounds([right_xbins],[right_ybins]);right.center.set_xlim(xl);right.top.set_xlim(xl)
                for ax in [left.center,left.top,right.center,right.top]:ax.set_xticks(np.linspace(xl.round(1).min(),xl.round(1).max(),3)[-2:])
                for ax in [left.side,right.side]:ax.set_xlim(left=-np.diff(ax.get_xlim())[0]*xscl)
                if not ratio:
                    for ax in [left.top,right.top]:ax.set_ylim(bottom=-np.diff(ax.get_ylim())[0]*xscl)
                    # for ax in [left.side,left.center,right.center,right.side,ax_right]:ax.grid(visible=True,which='major',alpha=1)
                if ratio:
                    for ax in [left.side,left.center,right.center,right.side,ax_right]:ax.set_ylim(yl)
                xl=np.array([ax.get_xlim() for ax in [left.center,right.center,left.top,right.top]]);xl=np.array([xl.min(),xl.max()])
                yl=np.array([ax.get_ylim() for ax in [left.center,right.center,left.side,right.side,ax_right]]);yl=np.array([yl.min(),yl.max()])
                for ax in [left.center,right.center,left.side,right.side,ax_right]:ax.set_ylim(yl)
                for ax in [left.center,right.center,left.top,right.top]:ax.set_xlim(xl)
                # ---------
                if symmetry:
                    for ax in [left.side,right.center]:ax.invert_xaxis() if (not ax.xaxis_inverted()) else None
                if ratio:
                    left.side.set_ylabel(f'{p} SNR RATIO using TF, log10')
                    ax_right.set_ylabel(f'{p} SNR RATIO using HPS, log10',rotation=-90,labelpad=15)
                else:
                    left.side.set_ylabel(f'{p} SNR using TF, log10')
                    ax_right.set_ylabel(f'{p} SNR using HPS, log10',rotation=-90,labelpad=15)
                assert np.all(np.array(ysides[0].get_ylim())==np.array([ax.get_ylim() for ax in ysides])), 'failed ylim update'
                assert np.all(np.sort(np.array(xsides[0].get_xlim()))==np.array([np.sort(ax.get_xlim()) for ax in xsides])), 'failed xlim update'
            bandname=f'{'_'.join([str(int(i)).zfill(2) for i in b.split('_')])}s'
            phasename=f'{p}.wave'
            keyname=key

            file=f'{'_'.join([phasename,keyname,bandname])}.{save_format}';file=file.replace('StaDepth','WaterDepth')
            file=f'{str(pi).zfill(2)}.{str(bi).zfill(2)}_{file}'
            if ratio:file=f'Ratio.w.Orig_{file}'
            if key=='StaDepth':file=f'_{file}'
            if scen:file=f'Scenario.{scen}_{file}' 
            save_tight(plotfolder/file,dpi=800)
