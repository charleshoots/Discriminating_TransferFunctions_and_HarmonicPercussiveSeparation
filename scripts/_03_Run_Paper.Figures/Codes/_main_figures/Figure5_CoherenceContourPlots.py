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

from scipy.stats import spearmanr, pearsonr, norm
# Coherence contour plot
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
# function dataset averaged coherence plot
def dataset_averaged_coherence_plot(f,z,coh,ms=35,zconnect=False,figsize=[6,6],cmap='viridis',checkerboard=(40,9),
    # title value
    title='Station Averaged Coherence',fontsize=4,levels=None,fig=None,ax=None,octav=True,fmin=1/100,fnlinewidth=0.4):
    # font value
    font = {'weight':'normal','size':fontsize};matplotlib.rc('font', **font)
    z = np.round(z)
    i = np.argsort(z)
    z,coh = z[i],coh[i,:]
    if levels is None:levels=np.linspace(np.min(coh),np.max(coh),20)
    if ax is None:fig,ax=plt.subplots(figsize=figsize)
    # cp value
    cp = coh
    # cp = np.array([smooth(d,k=3) for d in coh]);cp = gaussian_filter(coh,.5)

    checkerdensity,brightness=checkerboard
    # checkerboard value
    checkerboard = np.indices((checkerdensity, checkerdensity*2)).sum(axis=0) % 2  # 0/1 checker pattern
    from matplotlib.colors import LinearSegmentedColormap
    # dark greys value
    dark_greys = LinearSegmentedColormap.from_list('dark_greys', ['#777', f'#{str(brightness)*3}'])
    ax.imshow(checkerboard,cmap=dark_greys,interpolation='nearest',aspect='auto',extent=[0,1,0,1],origin='lower',
    # zorder value
    zorder=-1e3,transform=ax.transAxes,
    # alpha value
    alpha=1.0)  # Adjust transparency as needed)

    # s value
    s=ms
    # faxis value
    faxis=(f>=(fmin))&(f<=1)
    # F value
    F=f.copy()[faxis]
    # C value
    C=np.asarray(cp[:,faxis])            # same shape

    if octav:F,C = octavg(C,F)
    # C value
    C=np.array([C[z==zi,:].mean(axis=0) for zi in np.unique(z)])
    if zconnect:F, Z = np.meshgrid(F, np.arange(0,len(np.unique(z))))
    else:F, Z = np.meshgrid(F, np.unique(z)) 

    # norm value
    norm = mpl.colors.Normalize(vmin=0, vmax=1.0)
    # F,Z,C=np.flip(F,axis=0),np.flip(Z,axis=0),np.flip(C,axis=0) #Looks better with this off
    sc = ax.scatter(F.ravel(), Z.ravel(), c=C.ravel(), marker='s', s=s, cmap=cmap, norm=norm, edgecolors='none')
    ax.set_xscale('log')

    # fn z value
    fn_z = np.linspace(z.min(),z.max(),len(z)*10)
    # fn value
    fn = [fnotch(i) for i in fn_z]
    ax.plot(fn,fn_z,linestyle=':',color='k',linewidth=fnlinewidth)
    ax.set_xlim(fmin,1);ax.set_xscale('log')
    # n value
    n=1 if zconnect else 70;ax.set_ylim(np.min(Z)-n,np.max(Z)+n)
    # fticks value
    fticks = np.array([1/100,1/50,1/30,1/10,1])
    ax.set_xticks(fticks)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticklabels(np.array(1/fticks,dtype=int))
    if zconnect:
        # zy value
        zy=np.arange(0,6500,1500)
        # yt value
        yt=np.round(np.interp(zy,np.unique(z),np.arange(0,len(np.unique(z)),1)))
        ax.set_yticks(yt)
        ax.set_yticklabels(zy)
    ax.set_ylabel('Water depth, m',fontweight='bold')
    ax.set_xlabel('Period, s',fontweight='bold')
    # ax.set_title(f'{coh.shape[0]} {title}')
    ax.set_facecolor('k')
    if fig is not None:fig.suptitle(title)
    # plt.colorbar(cnt)
    plt.tight_layout()
    # norm value
    norm = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)
    # cbar value
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax,extend=None)
    cbar.set_ticks(np.arange(0, 1.01, 0.25))
    if fig is not None:return fig
    if ax is not None:return ax

# make bands
def make_bands(N, width=None,line=np.linspace, lo=1.0, hi=100.0):
    if width is None:width=(hi-lo)/N
    # if width<=0 or width>(hi-lo): raise ValueError("width must be in (0, hi-lo]")
    s=line(lo, hi-width, N)
    return np.c_[s, s+width]
# function local corr
def local_corr(x,y,N=200,bins=500,quantile=True,kind='pearson',nmin=1):
    x,y=np.asarray(x),np.asarray(y); m=np.isfinite(x)&np.isfinite(y); x,y=x[m],y[m]
    # edges=np.quantile(x,np.linspace(0,1,bins+1)) if quantile else np.linspace(np.min(x),np.max(x),bins+1)
    # dbin(np.unique(x))
    # edges=np.sort(np.unique(x))[::3]
    c,lo,hi,xc=[],[],[],[]
    # aa,bb=edges[:-1],edges[1:]
    edges=make_bands(N,bins,lo=min(x)-bins,hi=max(x)+bins).T
    aa,bb=edges
    for _ in range(50):
        # idx value
        idx=np.array([np.sum((x>=a)&(x<=b))>0 for a,b in zip(aa,bb)])
        # edges value
        edges=edges[:,idx]
        aa,bb=edges

    for a,b in zip(aa,bb):
        i=(x>=a)&(x<=b); n=i.sum()
        if n<nmin:
            c.append(np.nan); lo.append(np.nan); hi.append(np.nan); xc.append(.5*(a+b)); continue
        # r value
        r = spearmanr(x[i],y[i]).correlation if kind=='spearman' else pearsonr(x[i],y[i])[0]
        if np.isnan(r):
            j=0
        z=.5*np.log((1+r)/(1-r)) if abs(r)<1 else np.sign(r)*np.inf
        # se value
        se=1/np.sqrt(max(n-3,1)); zlo,zhi=z-1.96*se,z+1.96*se
        # rlo value
        rlo=(np.exp(2*zlo)-1)/(np.exp(2*zlo)+1); rhi=(np.exp(2*zhi)-1)/(np.exp(2*zhi)+1)
        c.append(r); lo.append(rlo); hi.append(rhi); xc.append(.5*(a+b))
    return np.array(xc),np.array(c),np.array(lo),np.array(hi)
# plot local corr
def plot_local_corr(ax,x,y,clr='blue',ms=3,alpha=.3,**kw):
    xc,c,lo,hi=local_corr(x,y,**kw)
    ax.fill_betweenx(xc,lo,hi,alpha=alpha,linewidth=0,color=clr)
    ax.plot(c,xc,'-o',ms=ms, lw=1,color=clr,alpha=alpha)
    ax.axvline(0,color='k',ls=':',lw=1)
    return c




# plotfolder value
plotfolder=dirs.Ch1/'_main_figures'/'Figure5_CoherenceContour';plotfolder.mkdir(parents=True,exist_ok=True)

mpl.rcParams.update({
"font.size": 4,              # base text size (fallback for everything)
"axes.titlesize": 4,         # axes titles
"axes.labelsize": 4,         # x/y labels (also used by colorbar label)
"xtick.labelsize": 4,        # x tick labels (affects horizontal colorbar ticks)
"ytick.labelsize": 4,        # y tick labels (affects vertical colorbar ticks)
"legend.fontsize": 4,        # legend text
"legend.title_fontsize": 4,  # legend title
"figure.titlesize":4,    # suptitle
'ytick.major.width':0.5,
'xtick.major.width':0.5,
'axes.edgecolor':'k',
'axes.linewidth':0.5})

# COHERENCE CONTOURS

icat = catalog.sr.copy()
usnr = unpack_metrics(icat)

# Coherence contours
z=np.array(list(icat.StaDepth))
zi=np.unique(z)[0];m='TF_Z';z==zi
stat = np.median
coh=usnr.coh.__dict__[m].D.copy()
f=1/usnr.coh.bands
foct,coh_oct=octave_average(coh,f,domain='geo')
coh_oct=np.array([stat(coh_oct[z==zi,:],axis=0) for zi in np.unique(z)])

# idx = (icat.Magnitude<=(6.5))
# idx = (icat.Magnitude>=(7.0))
# idx = (icat.Magnitude>=(7.0))&(icat.Magnitude<=(7.5))
# icat=icat[idx]
# ms=33 #markersize
ms=25 #markersize, default is 25
# fig,axes=figs(2,2,f=[4,6])


figsize=(6,4)

fig=plt.figure(figsize=figsize, constrained_layout=False)
# gs = gridspec.GridSpec(2, 4, width_ratios=[3,.9,3,.9], height_ratios=[1,1], wspace=.1, hspace=.25)
# contour1=fig.add_subplot(gs[0, 0])
# sc1=fig.add_subplot(gs[0, 1])
# contour2=fig.add_subplot(gs[0, 2],sharex=contour1)
# sc2=fig.add_subplot(gs[0, 3])
# contour3=fig.add_subplot(gs[1, 0])
# sc3=fig.add_subplot(gs[1, 1])
# contour4=fig.add_subplot(gs[1, 2],sharex=contour1)
# sc4=fig.add_subplot(gs[1, 3])

# Inner wspace (0.06) controls the tiny gap between a wide contour and its thin scatter panel.
# Outer wspace (0.28) keeps the two pairs visually separate.
# Outer 2x2 of "pairs" (each pair = [wide, thin])
gapbetweenpairs=-.2
gapbetweenquads=0.15
figuregap=0.01

outer = gridspec.GridSpec(2, 2, figure=fig, wspace=gapbetweenquads, hspace=0.25)  # spacing between pairs

# Top-left pair
TL=gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0,0],width_ratios=[1.0, 0.28], wspace=gapbetweenpairs)
contour1=fig.add_subplot(TL[0,0]);sc1=fig.add_subplot(TL[0,1])

# Top-right pair
TR=gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0,1],width_ratios=[1.0, 0.28], wspace=gapbetweenpairs)
contour2=fig.add_subplot(TR[0,0], sharex=contour1);sc2=fig.add_subplot(TR[0,1])

# Bottom-left pair
BL=gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1,0],width_ratios=[1.0, 0.28], wspace=gapbetweenpairs)
contour3=fig.add_subplot(BL[0,0]);sc3=fig.add_subplot(BL[0,1])

# Bottom-right pair
BR=gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1,1],width_ratios=[1.0, 0.28], wspace=gapbetweenpairs)
contour4=fig.add_subplot(BR[0,0], sharex=contour1);sc4=fig.add_subplot(BR[0,1])
# fig.subplots_adjust(left=0.07, right=0.98, top=0.96, bottom=0.08,wspace=figuregap)



axes=np.array([[contour1,sc1,contour2,sc2],[contour3,sc3,contour4,sc4]])
cohaxes=np.array([axes[0,0],axes[0,2],axes[1,0],axes[1,2]])
scatteraxes=np.array([axes[0,1],axes[0,3],axes[1,1],axes[1,3]])

yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
methods=['TF_Z','HPS_Z','HPS_1','HPS_2']
for m,ax in zip(methods,cohaxes.reshape(-1)):
    z=np.array(list(icat.StaDepth))
    stat = np.median
    coh=usnr.coh.__dict__[m].D.copy()
    f=1/usnr.coh.bands
    pz,pf,pcoh=z,f,coh
    foct,coh_oct=octave_average(coh,f,domain='geo')
    pz,pf,pcoh=z,foct,coh_oct
    cmap=darken(mpl.colormaps.get_cmap('viridis'),.95)
    cmap=cm.lipari
    pz,pcoh=np.unique(z),np.array([stat(pcoh[pz==zi,:],axis=0) for zi in np.unique(pz)])
    dataset_averaged_coherence_plot(pf,pz,pcoh,ax=ax,ms=ms,figsize=[6,6],octav=False,
    cmap=cmap)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
[ax.set_xlabel(None) for ax in cohaxes]
[ax.set_ylabel(None) for ax in [cohaxes[1],cohaxes[3]]]
[ax.set_ylim(min(z),max(z)) for ax in axes.reshape(-1)]
[ax.set_yticks([]) for ax in [cohaxes[1],cohaxes[3]]]
cbar_axes=[a for a in fig.axes if a.get_label() == '<colorbar>']
[a.remove() for a in cbar_axes]
ax=cohaxes[0]
norm = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)
cax=inset_axes(ax, width='260%', height='3%',
loc="lower left", bbox_to_anchor=(.07,1.06, .97, 1.0), #[x0, y0, width, height] in axes coordinated
bbox_transform=ax.transAxes, borderpad=0)
cbar=plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),cax=cax,extend=None,orientation='horizontal')
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.xaxis.set_label_position('top')
cbar.ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
cbar.set_label('Coherence')
# _=[ax.yaxis.set_label_position('right') for ax in scatteraxes]
_=[ax.yaxis.set_ticks_position('none') for ax in scatteraxes]
[ax.tick_params(left=False, right=False, labelright=False, labelleft=False) for ax in scatteraxes]
# [ax.set_yticks([]) for ax in scatteraxes]
cbar.set_ticks(np.arange(0, 1.01, 0.25))

alpha=0.5
fns=['IG','MS',None]
# fns=['IG']
for fn in fns:
    ms=1
    print(fn if (fn is not None) else 'Full')
    x=np.array(list(icat.StaDepth))
    band=(1,100)
    N=len(np.unique(x)) #98 or higher is good
    bins=np.diff(np.sort(np.unique(x))).max()*1.3
    for axi,(mthd,ax)in enumerate(zip(methods,scatteraxes)):
        ax.cla()
        # y1=usnr.coh.TF_Z.Average(band,fn=fn,agg='median')
        y1=usnr.snr.__dict__[mthd].Average(band,fn=fn,agg='median').mean(axis=1)
        y2=usnr.coh.__dict__[mthd].Average(band,fn=fn,agg='median')
        c1=plot_local_corr(ax,x,y1,clr='orangered',bins=bins,quantile=True,kind='pearson',N=N,ms=ms,alpha=alpha)
        c2=plot_local_corr(ax,x,y2,clr='skyblue',bins=bins,quantile=True,kind='pearson',N=N,ms=ms,alpha=alpha)
        ax.set_ylim(x.min(),x.max())

    file=f'Contours.{fn if (fn is not None) else 'Full'}.png'
    save_tight(plotfolder/file,fig,dpi=800)