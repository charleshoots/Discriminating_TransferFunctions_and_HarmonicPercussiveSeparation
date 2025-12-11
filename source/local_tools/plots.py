# from imports import *
from modules import *
from local_tools.math import avg_meter
import local_tools as lt
from local_tools import ObsQA
import scipy.stats as stats
from matplotlib.gridspec import GridSpec

octavg=lt.math.octave_average
from scipy.stats import iqr
from local_tools.math import cohstats
from scipy.ndimage import gaussian_filter
from local_tools.math import fnotch
from local_tools.LabeledMatrix import AggregateMeasurements as LM

def argsort_luminence(cc,cmap):
    cc=np.array(cc).ravel()
    luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098
    color=[cmap(c) for c in enumerate(cc)]
    zorder=np.argsort(np.array([luminance(c) for c in color]))
    return zorder
def color_master(key,sets=np.nan,N=100):
    if isinstance(sets,int):sets=np.arange(0,sets,1)
    elif (np.shape(np.atleast_1d(np.array(sets)))[0]-1)==0:sets=np.arange(0,N,1)
    #reverses the colormap of any key in this variable, regardless of other conditions.
    color_overrides=['Pressure_Gauge','Magnitude','Sediment_Thickness_m']
    #Luminnance model, on a scale from 0 (black) to 1 (white)
    darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)],name=cmap.name).resampled(N)
    luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
    cmap = darken(cm.cmaps['bamako'].copy(),frac=.8)
    if key=='Instrument_Design':cmap=ListedColormap([ColorStandard.instrument[s] for s in sets], name='custom_cmap')
    elif key=='Network':cmap=ListedColormap([ColorStandard.network[s] for s in sets], name='custom_cmap')
    elif key=='Magnitude':cmap=darken(cm.cmaps['lipari'],.8)
    elif key=='StaDepth':cmap=cm.cmaps['glasgow'].reversed()
    elif key=='Distance':cmap=darken(cm.cmaps['devon'])
    elif key=='Sediment_Thickness_m':cmap=darken(cm.cmaps['lapaz'].copy(),.9).reversed()
    elif key=='Coherence':cmap=darken(cm.cmaps['lajolla'].copy(),.9)
    else:cmap=darken(darken(cm.cmaps['glasgow'].resampled(N),0.7).reversed(),0.9)
    cmap=cmap.resampled(len(sets))
    if (luminance(cmap(0))<luminance(cmap(1e3)))&(not isinstance(sets[0],str)):cmap=cmap.reversed()
    if key in color_overrides:cmap=cmap.reversed()
    color=[cmap(si/len(sets)) for si,_ in enumerate(sets)]
    #suggested zorder based on luminance, will put the darker colors in the back
    zorder = np.array([luminance(c) for c in color])
    return color,zorder,cmap
# ########################################################################################
# ########################################################################################

def meancoh(icat,methods=['TF','HPS_Z','HPS_1','HPS_2'],bands=['1_10','10_30','30_100'],notched=False,f=[None]):
    if f[0]==None:f=icat.iloc[0].Data.Coherence().f
    #Calculates average coherence for every source-receiver pair in every method and every band defined in bands.
    aggregate=lambda b,method,s,agg=np.mean,fn=True:agg(s.Coherence[method].reshape(-1)[((1/f)>=min(b))&((1/f)<=max(b)) & fn])
    coh=AttribDict({b:AttribDict({m:pd.Series([aggregate([int(i) for i in b.split('_')],m,s,fn=f<=fnotch(s.StaDepth) if notched else True)
    for s in icat.iloc]) for m in methods}) for b in bands})
    return coh
def reduce_snr(icat,ratio=False,mergephases=True,log=False,win=None):
    sr=icat.iloc[0]
    bands=list(sr.SNR.keys());b=bands[0]
    phases=np.array(list(sr.SNR[bands[0]].keys()));phases[(~(phases=='LT'))&(~(phases=='stnm'))];p=phases[0]
    phases = ['P','Pdiff','S','Sdiff','Rg']
    if win=='LT':mergephases=False
    methods=list(sr.SNR[b][p].keys());m=methods[0]
    if ratio:
        # methods=['TF.HZ','HPS.HZ',['HPS.H1','HPS.H2',],'Original.HZ','Original.HZ',['Original.H1','Original.H2',]]
        methods=['TF.HZ','HPS.HZ','HPS.H1','HPS.H2']
    else:
        # methods=['TF.HZ','HPS.HZ',['HPS.H1','HPS.H2',],'Original.HZ',['Original.H1','Original.H2',]]
        methods=['TF.HZ','HPS.HZ','HPS.H1','HPS.H2','Original.HZ','Original.H1','Original.H2']
        # methods=['TF.HZ','HPS.HZ',['HPS.H1','HPS.H2']]
    # if win not in ['ST','LT']:
    #     D=[[[[(sr.SNR[b][p][m] if p in sr.SNR[b].keys() else None) if not isinstance(m,list) else (np.mean([sr.SNR[b][p][m[0]],sr.SNR[b][p][m[1]]]) if p in sr.SNR[b].keys() else None) for b in bands] for p in phases] for m in methods] for sr in icat.iloc]
    # elif win=='ST':
    #     D=[[[[(sr.SNR[b][win][p][m] if p in sr.SNR[b][win].keys() else None) if not isinstance(m,list) else (np.mean([sr.SNR[b][win][p][m[0]],sr.SNR[b][win][p][m[1]]]) if p in sr.SNR[b][win].keys() else None) for b in bands] for p in phases] for m in methods] for sr in icat.iloc]
    # elif win=='LT':
    #     D=[[[[(sr.SNR[b][win][m] ) if not isinstance(m,list) else (np.mean([sr.SNR[b][win][m[0]],sr.SNR[b][win][m[1]]]) ) for b in bands] for p in phases] for m in methods] for sr in icat.iloc]
    if win not in ['ST','LT']:
        D=[[[[sr.SNR[b][p][m] if p in sr.SNR[b].keys() else None for b in bands] for p in phases] for m in methods] for sr in icat.iloc]
    elif win=='ST':
        D=[[[[sr.SNR[b][win][p][m] if p in sr.SNR[b][win].keys() else None for b in bands] for p in phases] for m in methods] for sr in icat.iloc]
    elif win=='LT':
        D=[[[[sr.SNR[b][win][m] for b in bands] for p in phases] for m in methods] for sr in icat.iloc]
    D=np.array(D) # sr X method X phase X band # (4307, 3, 5, 10)
    D[D==None]=np.nan
    if ratio:D=D[:,:3,:,:]/D[:,3:,:,:]
    if mergephases:
        #Moves axis to the shape: method, sr, bands, ph = 100, 50, 3
        # method X sr X band X phase
        Merged=D.copy()
        Merged = np.array([np.nanmean(Merged[:,:,0:2,:],axis=2),np.nanmean(Merged[:,:,2:4,:],axis=2),Merged[:,:,-1,:]])
        phases=['P','S','Rg'];Shaped=np.moveaxis((np.moveaxis(Merged,0,-1)),1,0)
        D=Shaped
    else:D=np.moveaxis(np.moveaxis(D,1,0),2,3)
    Nmethods,Nsr,Nbands,Nphases=len(methods),len(icat),len(bands),len(phases)
    assert (D.shape[0]==Nmethods)
    assert (D.shape[1]==Nsr)
    assert (D.shape[2]==Nbands)
    assert (D.shape[3]==Nphases)
    if log:D=np.log10(D.tolist())
    # Outputs should always be this shape:# method X sr X band X phase
    return D
def unpack_coh(icat,flim=0.01):
    f=icat.iloc[0].Data.Coherence().f
    faxis=(f>=flim)&(f<=1)#;f=f[faxis]
    mthds=['TF_Z','HPS_Z','HPS_1','HPS_2']
    bands=1/f[faxis]
    D=np.array([np.array([sr.Coherence[m].reshape(-1)[faxis] for sr in icat.iloc]) for m in mthds])
    return bands,D
def unpack_metrics(icat=None,new=False,ratio=False,mergephases=True,log=False,flim=1/200):
    if not new:
        icatpairs=np.array(list(icat.StaName + '-' + icat.Name))
        depth = np.array(list(icat.StaDepth))
        dirs = lt.io.dir_libraries()
        D=pd.read_pickle(dirs.Data/'SNR_Models'/'BulkHold.pkl')
        dpairs=np.array([f'{a}-{b}' for a,b in zip(D['stnm'],D['evname'])])
        idx=np.array([np.where(dpairs==i)[0][0] for i in icatpairs])
        keys=[['coh','TF_Z', 'HPS_Z','HPS_1','HPS_2'], ['snr','TF_Z','HPS_Z','HPS_1','HPS_2', 'Original_Z', 'Original_1','Original_2'], ['ST','TF_Z', 'HPS_Z','HPS_1','HPS_2', 'Original_Z', 'Original_1','Original_2'], ['LT','TF_Z','HPS_Z','HPS_1','HPS_2','Original_Z','Original_1','Original_2']]
        Q={}
        for k in keys:
            M={}
            phases=None
            bands=None
            k0=k[0]
            D0=D[k0]
            bands=D0['bands'].copy()
            phaseinvariant=True if k0 in ['coh','LT'] else False
            if k0 in ['snr','ST']:phases=D0['phases'].copy()
            for j in k[1:]:
                data=D0[j].copy()
                inputs={'D':data[idx,:] if phaseinvariant else data[idx,:,:],'phases':phases,'bands':bands,'depth':depth,'comp':f'Original_{j.split('_')[-1]}'}
                A=LM(**inputs)
                M.update({j:A})
                M.update({'phases':phases,'depth':depth,'bands':bands})
            M=LM(**M)
            Q.update({k0:M})
        Q.update({'phases':None,'bands':None})
        output=LM(**Q)
        return output
    reray=lambda w:np.array(w.tolist())
    D=reduce_snr(icat,ratio=ratio,mergephases=mergephases,log=log) #['TF.HZ','HPS.HZ','HPS.H1','HPS.H2','Original.HZ','Original.H1','Original.H2']
    sr=icat.iloc[0]
    bands=list(sr.SNR.keys());b=bands[0]
    bands_centers=np.array([[float(j) for j in i.split('_')] for i in bands])[:,1]
    phases=["P", "S", "Rg"]
    TF=LM(D=reray(D[0]), bands=bands_centers, phases=phases,comp='Original_Z')
    HPS_Z=LM(D=reray(D[1]),  bands=bands_centers, phases=phases,comp='Original_Z')
    HPS_1=LM(D=reray(D[2]),  bands=bands_centers, phases=phases,comp='Original_1')
    HPS_2=LM(D=reray(D[3]),  bands=bands_centers, phases=phases,comp='Original_2')
    Original_Z=LM(D=reray(D[4]),  bands=bands_centers, phases=phases,comp='Original_Z')
    Original_1=LM(D=reray(D[5]),  bands=bands_centers, phases=phases,comp='Original_1')
    Original_2=LM(D=reray(D[6]),  bands=bands_centers, phases=phases,comp='Original_2')
    SNR=LM(TF_Z=TF,HPS_Z=HPS_Z,HPS_1=HPS_1,HPS_2=HPS_2,Original_Z=Original_Z,Original_1=Original_1,Original_2=Original_2,bands=bands_centers,comp=None)  # top-level container
    # All phases reduced over bands -> (sr, ph)
    all_phases = SNR.TF_Z.Average(band=(1, 10))
    # Single-phase as a child LM -> (sr,)
    p_line = SNR.TF_Z.P.Average(band=(1, 10))
    s_line = SNR.TF_Z.S.Average(band=(1, 10))
    # Also works via dict style if you prefer:
    rg_line = SNR.TF_Z["Rg"].Average(band=(1, 10))
    bands,D=unpack_coh(icat,flim)
    TF=LM(D=D[0],bands=bands)
    HPS_Z=LM(D=reray(D[1]),bands=bands,comp='Original_Z')
    HPS_1=LM(D=reray(D[2]),bands=bands,comp='Original_1')
    HPS_2=LM(D=reray(D[3]),bands=bands,comp='Original_2')
    COH=LM(TF_Z=TF,HPS_Z=HPS_Z,HPS_1=HPS_1,HPS_2=HPS_2, bands=bands,comp=None)

    QLT=reduce_snr(icat,win='LT')#['TF.HZ','HPS.HZ','HPS.H1','HPS.H2','Original.HZ','Original.H1','Original.H2']
    TF_Z=LM(D=reray(QLT[0,:,:,0]),bands=bands_centers,comp='Original_Z')
    HPS_Z=LM(D=reray(QLT[1,:,:,0]),bands=bands_centers,comp='Original_Z')
    HPS_1=LM(D=reray(QLT[2,:,:,0]),bands=bands_centers,comp='Original_1')
    HPS_2=LM(D=reray(QLT[3,:,:,0]),bands=bands_centers,comp='Original_2')
    Original_Z=LM(D=reray(QLT[4,:,:,0]),bands=bands_centers,comp='Original_Z')
    Original_1=LM(D=reray(QLT[5,:,:,0]),bands=bands_centers,comp='Original_1')
    Original_2=LM(D=reray(QLT[6,:,:,0]),bands=bands_centers,comp='Original_2')
    LT=LM(TF_Z=TF_Z,HPS_Z=HPS_Z,HPS_1=HPS_1,HPS_2=HPS_2,Original_Z=Original_Z,Original_1=Original_1,Original_2=Original_2,bands=bands_centers,comp=None)
    QST=reduce_snr(icat,win='ST')#['TF.HZ','HPS.HZ','HPS.H1','HPS.H2','Original.HZ','Original.H1','Original.H2']
    TF_Z=LM(D=reray(QST[0]),bands=bands_centers,phases=phases,comp='Original_Z')
    HPS_Z=LM(D=reray(QST[1]),bands=bands_centers,phases=phases,comp='Original_Z')
    HPS_1=LM(D=reray(QST[2]),bands=bands_centers,phases=phases,comp='Original_1')
    HPS_2=LM(D=reray(QST[3]),bands=bands_centers,phases=phases,comp='Original_2')
    Original_Z=LM(D=reray(QST[4]),bands=bands_centers,phases=phases,comp='Original_Z')
    Original_1=LM(D=reray(QST[5]),bands=bands_centers,phases=phases,comp='Original_1')
    Original_2=LM(D=reray(QST[6]),bands=bands_centers,phases=phases,comp='Original_2')
    ST=LM(TF_Z=TF_Z,HPS_Z=HPS_Z,HPS_1=HPS_1,HPS_2=HPS_2,Original_Z=Original_Z,Original_1=Original_1,Original_2=Original_2,bands=bands_centers,comp=None)

    output=LM(coh=COH,snr=SNR,ST=ST,LT=LT)
    return output


# ########################################################################################
# ########################################################################################


def fig_setup(figsize=(6, 4),width_ratios=[.4, .5,.5, .4],height_ratios=[.2, .5],debug=False):
    fig=plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(nrows=2, ncols=4, width_ratios=width_ratios, height_ratios=height_ratios)
    ax_top_center_left = fig.add_subplot(gs[0, 1])
    ax_top_center_right = fig.add_subplot(gs[0, 2],sharey=ax_top_center_left)
    ax_bottom_center_left = fig.add_subplot(gs[1, 1],sharex=ax_top_center_left)
    ax_bottom_center_right = fig.add_subplot(gs[1, 2],sharex=ax_top_center_right)#,sharey=ax_bottom_center_left)
    ax_bottom_left = fig.add_subplot(gs[1, 0])#,sharey=ax_bottom_center_left)
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
    left=AttribDict({});left.top=ax_top_center_left;left.side=ax_bottom_left;left.center=ax_bottom_center_left
    right=AttribDict({});right.top=ax_top_center_right;right.side=ax_bottom_right;right.center=ax_bottom_center_right
    left.name='left';right.name='right'
    fig.subplots_adjust(wspace=.35,hspace=0.1)
    return fig,left,right

def simple_hist(isnr,usnr,key='StaDepth',sets=np.arange(0,6500,500),meter='snr',p='Rg',y_nbins=20,ratio=True,filledstack=False,figsize=(4,4),
    log=False,cumulative=False,density=False,stacked=True,methods=['TF.Z', 'HPS.Z', 'HPS.H'],xl=[None],orientation='vertical'):
    # snrmin,snrmax=[-1.0,2.0]
    snrmin,snrmax=[-1.0,3.0]
    # snrmin,snrmax=[-1.0,6.0]
    cohmin,cohmax=[0,1.0]
    if meter=='coh':xmin,xmax=cohmin,cohmax
    else:xmin,xmax=snrmin,snrmax
    bins_collect=[]
    figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
    yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
    metername={'snr':'SNR','coh':'Coherence','ST':'Signal','LT':'Noise',}[meter]
    if meter=='snr':title=lambda:f'{'Relative' if ratio else ''} SNR for {p.replace('Rg','Rayleigh')} waves using {mtitle}'
    # if meter=='coh':title=lambda:f'{metername} {m}'
    if meter=='coh':title=lambda:f'{yttl(f'{mtitle}')}'
    if meter=='ST':title=lambda:f'{'Relative' if ratio else ''} Signal for {p.replace('Rg','Rayleigh')} waves using {mtitle}'
    if meter=='LT':title=lambda:f'{'Relative' if ratio else ''} Pre-arrival noise using {mtitle}'
    if meter=='coh':p='Rg';ratio=False

    lw=0.4 if filledstack else 1.8
    # lw=.7
    color,zorder,cmap=color_master(key,sets)
    fig,axes=figs(3,f=figsize,y='all')
    if (not ratio)&(not filledstack)&(not meter=='coh'):
        ax=axes[0];m='Original'
        d=np.array([usnr[meter][m][i] for i in ['P','S','Rg']]) if p=='all' else np.atleast_2d(np.array(usnr[meter][m][p]))
        d=[d[:,(isnr[key]>=sets[si])&(isnr[key]<=sets[si+1])].reshape(-1) for si in range(len(sets)-1)]
        bins=[np.linspace(di.min(),di.max(),y_nbins) for di in d]
        _=[ax.hist(d[si],bins=bins[si],orientation=orientation,density=density,cumulative=cumulative,stacked=stacked,histtype='step',
        alpha=.7,linewidth=lw,linestyle='--',color=color[si],rwidth=1,zorder=zorder[si],log=log) for si in range(len(sets)-1)]
        bins_collect.append(bins)
    ax=axes[0];m=methods[0];mtitle=m.replace('HZ','Z').replace('_','.').replace('.',' ')
    ax.set_title(title())
    d=np.array([usnr[meter][m][i] for i in ['P','S','Rg']]) if p=='all' else np.atleast_2d(np.array(usnr[meter][m][p]))
    d=[d[:,(isnr[key]>=sets[si])&(isnr[key]<=sets[si+1])].reshape(-1) for si in range(len(sets)-1)]
    
    # if meter=='coh':bins=np.linspace(0,1,y_nbins)
    if meter=='coh':bins=[np.linspace(di.min(),di.max(),y_nbins) for di in d]
    else:
        bins=[np.linspace(max([snrmin,di.min()]),min([snrmax,di.max()]),y_nbins) for di in d]
        bins=[np.array(np.sort(np.hstack([np.linspace(min(b),0,y_nbins),np.linspace(0,max(b),y_nbins)]))) for b in bins]
    # --------------

    if filledstack:
        bins=np.linspace(np.min(bins),np.max(bins),y_nbins)
        ax.hist(d,bins=bins,orientation=orientation,density=density,cumulative=cumulative,stacked=stacked,histtype='barstacked',
        alpha=1,linewidth=lw,color=color[:len(d)],rwidth=1,edgecolor='w',log=log)
    else:
        _=[ax.hist(d[si],bins=bins[si],orientation=orientation,density=density,cumulative=cumulative,
        stacked=stacked,histtype='step',alpha=1.0,linewidth=lw,color=color[si],rwidth=1,zorder=zorder[si],log=log) for si in range(len(sets)-1)]
    bins_collect.append(bins)
    # --------------

    ax=axes[1];m=methods[1];mtitle=m.replace('HZ','Z').replace('_','.').replace('.',' ')
    ax.set_title(title())
    d=np.array([usnr[meter][m][i] for i in ['P','S','Rg']]) if p=='all' else np.atleast_2d(np.array(usnr[meter][m][p]))
    d=[d[:,(isnr[key]>=sets[si])&(isnr[key]<=sets[si+1])].reshape(-1) for si in range(len(sets)-1)]
    if filledstack:
        bins=np.linspace(np.min(bins),np.max(bins),y_nbins)
        ax.hist(d,bins=bins,orientation=orientation,density=density,cumulative=cumulative,stacked=stacked,histtype='barstacked',
        alpha=1,linewidth=lw,color=color[:len(d)],rwidth=1,edgecolor='w',log=log)
    else:
        _=[ax.hist(d[si],bins=bins[si],orientation=orientation,density=density,cumulative=cumulative,
        stacked=stacked,histtype='step',alpha=1.0,linewidth=lw,color=color[si],rwidth=1,zorder=zorder[si],log=log) for si in range(len(sets)-1)]
    bins_collect.append(bins)
    # --------------
    ax=axes[2];m=methods[2];mtitle=m.replace('HZ','Z').replace('_','.').replace('.',' ') if meter=='coh' else 'HPS(H)'
    ax.set_title(title())
    d=np.array([usnr[meter][m][i] for i in ['P','S','Rg']]) if p=='all' else np.atleast_2d(np.array(usnr[meter][m][p]))
    d=[d[:,(isnr[key]>=sets[si])&(isnr[key]<=sets[si+1])].reshape(-1) for si in range(len(sets)-1)]
    if filledstack:
        bins=np.linspace(np.min(bins),np.max(bins),y_nbins)
        ax.hist(d,bins=bins,orientation=orientation,density=density,cumulative=cumulative,stacked=stacked,histtype='barstacked',
        alpha=1,linewidth=lw,color=color[:len(d)],rwidth=1,edgecolor='w',log=log)
    else:
        _=[ax.hist(d[si],bins=bins[si],orientation=orientation,density=density,cumulative=cumulative,
        stacked=stacked,histtype='step',alpha=1.0,linewidth=lw,color=color[si],rwidth=1,zorder=zorder[si],log=log) for si in range(len(sets)-1)]
    bins_collect.append(bins)
    cbar_ax=fig.add_axes([1.06, 0.1, 0.02,0.8])
    if key in ['Instrument_Design', 'Network', 'Pressure_Gauge', 'Environment', 'Seismometer']:
        ncat=len(sets)
        boundaries=np.arange(ncat + 1);norm=mpl.colors.BoundaryNorm(boundaries, ncat)
        cbar_ticks=boundaries[:-1]+0.5 # centers
        cbar_ticklabels=[str(s) for s in sets]
        label=key.replace('_', ' ')
    else:
        # For continuous or binned variables (e.g., StaDepth, Magnitude)
        boundaries=np.sort(np.unique(sets));norm=mpl.colors.Normalize(vmin=boundaries.min(), vmax=boundaries.max())
        cbar_ticks = boundaries
        cbar_ticklabels = [str(s) for s in boundaries]
        label = {'StaDepth':'Water depth, m','Magnitude':'Magnitude, Mw',
        'Distance':'Distance, °',
        'Sediment_Thickness_m':'Sediment thickness, m'}[key]
    if ratio & (not meter=='coh'):_=[ax.axvline(0.0,ls='--',lw=0.8) for ax in axes]
    if (not cumulative==False):_=[ax.axhline(len(isnr)/2,ls='--',lw=0.8) for ax in axes]
    if xl[0]==None:
        dx=.05
        xl=[xmin-dx,xmax+dx] if meter=='snr' else [-.05,1.05]
    else:
        dx = min(np.diff(np.sort(np.unique(bins_collect))))/2
        xmin,xmax=np.min(bins_collect),np.max(bins_collect)
        xl=[xmin-dx,xmax+dx]
    _=[ax.set_xlim(xl) for ax in axes]
    if meter=='snr':
        if log:
            _=[ax.set_ylim([1,3000]) for ax in axes];_=[ax.set_yticks([10,100,1000])]
    sm=mpl.cm.ScalarMappable(cmap=cmap, norm=norm);sm.set_array([])
    cbar_ticks=cbar_ticks if not (key in ['StaDepth','Sediment_Thickness_m']) else cbar_ticks[::2]
    cbar_ticklabels=cbar_ticklabels if not (key in ['StaDepth','Sediment_Thickness_m']) else cbar_ticklabels[::2]
    cbar=fig.colorbar(sm, cax=cbar_ax, boundaries=boundaries, orientation='vertical', label=label,shrink=0.7, aspect=30)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(cbar_ticklabels)


def alpha_empties(y,bps,alpha=0.0,minvals=0,sp=False):
    for i,bp in zip(y,bps):
        if len(i)<=(minvals):[[j.set_alpha(alpha) for j in bp[k]] for k in bp.keys()]
    for i,bp in zip(y,bps):
        for mean in bp['means']:
            if len(i)<=(minvals):mean.set_alpha(alpha)
    for i,bp in zip(y,bps):
        for median in bp['medians']:
            if len(i)<=(minvals):median.set_alpha(0.0)
    for i,bp in zip(y,bps):
        for whiskers in bp['caps']:
            if len(i)<=(minvals):whiskers.set_alpha(0.0)
def dfilt(xx,y,bins,igp_parse_for_tf=True):
    y=np.array(y)
    gdaxis=lambda x,y,bins:[np.array(x>=bins[i])&np.array(x<=bins[i+1])&np.array((~np.isnan(y))) for i in range(len(bins)-1)]
    gdaxis_str=lambda x,y,bins:[np.array(x==bins[i])&np.array((~np.isnan(y))) for i in range(len(bins))]
    if isinstance(bins[0],str):
        daxis=gdaxis_str(xx,y,bins)
    else:
        y[xx<np.min(bins)]=np.nan;y[xx>np.max(bins)]=np.nan
        xx[xx<np.min(bins)]=np.min(bins);xx[xx>np.max(bins)]=np.max(bins)
        daxis=gdaxis(xx,y,bins)
    daxis=[d&igp_parse_for_tf for d in daxis];y=[np.array(y)[d] for d in daxis]
    return daxis,y
def dfilt_band(mthds,mtr,usnr,isnr,p,comp='Original_Z'):
    bands=list(usnr[mtr].bands);bands.append(1);bands=np.sort(np.array(bands))
    if mtr=='snr':y0=np.array([[usnr[mtr].__dict__[m].R(comp)[p].Average((bands[bi],bands[bi+1])) for bi in range(len(bands[1:]))] for m in mthds]).squeeze()
    else:y0=np.array([[usnr[mtr].__dict__[m].Average((bands[bi],bands[bi+1])) for bi in range(len(bands[1:]))] for m in mthds]).squeeze()
    y=y0
    daxis=np.array(~np.isnan(y),dtype=int)
    positions=[np.mean([bands[bi],bands[bi+1]]) for bi in range(len(bands[1:]))]
    return y,positions,daxis


def hline(ax,y,x=None):
    lw=3;ym=pd.Series([np.median(pd.Series(i)) for i in y]).median()
    ax.axhline(ym,c='k',linewidth=lw+2,zorder=-1e10,alpha=0.2)
    ax.axhline(ym,c='gainsboro',linewidth=lw,zorder=-1e10,alpha=0.2)
def BoxPlots(xcats,usnr,isnr,mtrs=['snr','coh'],minvals=10,share_pairs=True,showfliers=False,zbins=[[0,6000]],ratio=True):
    figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
    #BOX PLOTS
    boxprops = dict(linestyle='-', linewidth=1,facecolor='gainsboro',color='k')
    flierprops = dict(marker='o', markerfacecolor='k', markersize=1,markeredgecolor='none')
    medianprops = dict(linestyle='-', linewidth=1.2, color='k',alpha=1)
    meanprops = dict(marker='D', markeredgecolor='k',markerfacecolor='gainsboro',linestyle='-', linewidth=.1,markersize=3,alpha=1.0)
    prd_boxprops = dict(linestyle='--', linewidth=0.8,facecolor=None,alpha=0.3,color='k')
    prd_flierprops = dict(marker='o', markerfacecolor='k', markersize=1,markeredgecolor='none',alpha=0.3)
    prd_medianprops = dict(linestyle='-', linewidth=0.3, color='k',alpha=0.3)
    prd_meanprops = dict(marker='D', markeredgecolor='k',markerfacecolor='k',linestyle='-',linewidth=0.0,markersize=1.8,alpha=0.3)
    gargs=lambda:dict(widths=dx*.9,showfliers=showfliers,boxprops=boxprops,medianprops=medianprops,meanprops=meanprops,showmeans=True,flierprops=flierprops)
    prd_gargs=lambda:dict(widths=dx*.9,zorder=-1e5,whiskerprops=dict(alpha=0.0),showcaps=False,showfliers=False,boxprops=prd_boxprops,medianprops=prd_medianprops,meanprops=prd_meanprops,showmeans=True,flierprops=prd_flierprops)
    status=lambda:f'{xci+1}/{len(xcats)} | {mi+1}/{len(mtrs)} | {zi+1}/{len(zbins)} || {xcat} | {mtr}'
    preferred_pbands={'P':'1_10','S':'10_30','Rg':'30_100'}
    yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
    nstd_pctile=lambda i,n=1:[sum(i<=(np.std(i)-np.mean(i)))/len(i),sum(i<=(np.std(i)+np.mean(i)))/len(i)]
    isnr['Distance']=np.array([np.round(i) for i in isnr['Distance']])
    for xci,xcat in enumerate(xcats):
        for mi,mtr in enumerate(mtrs):
            # if (xcat=='band')&(not mtr=='coh'):continue
            if xcat==mtr:continue
            for zi,zb in enumerate(zbins):
                yco={}
                zfile=f'{zb[0]}_{zb[1]}'
                if len(zbins)>1:
                    # isnr=catalog.sr.copy()
                    isnr=isnr[(isnr.StaDepth>=zb[0])&(isnr.StaDepth<=zb[1])]
                    if len(isnr)==0:continue
                    ratio=True;usnr=unpack_metrics(isnr,ratio=ratio,log=True)
                print(status())

                if xcat=='band':fig,axes=figs(1,3,f=(6,2),y='row');axes=np.atleast_2d(axes) # if mtr=='snr' else 'all'
                else:fig,axes=figs(3,3,f=(6.5,5),y='row') 
                for ki,(k,_) in enumerate(zip(['P','S','Rg',],axes[:,0])):

                    if xcat=='NoiseAverage':x=isnr[xcat];bins=np.arange(-280,-70+20,20);xl=[-270,-75]
                    if xcat=='Seismometer':x=isnr[xcat];bins=x.unique();xl=[-1,len(bins)]
                    if xcat=='Instrument_Design':x=isnr[xcat];bins=x.unique();xl=[-1,len(bins)]
                    if xcat=='StaDepth':x=isnr[xcat];bins=np.arange(0,6000+500,500);xl=[-100,6100]
                    if xcat=='Magnitude':x=isnr[xcat];bins=np.arange(6.0,8.0+.5,.5);xl=[5.9,8.1]
                    if xcat=='Distance':x=isnr[xcat];bins=np.array([30,60,90,120]);xl=[28,122]
                    if xcat=='snr':bins=np.arange(-3,8,.5);xl=[-3.1,7.6]
                    if xcat=='coh':bins=np.arange(0,1.05,.05);xl=[-0.1,1.1]
                    if xcat=='band':bins=np.arange(0,105,5);bins[0]=1;bins=bins[1:];x=bins;xl=[-1,102]
                    if isinstance(bins[0],str):dx=1
                    else:dx=max(abs(np.diff(bins)))

                    if not (xcat in ['snr','coh']):xx=x.copy()
                    if xcat=='NoiseAverage':xx=np.array([i[preferred_pbands[k]] for i in x.copy()])
                    band=f'{'-'.join(preferred_pbands[k].split('_'))}s'

                    igp_parse_for_tf=fnotch(isnr.StaDepth)>(1/int(preferred_pbands[k].split('_')[1]))
                    # igp_parse_for_tf=True

                    raxes=axes[ki,:]
                    
                    if xcat in ['StaDepth','Magnitude','NoiseAverage','coh','snr','Distance']:positions=bins[:-1]+.5*dx
                    else:positions=np.arange(0,len(bins),1)

                    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    mthds = ['TF.Z','HPS.Z','HPS.H']
                    colors=['thistle','turquoise','salmon']
                    methodcolors={m:c for m,c in zip(mthds,colors)}
                    for ri,(ax,mthd) in enumerate(zip(raxes,mthds)):
                        color = methodcolors[mthd]
                        if not xcat=='band':
                            y=usnr[mtr][mthd][k]
                            if xcat in ['snr','coh']:x=xx=usnr[xcat][mthd][k]
                            daxis,y=dfilt(xx,y,bins)
                        else:
                            y,positions,daxis=dfilt_band([mthd],isnr)
                        yco['.'.join([mthd,k,mtr])]=y
                        boxprops['facecolor']=color
                        args=gargs().copy();percents=y
                        bps=[ax.boxplot(i if len(i)>minvals else [0 if mtr=='snr' else 1],
                        positions=[positions[bi]],
                        patch_artist=True,**args) for bi,(p,i) in enumerate(zip(percents,y))]
                        alpha_empties(y,bps,0.0,minvals=minvals)
                        ax.text(0.98, 0.00,f'n={sum(np.nansum(daxis,axis=0)>0)}',transform=ax.transAxes,ha='right', va='bottom')
                        if share_pairs:
                            if not xcat=='band':
                                y=usnr[mtr][mthd][k]
                                if xcat in ['snr','coh']:x=xx=usnr[xcat][mthd][k]
                                _,y=dfilt(xx,y,bins,True if (not share_pairs) else igp_parse_for_tf)
                            else:
                                y,positions,daxis=dfilt_band([mthd],isnr)
                            prd_boxprops['facecolor']=color
                            shared_pairs_args = prd_gargs().copy();percents=y
                            bps=[ax.boxplot(i if len(i)>minvals else [0 if mtr=='snr' else 1],
                            positions=[positions[bi]],
                            patch_artist=True,**shared_pairs_args) for bi,(p,i) in enumerate(zip(percents,y))]
                            alpha_empties(y,bps,0.0,minvals=minvals)
                            xtrabox=shared_pairs_args.copy()
                            b=xtrabox['boxprops'].copy()
                            b.update({'facecolor':'None'});b.update({'alpha':1.0});b.update({'color':'k'})
                            b.update({'linestyle':':'});b.update({'linewidth':.7})
                            xtrabox.update({'boxprops':b.copy()})
                            xtrabox.update({'zorder':1e5})
                            bps=[ax.boxplot(i if len(i)>minvals else [0 if mtr=='snr' else 1],
                            positions=[positions[bi]],
                            patch_artist=True,**xtrabox) for bi,(p,i) in enumerate(zip(percents,y))]
                            alpha_empties(y,bps,minvals=minvals,sp=True)

                        if (ki==1)&(ri==0):ax.set_ylabel(f'{yttl('')}' if mtr=='coh' else r'$R_{SNR}$',fontsize=12)
                        if (ki==0):ax.set_title(f'{mthd.split('.')[0]} ({mthd.split('.')[1]})')
                        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

                    ax2=ax.twinx();ax2.set_yticks([])
                    if mtr=='coh':ax2.set_ylabel(f'{band}',rotation=-90,labelpad=15)
                    else:ax2.set_ylabel(f'{k.replace('Rg','Rayleigh')} ({band})',rotation=-90,labelpad=15)
                    if ki==0:ax.set_title(f'{mthd.split('.')[0]} ({mthd.split('.')[1]})')
                    for ax in raxes:
                        if xcat=='snr':ax.axvline(0.0,linewidth=0.5,ls=':',c='k')
                    for ax in raxes:
                        ax.set_xlim(xl)
                        if isinstance(bins[0],str):xt=np.arange(0,len(bins),1);xtlabels=bins;rot=True
                        else:
                            # xt=bins if len(bins)<6 else bins[::3]
                            xt=bins if len(bins)<6 else np.linspace(bins.min(),bins[:-2].max(),5)
                            if xcat=='snr':xt=np.arange(-3,8,2.5)
                            if xcat=='band':xt=np.array([1,20,40,60,80,100])
                            if xcat=='StaDepth':xt=np.array([0,2500,5000])
                            xtlabels=[f'{i:.1f}' for i in xt] if xcat in ['snr','coh'] else [f'{int(i)}' for i in xt];rot=False
                        if xcat in ['Seismometer']:xtlabels=['\n'.join(i.split(' ')) for i in xtlabels]
                        _=ax.set_xticks(xt)
                        # _=ax.set_xticklabels(xtlabels,rotation=90 if rot else 0)
                        _=ax.set_xticklabels(xtlabels,rotation=90 if rot else 0)
                        ax.set_xlim(xl)
                    if (ki==2):
                        if xcat=='snr':raxes[1].set_xlabel('Improvements to SNR, log10' if ratio else 'SNR, log10')
                        else:raxes[1].set_xlabel(xcat.replace('StaDepth','Water depth, m').replace('Magnitude','Magnitude, Mw').replace('NoiseAverage','Station mean noise, dB'))
                    if xcat=='band':raxes[1].set_xlabel('Period, s')
                    if showfliers:
                        if (ki==0)&(mtr=='coh'):ax.set_ylim(0.90,1.01)
                        elif mtr=='coh':ax.set_ylim(-.15,1.05)
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


def analyze_distribution(f,fn,data,method):
    """Analyze and visualize the distribution of each column in a numpy array in a single figure."""
    num_columns = data.shape[1]
    
    # Compute statistics for each column
    means = np.mean(data, axis=0)
    std_devs = np.std(data, axis=0)
    skewnesses = stats.skew(data, axis=0)
    kurtoses = stats.kurtosis(data, axis=0)
    fig = plt.figure(figsize=(10, 9))
    gs = GridSpec(3, 2, height_ratios=[1.5, 1, 1], width_ratios=[1, 1])

    ax=fig.add_subplot(gs[0, :])
    _=[ax.scatter(f,y, alpha=0.1, s=5,c='k') for y in data]
    # ax.set_title("Scatter Plot of Data")
    ax.set_xlabel("Column Index")
    ax.set_ylabel("Value")
    ax.axvline(fn, color='k', linestyle='--',linewidth=1,alpha=.3)
    ax.set_xlim([f.min(),f.max()])
    ax.set_xscale('log')

    # Plot means
    ax=fig.add_subplot(gs[1, 0])
    ax.plot(f,means, marker='o', linestyle='-', color='b', alpha=0.2,markersize=3)
    ax.set_title("Mean of Each Frequency")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Mean Value")
    ax.axvline(fn, color='k', linestyle='--',linewidth=1,alpha=.3)
    ax.set_xlim([f.min(),f.max()])
    ax.set_xscale('log')
    
    # Plot standard deviations
    ax=fig.add_subplot(gs[1, 1])
    ax.plot(f,std_devs, marker='o', linestyle='-', color='r', alpha=0.2,markersize=3)
    ax.set_title("Standard Deviation of Each Frequency")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Standard Deviation")
    ax.axvline(fn, color='k', linestyle='--',linewidth=1,alpha=.3)
    ax.set_xlim([f.min(),f.max()])
    ax.set_xscale('log')
    
    # Plot skewness
    ax=fig.add_subplot(gs[2, 0])
    ax.plot(f,skewnesses, marker='o', linestyle='-', color='g', alpha=0.2,markersize=3)
    ax.set_title("Skewness of Each Frequency")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Skewness")
    ax.axvline(fn, color='k', linestyle='--',linewidth=1,alpha=.3)
    ax.set_xlim([f.min(),f.max()])
    ax.set_xscale('log')
    # Plot kurtosis
    ax=fig.add_subplot(gs[2, 1])
    ax.plot(f,kurtoses, marker='o', linestyle='-', color='purple', alpha=0.2,markersize=3)
    ax.set_title("Kurtosis of Each Frequency")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Kurtosis")
    ax.axvline(fn, color='k', linestyle='--',linewidth=1,alpha=.3)
    ax.set_xlim([f.min(),f.max()])
    ax.set_xscale('log')
    fig.suptitle(method)
    plt.tight_layout()
    plt.show()

def get_gridplot():
    fig = plt.figure(figsize=(20,15),layout="constrained")
    height_ratios = [1,1,0.6,0.6,0.6,0.6]
    gs = gridspec(6, 6, figure=fig,height_ratios=height_ratios)
    ax_23_0=fig.add_subplot(gs[0,0:2]);ax_23_1=fig.add_subplot(gs[0,2:4]);ax_23_2=fig.add_subplot(gs[0,4:])
    ax_23_3=fig.add_subplot(gs[1,0:2]);ax_23_4=fig.add_subplot(gs[1,2:4]);ax_23_5=fig.add_subplot(gs[1,4:])
    ax23 = np.array([[ax_23_0,ax_23_1,ax_23_2],[ax_23_3,ax_23_4,ax_23_5]])
    ax_42_0=fig.add_subplot(gs[2,:3]);ax_42_1=fig.add_subplot(gs[2,3:])
    ax_42_2=fig.add_subplot(gs[3,:3]);ax_42_3=fig.add_subplot(gs[3,3:])
    ax_42_4=fig.add_subplot(gs[4,:3]);ax_42_5=fig.add_subplot(gs[4,3:])
    ax_42_6=fig.add_subplot(gs[5,:3]);ax_42_7=fig.add_subplot(gs[5,3:])
    ax42 = np.array([[ax_42_0,ax_42_1],[ax_42_2,ax_42_3],[ax_42_4,ax_42_5],[ax_42_6,ax_42_7]])
    return fig,ax23,ax42

def station_event_page(st_hold,sta,evmeta,method,type='stream',**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[1,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')];defargs.Noise=True
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    if args.csd_pairs[0][0][0]==args.csd_pairs[0][0][1]:args.Noise=False
    # ------------
    # ------------
    if type.lower()=='stream':
        note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:
        note = 'Noise (gray) | Raw (red)'
    staname = st_hold[0].stats.network+'.'+st_hold[0].stats.station
    stastr = ' | '.join([method.upper(),staname,sta.Experiment,
    'Depth: '+str(int(1000*abs(st_hold[0].stats.sac.stel)))+'m',
    'F-Notch: '+str(int(1/fnotch(1000*abs(st_hold[0].stats.sac.stel))))+'s',
    note])
    if not args.nev:args.nev = len(st_hold.select(location='*Raw*'))
    if type.lower()=='stream':columns = args.bands
    elif type.lower()=='metrics':columns=['Coherence','Phase']
    nev_per_plot = args.nev
    nplots = int(np.ceil(nev_per_plot/len(evmeta)))
    for ploti in range(nplots):
        fig,axes = plt.subplots(nrows=nev_per_plot,ncols=len(columns),layout='constrained',sharex='all',squeeze=True,figsize=(args.figwidth*args.figaspect[0],nev_per_plot*args.figaspect[1]))
        axes = np.atleast_2d(axes).reshape((nev_per_plot,len(columns)))
        fig.suptitle(stastr)
        fn = 1/fnotch(1000*abs(st_hold[0].stats.sac.stel))
        for bi,b in enumerate(columns):
            # print(method+'-'+type + '| Column:'+str(bi+1)+'/'+str(len(columns)))
            band_ax = axes[:,bi]
            st_band = st_hold.copy()
            correct_hold = st_band.select(location='*Correct*').copy()
            raw_hold = st_band.select(location='*Raw*').copy()
            if type.lower()=='stream':
                st_band.filter('bandpass',freqmin=1/b[1],freqmax=1/b[0],zerophase=True,corners=4)
                correct_hold = st_band.select(location='*Correct*').copy()
                raw_hold = st_band.select(location='*Raw*').copy()
                # --------------
                # This handles large scale amplitude differences between the two traces.
                # These are traces where the raw amplitudes exceed the vertical_scale threshold wrt the corrected trace amplitudes.
                ylim = np.array([np.max(np.abs(c.data))*args.vertical_scale for c in correct_hold])
                out_scaled = np.array([np.max(np.abs(c.data)) for c in raw_hold]) > ylim
                if np.any(out_scaled):
                    print('Large amplitude scale differences detected in '+staname)
                    for tr_ind,tr in enumerate(raw_hold):
                        if out_scaled[tr_ind]:tr.data = (tr.data/np.max(np.abs(tr.data)))*ylim[tr_ind]
                # --------------
            for si,s in enumerate(['Raw','Corrected']):
                st = {'Raw':raw_hold,'Corrected':correct_hold}[s]
                for tri in range(nev_per_plot):
                    tr_ind = tri*(ploti+1)
                    tr = st[tr_ind].copy()
                    tr.trim(evmeta[tr_ind].origins[0].time,evmeta[tr_ind].origins[0].time+7200,pad=True,fill_value=0).taper(0.0001)
                    colors = args.phasecolors
                    ax = band_ax[tr_ind]
                    if tr_ind==0:
                        if type.lower()=='stream':ax.set_title(''.join([str(b[0]),'s-',str(b[1]),'s']))
                        else:ax.set_title(b+':'+args.csd_pairs[0][0])
                    if type.lower()=='stream':
                        x = tr.times('relative');y = tr.data;ax.plot(x,y,color=args.linecolor[si],alpha=args.alpha[si],linewidth=args.linewidth[si])
                    elif type.lower()=='metrics':
                        x = tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[0];ind=x<=1
                        y=tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[1][ind]
                        if s=='Raw':
                            if args.Noise:[ax.scatter(xy[0],xy[1],c='darkgrey',s=0.1,label=args.csd_pairs[pi][0]+':Noise') 
                            for pi,xy in enumerate([avg_meter(st_hold.Noise,b,p[0]) for p in args.csd_pairs])]
                            [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.8,color='r') for p in args.csd_pairs]
                            [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linewidth=0.2,alpha=0.4,color='r') for p in args.csd_pairs]
                        if s=='Corrected':
                            [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.4,color=p[1]) for p in args.csd_pairs]
                            [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linestyle=':',linewidth=0.45,alpha=0.4,color=p[1]) for p in args.csd_pairs]
                    if type.lower()=='metrics':
                        ax.set_xlim(1/500,1);ax.set_xscale('log')
                        if b.lower()=='phase':ax.set_ylim(-180,180)
                        if b.lower()=='coherence':ax.set_ylim(0,1.1)
                    else:ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim[tr_ind],ylim[tr_ind]);ax.set_yticklabels('')
                    if s=='Corrected':
                        stallaz=[tr.stats.sac.stla,tr.stats.sac.stlo,tr.stats.sac.stel]
                        evllaz=[evmeta[tr_ind].origins[0].latitude,evmeta[tr_ind].origins[0].longitude,evmeta[tr_ind].origins[0].depth/1000]
                        tr.stats.sac.gcarc = lt.math.distance(sta,evmeta[tr_ind])
                        evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name,'M'+str(evmeta[tr_ind].magnitudes[0].mag),str(np.round(tr.stats.sac.gcarc,2))+'°'])
                        ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='center',fontsize=10)
                        if type.lower()=='stream':
                            if tr.stats.sac.gcarc<=100:phases=args.phases
                            else:phases=args.shadow_phases;[colors.update({p:'k'}) for p in phases]
                            ar = get_arrivals(stallaz,evllaz,model = 'iasp91',phases=phases)
                            [ax.axvline(a[1],linewidth=0.1,color=colors[a[0]]) for a in ar]
                            [ax.text(a[1],y.max(),a[0],fontsize=6,color='k',verticalalignment='bottom',horizontalalignment='center') for a in ar]
            if type.lower()=='metrics':ax.set_xlabel('frequency (hz)')
            else:ax.set_xlabel('seconds')
    return fig
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# ----------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def station_event_page_averages(st_hold,sta,evmeta,type='Metrics',raw_reference=None,**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = {'ZP':'#0c51a6','ZZ':'#2a7e93','Z1':'#7370cb','Z2':'#4f86c5'}
    # ------------
    if type.lower()=='stream':note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:note = '\nVariance (shaded)'
    staname = st_hold[0].stats.network+'.'+st_hold[0].stats.station
    # if method.lower()=='atacr':tf='('+st_hold.select(location='*Correct*')[0].stats.location.split('.')[1]+')'
    tf=''
    stastr = ' | '.join([staname+tf,sta.Experiment,'Depth: '+str(int(1000*abs(st_hold[0].stats.sac.stel)))+'m, Notch: '+str(int(1/fnotch(1000*abs(st_hold[0].stats.sac.stel))))+'s',note])
    if not args.nev:args.nev = len(st_hold.select(location='*Raw*'))
    nrows = 1;ncols=2
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # if type.lower()=='metrics':columns=['Coherence']
    # if type.lower()=='metrics':columns=['Coherence','Phase','Admittance']
    fig,axes = plt.subplots(nrows=ncols,ncols=nrows,layout='constrained',sharex='all',squeeze=True,figsize=(8,7))
    axes = axes.reshape(-1)
    fig.suptitle(stastr)
    fn = 1/fnotch(1000*abs(st_hold[0].stats.sac.stel))
    methods=['ATaCR','NoiseCut']
    for bi,b in enumerate(columns):
        method = methods[bi]
        metric = b[1][0]
        pair = b[1][1]
        # for pi,pair in enumerate(args.csd_pairs):
        channels = pair
        if channels[0]==channels[1]:args.Noise=False
        else: args.Noise=True
        # print(method+'-'+type + '| Column:'+str(bi+1)+'/'+str(len(columns)))
        ax = axes[bi]
        st_band = st_hold.copy()
        correct_hold = st_band.select(location=f'*{method}*').copy()
        raw_hold = st_band.select(location='*Raw*').copy()
        if raw_reference:raw_hold = raw_reference
        for si,s in enumerate(['Corrected']):
            st = {'Raw':raw_hold,'Corrected':correct_hold}[s].copy()
            for tri in range(len(correct_hold)):
                tr_ind = tri
                tr = st[tr_ind].copy()
                if tr_ind==0:
                    if type.lower()=='stream':ax.set_title(''.join([str(metric),'s-',str(pair),'s']))
                    else:ax.set_title(f'{method} {pair} {metric}')
                if type.lower()=='metrics':
                    x = tr.Metrics.__getattribute__(metric)(pair)[0];ind=x<=1
                    y=tr.Metrics.__getattribute__(metric)(pair)[1][ind]
                    if s=='Raw':
                        if args.Noise:[
                        ax.scatter(xy[0],np.abs(xy[1]),c='darkgrey',s=0.1,label=pair+':Noise')
                        for ni,xy in enumerate([avg_meter(st_hold.Noise,metric,pair) for p in [1]])]
                        [ax.scatter(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),s=0.4,color='r',alpha=0.1) for p in [1]]
                        [ax.plot(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),linewidth=0.05,alpha=0.05,color='r') for p in [1]]
                    if s=='Corrected':
                        [ax.scatter(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),s=0.2,color=args.csd_pairs[pair],alpha=0.1) for p in [1]]
                        [ax.plot(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),linestyle=':',linewidth=0.05,alpha=0.05,color=args.csd_pairs[pair]) for p in [1]]
                if type.lower()=='metrics':
                    ax.set_xlim(1/500,1);ax.set_xscale('log')
                    # if metric.lower()=='phase':ax.set_ylim(-180,180)
                    if metric.lower()=='phase':ax.set_ylim(0,180)
                    if metric.lower()=='coherence':ax.set_ylim(0,1.01)
                else:ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim[tr_ind],ylim[tr_ind]);ax.set_yticklabels('')
                # ------------------------------------------------------------------------------------------
            y_avg = np.mean(([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)
            y_var = np.std(([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)**2

            # Should I set Phase domain to [0,180] instead of [-180,180]? Would be easier to read...
            # y_avg = np.mean(np.abs([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)
            # y_var = np.std(np.abs([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)**2
            if s=='Corrected':
                stallaz=[tr.stats.sac.stla,tr.stats.sac.stlo,tr.stats.sac.stel]
                evllaz=[evmeta[tr_ind].origins[0].latitude,evmeta[tr_ind].origins[0].longitude,evmeta[tr_ind].origins[0].depth/1000]
                tr.stats.sac.gcarc = lt.math.distance(sta,evmeta[tr_ind])
                # if type.lower()=='metrics':
                    # evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name])
                # el:
                if not type.lower()=='metrics':
                    evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name,'M'+str(evmeta[tr_ind].magnitudes[0].mag),str(np.round(tr.stats.sac.gcarc,2))+'°'])
                    ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='center',fontsize=10)
                # ------
                ax.scatter(x[ind],y_avg,label=':'.join([s,pair]),s=0.5,color=args.csd_pairs[pair])
                ax.plot(x[ind],y_avg,label=':'.join([s,pair]),linewidth=0.8,alpha=0.8,color=args.csd_pairs[pair])
                ax.fill_between(x[ind],y_avg-y_var,y_avg+y_var, alpha=0.2,color=args.csd_pairs[pair])
            elif s=='Raw':
                ax.scatter(x[ind],y_avg,label=':'.join([s+'-Average',pair]),s=0.5,color='r')
                ax.plot(x[ind],y_avg,label=':'.join([s+'-Average',pair]),linewidth=0.8,alpha=0.5,color='r')
        if type.lower()=='metrics':ax.set_xlabel('frequency (hz)')
        else:ax.set_xlabel('seconds')
        ax.axvline(1/fn,alpha=0.4,linewidth=1,color='k',linestyle='-.')
        ax.text(1/fn,1,str(int(fn))+'s',verticalalignment='top',horizontalalignment='right')
    return fig


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# ----------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def station_event_single(st_hold,sta,evmeta,type='stream',**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[1,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')];defargs.Noise=True
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    if args.csd_pairs[0][0][0]==args.csd_pairs[0][0][1]:args.Noise=False
    # ------------# ------------
    method=''
    if type.lower()=='stream':
        note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:
        note = 'Noise (gray) | Raw (red)'
    staname = st_hold.Raw.stats.network+'.'+st_hold.Raw.stats.station
    stastr = ' | '.join([method.upper(),staname,sta.Experiment,
    'Depth: '+str(int(1000*abs(st_hold.Raw.stats.sac.stel)))+'m',
    'F-Notch: '+str(int(1/fnotch(1000*abs(st_hold.Raw.stats.sac.stel))))+'s',note])
    if not args.nev:args.nev = len([st_hold.Raw])
    if type.lower()=='stream':columns = args.bands
    elif type.lower()=='metrics':columns=['Coherence','Phase']
    nev_per_plot = args.nev
    nplots = int(np.ceil(nev_per_plot/len(evmeta)))

    channel_color=AttribDict()
    channel_color.ZP='#0c51a6';channel_color.ZZ='#2a7e93'
    channel_color.Z1='#7370cb';channel_color.Z2='#4f86c5'

    # XX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XX
    fig,ax23,ax42 = get_gridplot()
    nax = ax42.size + ax23.size
    plt.tight_layout()
    axis = [0,1,2,3,4,5,0,1,2,3,4,5,6,7]
    fig.suptitle(stastr,y=1.025)
    fn = 1/fnotch(1000*abs(st_hold.Raw.stats.sac.stel))
    meters = ['Coherence','Phase','Coherence','Phase','Coherence','Phase','Coherence','Phase']
    channels = ['ZP','ZP','ZZ','ZZ','Z1','Z1','Z2','Z2']
    bands = [[1,10],[10,30],[30,100],[1,10],[10,30],[30,100]]
    methods = ['ATaCR','ATaCR','ATaCR','Noisecut','Noisecut','Noisecut']
    # XX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XX

    # x = tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[0];ind=x<=1
    # y=tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[1][ind]
    # if s=='Raw':
    #     if args.Noise:[ax.scatter(xy[0],xy[1],c='darkgrey',s=0.1,label=args.csd_pairs[pi][0]+':Noise') 
    #     for pi,xy in enumerate([avg_meter(st_hold.Noise,b,p[0]) for p in args.csd_pairs])]
    #     [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.8,color='r') for p in args.csd_pairs]
    #     [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linewidth=0.2,alpha=0.4,color='r') for p in args.csd_pairs]
    # if s=='Corrected':
    #     [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.4,color=p[1]) for p in args.csd_pairs]
    #     [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linestyle=':',linewidth=0.45,alpha=0.4,color=p[1]) for p in args.csd_pairs]
    # if type.lower()=='metrics':
    # ax.set_xlim(1/500,1);ax.set_xscale('log')
    # if b.lower()=='phase':ax.set_ylim(-180,180)
    # if b.lower()=='coherence':ax.set_ylim(0,1.1)

    # if args.Noise:
    #     [ax.scatter(xy[0],np.abs(xy[1]),c='darkgrey',s=0.1,label=channels+':Noise') 
    #     for ni,xy in enumerate([avg_meter(st_hold.Noise,b,channels) for p in [pair]])]
    for axi in range(nax):
        print(str(evmeta.indi)+'-'+staname+'- Single event page- '+str(axi+1)+'/'+str(nax))
        if axi<ax23.size:
            ax=ax23.reshape(-1)[axis[axi]];c_band=bands[axi];c_method=methods[axi]
            # ------------
            # Trace section
            # ------------
            raw = st_hold.Raw.copy()
            atacr = st_hold.ATaCR.copy()
            hps = st_hold.Noisecut.copy()

            raw.taper(0.0001)
            raw.filter('bandpass',freqmin=1/c_band[1],freqmax=1/c_band[0],zerophase=True,corners=4)
            raw.trim(evmeta[0].origins[0].time,evmeta[0].origins[0].time+7200,pad=True,fill_value=0)

            hps.taper(0.0001)
            hps.filter('bandpass',freqmin=1/c_band[1],freqmax=1/c_band[0],zerophase=True,corners=4)
            hps.trim(evmeta[0].origins[0].time,evmeta[0].origins[0].time+7200,pad=True,fill_value=0)

            atacr.taper(0.0001)
            atacr.filter('bandpass',freqmin=1/c_band[1],freqmax=1/c_band[0],zerophase=True,corners=4)
            atacr.trim(evmeta[0].origins[0].time,evmeta[0].origins[0].time+7200,pad=True,fill_value=0)

            if c_method=='ATaCR':
                method_tr=atacr
            else:
                method_tr=hps
            ylim = np.max(np.array(np.max([np.max(np.abs(c.data))*args.vertical_scale for c in [atacr]])))
            out_scaled = np.max(np.abs(raw.data)) > ylim
            if np.any(out_scaled):
                print('Large amplitude scale differences detected in '+staname)
                raw.data = (raw.data/np.max(np.abs(raw.data)))*ylim
            x = raw.times()
            ax.plot(raw.times(),raw.data,color='red',alpha=1.0,linewidth=0.1)
            ax.plot(method_tr.times(),method_tr.data,color='black',alpha=1.0,linewidth=0.2)
            ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim,ylim);ax.set_yticklabels('')
            ax.set_title(''.join([str(c_band[0]),'s-',str(c_band[1]),'s']),bbox=dict(facecolor='white', alpha=1),pad=0,verticalalignment='baseline',loc='left',y=0)
            if axi>2:
                ax.set_xlabel('seconds')
            else:
                ax.set_xticklabels('')

            # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            stallaz=[raw.stats.sac.stla,raw.stats.sac.stlo,raw.stats.sac.stel]
            evllaz=[evmeta[0].origins[0].latitude,evmeta[0].origins[0].longitude,evmeta[0].origins[0].depth/1000]
            raw.stats.sac.gcarc = locations2degrees(sta,evmeta)
            evstr = '|'.join(['['+str(evmeta.indi)+'] '+c_method,evmeta[0].Name,'M'+str(evmeta[0].magnitudes[0].mag),str(np.round(raw.stats.sac.gcarc,2))+'°'])
            ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='bottom',fontsize=10)
            if type.lower()=='stream':
                colors = args.phasecolors
                if raw.stats.sac.gcarc<=100:phases=args.phases
                else:phases=args.shadow_phases;[colors.update({p:'k'}) for p in phases]
                ar = get_arrivals(stallaz,evllaz,model = 'iasp91',phases=phases)
                [ax.axvline(a[1],linewidth=0.1,color=colors[a[0]]) for a in ar]
                [ax.text(a[1],ylim,a[0],fontsize=6,color='k',verticalalignment='bottom',horizontalalignment='center') for a in ar]
            # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        else:
            ax=ax42.reshape(-1)[axis[axi]];c_chan=channels[axis[axi]];c_meter=meters[axis[axi]]
            # ------------
            # Metrics section
            # ------------
            method_color=dict()
            method_color['ATaCR']='magenta'
            method_color['Noisecut']='blue'
            for ci,c_method in enumerate(['ATaCR','Noisecut']):
                tr = st_hold[c_method]
                x = tr.Metrics.__getattribute__(c_meter)(c_chan)[0]
                ind=x<=1
                y=tr.Metrics.__getattribute__(c_meter)(c_chan)[1][ind]
                if ci==0:
                    if not c_chan=='ZZ':
                        if args.Noise:
                            [ax.scatter(xy[0],xy[1],c='darkgrey',edgecolor='darkgrey',facecolor='darkgrey',s=0.1,label='Noise') 
                            for pi,xy in enumerate([avg_meter(st_hold.Noise,c_meter,p) for p in [c_chan]])]    
                            # [ax.plot(xy[0],xy[1],c='darkgrey',linewidth=0.3,label='Noise') 
                            # for pi,xy in enumerate([avg_meter(st_hold.Noise,c_meter,p) for p in [c_chan]])]
                    [ax.scatter(x[ind],st_hold.Raw.Metrics.__getattribute__(c_meter)(p)[1][ind],s=0.8,color='r') for p in [c_chan]]
                    [ax.plot(x[ind],st_hold.Raw.Metrics.__getattribute__(c_meter)(p)[1][ind],label='Raw',linewidth=0.8,alpha=0.4,color='r') for p in [c_chan]]
                [ax.scatter(x[ind],tr.Metrics.__getattribute__(c_meter)(p)[1][ind],s=2,facecolor=channel_color[c_chan],edgecolors=method_color[c_method]) for p in [c_chan]]
                [ax.plot(x[ind],tr.Metrics.__getattribute__(c_meter)(p)[1][ind],label=c_method,linewidth=0.8,alpha=0.4,color=method_color[c_method]) for p in [c_chan]]
                ax.set_xlim(1/500,1);ax.set_xscale('log')
                if c_meter.lower()=='phase':ax.set_ylim(-180,180)
                if c_meter.lower()=='coherence':ax.set_ylim(0,1.1)
                ax.set_title(c_meter+':'+c_chan,bbox=dict(facecolor='white', alpha=1),pad=0,verticalalignment='baseline',loc='left',y=0)
                if axi>=12:
                    ax.set_xlabel('frequency (hz)')
                if axi==6:
                    lg = ax.legend(ncols=2,loc='upper left')
                    [l.set_linewidth(4) for l in lg.legendHandles]
                    [l.set_sizes([5,5,5]) for l in [lg.legendHandles[0]]]
                if ci==0:
                    ax.axvline(1/fn,alpha=0.4,linewidth=1,color='k',linestyle='-.')

    return fig

# def dataset_averaged_coherence_plot(cat,dirs,
#     title='Station Averaged Noise',
#     meter='Coherence',
#     chans='ZP',
#     figsize=[10,4],fontsize=12,lw=None):
#     font = {'weight':'bold','size':fontsize};matplotlib.rc('font', **font)
#     f = avg_meter(get_noise(dirs,cat.StaName[0]),meter,chans)[0]
#     coh = np.array([avg_meter(get_noise(dirs,stanm),meter,chans)[1] 
#     for stanm in cat.StaName[np.argsort(cat.StaDepth)]])
#     x = np.round(np.sort(cat.StaDepth))
#     # ---------------------------------------------------------------
#     # More spatially accurate approach but matrix is sparse.
#     # xx = np.arange(int(x.min()),int(x.max())+1,int(np.diff(x).min()))
#     # zz = np.empty((len(xx),len(f)));zz[:] = 0
#     # for ii,xi in enumerate(x):zz[np.where(xi==xx)[0][0],:] = cp[ii,:]
#     # ---------------------------------------------------------------
#     fig,ax=plt.subplots(figsize=figsize);x = np.round(np.sort(cat.StaDepth))
#     # cp = np.array([smooth(d,k=3) for d in coh]);cp = gaussian_filter(coh,.5)
#     cp = coh
#     ax.contourf(f,x,cp,linewidth=2,levels=10,extend='max',vmin=coh.min(),vmax=coh.max(),linewidths=lw)
#     fn = [fnotch(fq,n=1) for fq in x]
#     ax.plot(fn,x,linestyle='dashed',color='w',linewidth=0.4*fontsize)
#     ax.set_xlim(1/500,1);ax.set_xscale('log')
#     fticks = np.array([1/500,1/300,1/100,1/50,1/10,1])
#     ax.set_xticks(fticks)
#     ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#     ax.set_xticklabels(np.array(1/fticks,dtype=int))
#     ax.set_ylabel('Station depth (m)',fontweight='bold')
#     ax.set_xlabel('Period (s)',fontweight='bold')
#     ax.set_title(f'{coh.shape[0]} {title}')
#     plt.tight_layout()
#     return fig
def catalog_map_plot(sta_inv,ev_cat,
    mode='Depth',
    water_fill_color='lightblue',
    continent_fill_color='lightgrey',
    projection = 'ortho',
    figsize = (15,8),
    cmap=None):
    # projection = ['ortho','global']
    inv = sta_inv.copy()
    evm_plottable = ev_cat.copy()
    if cmap==None:cmap={'mag':'bwr','depth':'tab20c'}[mode.lower()]
    if mode.lower()=='mag':
        for e in evm_plottable:e.origins[0].depth = e.magnitudes[0].mag*1000
    clear_output(wait=False)
    nev = len(evm_plottable)
    nsta = len(np.unique(list(itertools.chain.from_iterable([e.Stations for e in evm_plottable]))))
    fig = evm_plottable.plot(resolution='f',
    water_fill_color=water_fill_color,continent_fill_color=continent_fill_color,
    projection=projection,
    color='depth',label=None)
    settings = AttribDict();settings.map_ortho = AttribDict();settings.map_global=AttribDict()
    settings.map_ortho.shrink=0.4;settings.map_global.shrink=0.7
    title=' | '.join([str(nev)+' Events',str(nsta)+ ' Stations'])
    clear_output(wait=False)
    events_plot = fig.get_axes()[0].get_children()[-14]
    events_plot.set_cmap(cmap)
    events_plot.set_edgecolor('k');events_plot.set_linewidth(0.3)
    sizes = events_plot.get_sizes();events_plot.set_sizes(sizes/3)
    fig.set_size_inches(figsize);fig.set_tight_layout('tight')
    norm = fig.get_axes()[1]._colorbar.norm
    if mode.lower()=='mag':norm._vmin=6.0;norm._vmax=8.0
    fig.get_axes()[1].remove()
    if len(np.unique([e.magnitudes[0].magnitude_type for e in evm_plottable]))>1:mtype='M'
    else:mtype=np.unique([e.magnitudes[0].magnitude_type for e in evm_plottable])[0]
    if mode.lower()=='mag':colorbarlabel = 'Magnitude ('+mtype+')'
    else:colorbarlabel = 'Depth (km)'
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    ax=fig.get_axes()[0], orientation='horizontal',
    label=colorbarlabel,
    shrink=settings['map_'+projection].shrink,aspect=75,pad=0.01)
    plt.draw()
    fig.suptitle(title,y=0.93)
    # ------------
    # -- Station plot
    fig = inv.plot(size=50,label=False,fig=fig,color_per_network=True,marker='v')
    stations_plot = fig.get_axes()[0].get_children()[-13]
    stations_plot.set_edgecolor('k')
    stations_plot.set_linewidth(0.5)
    fig.show()
    lg = fig.get_axes()[0].get_legend()
    labels=lg.texts;handles=lg.legend_handles;lg.remove()
    labels = [l._text for l in labels]
    lg = fig.get_axes()[0].legend(handles, labels, ncols=11,loc='lower left',
    labelspacing=0.0,columnspacing=0.1,
    handletextpad=-.5,borderpad=0.2,edgecolor='k',fontsize=11)
    lg.set_zorder(1000)
    fig.show()
    return fig


def plot_spec_coh_adm_ph(Metrics):
    pairs = ['ZP']
    meters = ['psd','Coherence','Admittance','Phase']
    fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(8,10),layout='constrained',squeeze=False,sharey='row',sharex='all')
    axes = axes.reshape(-1)
    stam = Metrics['ATaCR'].traces[0].stats.network + '.' + Metrics['ATaCR'].traces[0].stats.station
    label = Metrics['ATaCR'].traces.select(channel='*Z')[0].stats.location
    Pre = Metrics['Raw']
    Post = Metrics['ATaCR']
    Noise = Metrics['Noise']
    fn = fnotch(abs(Post.traces[0].stats.sac.stel*1000))
    tstamp = Pre.traces[0].stats.starttime.strftime('%Y.%j.%H.%M')
    for pi,(ax,m) in enumerate(zip(axes,meters)):
        if m=='psd':
            p = 'Z'
        else:
            p = pairs[0]
        evf,prey = Pre.__getattribute__(m)(p)
        evf,posty = Post.__getattribute__(m)(p)
        if m=='psd':
            noisef,noisey = Noise.f,Noise.StaNoise.power.__dict__['cZZ']
        else:
            noisef,noisey = Noise.__getattribute__(m)(p)
        noisey = noisey[noisef>0]
        noisef = noisef[noisef>0]
        if m=='psd':
            noisey = 10*np.log10(noisey)
            prey = 10*np.log10(prey)
            posty = 10*np.log10(posty)
        ax.scatter(noisef,noisey,s=0.5,c='gray',label='Noise')
        if pi==0:
            lbl = stam + ' | ' + label + ' ' + tstamp + ' | '
        else:
            lbl = ''
        ax.set_xlabel('Frequency')
        ax.set_ylabel(m.replace('psd','Power Density'))
        ax.set_title(lbl + p + '-' + m.replace('psd','PSD'),fontweight='bold')
        ax.scatter(evf,prey,c='k',label='PRE',marker='o',s=1)
        ax.scatter(evf,posty,c='m',label='POST',marker='o',s=0.5)
        ax.axvline(fn,linewidth=0.2,color='k')
        if pi==0:
            ax.text(fn*1.05,0.99*min(ax.get_ylim()),'Fn:' + str(round(1/fn*100)/100) + 's',alpha=0.4)
            ax.set_xscale('log')
            ax.set_xlim(evf[1],evf[-1])
            ax.legend(markerscale=10,ncols=len(meters))
    plt.tight_layout()
    return fig
def dataset_averaged_coherence_plot(f,z,coh,ms=35,zconnect=False,figsize=[6,6],cmap='viridis',checkerboard=(40,9),
    title='Station Averaged Coherence',fontsize=5,levels=None,fig=None,ax=None,octav=True,fmin=1/100):
    font = {'weight':'bold','size':fontsize};matplotlib.rc('font', **font)
    z = np.round(z)
    i = np.argsort(z)
    z,coh = z[i],coh[i,:]
    if levels is None:levels=np.linspace(np.min(coh),np.max(coh),20)
    if ax is None:fig,ax=plt.subplots(figsize=figsize)
    cp = coh
    # cp = np.array([smooth(d,k=3) for d in coh]);cp = gaussian_filter(coh,.5)

    checkerdensity,brightness=checkerboard
    checkerboard = np.indices((checkerdensity, checkerdensity*2)).sum(axis=0) % 2  # 0/1 checker pattern
    from matplotlib.colors import LinearSegmentedColormap
    dark_greys = LinearSegmentedColormap.from_list('dark_greys', ['#777', f'#{str(brightness)*3}'])
    ax.imshow(checkerboard,cmap=dark_greys,interpolation='nearest',aspect='auto',extent=[0,1,0,1],origin='lower',
    zorder=-1e3,transform=ax.transAxes,
    alpha=1.0)  # Adjust transparency as needed)

    s=ms
    faxis=(f>=(fmin))&(f<=1)
    F=f.copy()[faxis]
    C=np.asarray(cp[:,faxis])            # same shape

    if octav:F,C = octavg(C,F)
    C=np.array([C[z==zi,:].mean(axis=0) for zi in np.unique(z)])
    if zconnect:F, Z = np.meshgrid(F, np.arange(0,len(np.unique(z))))
    else:F, Z = np.meshgrid(F, np.unique(z)) 

    norm = mpl.colors.Normalize(vmin=0, vmax=1.0)
    # F,Z,C=np.flip(F,axis=0),np.flip(Z,axis=0),np.flip(C,axis=0) #Looks better with this off
    sc = ax.scatter(F.ravel(), Z.ravel(), c=C.ravel(), marker='s', s=s, cmap=cmap, norm=norm, edgecolors='none',linewidths=0)
    ax.set_xscale('log')

    fn_z = np.linspace(z.min(),z.max(),len(z)*10)
    fn = [fnotch(i) for i in fn_z]
    ax.plot(fn,fn_z,linestyle=':',color='k',linewidth=0.4*fontsize)
    ax.set_xlim(fmin,1);ax.set_xscale('log')
    n=1 if zconnect else 70;ax.set_ylim(np.min(Z)-n,np.max(Z)+n)
    fticks = np.array([1/100,1/50,1/30,1/10,1])
    ax.set_xticks(fticks)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticklabels(np.array(1/fticks,dtype=int))
    if zconnect:
        zy=np.arange(0,6500,1500)
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
    norm = mcolors.TwoSlopeNorm(vmin=0, vcenter=0.5, vmax=1)
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax,extend=None)
    cbar.set_ticks(np.arange(0, 1.01, 0.25))
    if fig is not None:return fig
    if ax is not None:return ax

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| HPS Spectrogram Plots (Replicates Zali) |||||||||||||||||||||||||
# _________________________________________________________________________________________________________________

def hps_spectrograms(S_full, S_background, S_hps, frequencies, times,figsize=(11,12)):
        t = times
        f = frequencies
        s = S_full
        cmap = 'magma'
        xlabel = True
        yscale = 'log'
        ax = None
        vmin,vmax = None,None
        fig, axes = plt.subplots(nrows=3, ncols=1,figsize=figsize,height_ratios=[1,1,1],width_ratios=[1],layout='constrained',squeeze=False,sharey='col',sharex='row')
        titles = ['Raw','Noise','Noise Removed']
        for r,s in enumerate([S_full, S_background, S_hps]):
                gax = axes[r,0]
                pc = gax.pcolormesh(t,f, 10*np.log10(s), cmap = cmap, shading= 'auto')
                if (vmin is not None) and (vmax is not None):
                        pc.set_clim(vmin, vmax)
                else:
                        vmin,vmax = pc.get_clim()
                if ylabel:
                        gax.set_ylabel('Frequency (Hz)',fontweight='bold')
                if xlabel:
                        gax.set_xlabel('Time (s)',fontweight='bold')
                gax.set_yscale(yscale)
                gax.set_ylim(f[1],f[-1])
                gax.set_xlim(t[0],t[-1])
                gax.set_title(titles[r],fontweight='bold')
                if fig is not None:
                        fig.colorbar(pc, ax=gax, pad=0.01, label='dB')

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| EVENT RECORD FUNCTION |||||||||||||||||||||||||||||||||||||||||||
# _________________________________________________________________________________________________________________
def event_record_plot(evstream,evstream_back=None,linewidth=0.2,trim = (None,None),scales = [1,1],band = None,facecolor=('b','r'),norm = 'trace',figsize = (20,13),phases = ('P','S','SKS','PKiKP','SKiKS','SKSSKS',),evdepth=None,title='',sortindex=None,ax=None,normscale=1.0,residual_fraction=1.0):
        # linewidth=0.2
        # prepare traces
        if evstream_back:
            sets = [evstream_back,evstream] # [UNCORRECT_SET , CORRECT_SET] 
        else:
            sets = [evstream]
        sets = [preparetraces(stream,trim=trim,band=band,sortindex=sortindex) for stream in sets]
        if len(sets)==2:
            residuals = [ev0.data - ev.data for ev0,ev in zip(sets[0],sets[1])]
            residuals = [normscale * residual_fraction * np.array(res) / np.max(np.abs(np.array(res))) for res in residuals]
        fig = None
        if ax is None:
            fig, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
            ax = axes[0,0]
        # norming for plots
        postsetmax = [np.max(ev.data) for ev in  sets[-1]]
        for si in range(len(sets)):
            for i in range(len(sets[0])): #norm the uncorrected set to the max in the corrected set
                # -----------
                if isinstance(norm,list):
                        sets[si][i].data = sets[si][i].data/norm[i]
                elif norm.lower()=='postset':
                        sets[si][i].data = sets[si][i].data/postsetmax[i]
                elif norm.lower()=='trace':
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                elif norm.lower()=='col':
                        # print('Trace norm scaling by r')
                        dist = [s.stats.sac.dist for s in sets[si]]
                        norms = [d/np.max(dist) for d in dist]
                        norms = [d**(1) for d in norms]
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                        sets[si][i].data = sets[si][i].data / norms[i]
                        colmax = np.max([abs(d.data.max()) for d in sets[si]])
                        sets[si][i].data = sets[si][i].data / colmax
                # -----------
        correctset = sets[-1]
        [ax.plot(tr.times(),tr.data + ysep,linewidth=linewidth,color='k') for ysep,tr in enumerate(correctset)]
        if len(sets)>1:
            for si,s in enumerate(sets):
                [ax.fill_between(tr.times(),tr.data*scales[si] + ysep,tr.data*0 + ysep, where=np.abs(tr.data)>=0, facecolor=facecolor[si]) for ysep,tr in enumerate(s)]
        [ax.plot(tr.times(),tr.data*0 + ysep,linewidth=0.4,color='k') for ysep,tr in enumerate(correctset)]
        if evdepth is not None:
            arrivals = [ObsQA.io.get_arrivals(sta_llaz=(sta.stats.sac.stla,sta.stats.sac.stlo,sta.stats.sac.stel),ev_llaz=(sta.stats.sac.evla,sta.stats.sac.evlo,evdepth),phases=phases) for sta in correctset]
            ardict = dict()
            corephase_dy = 0.2
            direcphase_dy = 0.03
            [[ardict.update({ph[0]:[]}) for ph in a] for a in arrivals]
            [[ardict[ph[0]].extend([ph[1]]) for ph in a] for a in arrivals]
            [ardict.update({k:np.max(ardict[k])}) for k in list(ardict.keys())]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='b',linewidth=0.2) for ph in a] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='P'] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='S'] for ysep,a in enumerate(arrivals)]
            [ax.text(ardict[k], len(correctset) - corephase_dy, k, color='b',horizontalalignment='center') for k in list(ardict.keys()) if (k!='P') and (k!='S')]
            [ax.text(ardict[k], len(correctset) + direcphase_dy, k, color='r',horizontalalignment='center') for k in list(ardict.keys()) if ((k=='P') or (k=='S'))]
        yl = (-1,len(correctset))
        ax.set_ylim(yl)
        ax.set_xlim(correctset[0].times()[0],correctset[0].times()[-1])
        ax.set_yticks([i for i in range(len(correctset))])
        labels = [str(int(ev.stats.sac.dist)) +'km' + ' [' + ev.stats.network + '] ' + ev.stats.station  + '\n depth:' + str(int(abs(ev.stats.sac.stel*1000))) + 'm' for ev in correctset]
        ax.set_yticklabels(labels)
        ax.set_xlabel('Time(s)')
        if fig is not None:
            fig.suptitle(title,fontweight='bold',fontsize=15)
        return ax


def simple_ax_sta_metrics(ax,sr,sta,method,event_name=None,flim=[1/500,1],octave_av=False,y_faxis=True,metric='coh',tf='ZP-21'):
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = {'Z1':'#0c51a6','ZP':'#2a7e93','ZZ':'#7370cb','Z2':'#4f86c5'}
    note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    stanm = sr.StaName
    tf=tf.replace('-','_').lower()
    # stastr = ' | '.join([stanm+tf,sta.Experiment,'Depth: '+str(int(abs(sta.StaDepth)))+'m, Notch: '+str(int(1/fnotch(1000*abs(sta.StaDepth))))+'s',note])
    alpha=0.4
    outlierprops={'color':'r','s':10,'alpha':0.09}
    inlierprops={'color':'dodgerblue','s':13,'alpha':alpha} #'#7370cb' royalblue dodgerblue
    whiskerprops={'linewidth':0.25,'color':'k','alpha':alpha}
    whiskerwidth=0.002;margin=1.02
    midlineprops={'linewidth':0.5,'color':'k','alpha':0.8}
    event_M_lnie={'linewidth':0.75,'color':'r','alpha':0.8}
    midscatterprops={'s':1,'color':'k'} #args.csd_pairs[pair]
    notchprops = {'alpha':0.4,'linewidth':0.5,'color':'k','linestyle':'-.'}
    labelprops = {'fontsize':5} #'fontweight':'bold'
    methodlabelprops = {'fontweight':'bold','fontsize':5,'verticalalignment':'bottom','horizontalalignment':'right'}
    if metric=='coh':metric_label = 'Coherence'
    if metric=='adm':metric_label = 'Admittance'
    stastr=f'"Figure: ZZ {metric_label} for {stanm} | {sr.Experiment} | {str(int(abs(sr.StaDepth)))}m"'

    keys=['TF_ZP-21']
    fn = 1/fnotch(sr.StaDepth)
    methods=[method]
    metric=metric.lower()
    COH = sta.Data[0].Coherence()
    f = COH.f
    for bi,b in enumerate(methods):
        method = methods[bi]
        pair = 'ZZ'
        # ax.set_title(f'i.',**labelprops)
        # ax.set_title(f'i. {pair} {metric}',**labelprops)
        ind=(f>0)&(f<=1)
        f=f[ind]
        if method=='NoiseCut':
            M=COH.HPS.zz.coh[:,ind]
            if event_name:
                event_M=sr.Coherence.HPS_Z[:,ind]
                if octave_av:event_M = octavg(event_M,f)[1]
            if octave_av:M = octavg(M,f)[1]
        else:
            M=COH.ATaCR.zp_21.coh[:,ind]
            if event_name:
                event_M=sr.Coherence.TF[:,ind]
                if octave_av:event_M = octavg(event_M,f)[1]
            if octave_av:M = octavg(M,f)[1]
        if octave_av:
            M = octavg(M,f)[1]
            event_M = octavg(event_M,f)[1]
            f = octavg(M[0],f)[0]

        upper,lower,midline,outliers,inliers=cohstats(M,margin=margin) 
        #midline is currently defined as the MEAN of the statistical inliers within the Q1 to Q3 margin (50% of all observations for a given frequency).
        y_mid = midline
        # y_mid = M.mean(axis=0)

        if metric_label.lower()=='phase':ax.set_ylim(0,180)
        if metric_label.lower()=='coherence':ax.set_ylim(0,1.02)
        if y_faxis:
            # [ax.scatter(x,yy,label=':'.join([pair]),s=5,facecolor=args.csd_pairs[pair],alpha=0.01) for yy in M]
            # ax.set_xlim(f[0],1)
            # y_mid = gaussian_filter(midline, sigma=1)
            # y_mid = midline
            ax.scatter(y_mid,f,label=':'.join([pair]),**midscatterprops)
            ax.hlines(f,xmin=y_mid,xmax=upper,**{'linewidth':0.3,'color':'k','alpha':alpha})
            ax.hlines(f,xmin=lower,xmax=y_mid,**{'linewidth':0.3,'color':'k','alpha':alpha})
            ax.vlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.vlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.plot(y_mid,f,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]

            ax.plot(event_M.reshape(-1),f,**event_M_lnie)

            ax.set_ylabel('Frequency, Hz',**labelprops)
            ax.set_xlabel(metric_label,**labelprops)
            ax.axhline(1/fn,**notchprops)
            # ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            # if bi==0:ax.set_yticklabels([]);ax.set_ylabel('')
            ax.set_ylim(flim[0],flim[1])
            if metric=='coh':ax.set_xlim(0,1)
            _ = ax.set_yscale('log')
        else:
            ax.set_xlim(f[0],1)
            y_mid = M.mean(axis=0)
            # y_mid = gaussian_filter(midline, sigma=1)
            # y_mid = midline
            ax.scatter(f,y_mid,label=':'.join([pair]),**midscatterprops)
            ax.vlines(f,ymin=y_mid,ymax=upper,**{'linewidth':0.3,'color':'k'})
            ax.vlines(f,ymin=lower,ymax=y_mid,**{'linewidth':0.3,'color':'k'})
            ax.hlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.hlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.plot(f,y_mid,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]
            ax.set_xlabel('Frequency, Hz',**labelprops)
            ax.set_ylabel(metric_label,**labelprops)
            ax.axvline(1/fn,**notchprops)

            _ = [ax.set_xscale('log') for ax in axes]
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            # ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
    # ax.text(0.001+1/fn,0,str(int(np.round(fn)))+'s',verticalalignment='bottom',horizontalalignment='left',fontweight='bold',fontsize=11)
    return ax
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # ++++++++++++++++++++++++++ CONSTRUCTOR AREA ++++++++++++++++++++++++++++++
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def ax_sta_metrics(ax,report,sta,method,event_name=None,flim=[1/500,1],octave_av=False,y_faxis=True,metric='coh',tf='ZP-21'):
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = {'Z1':'#0c51a6','ZP':'#2a7e93','ZZ':'#7370cb','Z2':'#4f86c5'}
    note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    stanm = sta.StaName
    tf=tf.replace('-','_').lower()
    # stastr = ' | '.join([stanm+tf,sta.Experiment,'Depth: '+str(int(abs(sta.StaDepth)))+'m, Notch: '+str(int(1/fnotch(1000*abs(sta.StaDepth))))+'s',note])
    alpha=0.4
    outlierprops={'color':'r','s':10,'alpha':0.09}
    inlierprops={'color':'dodgerblue','s':13,'alpha':alpha} #'#7370cb' royalblue dodgerblue
    whiskerprops={'linewidth':0.25,'color':'k','alpha':alpha}
    whiskerwidth=0.002;margin=1.02
    midlineprops={'linewidth':0.5,'color':'k','alpha':0.8}
    event_M_lnie={'linewidth':0.75,'color':'r','alpha':0.8}
    midscatterprops={'s':1,'color':'k'} #args.csd_pairs[pair]
    notchprops = {'alpha':0.4,'linewidth':0.5,'color':'k','linestyle':'-.'}
    labelprops = {'fontsize':5} #'fontweight':'bold'
    methodlabelprops = {'fontweight':'bold','fontsize':5,'verticalalignment':'bottom','horizontalalignment':'right'}
    if metric=='coh':metric_label = 'Coherence'
    if metric=='adm':metric_label = 'Admittance'
    stastr=f'"Figure: ZZ {metric_label} for {stanm} | {sta.Experiment} | {str(int(abs(sta.StaDepth)))}m"'

    keys=['TF_ZP-21']
    fn = 1/fnotch(sta.StaDepth)
    methods=[method]
    metric=metric.lower()
    for bi,b in enumerate(methods):
        method = methods[bi]
        pair = 'ZZ'
        f=report.f
        # ax.set_title(f'i.',**labelprops)
        # ax.set_title(f'i. {pair} {metric}',**labelprops)
        ind=(f>0)&(f<=1)
        f=f[ind]
        if method.lower()=='noisecut':
            M=report[method][stanm][metric][:,ind]
            events=report[method][stanm].events[~np.any(np.isnan(M),axis=1)]
        else:
            M=report[method][tf][stanm][metric][:,ind]
            events=report[method][tf][stanm].events[~np.any(np.isnan(M),axis=1)]
             
        M = M[~np.any(np.isnan(M),axis=1)]
        if event_name:event_M=M[np.where(events==event_name)[0],:].reshape(-1)
        if octave_av:
            avged=np.atleast_2d(np.squeeze(np.array([octavg(c,f)[1] for c in M])))
            event_M = octavg(event_M,f)[1]
            f = octavg(M[0],f)[0]
            M=avged

        upper,lower,midline,outliers,inliers=cohstats(M,margin=margin) 
        #midline is currently defined as the MEAN of the statistical inliers within the Q1 to Q3 margin (50% of all observations for a given frequency).
        y_mid = midline
        # y_mid = M.mean(axis=0)

        if metric_label.lower()=='phase':ax.set_ylim(0,180)
        if metric_label.lower()=='coherence':ax.set_ylim(0,1.02)
        if y_faxis:
            # [ax.scatter(x,yy,label=':'.join([pair]),s=5,facecolor=args.csd_pairs[pair],alpha=0.01) for yy in M]
            # ax.set_xlim(f[0],1)
            # y_mid = gaussian_filter(midline, sigma=1)
            # y_mid = midline
            ax.scatter(y_mid,f,label=':'.join([pair]),**midscatterprops)
            ax.hlines(f,xmin=y_mid,xmax=upper,**{'linewidth':0.3,'color':'k','alpha':alpha})
            ax.hlines(f,xmin=lower,xmax=y_mid,**{'linewidth':0.3,'color':'k','alpha':alpha})
            ax.vlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.vlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.plot(y_mid,f,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]

            ax.plot(event_M,f,**event_M_lnie)

            ax.set_ylabel('Frequency, Hz',**labelprops)
            ax.set_xlabel(metric_label,**labelprops)
            ax.axhline(1/fn,**notchprops)
            # ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            # if bi==0:ax.set_yticklabels([]);ax.set_ylabel('')
            ax.set_ylim(flim[0],flim[1])
            if metric=='coh':ax.set_xlim(0,1)
            _ = ax.set_yscale('log')
        else:
            ax.set_xlim(f[0],1)
            y_mid = M.mean(axis=0)
            # y_mid = gaussian_filter(midline, sigma=1)
            # y_mid = midline
            ax.scatter(f,y_mid,label=':'.join([pair]),**midscatterprops)
            ax.vlines(f,ymin=y_mid,ymax=upper,**{'linewidth':0.3,'color':'k'})
            ax.vlines(f,ymin=lower,ymax=y_mid,**{'linewidth':0.3,'color':'k'})
            ax.hlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.hlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.plot(f,y_mid,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]
            ax.set_xlabel('Frequency, Hz',**labelprops)
            ax.set_ylabel(metric_label,**labelprops)
            ax.axvline(1/fn,**notchprops)

            _ = [ax.set_xscale('log') for ax in axes]
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            # ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
    # ax.text(0.001+1/fn,0,str(int(np.round(fn)))+'s',verticalalignment='bottom',horizontalalignment='left',fontweight='bold',fontsize=11)
    return ax
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # ++++++++++++++++++++++++++ CONSTRUCTOR AREA ++++++++++++++++++++++++++++++
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++