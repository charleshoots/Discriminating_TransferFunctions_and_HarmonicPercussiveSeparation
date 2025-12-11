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
from scipy.stats import iqr
from local_tools.quick_class import *
from local_tools.math import spectra
from obspy.geodetics import kilometers2degrees
# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
import statsmodels.api as sm
import local_tools.dataspace as ds

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
# rmse value
rmse=lambda y:( (  ( abs(y)-abs(y).mean() )**2  ).mean())**.5 
# rms value
rms=lambda y:np.mean(y**2)**0.5
# s value
s=cat.r.iloc[0];faxis=(s.Noise.f>(1/100) )& (s.Noise.f<=1)
# f value
f=s.Noise.f[faxis]
cat.r['NoiseAverage']=[{f'{b[0]}_{b[1]}':-rms(s.Noise.Z[faxis][(f<=(1/b[0]))&(f>=(1/b[1]))]) for b in [[1,10],[10,30],[30,100]]} for s in cat.r.iloc]
cat.sr['NoiseAverage']=[cat.r.loc[sr.StaName].NoiseAverage[0] for sr in cat.sr.iloc]

# figs value
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',w=None,layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout,width_ratios=np.ones(c) if w is None else w)
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
# baz value
baz=lambda s:obspy.geodetics.base.gps2dist_azimuth(s.Latitude,s.Longitude,s.LaLo[0],s.LaLo[1])[1]
# bootstrap value
bootstrap = lambda y,nruns=10000,nchoose=100,aggregate=np.mean: np.mean([aggregate(np.random.choice(y[~np.isnan(y)],nchoose)) for _ in range(nruns)])
# centers value
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
# dbin value
dbin = lambda x:np.array([x[:-1],x[1:]]).T
# dirs value
dirs=io.dir_libraries()
# CSV value
CSV=dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx'
df=pd.read_excel(CSV)
# theta deg value
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in cat.sr.StaName])
cat.sr['TiltDirection']=theta_deg
# oriencohere value
oriencohere=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].oriencohere for stnm in cat.sr.StaName])
cat.sr['TiltCoherence']=oriencohere
# make bands
def make_bands(N, width=None,line=np.linspace, lo=1.0, hi=100.0):
    if width is None:width=(hi-lo)/N
    if width<=0 or width>(hi-lo): raise ValueError("width must be in (0, hi-lo]")
    # s value
    s=line(lo, hi-width, N)
    return np.c_[s, s+width]



from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl


#An HPS example spectrogram
mpl.rcParams['axes.linewidth']= 0.5     # subplot borders
mpl.rcParams['axes.edgecolor'] = 'k'
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['ytick.major.width'] = 0.5


# plot spec
def plot_spec(fig, ax, S, f, t, cmap='magma', cax=None,fontsize=7):
    # units value
    units = ['seconds after origin','minutes after origin','hours after origin']
    # t value
    t = t/60; units.pop(0)
    # pcm value
    pcm = ax.pcolormesh(t, f, librosa.power_to_db(np.abs(S)), cmap=cmap, shading='auto')

    # Colorbar that does NOT change the layout:
    if cax is None:
        # cax value
        cax = inset_axes(ax, width="2%", height="100%", loc="lower left",
        # bbox to anchor value
        bbox_to_anchor=(1.02, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=0)
    # cbar value
    cbar = fig.colorbar(pcm, cax=cax,label='dB')
    # cbar.# ax.tick_params(size=fontsize,labelsize=fontsize)
    ax.set_yscale('log')
    ax.set_ylim(0, 1)
    ax.set_xlabel(units[0])
    return cbar,pcm



# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure2_Methods.Example' #Folder maps are saved to.
plotfolder.mkdir(parents=True,exist_ok=True)

mpl.rcParams.update({
'lines.antialiased': False,
# 'lines.solid_joinstyle': 'miter',
'lines.solid_capstyle': 'butt'
})

# sta='7D.FS42D'
sta='7D.FN12C'
# ev='2015.125.01.44'
# ev='2015.208.21.41'
# ev='2015.150.11.23'
# ev='2015.140.22.48'
# ev='2015.088.23.48'
ev='2013.298.17.10'
# files value
files = [[sta,ev]]




tfzcolor='orangered'
hpszcolor='skyblue'

# stas=['7D.FS42D']
# 7D.FS42D.2015.088.23.48.pkl
status=lambda:print(f'{si+1}/{len(stas)} {ei+1}/{len(evs)}')

# fontsize=4

# figsize=(4,3)
figsize=(6,3.5)

for fi in files:
    stas=[fi[0]]
    evs=[fi[1]]
    for si,sta in enumerate(stas):
        # evs=list(cat.sr.aloc[sta].Name)
        # evs=['2015.088.23.48']
        for ei,ev in enumerate(evs):
            status()
            fold=dirs.Events_HPS/'corrected'/sta/'Spectrograms'
            file=f'{sta}.{ev}.pkl'
            d=load_pickle(fold/file)
            d=pd.DataFrame(d).Spectrograms

            data=cat.sr.aloc[ev].aloc[sta].iloc[0]


            # fig,axes = figs(2,2,f=(6,4))
            fig=plt.figure(figsize=figsize, constrained_layout=False)
            gs = gridspec.GridSpec(2, 2, width_ratios=[.1,1], height_ratios=[.4,1], wspace=.1, hspace=.05)
            blank_ax=fig.add_subplot(gs[0, 0])
            tr_ax=fig.add_subplot(gs[0, 1])
            coh_ax=fig.add_subplot(gs[1, 0])
            spec_ax=fig.add_subplot(gs[1,1])


            lw=0.4

            blank_ax.set_visible(False)
            tr_ax.set_yticks([])

            xmax=100

            # -------Coherences
            ax=coh_ax
            tff=data.Data.TF().f;faxis=(tff>0)&(tff<=1);tff=tff[faxis]
            tf=data.Data.TF().transfunc['ZP-21']['TF_ZP-21'][faxis]
            # tf=abs(tf);
            tf=tf/tf.max()
            tf=tf+abs(tf.min())
            tf=tf/tf.max()
            coh_f=data.Data.Coherence().f
            faxis=(coh_f>(1/200))&(coh_f<=1)
            coh_f=coh_f[faxis]
            tfz_coh=data.Coherence['TF_Z'].reshape(-1)[faxis]
            hpsz_coh=data.Coherence['HPS_Z'].reshape(-1)[faxis]
            fbands=make_bands(1000,np.max(np.diff(coh_f)),lo=min(coh_f),hi=1)
            tfz_coh=np.array([tfz_coh[(coh_f>=fb[0])&(coh_f<=fb[1])].mean() for fb in fbands])
            hpsz_coh=np.array([hpsz_coh[(coh_f>=fb[0])&(coh_f<=fb[1])].mean() for fb in fbands])
            coh_f=fbands.mean(axis=1)
            ax.plot(tfz_coh,coh_f,color=tfzcolor,lw=lw*2)
            # ax.plot(tf,tff,color='k',ls=':')
            ax.plot(hpsz_coh,coh_f,color=hpszcolor,lw=lw*2)
            ax.set_yscale('log')
            ax.set_ylim(1/100,1)
            ax.set_xlabel(r'$\gamma$')
            ax.axhline(fnotch(data.StaDepth),ls=':',lw=0.4,c='k')
            ax.set_ylabel('frequency, Hz')
            # ax.tick_params(size=fontsize,labelsize=fontsize)

            noise=data.Data.Noise.Averaged()
            cf,czh=noise.coherence('ZH',return_f=True)
            faxis=(cf>(1/200))&(cf<=1);cf=cf[faxis];czh=czh[faxis]
            czp=noise.coherence('ZP')[faxis]
            zpcolor='black'
            zhcolor='gray'
            rs=20
            ax.plot(czp[::rs],cf[::rs],color=zpcolor,lw=lw,ls='-.',zorder=-1e4)
            ax.plot(czh[::rs],cf[::rs],color=zhcolor,lw=lw,ls='-.',zorder=-1e4)



            ax=spec_ax
            # -------Spectrogram
            faxis=(d[0]['f']>0)&(d[0]['f']<=1)
            S_full,S_background,S_hps,f,t=d[0]['Full'][faxis,:],d[0]['Background'][faxis,:],d[0]['HPS'][faxis,:],d[0]['f'][faxis],d[0]['t']
            t=t-t.min()
            S=S_full
            cbar,pcm=plot_spec(fig,ax,S,f,t,cmap='magma')
            # pcm.set_rasterized(True)
            # for c in pcm.collections:
            pcm.set_rasterized(True)
            pcm.set_edgecolor("face")
            pcm.set_linewidth(0)

            cbar.set_label('dB')
            ax.set_ylim(1/100,1)
            ax.set_yticks([])
            ax.set_xticks([30,60,90,120])
            ax.set_xlim(0,xmax)
            # ax.tick_params(size=fontsize,labelsize=fontsize)
            ax.set_xlabel('minutes after origin')

            # -------Traces
            ax=tr_ax
            orig_c='black'
            band=(1/fnotch(data.StaDepth),60)
            band=(8,60)
            traces=data.Traces()
            traces.taper(0.01).filter('bandpass',freqmin=1/max(band),freqmax=1/min(band))
            t=traces[0].times()/60
            orig_tr=traces.select(location='Original')[0].data
            maxamp=np.max(abs(orig_tr))
            hps_tr=traces.select(location='NoiseCut')[0].data
            tf_tr=traces.select(location='ATaCR')[0].data
            orig_tr=orig_tr/maxamp;hps_tr=hps_tr/maxamp;tf_tr=tf_tr/maxamp
            offset=0
            ax.plot(t,orig_tr+offset,lw=lw,color=orig_c,alpha=1);ax.plot(t,hps_tr+offset,lw=lw,color=hpszcolor,alpha=1)
            offset=2
            ax.plot(t,orig_tr+offset,lw=lw,color=orig_c,alpha=1);ax.plot(t,tf_tr+offset,lw=lw,color=tfzcolor,alpha=1)
            ax.set_xlim(0,xmax)
            ax.set_xticks([])
            # ax.tick_params(size=fontsize,labelsize=fontsize)
            phases=event_stream_arrivals(traces[0],data.Event)
            if (('S' in phases.keys())&('P' in phases.keys())):
                PS=[phases['P'][0]/60,phases['S'][0]/60]
                ax.axvline(PS[0],lw=.4,zorder=-1e3,ls=':',c='k');ax.text(PS[0],0.99*max(ax.get_ylim()),'P',va='top',ha='right')
                ax.axvline(PS[1],lw=.4,zorder=-1e3,ls=':',c='k');ax.text(PS[1],0.99*max(ax.get_ylim()),'S',va='top',ha='right')
                spec_ax.axvline(PS[0],lw=.4,zorder=1e3,ls=':',c='k')
                spec_ax.axvline(PS[1],lw=.4,zorder=1e3,ls=':',c='k')
            file=f'{sta}.{ev}.Methods.Example.{save_format}'
            fig.savefig(str(plotfolder/('S03.Figure02_'+file)),pad_inches=0.05,dpi=3000,format=save_format)
            plt.close('all')