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
# function custom cmap
def custom_cmap(ind=0,nbins=5):
    if ind==0:cmap = cm.cmaps['glasgow'].reversed().resampled(nbins)
    if ind==1:cmap = cm.cmaps['batlow'].reversed().resampled(nbins)
    return cmap
# figs value
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',w=None,layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout,width_ratios=np.ones(c) if w is None else w)
from obspy.signal.trigger import classic_sta_lta,carl_sta_trig,recursive_sta_lta
# stalta methods value
stalta_methods={'classic':classic_sta_lta,'carl':carl_sta_trig,'recurssive':recursive_sta_lta}
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
# suspect stations value
suspect_stations=np.array(['ZA.B02','YL.C09W','7D.G25B','7D.FS08D','7D.G17B','YL.A14W'])
# baz value
baz=lambda s:obspy.geodetics.base.gps2dist_azimuth(s.Latitude,s.Longitude,s.LaLo[0],s.LaLo[1])[1]
# bootstrap value
bootstrap = lambda y,nruns=10000,nchoose=100,aggregate=np.mean: np.mean([aggregate(np.random.choice(y[~np.isnan(y)],nchoose)) for _ in range(nruns)])
# centers value
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
# dbin value
dbin = lambda x:np.array([x[:-1],x[1:]]).T
# phases value
phases=['P','S','Rg'];preferred_pbands={'P':'1_10','S':'10_30','Pdiff':'1_10','Sdiff':'10_30','Rg':'30_100'}
# methods value
methods=['NoiseCut','ATaCR'];mnames={'NoiseCut':'HPS','ATaCR':'TF'}
# mnames r value
mnames_r={mnames[k]:k for k in mnames.keys()};mname_comp={f'HPS_Z':'NoiseCut','TF':'ATaCR','Original':'Original'}
# mnames comp r value
mnames_comp_r={mname_comp[k]:k for k in mname_comp.keys()}
# cohnames2snrnames value
cohnames2snrnames=c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_H':'HPS.H'}
# c2s value
c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_H':'HPS.H'}
# c2s r value
c2s_r={c2s[k]:k for k in c2s.keys()}
# dirs value
dirs= io.dir_libraries()
# OUT CSV value
OUT_CSV=dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx'
df=pd.read_excel(OUT_CSV)
# theta deg value
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in cat.sr.StaName])
cat.sr['TiltDirection']=theta_deg
# oriencohere value
oriencohere=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].oriencohere for stnm in cat.sr.StaName])
cat.sr['TiltCoherence']=oriencohere



# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure3_DirectMetricComparison';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'







# icat value
icat=cat.sr.copy()
# usnr value
usnr=unpack_metrics(icat)


# function cbarlam
def cbarlam(axes,fig,cmap=cm.cmaps['glasgow'],vmin=0,vmax=1,orientation='vertical',fraction=0.025,pad=0.04,frameon=True):
    norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar=fig.colorbar(sm, ax=axes, orientation=orientation, fraction=fraction, pad=pad)
    return cbar
def cscale(c_raw,n=None,smin=1,smax=150,sgamma=1):
    if n is None:n=len(c_raw)
    c=np.nan_to_num(c_raw,nan=np.nanmin(c_raw));a,b=np.nanmin(c),np.nanmax(c)
    if b==a:b=a+1
    cn=((c-a)/(b-a)).clip(0,1)**sgamma
    return np.tile((smin+cn*(smax-smin))[:,None],n)
def band_xy(mtr,band,ph,fn=None,octave=True):
    if mtr=='coh':
        x=usnr.coh.HPS_Z.Average(band,fn=fn,octave=octave); y=usnr.coh.TF_Z.Average(band,fn=fn,octave=octave)
        # g=np.ones(len(x),dtype=bool)
    else:
        # x=usnr.snr.HPS_Z.R().__dict__[ph].Average(band,fn=fn,octave=octave)
        # y=usnr.snr.TF_Z.R().__dict__[ph].Average(band,fn=fn,octave=octave)
        x=usnr.snr.HPS_Z.R().Average(band,fn=fn,octave=octave)
        y=usnr.snr.TF_Z.R().Average(band,fn=fn,octave=octave)
        # g=np.ones(len(x),dtype=bool)
        bar = np.log10(.9) #SNR loss of 50%
        g=(x<0)&(x>bar); x[g]=np.nan;y[g]=np.nan
        g=(y<0)&(y>bar); x[g]=np.nan;y[g]=np.nan
    return x,y


icat=cat.sr.copy()
usnr=unpack_metrics(icat)


inst_color=ColorStandard.instrument
seis_marker=ColorStandard.seismometer_marker
u=np.unique(np.array(list(icat.StaName)))
cols=[inst_color[cat.r.loc[s].iloc[0].Instrument_Design] for s in u]
mks =[seis_marker[cat.r.loc[s].iloc[0].Seismometer]       for s in u]

phases=np.array(['P','S','Rg'])
# pref={tuple(b):ph for b,ph in zip(bands,phases)}
pcol={'P':'red','S':'royalblue','Rg':'violet'}
# pcol={'P':'lightsteelblue','S':'slategrey','Rg':'black'}
# bands=[tuple(b) for b in bands]
bands=np.array([(1,10),(10,30),(30,100)])
phaselabels={'P':'1-10s','S':'10-30s','Rg':'30-100s'}
phaselabels = {k:f'{b[0]}-{b[1]}s' for k,b in zip(phases,bands)}
phases = list(phaselabels.keys())


yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# msize=10
msize=25
msizes=cscale(np.array(bands).mean(axis=1),smin=msize*.4,smax=msize*.7,sgamma=2)[:,0]
msizes = msizes*0+msize
figsize=(5,2)
cmap=mpl.colormaps.get_cmap('viridis')
xcat='Magnitude';dx=.2;u=dbin(np.arange(6,8+dx,dx))
xcat='StaDepth';dx=500;u=dbin(np.arange(0,6000+dx,dx));cmap=cm.cmaps['glasgow']
# xcat='Magnitude';dx=2;u=dbin(np.arange(6,8+dx,dx))
norm=mpl.colors.Normalize(vmin=np.min(u), vmax=np.max(u))
cmap = cmap.resampled(u.shape[0])
octave=True

pref = ['P','S','Rg']

stat=np.nanmean

scat_alpha=.6;stemlw=1.1;stemalpha=0.5;cap_ratio=1.1  # ~±5% in log-x
colors={'HPS_1':'olive','HPS_2':'goldenrod','TF_Z':'orangered','HPS_Z':'skyblue',}

def getkw(mg,ml,ph,msize,scat_alpha,lw=None,ec=None,marker=None,c=None,cmap=None,norm=None):
    colors=pcol.copy();colors.update({'blank':'w'})
    if ec is None:ec=colors[ph]
    if marker is None:marker=ml
    s=msize if max(mg)==8.0 else msize/3
    alpha=scat_alpha if max(mg)==8.0 else 1.0
    if c is None:c='k' if ph=='blank' else 'grey'
    if lw is None:lw=0.35 if ph=='blank' else 0.35
    zorder=-100 if max(mg)==8.0 else 100
    kw=dict(s=s,marker=marker,alpha=alpha,c=c,ec=ec,lw=lw,zorder=zorder,cmap=cmap,norm=norm)
    return kw

snravg = np.nanmax;snravgname = 'snrmax'
# snravg = np.nanmean;snravgname = 'snrmean'

magwins=[[6,8.0]];mg=magwins[0]
cbar = True

mtrlims={'coh':[],'snr':[]}
fig,axes=figs(1,2,f=figsize,x=False,y=False)

mlwins,fnwins=['<','>'],['IG',None]
for ml,fn in zip(mlwins,fnwins):

    data={'coh':{},'snr':{}}
    for bi,b in enumerate(bands):
        # ph=pref[tuple(b)]
        ph = pref[bi]
        x,y=band_xy('coh',b,ph,fn=fn,octave=octave)
        data['coh'][tuple(b)]=x,y
        x,y=band_xy('snr',b,ph,fn=fn,octave=octave)
        x,y = snravg(x,axis=1),snravg(y,axis=1)
        data['snr'][tuple(b)]=x,y
    for ax,mtr in zip(axes,['coh','snr']):
        # base limits / refs
        if mtr=='coh':
            lim=[.07,1.03]; ax.set_xlim(lim); ax.set_ylim(lim)
            # lim=[0,1]
            ax.set_ylabel(yttl('TF Z')); ax.set_xlabel(yttl('HPS Z'))
        else:
            lim=[1e-6,.8]
            lim=[-.05,.8]
            ax.set_xlim(lim[::-1]); ax.set_ylim(lim[::-1]); 
            # ax.set_xscale('symlog'); ax.set_yscale('symlog')
            ax.yaxis.set_label_position('right'); ax.yaxis.tick_right()
            ax.set_ylabel(r'$\eta_{\;TF Z}$'); ax.set_xlabel(r'$\eta_{\;HPS Z}$')
        ax.plot(lim,lim,c='k',zorder=1e3,lw=1,ls=':')
        for bi,b in enumerate(bands[::-1]):
            alpha=1.0 if bi<2 else 1.0
            # ph=pref[tuple(b)]
            ph = pref[bi]
            ec=pcol[ph]; size=msizes[bi] if bi<2 else msizes[bi]*.5
            x,y=data[mtr][tuple(b)]
            g=np.ones(len(x),dtype=bool)
            jj=np.array(list(icat[xcat]))[g]
            xn=np.array([stat(x[(jj>=s[0])&(jj<=s[1])]) for s in u])
            yn=np.array([stat(y[(jj>=s[0])&(jj<=s[1])]) for s in u])
            kw=getkw(mg,ml,ph,msize,scat_alpha,c=u.max(axis=1),cmap=cmap,norm=norm)
            ax.scatter(xn,yn,**kw)
            mtrlims[mtr].append([xn,yn])


        if (mtr=='coh')&(ml==mlwins[-1]):
            hdls=[ax.scatter(np.nan,np.nan,label=phaselabels[ph],**getkw(mg,ml,ph,msize,scat_alpha,c='None',lw=.5)) for ph in phases]

            sch=[]
            ml,fn=mlwins[0],fnwins[0];m='blank';mg=magwins[-1]
            kw=getkw(mg,ml,'blank',msize,scat_alpha);kw.update({'alpha':1.0,'c':['k']})
            sc=ax.scatter(np.nan,np.nan,label=r'$f$  $\leq$  $f_{n}$', **kw)
            sch.append(sc)
            ml,fn=mlwins[1],fnwins[1];m='blank';mg=magwins[-1]
            kw=getkw(mg,ml,'blank',msize,scat_alpha);kw.update({'alpha':1.0,'c':['k']})
            sc=ax.scatter(np.nan,np.nan,label=r'$f$  >  $f_{n}$' if not np.isin(None,fnwins) else r'$f$ $\leq$ $f_{n}$ $\leq$ $f$', **kw)
            sch.append(sc)
            hdls.append(sch[0]);hdls.append(sch[1])

            leg=ax.legend(handles=hdls,facecolor='None',edgecolor='k',frameon=True,fancybox=False,framealpha=1.0, ncols=1,loc='upper left',markerscale=1)
            # leg.set_linewidth(0.2) 
            frame = leg.get_frame()
            frame.set_linewidth(0.2) 
            # for h in leg.legend_handles:h.set_linewidth(.2)
        if (mtr=='snr')&(ml==mlwins[-1]):
            if cbar:
                cbarlam(ax,fig,cmap=cmap,vmin=np.min(u),vmax=np.max(u),orientation='vertical',fraction=0.025,pad=0.04)

mtr='coh';ax=axes[0]
lims=[np.floor(np.nanmin(mtrlims[mtr])*10)/10,1.03]
ax.set_xlim(lims);ax.set_ylim(lims)

mtr='snr';ax=axes[1]
lims=[np.ceil(np.nanmax(mtrlims[mtr])*10)/10,-.02]
ax.set_xlim(lims);ax.set_ylim(lims)


file=f"{xcat}.SNR.and.COH.Comparisons{f'.{fn}' if fn is not None else ''}.FishPlot.{save_format}"
_=save_tight(plotfolder/('S03.Figure03_'+file),fig,dpi=900)
k=0
