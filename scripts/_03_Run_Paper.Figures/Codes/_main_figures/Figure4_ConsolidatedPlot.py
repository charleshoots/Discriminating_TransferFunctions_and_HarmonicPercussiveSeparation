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

from local_tools.quick_class import *

# cat value
cat = catalog.copy()
# octavg value
octavg = lt.math.octave_average

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
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
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
cohnames2snrnames=c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_1':'HPS.1','HPS_2':'HPS.2'}
# c2s value
c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_1':'HPS.1','HPS_2':'HPS.2'}
# c2s r value
c2s_r={c2s[k]:k for k in c2s.keys()}
# make bands
def make_bands(N, width=None,line=np.linspace, lo=1.0, hi=100.0):
    if width is None:width=(hi-lo)/N
    if width<=0 or width>(hi-lo): raise ValueError("width must be in (0, hi-lo]")
    # s value
    s=line(lo, hi-width, N)
    return np.c_[s, s+width]


# icat value
icat=cat.sr.copy()
# usnr value
usnr=unpack_metrics(icat)

# function argsort luminence
def argsort_luminence(cc,cmap):
    # cc value
    cc=np.array(cc).ravel()
    luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098
    color=[cmap(c) for c in enumerate(cc)]
    zorder=np.argsort(np.array([luminance(c) for c in color]))
    return zorder
from matplotlib.colors import LinearSegmentedColormap

dirs= io.dir_libraries()
CSV=dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx'
df=pd.read_excel(CSV)
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in icat.StaName])
cat.sr['Tilt']=theta_deg


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



# --- Plot
avged=False
avged=True

mtrsets='R.vs.SNR'
# mtrsets='ST.vs.LT'

# -------Options--------
file='ConsolidatedFigure.png'
plotfolder=dirs.Ch1/'_main_figures'/'Figure4_ConsolidatedPlot';plotfolder.mkdir(parents=True,exist_ok=True)
phase_letters=['P','S','R']  # order matches third dim of snr
msize=25; scat_alpha=.6
colors={'HPS_1':'olive','HPS_2':'goldenrod','TF_Z':'orangered','HPS_Z':'skyblue',}
mthds=list(colors.keys())
# mthds=np.flip(mthds)
prefbands={'P':(1.0,10.0),'S':(10.0,30.0),'R':(30.0,100.0)}

# -------Plot code:--------
# --------------------------------------------------------------------------------
fig,axes=figs(2,2,f=(5.5,4.5),x=False,y=False)
sectors=[['all_coh','1_10'],['10_30','30_100']]

prefph = {'1_10':'(P/Pdiff)','10_30':'(S/Sdiff)','30_100':'(Rayleigh)'}
prefph = {'1_10':'','10_30':'','30_100':''}


stnm=np.array(list(icat.StaName))
fn=fnotch(np.array(list(icat.StaDepth)))
mags=np.array(list(icat.Magnitude))
quartile=lambda A,p:np.nanpercentile(A,p,0)
yttl = lambda c: fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
loc=lambda s,mg=(6,8):(stnm==s)&(mags>=mg[0])&(mags<=mg[1])
savefile=file
def getkw(mg,ml,msize,colors,m,scat_alpha):
    colors=colors.copy();colors.update({'blank':'w'})
    marker=ml
    s=msize if max(mg)==8.0 else msize/3
    marker=marker
    c=[colors[m]]
    alpha=scat_alpha if max(mg)==8.0 else 1.0
    ec='grey' if max(mg)==8.0 else 'k'
    lw=0.15
    zorder=-100 if max(mg)==8.0 else 100
    kw=dict(s=s,marker=marker,c=c,alpha=alpha,ec=ec,lw=lw,zorder=zorder)
    return kw
# argsort_luminence(np.flip(np.arange(1,5,1)/len(mthds)),LinearSegmentedColormap.from_list('dark_greys', list(colors.values())).resampled(len(mthds)))

stemlw=1.1;stemalpha=0.5;  # ~±5% in log-x
cap_dec = 0.007 # ← this is the “~±5% in log-x”
s = 10**cap_dec # ≈ 1.122 (nice & small)

snravg = np.nanmax;snravgname = 'snrmax'
# snravg = np.nanmean;snravgname = 'snrmean'

xhold = {'TF_Z':[],'HPS_Z':[],'HPS_1':[],'HPS_2':[]}
for ri,(raxes,rsect) in enumerate(zip(axes,sectors)):
    for axi,(ax,sect) in enumerate(zip(raxes,rsect)):
        if sect=='all_coh':
            f=1/usnr.coh.bands; msk=(f>=1/100)&(f<=1); f=f[msk]
            coh=np.array([usnr.coh.__dict__[m].D for m in mthds])
            coh=np.array([[quartile(A,p) for p in (25,50,75)] for A in coh])[:,:,msk]  # (method,3,F)
            P = 1/f
            xl, xr = P/s, P*s  # left/right cap x-ends around x0 = P
            # (optional: enforce order if you ever vary s by array)
            # xl, xr = np.minimum(xl, xr), np.maximum(xl, xr)
            hdls=[]
            for i,m in enumerate(mthds):
                zorder=100 if m in ['TF_Z','HPS_Z'] else -100
                # label=m.replace('_',' ')
                label=yttl(m.replace('_',' '))
                q1,med,q3 = coh[i,0],coh[i,1],coh[i,2]
                ax.plot(1/f, med, lw=1.3, color='k', label=label,alpha=stemalpha,zorder=zorder)
                h=ax.plot(1/f, med, lw=stemlw, color=colors[m], label=label,alpha=stemalpha,zorder=zorder)
                ax.vlines(1/f, q1, med, lw=stemlw/2, color=colors[m],alpha=stemalpha,zorder=zorder)
                ax.vlines(1/f, med, q3, lw=stemlw/2, color=colors[m],alpha=stemalpha,zorder=zorder)
                hdls.append(h[0])
                for ycap in (q1, q3):
                    ax.hlines(ycap, xl, xr, lw=stemlw, color=colors[m])  # boxplot-style caps
            ax.set_ylim(top=1.0412867859005928)
            # ax.set_xscale('log');
            ax.set_xlim(100,1)
            ax.set_xlabel('Period, s'); ax.set_ylabel('Coherence')
            leg=ax.legend(handles=hdls,frameon=False, ncols=1, fontsize=5,loc='lower right',markerscale=10)
            for h in leg.legend_handles:h.set_linewidth(4)
            # ax.grid(True,which='both',alpha=0.15,color='k')
        else:
            band=tuple(map(float,sect.split('_')))
            magwins=[[6.0,7.0],[7.0,8.0]]
            magwins=[[6.0,8.0]]
            for i,m in enumerate(mthds):
                for mgi,mg in enumerate(magwins):
                    idx=(mags>=mg[0])&(mags<=mg[1])
                    # mlwins,fnwins=['<','>'],['IG','MS']
                    mlwins,fnwins=['<','>'],['IG',None]
                    for ml,fn in zip(mlwins,fnwins):
                        if mtrsets=='R.vs.SNR':
                            x=np.array([usnr.snr.__dict__[m].R().Average(band,fn=fn) for m in mthds]) # (method, SR, phase)
                            x = snravg(x,axis=2)
                            y=np.array([usnr.coh.__dict__[m].Average(band,fn=fn) for m in mthds])     # (method, SR)
                        else:
                            x=np.array([usnr.LT.__dict__[m].R().Average(band,fn=fn) for m in mthds]) # (method, SR, phase)
                            # x=abs(x)
                            y=np.array([usnr.ST.__dict__[m].R().Average(band,fn=fn) for m in mthds])     # (method, SR)
                        if avged:
                            savefile=f'{file.split('.png')[0]}.Averaged';msize=25
                            
                            if mtrsets=='R.vs.SNR':
                                x=np.moveaxis(np.array([np.nanmean(x[:,loc(s,mg)],axis=1) for s in np.unique(stnm)]),0,1) #reduced to station averages
                                y=np.array([np.nanmean(y[:,loc(s,mg)],axis=1) for s in np.unique(stnm)]).T
                            else:
                                x=np.array([np.nanmean(x[:,loc(s,mg)],axis=1) for s in np.unique(stnm)]).T
                                y=np.moveaxis(np.array([np.nanmean(y[:,loc(s,mg),:],axis=1) for s in np.unique(stnm)]),0,1) #reduced to station averages
                            idx=np.bool_(np.ones(y.shape[1]))

                        for p,letter in enumerate(phase_letters):
                            if not band==prefbands[letter]:continue
                            if mtrsets=='R.vs.SNR':xx=x[i,idx,p] if x.ndim==3 else x[i,idx];yy=y[i,idx]
                            else:yy=y[i,idx,p];xx=x[i,idx];xx=10**xx;yy=10**yy
                            ax.scatter(xx,yy,**getkw(mg,ml,msize,colors,m,scat_alpha))
                        xhold[m].append(x)
            if mtrsets=='R.vs.SNR':ax.set_xlabel(rf'$\eta$ {prefph[sect]}');ax.set_ylabel(r'$\gamma$')
            else:ax.set_xlabel('Noise ratio, log10');ax.set_ylabel('Signal ratio, log10'); 
            ax.text(0.98,0.96,f'{sect.replace('_','-')}s',ha='right',va='center',transform=ax.transAxes,fontsize=8)
            # ax.grid(True,which='both',alpha=0.15,color='k')
            if avged:
                # if mtrsets=='R.vs.SNR':ax.set_xscale('log')
                if mtrsets=='R.vs.SNR':ax.set_ylim(top=1.0312867859005928 if axi==1 else 1.0412867859005928)
            else:pass
                # if mtrsets=='R.vs.SNR':ax.set_xscale('log')
                # else:ax.set_xscale('symlog');ax.set_yscale('symlog')
            # if not (mtrsets=='R.vs.SNR'):ax.set_yscale('log');ax.set_xscale('log')
            # if mtrsets=='R.vs.SNR':
            #     if (ri==1)&(avged):
            #         ax.set_ylim(bottom=0.1)
            #     else:
            #         ax.set_ylim(bottom=-.04,top=1.0412867859005928)
        if (axi==1)&(ri==0):
            sch=[]
            ml,fn=mlwins[0],fnwins[0];m='blank';mg=magwins[-1]
            kw=getkw(mg,ml,msize,colors,m,scat_alpha);kw.update({'alpha':1.0,'c':['k']})
            sc=ax.scatter(np.nan,np.nan,label=r'$f$  $\leq$  $f_{n}$', **kw)
            sch.append(sc)
            ml,fn=mlwins[1],fnwins[1];m='blank';mg=magwins[-1]
            kw=getkw(mg,ml,msize,colors,m,scat_alpha);kw.update({'alpha':1.0,'c':['k']})
            sc=ax.scatter(np.nan,np.nan,label=r'$f$  >  $f_{n}$' if not np.isin(None,fnwins) else r'$f$ $\leq$ $f_{n}$ $\leq$ $f$', **kw)
            sch.append(sc)
            ax.legend(frameon=False,handles=sch,loc='lower left',fontsize=5)



savefile = f'{savefile}.{snravgname}'

if not mtrsets=='R.vs.SNR':
    savefile = f'{mtrsets}.{savefile}'

_=save_tight(plotfolder/f'{savefile}.png',fig,dpi=2000)

print('Done.')