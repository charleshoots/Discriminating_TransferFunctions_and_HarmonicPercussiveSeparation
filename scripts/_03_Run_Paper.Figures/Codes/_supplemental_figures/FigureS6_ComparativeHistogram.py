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

import sys;from pathlib import Path
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
usnr=unpack_metrics(icat)

def argsort_luminence(cc,cmap):
    cc=np.array(cc).ravel()
    luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098
    color=[cmap(c) for c in enumerate(cc)]
    zorder=np.argsort(np.array([luminance(c) for c in color]))
    return zorder
from matplotlib.colors import LinearSegmentedColormap


dirs = io.dir_libraries()
OUT_CSV=Path(dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx')
df=pd.read_excel(OUT_CSV)
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in icat.StaName])
cat.sr['Tilt']=theta_deg

from scipy.stats import spearmanr, kendalltau, rankdata
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



mpl.rcParams.update({
"font.size": 4,              # base text size (fallback for everything)
"axes.titlesize": 4,         # axes titles
"axes.labelsize": 4,         # x/y labels (also used by colorbar label)
"xtick.labelsize": 4,        # x tick labels (affects horizontal colorbar ticks)
"ytick.labelsize": 4,        # y tick labels (affects vertical colorbar ticks)
"legend.fontsize": 6,        # legend text
"legend.title_fontsize": 6,  # legend title
"figure.titlesize":4,    # suptitle
'ytick.major.width':0.5,
'xtick.major.width':0.5,
'axes.edgecolor':'k',
'axes.linewidth':0.5})


# ---- Plot options
plotfolder=dirs.Ch1/'_supplemental_figures'/'FigureS6_ComparativeHistogram';plotfolder.mkdir(parents=True,exist_ok=True)

from scipy.stats import norm
mthds=['TF_Z','HPS_Z']
mtrs=['coh','snr']
prefbands={'P':(1.0,10.0),'S':(10.0,30.0),'Rg':(30.0,100.0)}
phases=['Rg','S','P']
cohmax=1.0;cohdx=0.08
snrmax=3.52;snrdx=0.3
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"
cmap=cm.grayC.resampled(100)
colors=[cmap(30),cmap(60),cmap(1000)]

fontsize = 7
filtband=(1,100) #Period band to do the averaging in
fn='IG'#IG = Only include values within a given band (ie (1,100)s) that are sensitive to the IG
# fn= None #None = Ignore all IG sensitivity when averaging over a given band.
# fn='Ms'#MS = Only include values within given band that are outside the sensitivity to the IG

fnsets=[fn]
fnsets=[None,'IG'] #Does both within and regardless of IG sensitivity



lims = {'snr':[],'coh':[]}
xcat = 'Instrument_Design'
cat.sr['All_data']='All data'
for fn in fnsets:
    for xcat in ['All_data','Instrument_Design','Pressure_Gauge','Seismometer']:
        xc = np.unique(cat.sr[xcat].to_numpy())
        for xci in xc:
            fig,axes=figs(2,2,f=(6,3.5),x=False,y=False)
            idx = (cat.sr[xcat]==xci).to_numpy()
            scatter_clr = np.array([ColorStandard.instrument[k] for k in cat.sr.Instrument_Design[idx]])
            scatter_mkr = np.array([ColorStandard.seismometer_marker[k] for k in cat.sr.Seismometer[idx]])
            for ri,raxes in enumerate(axes):
                for axi,(ax,mtr) in enumerate(zip(raxes,mtrs)):
                    if ri==0: #Scatter plots
                        x,y=[],[]
                        for ph,clr in zip(phases,colors):
                            band=prefbands[ph]
                            if mtr=='coh':
                                x.append(usnr.coh.__dict__[mthds[1]].Average(band,fn=fn)[idx])
                                y.append(usnr.coh.__dict__[mthds[0]].Average(band,fn=fn)[idx])
                            else:
                                for m in mthds:
                                    #Each phase's (P/S/Rg) band averaged SNR
                                    band_avg_snr = usnr.snr.__dict__[mthds[1]].R()[ph].Average(band,fn=fn) #Each pair's band-averaged SNR in a single phase
                                    x.append(band_avg_snr[idx]) #Keep a running list for each method examined
                                    band_avg_snr = usnr.snr.__dict__[mthds[0]].R()[ph].Average(band,fn=fn) #Each pair's band-averaged SNR in a single phase
                                    y.append(band_avg_snr[idx]) #Keep a running list for each method examined
                            xx,yy=x[-1],y[-1]
                            for ii in np.unique(scatter_mkr):
                                mkr=ii
                                mkr = 's'
                                jjx=np.array(scatter_mkr==ii)
                                # ax.scatter(xx[jjx],yy[jjx],c=scatter_clr[jjx],marker=ii,ec='k',s=20,lw=0.2)
                                # ax.scatter(xx[jjx],yy[jjx],c=clr,marker=ii,ec='k',s=20,lw=0.2)
                                ax.scatter(xx[jjx],yy[jjx],c=clr,marker=mkr,ec='k',s=10,lw=0.2)
                                lims[mtr].append([np.nanmin([xx,yy]),np.nanmax([xx,yy])])
                    if ri==1: #Histograms
                        yy=[]
                        for ph,clr in zip(phases,colors):
                            band=prefbands[ph]
                            if mtr=='coh':y=np.array([usnr.coh.__dict__[m].Average(band,fn=fn)[idx] for m in mthds])
                            else:y=np.array([usnr.snr.__dict__[m].R()[ph].Average(band,fn=fn)[idx] for m in mthds])
                            yy.append(y[0]-y[1]) #<---TFZ - HPSZ assuming mthds=['TF_Z','HPS_Z]

                        bins=np.arange(0,cohmax+cohdx,cohdx) if mtr=='coh' else np.arange(0,snrmax+snrdx,snrdx)
                        bins[-1]=cohmax if mtr=='coh' else snrmax
                        bins=list(bins);bins.extend(list(-np.array(bins)));bins=np.sort(np.unique(bins))
                        phname={'P':'P/Pdiff','S':'S/Sdiff','Rg':'Rayleigh'}
                        label=[f'{int(min(prefbands[ph]))}-{int(max(prefbands[ph]))}s ({phname[ph]})' if mtr=='snr' else f'{int(min(prefbands[ph]))}-{int(max(prefbands[ph]))}s' for ph in phases]

                        N_total = sum(len(i) for i in yy)
                        weights = [np.ones(len(i)) / N_total * 100 for i in yy]
                        # sorti=np.flip(np.argsort([np.mean(i[i>=0])-np.mean(i[i<=0]) for i in yy]))
                        sorti=[2,1,0] #Order of histogram plots
                        # if fn is None:sorti=np.flip(sorti)
                        ax.hist(list(np.array(yy)[sorti]),
                        bins=bins,weights=weights,color=np.array(colors)[sorti],
                        edgecolor='k',alpha=1.0 if (max(band) == 10) else 1.0,  # your logic
                        label=np.array(label)[sorti],stacked=True,zorder=1e3,)
                        from matplotlib.ticker import PercentFormatter
                        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100))
                        if (ri==1)&(axi==0):ax.set_ylabel(f'Percent of source-receivers ({f'{str(sum(idx))[0]},{str(sum(idx))[1:]}'})',ha='center',labelpad=2.5,fontsize=6)
                        ax.set_ylim((0.01528455995339294, 150.95276962368143))
                        ax.set_yscale('log')
                        yt=[1.e-01, 1.e+00, 1.e+01, 1.e+02];ax.set_yticks(yt);ax.set_yticklabels(yt)
                        if mtr=='coh':xlabel=[yttl('HPS Z'),yttl('TF Z')]
                        else:xlabel=[rf'$\eta_{{{'HPS Z'}}}$',rf'$\eta_{{{'TF Z'}}}$']
                        spaces=r'$\quad$'*3;arrowspaces=r'$\quad$'*14
                        xlabel = f'$\\leftarrow${arrowspaces}$\\rightarrow$ \n {xlabel[0]} is higher {spaces}{xlabel[1]} is higher'
                        ax.set_xlabel(xlabel,fontsize=6,labelpad=3)
                        ax.axvline(0,c='r',lw=1,ls='-.',zorder=1e5)

                    # if (ri==1)&(axi==1):
                    if (ri==1):
                        handles, labels = ax.get_legend_handles_labels()
                        ax.legend(handles=list(np.array(handles[:3])[np.argsort(sorti)]),frameon=True,fontsize=4,framealpha=1)
                    if ri==0:
                        ax.xaxis.set_ticks_position('top')      # ticks at top
                        ax.xaxis.set_label_position('top')      # x-label at top
                        ax.tick_params(labelbottom=False, labeltop=True)
                        if axi==1:
                            ax.yaxis.set_ticks_position('right')      # ticks at top
                            ax.yaxis.set_label_position('right')      # x-label at top
                            ax.tick_params(labelleft=False, labelright=True)
                    else:
                        if axi==1:
                            ax.yaxis.set_ticks_position('right')      # ticks at top
                            ax.yaxis.set_label_position('right')      # x-label at top
                            ax.tick_params(labelleft=False, labelright=True)


            for axi,ax in enumerate(axes[0,:]):
                if axi==0:
                    # ax.set(xlim=[-.05,1.05],ylim=[-.05,1.05])
                    line=[ax.get_xlim(),ax.get_ylim()];ax.set(xlim=[np.min(line),np.max(line)],ylim=[np.min(line),np.max(line)])
                    ax.plot([np.min([ax.get_ylim(),ax.get_xlim()]),np.max([ax.get_ylim(),ax.get_xlim()])],[np.min([ax.get_ylim(),ax.get_xlim()]),np.max([ax.get_ylim(),ax.get_xlim()])],ls=':',c='r')
                else:
                    line=[ax.get_xlim(),ax.get_ylim()]
                    ax.plot([np.min(line),np.max(line)],[np.min(line),np.max(line)],ls=':',c='r')
                    # ax.set(xlim=line[0],ylim=line[1])
                    ax.set(xlim=[np.min(line),np.max(line)],ylim=[np.min(line),np.max(line)])
            for ax in axes.reshape(-1):ax.grid(True,zorder=-1e10,alpha=0.3,color='k')
            for j,ax in zip([yttl,yttl_eta],axes[0,:]):ax.set_xlabel(j('HPS Z'),fontsize=fontsize,labelpad=7);ax.set_ylabel(j('TF Z'),fontsize=fontsize)
            sensitivity_text = {None:'regardless of infragravity sensitivity','IG':'within the infragravity sensitivity limit','MS':'outside the infragravity sensitivity limit'}
            center_text = f'{xci}\n{sensitivity_text[fn]}'
            ax=axes[1,0];ax.text(1.1,1.20,center_text,transform=ax.transAxes,ha='center',bbox=dict(facecolor='white', alpha=1.0,edgecolor='w',linewidth=0.1),fontsize=7,fontweight='bold')
            fig.tight_layout()
            fig.subplots_adjust(hspace=0.5,wspace=0.15)
            # for ax in axes[0,:]:ax.xaxis.set(labelpad=4)
            fold=plotfolder
            if not xcat=='All_data':fold=fold/xcat
            fold.mkdir(parents=True,exist_ok=True)
            file=f'{xci.replace(' ','_')}._{'RegardlessofIG' if fn is None else fn}_BasicComparativeHist.png'
            _=save_tight(fold/file,fig,dpi=900)
            print(f'{xcat} : {xci} - Done')
