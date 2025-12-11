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

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

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
figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
from obspy.signal.trigger import classic_sta_lta,carl_sta_trig,recursive_sta_lta
# stalta methods value
stalta_methods={'classic':classic_sta_lta,'carl':carl_sta_trig,'recurssive':recursive_sta_lta}
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
suspect_stations=np.array(['ZA.B02','YL.C09W','7D.G25B','7D.FS08D','7D.G17B','YL.A14W'])
baz=lambda s:obspy.geodetics.base.gps2dist_azimuth(s.Latitude,s.Longitude,s.LaLo[0],s.LaLo[1])[1]
bootstrap = lambda y,nruns=10000,nchoose=100,aggregate=np.mean: np.mean([aggregate(np.random.choice(y[~np.isnan(y)],nchoose)) for _ in range(nruns)])

phases=['P','S','Rg'];preferred_pbands={'P':'1_10','S':'10_30','Rg':'30_100'}
methods=['NoiseCut','ATaCR'];mnames={'NoiseCut':'HPS','ATaCR':'TF'}
mnames_r={mnames[k]:k for k in mnames.keys()};mname_comp={f'HPS_Z':'NoiseCut','TF':'ATaCR','Original':'Original'}
mnames_comp_r={mname_comp[k]:k for k in mname_comp.keys()}
cohnames2snrnames=c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_H':'HPS.H'}

yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"


cat=catalog.copy()
isnr=cat.sr.copy();ratio=True;usnr=unpack_metrics(isnr,ratio=ratio,log=True,notched=False)


fig,axes=figs(3,f=(6,6),x='all',y='all')

phases=['P','S','Rg']


# aggregate= lambda m:np.array([sr.Coherence[m].mean() for sr in isnr.iloc])
# aggregate= lambda m:np.array([np.median(sr.Coherence[m]) for sr in isnr.iloc])
aggregate= lambda m: np.nanmedian(np.array([usnr.coh[m][p] for p in phases]),axis=0)


ax=axes[0];m='TF.HZ'
x=np.nanmean(np.array([usnr.snr[m][p] for p in phases]),axis=0)
y=aggregate('TF.HZ')
_=[ax.scatter(xi,yi,edgecolors='k',color=ColorStandard.instrument[sr.Instrument_Design],marker=ColorStandard.seismometer_marker[sr.Seismometer]) for xi,yi,sr in zip(x,y,isnr.iloc)]
ax.set_ylabel(yttl(m.replace('.',' ')))

####################################################################################
ax=axes[1];m='HPS.HZ'
x=np.nanmedian(np.array([usnr.snr[m][p] for p in phases]),axis=0)
y=aggregate('HPS.HZ')
_=[ax.scatter(xi,yi,edgecolors='k',color=ColorStandard.instrument[sr.Instrument_Design],marker=ColorStandard.seismometer_marker[sr.Seismometer]) for xi,yi,sr in zip(x,y,isnr.iloc)]
ax.set_ylabel(yttl(m.replace('.',' ')))

####################################################################################
ax=axes[2];m='HPS.H'
x=np.nanmedian(np.array([usnr.snr[m][p] for p in phases]),axis=0)
y=aggregate('HPS.H')
_=[ax.scatter(xi,yi,edgecolors='k',color=ColorStandard.instrument[sr.Instrument_Design],marker=ColorStandard.seismometer_marker[sr.Seismometer]) for xi,yi,sr in zip(x,y,isnr.iloc)]
ax.set_ylabel(yttl(m.replace('.',' ')))

ax.set_xlabel(r'$R_{SNR}$')


ax.set_xlim([-1.05,2.05])
ax.set_ylim([-.05,1.05])


fold=dirs.SNR/'Plots'/'QuickScatter'
fold.mkdir(exist_ok=True,parents=True)

file='QuickScatter.png'
save_tight(fold/file,fig,dpi=800)