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

# 7D.J20D: 1 (2.9%)
# 2D.OBS23: 3 (5.4%)
# 7D.G34D: 8 (25.0%)
# 7D.J25A: 7 (15.9%)
# 7D.G26D: 9 (27.3%)
# 7D.J17D: 7 (28.0%)
# 7D.FS42D: 6 (24.0%)
# 7D.J41C: 6 (18.2%)
# 7D.J42C: 13 (44.8%)
# 7D.FN07A: 22 (48.9%)
# 7D.M07A: 8 (18.2%)
# 7D.FN14A: 20 (43.5%)
# 7D.FN12C: 11 (36.7%)
# 7D.J50A: 12 (27.3%)
# 7D.FS04D: 6 (24.0%)





# function station records
def station_records(s,b,fold,SR,SNR,srmax=15,methods=['TF.HZ','Original.HZ','HPS.HZ'],prenoise_buffer=.01,filtertype='causul',tmax=7200):
    # preferred bands value
    preferred_bands={'P':'1_10','S':'10_30','Rg':'30_100'}
    # preferred phases value
    preferred_phases={preferred_bands[i]:i for i in preferred_bands.keys()}
    # phases value
    phases=np.unique([preferred_phases[b],'Rg','Noise'])
    # phases=[]
    
    # msplit value
    msplit = lambda m:m.split('.')[0]
    # mdict value
    mdict={'TF':'ATaCR','HPS':'NoiseCut','Original':'Original'}
    # mdict r value
    mdict_r={mdict[k]:k for k in mdict.keys()}
    df=SR[SR.StaName==s].copy()
    df.sort_values(by=['Distance'],inplace=True)
    # itr df value
    itr_df=SNR[SNR.StaName==s].copy()
    itr_df.sort_values(by=['Distance'],inplace=True)
    assert sum((itr_df.Name==df.Name)&(itr_df.StaName==df.StaName))==len(df)

    # st value
    st=Stream(unravel([sta.Traces() for sta in df.iloc]))
    # st evs value
    st_evs=np.array(np.repeat(df.Name,3))
    # band value
    band=[int(i) for i in b.split('_')]
    freqmin,freqmax=1/max(band),1/min(band)
    ret_i=lambda tro,idf:(idf.index.str.startswith(tro.stats.starttime.strftime('%Y.%j.%H')))&(idf.StaName==f'{tro.stats.network}.{tro.stats.station}')
    for tri,(tr,ev) in enumerate(zip(st,st_evs)):
        sn=f'{st[tri].stats.network}.{st[tri].stats.station}'
        trs=df.aloc[ev].iloc[0]
        st[tri].stats.distance=degrees2kilometers(trs.Distance)
        st[tri].trim(trs.Origin,trs.Origin+7200)
    srmax=min([srmax,np.unique([tr.stats.location for tr in st],return_counts=True)[-1].min()])
    dfi=list(np.arange(0,len(df),srmax));dfi.append(len(df));dfi=np.array([[ii,dfi[i+1]] for i,ii in enumerate(dfi[:-1])]);ifo=0
    for fi in dfi:
        ifo+=1
        c_df=df.iloc[fi[0]:fi[1]].copy()
        c_itr_df=itr_df.iloc[fi[0]:fi[1]].copy()
        for method in methods:
            mst=mdict[msplit(method)]
            ret_sr=lambda tro,idf:idf.copy()[(idf.index.str.startswith(tro.stats.starttime.strftime('%Y.%j.%H')))&(idf.StaName==f'{tro.stats.network}.{tro.stats.station}')].copy()
            ret_ar=lambda s,b,p:(np.array(s.wins[b][p]) - s.Origin)/tmax

            label=lambda s,b,p:f'[{p}] SNR: {s.SNR[b][p.replace('Noise','LT')][method]:.2e} |S: {(s.SNR[b].ST[p][method]):.2e} |N: {s.SNR[b].LT[method]:.2e}'
            wintimes=lambda tr,b,p:np.array(ret_sr(tr,SNR).iloc[0].wins[b][p])- ret_sr(tr,SNR).iloc[0].Origin
            ist=st.select(location=mst).copy()[fi[0]:fi[1]].copy()
            x=np.array([i.stats.distance for i in ist])
            dist_degree=False
            if dist_degree:x=kilometers2degrees(x)

            itr=ist.copy()
            tmax=np.ceil(max(c_itr_df[c_itr_df.Distance==c_itr_df.Distance.max()].iloc[0].wins[b]['Rg']) - c_itr_df[c_itr_df.Distance==c_itr_df.Distance.max()].iloc[0].Origin)
            tmax=7200
            # wins=[ret_sr(tr,SNR).iloc[0].wins[b].copy() for tr in itr].copy()
            # wins=[{i:np.array(w[i])-ret_sr(tr,SNR).iloc[0].Origin for i in w.keys() if i in phases} for w,tr in zip(wins,itr)]
            itr.detrend('linear').detrend('demean').taper(prenoise_buffer)
            itr.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase={'causul':False,'acausul':True}[filtertype]).taper(prenoise_buffer)
            fig,axes=figs(r=len(itr),c=1,f=(10,len(itr)),x=True,y=False,layout=None)
            axes=np.atleast_1d(axes).reshape(-1)
            _=[ax.plot(tr.times(),tr.data,c='k',zorder=-1e5,linewidth=1) for tr,ax in zip(itr,axes)]
            _=[ax.set_xlim(0,tmax) for ax in axes]

            bbox=dict(facecolor='white',edgecolor='black',boxstyle='square,pad=0.2',linewidth=0.6)
            props=dict(fontsize=6,fontweight='bold',bbox=bbox)

            [ax.text(0.99,1.1,
            f'{int(s.Distance)}° | MW{s.Magnitude} |{s.Name}',transform=ax.transAxes,
            ha='left', va='center',fontsize=8, fontweight='bold',
            bbox=bbox) for tr,ax,s in zip(itr,axes,c_itr_df.iloc)]

            [[ax.text(1.01,yi,
            label(s,b,p),transform=ax.transAxes,
            ha='left', va='bottom',fontsize=8, fontweight='bold',
            bbox=dict(facecolor='white',edgecolor='black',boxstyle='square,pad=0.3',linewidth=0.4),zorder=1e10) for yi,p in zip([.65,.35,.05],phases) if (p in s.wins[b].keys())&(p in ['P','S','Rg'])] 
            for tr,ax,s in zip(itr,axes,c_itr_df.iloc)]
            phasecolors={'P':'lightblue','S':'red','Rg':'green','Noise':'lightgray','LT':'lightgray'};bn=b
            [[ax.axvspan(s.wins[b][k][0]-s.Origin, s.wins[b][k][1]-s.Origin,alpha=0.35,facecolor=phasecolors[k],edgecolor=phasecolors[k],linewidth=2 if b==bn else 1) for k in s.wins[b].keys() if k in phases] for tr,ax,s in zip(itr,axes,c_itr_df.iloc)]
            [[ax.axvspan(s.wins[b][k][0]-s.Origin, s.wins[b][k][1]-s.Origin,alpha=1,facecolor='None',edgecolor='darkblue',linewidth=1.5,zorder=1e4) for k in s.wins[b].keys() if (b==bn)&(k in phases)] for tr,ax,s in zip(itr,axes,c_itr_df.iloc)]
            fig.subplots_adjust(hspace=0.4)
            fig.tight_layout(rect=[0, 0, 1, 0.97])  # reserve top 3% for the title
            fig.suptitle(f'{method} | {s} ({int(1/fnotch(cat.r.aloc[s].iloc[0].StaDepth))}s) | {len(itr)} source-receivers | {int(1/freqmax)}-{int(1/freqmin)}s',y=1.0,transform=fig.transFigure)

            axes[-1].set_xlabel('Time after origins (s)')

            fold.mkdir(exist_ok=True,parents=True)
            file=f'{b}s_{s}__{str(ifo).zfill(2)}__{method}.png'
            save_tight(fold/file,fig=fig,dpi=400)
            plt.close(fig);plt.close('all')
            print(f'Folder: {str(fold)}')
            print(f'File: {file}')





# # -------------------------
filtertype='acausul';note='V04';fold=dirs.SNR/'SNR.Models';file =f'SNR_{filtertype}.filter_{note}.pkl'
SNR=load_pickle(fold/file)
SR=catalog.sr.copy();SR.sort_values(by=['Name','StaName'],inplace=True)
SNR.sort_values(by=['Name','StaName'],inplace=True)
assert sum((SR.Name==SNR.Name)&(SR.StaName==SNR.StaName))==len(SNR), 'failed start sets'
# snr=unpacksnr(SNR.copy(),methods=['TF.HZ','HPS.HZ','Original'])
# # -------------------------
suspect_stations=cat.r.StaName
suspect_stations=np.array(['ZA.B02','YL.C09W','7D.G25B','7D.FS08D','7D.G17B','YL.A14W'])

prenoise_buffer=.01
methods=['TF.HZ','Original.HZ','HPS.HZ']
bands=['30_100','10_30','1_10']

for b in bands:
    for si,s in enumerate(suspect_stations):
        print(f'{si+1}/{len(suspect_stations)} - {b} - {s}')
        fold=dirs.SNR/'Plots'
        fold=fold/'station_records'/s
        # if fold.exists():continue
        station_records(s,b,fold,SR,SNR,srmax=15,methods=methods,prenoise_buffer=.01,filtertype=filtertype,tmax=7200)
        plt.close('all')