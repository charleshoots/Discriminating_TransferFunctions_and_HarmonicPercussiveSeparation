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
import time;start=time.time()
import inspect
from textwrap import fill
# runtime value
runtime=lambda:int(time.time()-start)
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
figs = lambda r=3,c=1,f=(5,6),x='all',y='all':plt.subplots(r,c,figsize=f,sharex='all',sharey='all',layout='constrained')




# methods value
methods = ['Original','NoiseCut','ATaCR']
# SR value
SR = cat.sr.copy();bands = ['1_10','10_30','30_100']
# noisewlen value
noisewlen=300;bleed_buffer = 24
# prenoise buffer value
prenoise_buffer=.01;p_lead=20;s_lead=20

# rmse value
rmse=lambda noise:( (  ( noise - noise.mean() )**2  ).mean())**.5
# snr equation value
snr_equation = f'RMSE(noise) = {inspect.getsource(rmse).split(':')[-1].replace('\n','')} | SNR = max(abs(signal)) / RMSE(abs(noise))'
# P snr wlen value
P_snr_wlen=AttribDict({'1_10':150,'10_30':150,'30_100':150,})
# S snr wlen value
S_snr_wlen=AttribDict({'1_10':150,'10_30':150,'30_100':150,})
# Rg snr wlen value
Rg_snr_wlen=AttribDict({'1_10':[2.0,4.2],'10_30':[2.0,4.2],'30_100':[2.0,4.2],})

# wlen value
wlen = {'P':P_snr_wlen,'S':S_snr_wlen,'Rg':Rg_snr_wlen}
# RgWin value
RgWin = lambda x,u:sorted([x/i for i in u])
SR.sort_values(by='Magnitude',inplace=True,ascending=False)
# evn value
evn = SR.Name.unique()
# evn=['2012.080.18.02','2010.129.05.59','2013.285.13.11','2010.199.05.56','2010.103.23.49','2010.064.16.07']

# evn value
evn=['2013.015.16.09','2009.344.02.30', '2009.358.00.23', '2010.094.22.40',
'2010.224.11.54', '2010.298.14.42', '2012.274.16.31',
'2013.143.17.19', '2013.247.02.32','2013.267.11.29',
'2014.103.12.36', '2014.268.17.51',
'2015.115.06.11', '2015.132.07.05']

# ==========================================================================================
# ------------------------------------------------------------------------------------------
# ==========================================================================================
filtertype='causul'
note='V04'
fold=dirs.SNR/'SNR.Models'
file =f'SNR_{filtertype}.filter_{note}.pkl'
SNR=load_pickle(fold/file)
# ==========================================================================================
# ------------------------------------------------------------------------------------------
# ==========================================================================================
bands = ['1_100','1_10','10_30','30_100']
status = lambda:print(f'{method} | {filtertype} | {ev} , Mw{iSR.Magnitude[0]} | {ei+1}/{len(evn)} | {bi+1}/{len(bands)}')
for bi,b in enumerate(bands):
    bn=b
    if bn=='1_100':bn='30_100'
    band=[int(i) for i in b.split('_')]
    freqmin,freqmax=1/max(band),1/min(band)
    for ei,ev in enumerate(evn):
        evSNR=SNR[SNR.Name==ev].copy()
        iSR=SR.loc[ev].copy()
        iSR.sort_values(by=['Distance'],inplace=True,ascending=True)
        st=Stream()
        _=[[st.append(tr) for tr in i.Traces()] for i in iSR.iloc]
        st.trim(iSR.Origin[0],iSR.Origin[0]+7200)
        st.detrend('linear').detrend('demean').taper(prenoise_buffer)
        # ==========================================================================================
        #A super dumb hack to force this plot type to sort by distance...it works though.
        dio=np.array([iSR[iSR.StaName==f'{st[tri].stats.network}.{st[tri].stats.station}'].iloc[0].Distance for tri,tr in enumerate(st)])
        tr_plot_sort=np.argsort(dio)
        for diind,dioi in enumerate(np.unique(dio[tr_plot_sort])):dio[dio==dioi]=diind+1
        tr_plot_sort=dio
        for tri,tr in enumerate(st):
            sn = f'{st[tri].stats.network}.{st[tri].stats.station}'
            trs = iSR[iSR.StaName==sn].iloc[0]
            st[tri].stats.distance=degrees2kilometers(trs.Distance)
            st[tri].stats.coordinates=AttribDict({'latitude':st[tri].stats.sac.stla,'longitude':st[tri].stats.sac.stlo})
            st[tri].stats.station=f'{st[tri].stats.network}.{st[tri].stats.station}'
            # st[tri].stats.network=int(st[tri].stats.distance)

            #very sensitive do not tweak this line
            override_sort=f'{str(int(tr_plot_sort[tri])).zfill(2)} | {str(np.round(kilometers2degrees(st[tri].stats.distance),2))}°' #very sensitive do not tweak this line
            #very sensitive do not tweak this line

            st[tri].stats.network=override_sort
        for method in methods:
            status()
            dist_degree=False
            plottend = (iSR.iloc[0].Traces()[0].stats.starttime + (1*(degrees2kilometers(iSR.Distance.max())/2.5))) 
            ist=st.select(location=method).copy()
            x=np.array([i.stats.distance for i in ist])
            if dist_degree:x=kilometers2degrees(x)
            itr=ist.copy()
            itr.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase={'causul':False,'acausul':True}[filtertype]).taper(prenoise_buffer)
            fig=itr.plot(equal_scale=False,type='relative',handle=True,ev_coord=iSR.iloc[0].LaLo,dist_degree=dist_degree,endtime=plottend)
            fig.set_figheight(10);fig.set_figwidth(17)
            axes=[i for i in fig.get_children() if isinstance(i,matplotlib.axes._axes.Axes)]
            ##### ysta=['.'.join(i.get_children()[1]._text.split('.')[:2]) for i in axes]
            # ysta=['.'.join(i.get_children()[1]._text.split('.')[2:-2]) for i in axes]
            ysta=['.'.join(i.get_children()[1]._text.split('°.')[-1].split('.')[:2]) for i in axes]
            phasecolors={'P':'lightblue','S':'red','Rg':'green','Noise':'lightgray'}
            for si,s in enumerate(ysta):
                s=iSR[iSR.StaName==s].iloc[0]
                Ph=s.Phases(phases=('P','S'))
                Ph={k:s.Origin+Ph[k] for k in Ph.keys()}
                sidx=np.array(icat.StaName==s)
                metrics = SNR[(SNR.StaName==s.StaName)&(SNR.Name==s.Name)].iloc[0].SNR
                wins = SNR[(SNR.StaName==s.StaName)&(SNR.Name==s.Name)].iloc[0].wins[bn].copy()
                if True:wins.update({k:[i-s.Origin for i in wins[k]] for k in wins.keys() if k in ['P','S','Rg','Noise']})
                # tr=s.Traces().copy()
                x=degrees2kilometers(s.Distance)
                ax = axes[si]
                # ax.get_children()[1].set_text(f'{ax.get_children()[1].get_text()}.{int(s.StaDepth)}m')
                ax.get_children()[1].set_text(f'{ax.get_children()[1].get_text()}.{int(s.StaDepth)}m ({(1/fnotch(s.StaDepth)):.1f}s)')
                ax.get_children()[1].set_ha('right');ax.get_children()[1].set_position((.99,.95));ax.get_children()[1].set_fontweight('bold')
                for bn in bands:
                    def gdat(usnr,icat,sta,mthd,b):
                        idx=np.array(icat.StaName==sta)
                        ST=usnr.__dict__['ST'].__dict__[mthd].R().Average((b))[idx]
                        LT=usnr.__dict__['LT'].__dict__[mthd].R().Average((b))[idx]
                        snr=usnr.__dict__['snr'].__dict__[mthd].R().Average((b))[idx]
                        coh=usnr.__dict__['coh'].__dict__[mthd].R().Average((b))[idx]
                        return ST,LT,snr,coh
                    ST,LT,snr,coh=gdat(usnr,icat,s,method,np.sort([int(i) for i in bn.split('_')]))
                    [ax.axvspan(wins[k][0], wins[k][1],alpha=0.35,facecolor=phasecolors[k],edgecolor=phasecolors[k],linewidth=3 if b==bn else 1) for k in wins.keys() if k in ['P','S','Rg','Noise']]
                    [ax.axvspan(wins[k][0], wins[k][1],alpha=1,facecolor='None',edgecolor='darkblue',linewidth=1.5,zorder=1e4) for k in wins.keys() if (b==bn)&(k in ['P','S','Rg','Noise'])]
                    if b==bn:
                        if bn=='1_100':bn='30_100'
                        # fontweight='normal'
                        fontweight='bold'
                        ym=np.min(ax.get_ylim());ymx = max(ax.get_ylim());xmx=max(ax.get_xlim())
                        texth=lambda k:{'P':'right' if PSwingap<240 else 'left','S':'left','Rg':'right'}[k]
                        xpos=lambda k:(10 if k in ['S','P'] else 0)+ {'P':np.min(wins[k]) if PSwingap<240 else np.max(wins[k]),'S':np.max(wins[k]),'Rg':min([0.99*np.max(wins[k]),xmx*.99])}[k]
                        text=lambda k:f'{k}: {metrics[bn][k][method]:.2f}\nST:{metrics[bn].ST[k][method]:.2e}\nLT:{metrics[bn].LT[method]:.2e})'
                        PSwingap = (min(wins['S'])-max(wins['P'])) if ('P' in wins.keys())&('S' in wins.keys()) else 0
                        [ax.text(xpos(k)/xmx,.1*ymx,text(k),zorder=1e5,ha=texth(k),fontweight=fontweight,
                        bbox=dict(facecolor='white', edgecolor='black',boxstyle='square,pad=0.2'),va='bottom',transform=ax.transAxes) for ki,k in enumerate(wins.keys()) if k in ['P','S','Rg']]
            fig.suptitle(fig.get_suptitle()+ f' |  {method} | Mw {iSR.iloc[0].Magnitude} ' + f'({int(iSR.Distance.min())}-{int(iSR.Distance.max())}°)' + f' | {1/freqmax} to {1/freqmin}s\n' + snr_equation,y=1.06)
            fig.set_figheight(10);fig.set_figwidth(17)
            fold = dirs.SNR/'Plots'/f'{note}_{filtertype}'
            fold = fold/'AllEvents'/f'Mw{iSR.Magnitude[0]}_{ev}'
            file = f'{ev}_Mw{iSR.Magnitude[0]}_{method}_{b}s_{filtertype}.filter_{note}.png'
            fold.mkdir(exist_ok=True,parents=True)
            save_tight(fold/file)

            plt.close()
    print(f"{b} Complete. Elapsed time: {(runtime())/60:.2f} minutes")
print(f"Run Complete. Elapsed time: {(runtime())/60:.2f} minutes")