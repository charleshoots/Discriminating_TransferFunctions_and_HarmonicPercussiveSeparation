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
from source.imports import *
from source.modules import *
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
# figs value
figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
from obspy.signal.trigger import classic_sta_lta,carl_sta_trig,recursive_sta_lta
# stalta methods value
stalta_methods={'classic':classic_sta_lta,'carl':carl_sta_trig,'recurssive':recursive_sta_lta}
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
# f value
f=cat.sr.iloc[0].Data.Coherence().f
# function meancoh
def meancoh(icat,methods=['TF','HPS_Z','HPS_H'],bands=['1_10','10_30','30_100']):
    #Calculates average coherence for every source-receiver pair in every method and every band defined in bands.
    aggregate=lambda b,method,s,agg=np.mean:agg(s.Coherence[method].reshape(-1)[((1/f)>=min(b))&((1/f)<=max(b))])
    # coh value
    coh=AttribDict({b:AttribDict({m:np.array([aggregate([int(i) for i in b.split('_')],m,s)
    for s in icat.iloc]) for m in methods}) for b in bands})
    return coh
# function unpacksnr
def unpacksnr(icat,bands=['1_10','10_30','30_100'],phases=['P','S','Rg'],methods=['ATaCR','NoiseCut'],ratio=False):
    # mnames value
    mnames={'ATaCR':'TF','NoiseCut':'HPS_Z','Original':'Original'}
    if ratio:
        # isnr value
        isnr=AttribDict({b:AttribDict({mnames[m]:AttribDict({p:np.array([s.SNR[b][p][m]/s.SNR[b][p]['Original'] if p in s.SNR[b].keys() else None
        for s in icat.iloc]) for p in phases}) for m in methods}) for b in bands})
    else:
        # isnr value
        isnr=AttribDict({b:AttribDict({mnames[m]:AttribDict({p:np.array([s.SNR[b][p][m] if p in s.SNR[b].keys() else None
        for s in icat.iloc]) for p in phases}) for m in methods}) for b in bands})
    return isnr

# function preproc
def preproc(sti,freqmin,freqmax,prenoise_buffer,filtertype,plottend):
    # nproc value
    nproc=2
    for _ in range(nproc):sti.detrend('simple').detrend('demean').detrend('simple')
    sti.taper(prenoise_buffer)
    sti.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase={'causul':False,'acausul':True}[filtertype]).taper(prenoise_buffer)
    for _ in range(nproc):sti.detrend('simple').detrend('demean').detrend('simple')
    sti.taper(prenoise_buffer)
    sti.trim(sti[0].stats.starttime,plottend)
    return sti

# plot cut record
def plot_cut_record(ev,icutdf,icutsnr,bands=['1_10','10_30','30_100'],methods=['NoiseCut.HZ','Original.HZ','ATaCR.HZ'],h_per_tr=1,width=17,note=''):
    # phasecolors value
    phasecolors={'P':'lightblue','S':'red','Rg':'green','Noise':'lightgray'}
    # prenoise buffer value
    prenoise_buffer=.01
    # evn value
    evn=[ev] if isinstance(ev,str) else ev
    assert len(evn)>0, 'not a string or list of string events'
    # figs value
    figs=[]
    # status value
    status=lambda:print(f'{b}s | {ev} | {method}')
    # forcedf value
    forcedf=lambda df:df if isinstance(df,pd.DataFrame) else pd.DataFrame([df])
    for bi,b in enumerate(bands):
        bn=b
        band=[int(i) for i in b.split('_')]
        freqmin,freqmax=1/max(band),1/min(band)
        for ei,ev in enumerate(evn):
            evidx = np.array(list(icutdf.Name==ev))
            stnm = np.array(list(icutdf.StaName))
            iSR=forcedf(icutdf.aloc[ev].copy())
            iSR.sort_values(by=['Distance'],inplace=True,ascending=True)
            st=Stream()
            _=[[st.append(tr) for tr in i.Traces()] for i in iSR.iloc]
            st.trim(iSR.Origin[0],iSR.Origin[0]+7200)

            # ==========================================================================================
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
                # override_sort=f'{str(np.round(kilometers2degrees(st[tri].stats.distance),2)).zfill(6)}°'
                override_sort=f'{str(int(tr_plot_sort[tri])).zfill(2)} | {str(np.round(kilometers2degrees(st[tri].stats.distance),2))}°' #very sensitive do not tweak this line
                st[tri].stats.network=override_sort
            for method in methods:

                metrics={'snr':lambda bn,s,ph:icutsnr.snr.__dict__[method].R().__dict__[ph].Average([float(bi) for bi in bn.split('_')])[evidx&(stnm==s)],
                'ST':lambda bn,s,ph:icutsnr.ST.__dict__[method].__dict__[ph].Average([float(bi) for bi in bn.split('_')])[evidx&(stnm==s)],
                'LT': lambda bn,s:icutsnr.LT.__dict__[method].Average([float(bi) for bi in bn.split('_')])[evidx&(stnm==s)],
                'coh': lambda bn,s:icutsnr.coh.__dict__[method].Average([float(bi) for bi in bn.split('_')])[evidx&(stnm==s)]}
                
                status()
                dist_degree=False
                plottend = (iSR.iloc[0].Traces()[0].stats.starttime + (1*(degrees2kilometers(iSR.Distance.max())/2.5))) 
                ist=st.select(location=method.split('.')[0].replace('TF','ATaCR').replace('HPS','NoiseCut'),channel=method.split('.')[-1]).copy()
                bst=st.select(location='Original.HZ'.split('.')[0].replace('TF','ATaCR').replace('HPS','NoiseCut'),channel='Original.HZ'.split('.')[-1]).copy()
                x=np.array([i.stats.distance for i in ist])
                if dist_degree:x=kilometers2degrees(x)


                bst=preproc(bst.copy(),freqmin,freqmax,prenoise_buffer,filtertype,plottend)
                itr=preproc(ist.copy(),freqmin,freqmax,prenoise_buffer,filtertype,plottend)


                yl=np.array([max([abs(i.data).max(),abs(j.data).max()]) for i,j in zip(itr,bst)])*1.12
                fig=itr.plot(equal_scale=False,type='relative',handle=True,ev_coord=iSR.iloc[0].LaLo,dist_degree=dist_degree,endtime=plottend)
                fig.set_figheight(h_per_tr*len(itr));fig.set_figwidth(width)
                axes=[i for i in fig.get_children() if isinstance(i,matplotlib.axes._axes.Axes)]
                [ax.set_ylim(-y,y) for y,ax in zip(yl,axes)]
                # ysta=['.'.join(i.get_children()[1]._text.split('.')[2:-2]) for i in axes]
                ysta=['.'.join(i.get_children()[1]._text.split('°.')[-1].split('.')[:2]) for i in axes]
                for si,s in enumerate(ysta):
                    s=iSR[iSR.StaName==s].iloc[0]
                    Ph=s.Phases(phases=('P','S'))
                    Ph={k:s.Origin+Ph[k][0] for k in Ph.keys()}
                    metrics = icutsnr[(icutsnr.StaName==s.StaName)&(icutsnr.Name==s.Name)].iloc[0].SNR.copy()
                    wins = icutsnr[(icutsnr.StaName==s.StaName)&(icutsnr.Name==s.Name)].iloc[0].wins[bn].copy()
                    if True:wins.update({k:[i-s.Origin for i in wins[k]] for k in wins.keys() if k in ['P','S','Rg','Noise']})
                    # tr=s.Traces().copy()
                    x=degrees2kilometers(s.Distance)
                    ax = axes[si]
                    # ax.get_children()[1].set_text(f'{ax.get_children()[1].get_text()}.{int(s.StaDepth)}m')
                    ax.get_children()[1].set_text(f'{ax.get_children()[1].get_text()}.{int(s.StaDepth)}m ({(1/fnotch(s.StaDepth)):.1f}s)')
                    ax.get_children()[1].set_ha('right');ax.get_children()[1].set_position((.99,.95));ax.get_children()[1].set_fontweight('bold')
                    for bn in bands:
                        [ax.axvspan(wins[k][0], wins[k][1],alpha=0.35,facecolor=phasecolors[k],edgecolor=phasecolors[k],linewidth=3 if b==bn else 1) for k in wins.keys() if k in ['P','S','Rg','Noise']]
                        [ax.axvspan(wins[k][0], wins[k][1],alpha=1,facecolor='None',edgecolor='darkblue',linewidth=1.5,zorder=1e4) for k in wins.keys() if (b==bn)&(k in ['P','S','Rg','Noise'])]
                        if b==bn:
                            if bn=='1_100':bn='30_100'
                            # fontweight='normal'
                            fontweight='bold'
                            ym=np.min(ax.get_ylim());ymx = max(ax.get_ylim());xmx=max(ax.get_xlim())
                            texth=lambda k:{'P':'right' if PSwingap<240 else 'left','S':'left','Rg':'right'}[k]
                            xpos=lambda k:(10 if k in ['S','P'] else 0)+ {'P':np.min(wins[k]) if PSwingap<240 else np.max(wins[k]),'S':np.max(wins[k]),'Rg':min([0.99*np.max(wins[k]),xmx*.99])}[k]
                            text=lambda k:f'{k}: {metrics['snr'](bn,s.StaName,k):.2f}\nST:{metrics['ST'](bn,s.StaName,k):.2e}\nLT:{metrics['LT'](bn,s.StaName):.2e})'
                            PSwingap = (min(wins['S'])-max(wins['P'])) if ('P' in wins.keys())&('S' in wins.keys()) else 0
                            [ax.text(xpos(k)/xmx,.1*ymx,text(k),zorder=1e5,ha=texth(k),fontweight=fontweight,
                            bbox=dict(facecolor='white', edgecolor='black',boxstyle='square,pad=0.2'),va='bottom',fontsize=7,transform=ax.transAxes) for ki,k in enumerate(wins.keys()) if k in ['P','S','Rg']]
                title=f'{fig.get_suptitle()} | {ev} |  {method} | Mw {iSR.iloc[0].Magnitude} ' + f'({int(iSR.Distance.min())}-{int(iSR.Distance.max())}°)' + f' | {1/freqmax} to {1/freqmin}s'
                if len(note)>0:title=f'{note}\n{title}'
                fig.suptitle(title,y=1.05,fontweight='bold')
                figs.append(fig)
    if len(figs)==1:figs=figs[0]
    return figs
def index(df):
    remcols=['Name_idx', 'StaName_idx', 'Station_idx', 'Network_idx',
    'Experiment_idx', 'Pressure_Gauge_idx', 'Seismometer_idx',
    'Instrument_Design_idx']
    cols=df.columns;cols=[c for c in cols if c not in remcols]
    df=df[cols]
    index_cols=['Name','StaName','Station','Network',
    'Experiment','Pressure_Gauge','Seismometer','Instrument_Design']
    idx_dummy=['_'*(ci+1) for ci,c in enumerate(index_cols)]
    for src, dst in zip(index_cols, idx_dummy):df[dst]=df[src]
    df.set_index(idx_dummy,drop=True,inplace=True)
    return df



icat=cat.sr.copy()
usnr=unpack_metrics(icat)

idx=(usnr.coh.HPS_Z.Average((1,10))<.94)&(icat.Magnitude>7.0);cutname='Mw7-8'
icat=icat[idx]
usnr=unpack_metrics(icat)

evn=icat.Name.unique()
# --------------------------------------------
# --------------------------------------------
# [5] PLOT THE ANSWER-------------------------
# --------------------------------------------
bands=['30_100','10_30','1_10',];phases=['Rg','P','S']
methods=['Original.HZ','TF.HZ','HPS.HZ']
mnames_r={'TF':'ATaCR','HPS_Z':'NoiseCut','HPS_H':'NoiseCut'}
plotfold=dirs.SNR/'Plots'/'CutSets'
(plotfold/cutname).mkdir(exist_ok=True,parents=True)
forcedf=lambda df:df if isinstance(df,pd.DataFrame) else pd.DataFrame([df])
for b in bands:
    fold=plotfold/cutname/b
    fold.mkdir(exist_ok=True,parents=True)
    # print(f'{ssi+1}/{len(scen)} . {sc} . {b}s')
    for evi,ev in enumerate(evn):
        print(f'{evi+1} / {len(evn)} events ')
        for mi,m in enumerate(methods):
            evfile = f'{ev}_Mw{forcedf(icat.loc[ev]).Magnitude[0]}'
            mfile=f'__{str(mi+1).zfill(2)}.{m}__{b}s'
            file=''.join([evfile,mfile,'.png'])
            if (fold/file).exists():continue
            fig=save_tight(fold/file,plot_cut_record(ev,icat,usnr,methods=[m],bands=[b]))
            plt.close(fig)