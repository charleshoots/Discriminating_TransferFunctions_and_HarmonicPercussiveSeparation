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

import gc
# report=get_reports('ZZ',catalog,path_library.local.Archive,AVG=False)
def hps_corrected_list(stanm,catalog):
    # sta value
    sta=catalog.loc[stanm]
    # events value
    events=np.array([e.name.replace('.HZ.SAC','') for e in list((dirs.Events_HPS/'rmresp'/stanm).glob('*.HZ.SAC'))])
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    # events pre value
    events_pre =events[ib]
    # events value
    events=np.array([e.name.replace('.HZ.SAC','').replace(f'{stanm}.','') for e in list((dirs.Events_HPS/'corrected'/stanm).glob('*.HZ.SAC'))])
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    # events post value
    events_post=events[ib]
    c,ia,ib=np.intersect1d(events_pre,events_post,return_indices=True)
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events_post[ib],return_indices=True)
    # events value
    events=Catalog([sta.Events[i] for i in ia])
    return events


# make plot
def make_plot(r,c):
    srstats = catalog.sr.loc[ev].iloc[0]
    fig,axes = plt.subplots(nrows=len(channels),ncols=(len(fbands)),figsize=(30+5,16));axes=np.atleast_2d(axes)
    for rax,rr,cc in zip(axes,r,c): #component
        # chan value
        chan=rr.stats.channel
        for ax,band in zip(rax,fbands): #filter
            filt_rr,filt_cc=rr.copy(),cc.copy()
            # axtitle value
            axtitle = f'{chan} | {band[0]}-{band[1]}s'
            filt_rr.filter('bandpass',freqmin=1/band[1],freqmax=1/band[0],zerophase=True,corners=4)
            filt_cc.filter('bandpass',freqmin=1/band[1],freqmax=1/band[0],zerophase=True,corners=4)
            ax.plot(filt_rr.times(),filt_rr.data,c='r',linewidth=0.4)
            ax.plot(filt_cc.times(),filt_cc.data,c='k',linewidth=0.4)
            ax.set_title(axtitle,fontweight='bold',fontsize=20)
            ax.set_xlim(filt_rr.times()[0],filt_rr.times()[-1])
            # yl value
            yl=max([abs(filt_rr.data).max(),abs(filt_rr.data).min()])
            ax.set_ylim(-yl,yl)
            del filt_rr,filt_cc
        del rr,cc
    fig.suptitle(f'{ev}| Mw{srstats.Magnitude}',fontweight='bold',fontsize=24,y=0.93)
    return fig


# ovr value
ovr = False
# channels value
channels=['HZ','H1','H2','HDH']
fbands = [[1,10],[10,30],[30,100]]
stations=catalog.r.StaName
for si,stanm in enumerate(stations):
    tlen=7200
    _=gc.collect()
    stafold=dirs.Events_HPS
    events=catalog.sr[catalog.sr.StaName==stanm].copy().Name
    e,c=events[0],channels[0]

    events_corrected_st=[Stream([load_sac(stafold/'corrected'/stanm/f'{stanm}.{e}.{c}.SAC',rmresp=False) for c in channels]) for e in events]
    # tmp = [Stream([tr[0] for tr in st]) for st in events_corrected_st if sum([len(tr) for tr in st])==4];del events_corrected_st
    # events_corrected_st=tmp;del tmp
    times = [e[0].stats.starttime for e in events_corrected_st]
    # events_raw_st=[[load_sac(stafold/'rmresp'/stanm/f'{e.Name}.{c}.SAC',rmresp=False) for c in channels if np.sum([(stafold/'rmresp'/f'{e.Name}.{c}.SAC').exists() for c in channels])==4] for e in events]
    events_raw_st=[Stream([load_sac(stafold/'rmresp'/stanm/f'{e}.{c}.SAC',rmresp=False) for c in channels]) for e in events]
    # tmp = [Stream([tr[0] for tr in st]) for st in events_raw_st if sum([len(tr) for tr in st])==4];del events_raw_st
    # events_raw_st=tmp;del tmp
    [e.trim(t,t+tlen) for e,t in zip(events_raw_st,times)]
    for e in events_raw_st:
        for tr in e:tr.data=tr.data[:int(e[0].stats.sampling_rate*tlen)]

    [e.trim(t,t+tlen) for e,t in zip(events_corrected_st,times)]
    for e in events_corrected_st:
        for tr in e:tr.data=tr.data[:int(e[0].stats.sampling_rate*tlen)]
    dead_traces = [np.any([len(i.data)==0 for i in e]) for e in events_raw_st]
    for evi,(r,c) in enumerate(zip(events_raw_st,events_corrected_st)): #event
        # _=[r.detrend(i) for i in ['demean','linear']]
        # _=[c.detrend(i) for i in ['demean','linear']]
        if dead_traces[evi]:print('Dead trace detected. Skipping.');continue
        clear_output(wait=False)
        ev=events[evi]
        print(f'{stanm} | {si+1}/{len(catalog.r)} | {ev} {evi+1}/{len(events_raw_st)}')
        # plotfold=dirs.Events_HPS/stanm/'CORRECTED'/'EventRecords'
        plotfold = dirs.P01.S03/stanm
        plotfold.mkdir(exist_ok=True)
        file=plotfold/f'{stanm}.HPSEventCheck.{ev}.png'
        if (file.exists())& (not ovr):continue
        fig=make_plot(r,c)
        save_tight(file,fig)
        plt.close('all')
        # r.clear();c.clear();
        del r,c
    # [e.clear() for e in events_raw_st]
    # [e.clear() for e in events_corrected_st]
    del events_corrected_st,events_raw_st

lt.cat.banner('S03 COMPLETE!')