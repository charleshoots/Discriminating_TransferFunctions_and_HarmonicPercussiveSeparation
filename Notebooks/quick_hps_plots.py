from imports import *
import gc
# report=get_reports('ZZ',catalog,path_library.local.Archive,AVG=False)
def hps_corrected_list(stanm,catalog):
    sta=catalog.loc[stanm]
    events=np.array([e.name.split(stanm+'.')[1].strip('.HZ.SAC') for e in list((dirs.Events_HPS/stanm/'CORRECTED').glob('*.HZ.SAC'))])
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    events =events[ib]
    events=events[np.where(np.array([np.sum([len(list((dirs.Events_HPS/stanm/'CORRECTED').glob(f'*{e}*{c}*.SAC'))) for c in ['Z','1','2','HDH']]) for e in events])==4)]
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    events=Catalog([sta.Events[i] for i in ia])
    return events


def make_plot(r,c):
    fig,axes = plt.subplots(nrows=len(channels),ncols=(len(fbands)),figsize=(30+5,16));axes=np.atleast_2d(axes)
    for rax,rr,cc in zip(axes,r,c): #component
        chan=rr.stats.channel
        for ax,band in zip(rax,fbands): #filter
            filt_rr,filt_cc=rr.copy(),cc.copy()
            axtitle = f'{chan} | {band[0]}-{band[1]}s'
            filt_rr.filter('bandpass',freqmin=1/band[1],freqmax=1/band[0],zerophase=True,corners=4)
            filt_cc.filter('bandpass',freqmin=1/band[1],freqmax=1/band[0],zerophase=True,corners=4)
            ax.plot(filt_rr.times(),filt_rr.data,c='r',linewidth=0.4)
            ax.plot(filt_cc.times(),filt_cc.data,c='k',linewidth=0.4)
            ax.set_title(axtitle,fontweight='bold',fontsize=20)
            ax.set_xlim(filt_rr.times()[0],filt_rr.times()[-1])
            yl=max([abs(filt_rr.data).max(),abs(filt_rr.data).min()])
            ax.set_ylim(-yl,yl)
            del filt_rr,filt_cc
        del rr,cc
    fig.suptitle(f'{ev.Name}| {ev.magnitudes[0].magnitude_type} {ev.magnitudes[0].mag}',fontweight='bold',fontsize=24,y=0.93)
    return fig

channels=['HZ','H1','H2','HDH']
fbands = [[1,10],[10,30],[30,100]]
for si,stanm in enumerate(catalog.StaName):
    _=gc.collect()
    stafold=dirs.Events_HPS/stanm
    events=hps_corrected_list(stanm,catalog)
    events_corrected_st=[[load_sac(stafold/'CORRECTED'/f'{stanm}.{e.Name}.{c}.SAC',rmresp=False) for c in channels if np.sum([(stafold/'rmresp'/f'{e.Name}.{c}.SAC').exists() for c in channels])==4] for e in events]
    tmp = [Stream([tr[0] for tr in st]) for st in events_corrected_st if sum([len(tr) for tr in st])==4];del events_corrected_st
    events_corrected_st=tmp;del tmp
    times = [e[0].stats.starttime for e in events_corrected_st]
    events_raw_st=[[load_sac(stafold/'rmresp'/f'{e.Name}.{c}.SAC',rmresp=False) for c in channels if np.sum([(stafold/'rmresp'/f'{e.Name}.{c}.SAC').exists() for c in channels])==4] for e in events]
    tmp = [Stream([tr[0] for tr in st]) for st in events_raw_st if sum([len(tr) for tr in st])==4];del events_raw_st
    events_raw_st=tmp;del tmp
    [e.trim(t,t+7200) for e,t in zip(events_raw_st,times)]
    for e in events_raw_st:
        for tr in e:tr.data=tr.data[:72000]

    for evi,(r,c) in enumerate(zip(events_raw_st,events_corrected_st)): #event
        r.detrend();c.detrend()
        clear_output(wait=False)
        ev=events[evi]
        print(f'{stanm} | {si+1}/{len(catalog)} | {ev.Name} {evi+1}/{len(events_raw_st)}')
        plotfold=dirs.Events_HPS/stanm/'CORRECTED'/'EventRecords'
        plotfold.mkdir(exist_ok=True)
        file=plotfold/f'{stanm}.HPSEventCheck.{ev.Name}.png'
        if file.exists():continue
        fig=make_plot(r,c)
        save_tight(file,fig)
        plt.close('all')
        # r.clear();c.clear();
        del r,c
    # [e.clear() for e in events_raw_st]
    # [e.clear() for e in events_corrected_st]
    del events_corrected_st,events_raw_st