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

# from imports import *
from modules import *
# from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from obspy import read
from obspy.geodetics import locations2degrees
# ========================================================================================================================================================
from IPython.display import clear_output
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from obspy import Stream
from modules import *
from functools import reduce
import local_tools
from local_tools import ObsQA
from local_tools.ObsQA.TOOLS.io import *
# from ObsQA.TOOLS.io import load_sac
k=1

def get_event_record_stream(event,dirs,method='atacr',tlen=7200,step='Corrected',trim_gap=0,P12_folder=None):
    if step.lower()=='raw':rmresp = True
    else:rmresp=False
    stas = event.Stations
    if step.lower()=='corrected':subfold='CORRECTED'
    else:subfold=''
    if method.lower()=='atacr':
        datafolder = dirs.Events;tf='ZP-21'
        files = [''.join([s+'.',event.Name,'.sta.',tf,'.HZ.SAC']) for s in stas]
    else:
        datafolder = dirs.Events_HPS;tf=''
        files = [''.join([event.Name,'.HZ.SAC']) for s in stas]
    if step.lower()=='raw':
        files = [f.replace('.sta.','').replace(tf,'').replace(s+'.','') for f,s in zip(files,stas)]
    st=Stream([ObsQA.TOOLS.io.load_sac(datafolder/s/subfold/('*'+f),rmresp=rmresp) for f,s in zip(files,stas)]);clear_output(wait=False)
    if P12_folder:
        st_P12=get_event_P12(stas,event.Name,P12_folder);st+=st_P12;clear_output(wait=False)
    st.taper(0.0001);st.trim(event.origins[0].time,event.origins[0].time+tlen,pad=True,fill_value=trim_gap)
    clear_output(wait=False)
    return st
def get_avg_noise(stanm,fold):return load_pickle(list((fold/stanm).glob('*sta.pkl'))[0])
# ----
def stream_metrics(st,evn,fold):
        stas = [s.stats.network+'.'+s.stats.station for s in st]
        for si,s in enumerate(stas):
            st[si].Noise = get_avg_noise(s,fold.parent/'AVG_STA')
            st_P12=get_event_P12([s],evn,fold);clear_output(wait=False) #s reqs a list
            st[si].Metrics = OBSM.Metrics(
            tr1=st_P12.select(channel='*1')[0].copy(),
            tr2=st_P12.select(channel='*2')[0].copy(),
            trP=st_P12.select(channel='*H')[0].copy(),
            trZ=Stream(st[si]).select(channel='*Z')[0].copy())
            del st_P12
        return st
def get_event_record(event,cat,dirs,method,P12_folder=None):
    # ________________________________________________________________________________
    for s in event.Stations:
        get_event_record

    st_corrected = get_event_record_stream(event,dirs,method=method,step='corrected',P12_folder=P12_folder)
    st_raw = get_event_record_stream(event,dirs,method=method,step='raw',P12_folder=P12_folder)
    for tri in range(len(st_raw)):
        st_raw[tri].stats.location='Raw'
        st_corrected[tri].stats.location='Corrected-'+method.replace('HPS','Noisecut')
    if P12_folder:
        for tri,tr in enumerate(st_corrected):st_corrected[tri].Metrics=st_raw[tri].Metrics.copy()/tr.Metrics.copy()
    st_hold = st_raw.copy()+st_corrected.copy()
    del st_raw,st_corrected
    # ________________________________________________________________________________
    return st_hold
def save_tight(filename,fig=None,dpi=200,format=None,margins=True):
        # Saves figure to PDF with no margins. Do not modify
        # plt.gca().set_axis_off()
        # plt.subplots_adjust(top = 1, bottom = 0.0, right = 1, left = 0,hspace = 0.07, wspace = 0.03)
        if margins:plt.margins(0.1,0.1)
        # plt.gca().xaxis.set_major_locator(plt.NullLocator())
        # plt.gca().yaxis.set_major_locator(plt.NullLocator())
        if fig is None:
                plt.savefig(str(filename), bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
                # print('Complete')
        else:
                fig.savefig(str(filename), bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
                # print('Complete')
        return plt.gcf()
def write_pickle(file,var):
    import pickle
    with open(str(file), 'wb') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # print('Saved to :' + str(file))
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        b = pickle.load(handle)
    return b
def get_noise(dirs,stanm):return load_pickle(list((dirs.SpectraAvg/stanm).glob('*sta.pkl'))[0])
def get_event_traces(stanm,evdir,ev,tf='sta.ZP-21.HZ.SAC'):
    stafold = evdir / stanm;label=tf.replace('sta.','').replace('.SAC','')
    raw = Stream([ObsQA.TOOLS.io.load_sac(stafold / (ev.Name + '.'+c+'.SAC'),rmresp=True) for c in ['H1','H2','HDH','HZ']])
    clear_output(wait=False)
    for i in range(len(raw)):raw[i].stats.location = 'Raw'
    corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
    corrected+=ObsQA.TOOLS.io.load_sac(stafold /  'CORRECTED' / '.'.join([stanm,ev.Name,tf]),rmresp=False)
    for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
# # ======================================================================================================================================================
def hps_corrected_list(stanm,catalog):
    sta=catalog.loc[stanm]
    events=np.array([e.name.split(stanm+'.')[1].strip('.HZ.SAC') for e in list((dirs.Events_HPS/stanm/'CORRECTED').glob('*.HZ.SAC'))])
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    events =events[ib]
    events=events[np.where(np.array([np.sum([len(list((dirs.Events_HPS/stanm/'CORRECTED').glob(f'*{e}*{c}*.SAC'))) for c in ['Z','1','2','HDH']]) for e in events])==4)]
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    events=Catalog([sta.Events[i] for i in ia])
    return events
def unravel(lst):return list(itertools.chain.from_iterable(lst))

def get_arrivals(sta_llaz=None,ev_llaz=None,zarc=None,model = 'iasp91',phases=('ttall',)):
        ''''
        Simple function pulls arrival times and phase names for a given event observed at a given station
        sta_llaz = List object containing [Lat,Lon] of the station
        ev_llaz = List object containing [Lat,Lon] of the event
        z must be in kilometers.
        phases = Tuple object containing a list of all desired phases. 'ttall'
        (default) will give every single phase available.
        -Charles Hoots,2022
        '''
        if zarc:z,arc=zarc
        else:z,arc=(ev_llaz[2],obspy.taup.taup_geo.calc_dist(ev_llaz[0],ev_llaz[1],sta_llaz[0],sta_llaz[1],6371,0))
        arrivals = obspy.taup.tau.TauPyModel(model=model).get_travel_times(z,arc,phase_list=phases)
        times = [[a.name,a.time] for a in arrivals]
        return times

def get_event_list(stanm,evdir,evmeta,tf = 'sta.ZP-21',mirror_fold=None):
    # sta,evdir,tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    stafold = evdir / stanm
    events = [f.name.replace(stanm + '.','').replace('.'+tf,'') for f in list((stafold /  'CORRECTED').glob('*.'+tf))]
    if mirror_fold:
        mirror = np.unique([g[:g.find('.Z')] for g in 
        [f.name.replace(stanm + '.','').replace('.sta','').replace('.day','').replace('.SA','').replace('.'+'HZ','') 
        for f in list((mirror_fold/stanm/'CORRECTED').glob('*'+'HZ'+'.'+'SAC'))]])
        events = list(np.intersect1d(events,mirror))
    eind = [[np.intersect1d(events,[e.Name for e in evmeta],return_indices=True)[1]][0],[np.intersect1d(events,[e.Name for e in evmeta],return_indices=True)[2]][0]]
    evmeta = Catalog([evmeta[e] for e in eind[1]]);events = [events[e] for e in eind[0]]
    return evmeta

def load_and_trim(hpsev_fold,stanm,ev,tlen,additional_chans,atacrevfold,label,tf=''):
    raw = Stream([ObsQA.TOOLS.io.load_sac(hpsev_fold /'rmresp'/stanm/(ev.Name + '.'+c+'.SAC'),rmresp=False) for c in ['HZ']])
    raw.trim(ev.origins[0].time,ev.origins[0].time+tlen,fill_value=0)
    clear_output(wait=False)
    raw.taper(.001)
    if len(additional_chans)>0:
        raw2 = Stream([ObsQA.TOOLS.io.load_sac(atacrevfold / 'rmresp'/stanm/(f'{ev.Name}.{c}.SAC'),rmresp=False) for c in additional_chans])
        clear_output(wait=False)
        raw2.trim(ev.origins[0].time,ev.origins[0].time+tlen,fill_value=0)
        trim_discrepancy = min([len(raw[0].data),min([len(r.data) for r in raw2])])
        for r in raw2:r.data=r.data[:trim_discrepancy]
        raw+=raw2
    for i in range(len(raw)):raw[i].stats.location='Raw'
    corrected=Stream(ObsQA.TOOLS.io.load_sac(hpsev_fold /  'corrected' / stanm/'.'.join([stanm,ev.Name,tf]),rmresp=False))
    corrected+= Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
    clear_output(wait=False)
    corrected.trim(raw[0].stats.starttime,raw[0].stats.endtime,fill_value=0)
    corrected.taper(.001)
    trim_discrepancy = min([min([len(ii.data) for ii in raw]),min([len(i.data) for i in corrected])])
    for c in corrected:c.data=c.data[:trim_discrepancy]
    for r in raw:r.data=r.data[:trim_discrepancy]
    for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
    # ++++++++++++++++++++++++++++++++++++++++++
    return corrected,raw

def add_metrics(corrected,raw):
    RtrZ=raw.select(channel='*Z')[0].copy()
    CtrZ=corrected.select(channel='*Z')[0].copy()
    RtrZ.Metrics=OBSM.Metrics(tr1=raw.select(channel='*1')[0].copy(),tr2=raw.select(channel='*2')[0].copy(),trP=raw.select(channel='*H')[0].copy(),trZ=raw.select(channel='*Z')[0].copy())
    CtrZ.Metrics = OBSM.Metrics(tr1=corrected.select(channel='*1')[0].copy(),tr2=corrected.select(channel='*2')[0].copy(),trP=corrected.select(channel='*H')[0].copy(),trZ=corrected.select(channel='*Z')[0].copy())
    CtrZ.Metrics = CtrZ.Metrics / RtrZ.Metrics
    RtrZ.Metrics = RtrZ.Metrics / CtrZ.Metrics
    return CtrZ,RtrZ



def reshapedict(c):
    n = np.array(list(c.keys()));n=n[n!='f']
    sn = unravel([[f'{ni[1:]}.{s}' for s in c[ni].keys()] for ni in n])
    r = AttribDict({f'{s}':c[f'n{s.split('.')[0]}'][s.split('.')[1]] for s in sn})
    r.f = c.f
    return r

def get_reports(sta=None,
    Archive=None,
    methods=['ATaCR','NoiseCut'],evna=None):
    methods = np.array(methods)
    Report=AttribDict()
    if Archive==None:dirs=dir_libraries();Archive=dirs.Data
    AnalysisFolder = Archive/'Analysis'/'Coherence'
    lp = lambda f:reshapedict(load_pickle(f))
    atacrfile = lambda:AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{tf}.zz.coherence.pkl'
    hpsfile = lambda :AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{c}.coherence.pkl'
    unionevents = lambda a:reduce(np.intersect1d,a)
    methods = [m.replace('NoiseCut','HPS') for m in methods]
    for m in methods:
        if m.lower()=='atacr':
            tfs = ['zp_21']
            d = AttribDict()
            for tf in tfs:d[tf]= lp(atacrfile())
            Report[m] = AttribDict()
            for tf in tfs:
                Report[m][tf] = AttribDict()
                stas = list(d[tf].keys()) if sta==None else [sta]
                for s in stas:
                    Report[m][tf][s] = d[tf][s]
                    if s=='f':continue
                    evs=unionevents([Report[m][tf][s].events for i in tfs])
                    ind = np.isin(Report[m][tf][s].events,evs)
                    Report[m][tf][s].coh=Report[m][tf][s].coh[ind,:]
                    Report[m][tf][s].events=np.array(Report[m][tf][s].events)[ind]
                    Report.f = d[tf].f
        elif m.lower()=='hps':
            Report[m]=AttribDict()
            d = AttribDict()
            comps = ['zz','11','22']
            for c in comps:d[c]= lp(hpsfile())
            for c in comps:Report[m][c]=AttribDict()
            Report.f = d[c].f
            stas = list(d[c].keys()) if sta==None else [sta]
            for c in comps:
                for s in stas:Report[m][c][s] = d[c][s]
            for s in stas:
                if s=='f':continue
                evs=unionevents([Report[m][c][s].events for c in comps])
                for c in comps:
                    ind = np.isin(Report[m][c][s].events,evs)
                    Report[m][c][s].coh=Report[m][c][s].coh[ind,:]
                    Report[m][c][s].events=np.array(Report[m][c][s].events)[ind]
    for m in methods:
        c = list(Report[m].keys())[0]
        if m=='ATaCR':
            # [del(Report['ATaCR'][k]) for k in list(Report['ATaCR'].keys()) if not (k=='zp_21')]
            for k in list(Report['ATaCR'].keys()):
                if not (k=='zp_21'):
                    del Report['ATaCR'][k]
        for s in stas:
            if s=='f':continue
            evs=unionevents([Report[m][list(Report[m].keys())[0]][s].events for m in methods])
            ind = np.isin(Report[m][c][s].events,evs)
            ikeys=list(Report[m].keys()) # if m=='HPS' else [c]
            for i in ikeys:
                Report[m][i][s].coh = Report[m][i][s].coh[ind,:]
                Report[m][i][s].events = np.array(Report[m][i][s].events)[ind]
    if len(stas)==1:
        for m in methods:
            for k in Report[m].keys():
                c = Report[m][k][stas[0]].copy()
                Report[m][k] = AttribDict()
                Report[m][k].coh = c.coh
                Report[m][k].events = c.events
    
    if not (evna==None):
        Report  = local_tools.cat.evcoh(Report,evna)
    return Report

def basic_preproc(st,tlen=7200,percent=None,lowpass=None):
        seconds = lambda A:A.hour*3600+A.minute*60+A.second
        tolerance=1 # seconds
        tend = np.min([tr.stats.endtime for tr in st])
        st.trim(tend-tlen,tend)
        trim_n = min([tr.data.shape[0] for tr in st])
        for tr in st:tr.data=tr.data[:trim_n]
        assert np.all([np.allclose(seconds(st[0].stats.endtime),seconds(tr.stats.endtime),tolerance) for tr in st]),'End times do not match'
        assert np.all([np.allclose(seconds(st[0].stats.starttime),seconds(tr.stats.starttime),tolerance) for tr in st]),'Start times do not match'
        st.detrend('linear');st.detrend('demean')
        if not percent==None:st.taper(percent)
        if not lowpass==None:st.filter('lowpass',freq=lowpass,zerophase=True)
        return st

def get_traces(stanm,event,channel='HZ',tf='sta.ZP-21',methods=['Original','NoiseCut','ATaCR'],preproc=True):
        if type(event)==obspy.core.event.event.Event:event=event.Name
        dirs=dir_libraries()
        rawdir=dirs.Events/'rmresp'/stanm
        atacr_correcteddir=dirs.Events/'corrected'/stanm
        hps_correcteddir=dirs.Events_HPS/'corrected'/stanm
        files=[]
        for m in methods:
                if m=='Original':rawfile=rawdir/f'{event}.{channel}.SAC';files.append(rawfile)
                elif m=='NoiseCut':hpsfile=hps_correcteddir/f'{stanm}.{event}.{channel}.SAC';files.append(hpsfile)
                elif m=='ATaCR':atacrfile=atacr_correcteddir/f'{stanm}.{event}.{tf}.{channel}.SAC';files.append(atacrfile)
        
        allexist=np.all([f.exists() for f in files])
        if allexist:
                st=Stream([ObsQA.TOOLS.io.load_sac(f) for f in files])
                for s,m in zip(st,methods):s.stats.location=m
                if preproc:st=basic_preproc(st)
                return st
        else:print(f'Not all files exist : {stanm} | {event}');return None




def get_station_events_hps(stanm,evdir,type='stream',tf='',evmeta=None,additional_chans=['H1','H2','HDH'],tlen=7200):
    # mirror_fold=None
    # sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    # --------------------------------------------- Event list
    # if not evmeta:evmeta = get_event_list(sta=sta,evdir=evdir[0],tf=tf,mirror_fold=mirror_fold)
    # --------------------------------------------- Load and store data
    label=tf.replace('sta.','').replace('.SAC','')
    if type=='stream':
        # For streams
        stafold = evdir[0]
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            if not (stafold / 'rmresp' / stanm / (ev.Name + '.HZ.SAC')).exists():
                continue
            hpsev_fold = evdir[0];atacrevfold = evdir[1]
            corrected,raw = load_and_trim(hpsev_fold,stanm,ev,tlen,additional_chans,atacrevfold,label,tf)
            st_hold+=raw;st_hold+=corrected
        st_hold.filter('bandpass',freqmin=1/1000,freqmax=st_hold[0].stats.sampling_rate/2,zerophase=True,corners=4)
    elif type=='metrics':
        # For metrics
        st_hold = Stream()
        hpsev_fold = evdir[0];atacrevfold = evdir[1]
        Noise = load_pickle(list((atacrevfold.parent/'AVG_STA'/stanm).glob('*sta.pkl'))[0])
        for evi,ev in enumerate(evmeta):
            if not (hpsev_fold /'rmresp'/stanm/(ev.Name +'.HZ.SAC')).exists():
                continue
            corrected,raw = load_and_trim(hpsev_fold,stanm,ev,tlen,additional_chans,atacrevfold,label,tf)
            assert len(corrected)==len(raw)
            if not (len(corrected)==4) or not (len(raw)==4):print('Data missing');continue
            corrected,raw = add_metrics(corrected,raw)
            corrected.stats.location = 'NoiseCut'
            corrected.Noise=Noise
            raw.Noise=Noise
            st_hold+=raw.copy();st_hold+=corrected.copy()
        st_hold.Noise=Noise
    return st_hold,evmeta


def get_station_events(stanm,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',evmeta=None):
    # mirror_fold=None
    # sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    stafold = evdir
    # --------------------------------------------- Event list
    # if evmeta is None:evmeta = get_event_list(stanm=stanm,evdir=evdir,tf=tf,mirror_fold=mirror_fold)
    # --------------------------------------------- Load and store data
    label=tf.replace('sta.','').replace('.SAC','')
    # For streams
    if type=='stream':
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            raw = Stream(ObsQA.TOOLS.io.load_sac(stafold / 'rmresp' / stanm/(ev.Name + '.HZ.SAC'),rmresp=False))
            clear_output(wait=False)
            corrected = Stream(ObsQA.TOOLS.io.load_sac(stafold /'corrected'/stanm/ '.'.join([stanm,ev.Name,tf]),rmresp=False))
            clear_output(wait=False)
            raw[0].stats.location = 'Raw'
            corrected[0].stats.location = 'ATaCR.'+label
            st = raw+corrected
            st.taper(0.001).filter('bandpass',freqmin=1/100,freqmax=1,zerophase=True,corners=4)
            st_hold = st_hold + st
    # For metrics
    elif type=='metrics':
        st_hold = Stream()
        Noise = load_pickle(list((evdir.parent/'AVG_STA'/stanm).glob('*sta.pkl'))[0])
        for evi,ev in enumerate(evmeta):
            raw = Stream([ObsQA.TOOLS.io.load_sac(stafold /'rmresp'/stanm/(ev.Name + '.'+c+'.SAC'),rmresp=False) for c in ['H1','H2','HDH','HZ']])
            clear_output(wait=False)
            for i in range(len(raw)):raw[i].stats.location = 'Raw'
            corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
            corrected+=ObsQA.TOOLS.io.load_sac(stafold /'corrected'/stanm/ '.'.join([stanm,ev.Name,tf]),rmresp=False)
            clear_output(wait=False)
            if not (len(corrected)==4) or not (len(raw)==4):print('Data missing');continue
            corrected,raw = add_metrics(corrected,raw)
            corrected.Noise=Noise;raw.Noise=Noise
            corrected.stats.location = 'ATaCR.'+label
            st_hold+=raw;st_hold+=corrected
        st_hold.Noise=Noise
    return st_hold,evmeta

