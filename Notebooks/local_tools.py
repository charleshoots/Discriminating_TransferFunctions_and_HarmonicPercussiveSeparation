from imports import *
from quick_class import *
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
def _calc_phase(ab,**args):
        ph = np.angle(ab,deg=True)
        return ph
def _calc_admittance(ab,bb,**args):
        ad = np.abs(ab)/bb
        return ad
def _calc_coherence(ab,aa,bb,**args):
        coh = ((np.abs(ab)**2)/(aa*bb))
        return coh
def avg_meter(avg,m,r):
    Meters={'Coherence':_calc_coherence,'Admittance':_calc_admittance,'Phase':_calc_phase}
    if not r=='ZP':r=''.join(sorted(r))
    AA=avg.power.__dict__['c'+r[0]+r[0]]
    BB=avg.power.__dict__['c'+r[1]+r[1]]
    AB=avg.cross.__dict__['c'+r[0]+r[1]]
    x=avg.f
    y=Meters[m](ab=AB,aa=AA,bb=BB)
    return x[x>=0],y[x>=0]
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
def get_station_events_hps(stanm,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None,evmeta=None):
    # sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    # --------------------------------------------- Event list
    if not evmeta:evmeta = get_event_list(sta=sta,evdir=evdir[0],tf=tf,mirror_fold=mirror_fold)
    # --------------------------------------------- Load and store data
    label=tf.replace('sta.','').replace('.SAC','')
    if type=='stream':
        # For streams
        stafold = evdir[0] / stanm
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            raw = load_sac(stafold / (ev.Name + '.HZ.SAC'),rmresp=True)[0]
            clear_output(wait=False)
            corrected = load_sac(stafold /  'CORRECTED' / '.'.join([stanm,ev.Name,tf]),rmresp=False)[0]
            clear_output(wait=False)
            raw[0].stats.location = 'Raw'
            corrected[0].stats.location = 'Corrected.'+label
            st = raw+corrected
            st.taper(0.001).filter('bandpass',freqmin=1/1000,freqmax=st[0].stats.sampling_rate/2,zerophase=True,corners=4)
            st_hold = st_hold + st
    elif type=='metrics':
        # For metrics
        st_hold = Stream()
        hpsev_fold = evdir[0];atacrevfold = evdir[1]
        # Noise = load_pickle(list((atacrevfold.parent.parent/'AVG_STA'/stanm).glob('*sta.pkl'))[0])
        for evi,ev in enumerate(evmeta):
            raw = Stream([load_sac(hpsev_fold /'rmresp'/stanm/(ev.Name + '.'+c+'.SAC'),rmresp=False)[0][0] for c in ['HZ']])
            clear_output(wait=False)
            raw2 = Stream([load_sac(atacrevfold / 'rmresp'/stanm/(f'{ev.Name}.{c}.SAC'),rmresp=False)[0][0] for c in ['H1','H2','HDH']])
            clear_output(wait=False)
            raw.taper(.001)
            raw.trim(raw2[0].stats.starttime,raw2[0].stats.endtime,pad=True,fill_value=0);raw+=raw2
            tlen = raw[0].stats.endtime-raw[0].stats.starttime
            for i in range(len(raw)):raw[i].stats.location='Raw'
            corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
            corrected+=load_sac(hpsev_fold /  'corrected' / stanm/'.'.join([stanm,ev.Name,tf]),rmresp=False)[0][0]
            clear_output(wait=False)
            corrected.taper(.001)
            corrected.trim(raw[0].stats.starttime,raw[0].stats.endtime,pad=True,fill_value=0)
            for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
            if not (len(corrected)==4) or not (len(raw)==4):print('Data missing');continue
            RtrZ=raw.select(channel='*Z')[0].copy()
            CtrZ=corrected.select(channel='*Z')[0].copy()
            # RtrZ.Metrics=OBSM.Metrics(tr1=raw.select(channel='*1')[0].copy(),tr2=raw.select(channel='*2')[0].copy(),trP=raw.select(channel='*H')[0].copy(),trZ=raw.select(channel='*Z')[0].copy())
            # CtrZ.Metrics = OBSM.Metrics(tr1=corrected.select(channel='*1')[0].copy(),tr2=corrected.select(channel='*2')[0].copy(),trP=corrected.select(channel='*H')[0].copy(),trZ=corrected.select(channel='*Z')[0].copy())
            # CtrZ.Metrics = CtrZ.Metrics / RtrZ.Metrics
            # RtrZ.Metrics = RtrZ.Metrics / CtrZ.Metrics
            # CtrZ.Noise=Noise;RtrZ.Noise=Noise
            st_hold+=RtrZ;st_hold+=CtrZ
        # st_hold.Noise=Noise
    return st_hold,evmeta
def get_station_events(stanm,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None,evmeta=None):
    # sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    stafold = evdir
    # --------------------------------------------- Event list
    if not evmeta:evmeta = get_event_list(sta=stanm,evdir=evdir,tf=tf,mirror_fold=mirror_fold)
    # --------------------------------------------- Load and store data
    label=tf.replace('sta.','').replace('.SAC','')
    # For streams
    if type=='stream':
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            raw = load_sac(stafold / (ev.Name + '.HZ.SAC'),rmresp=True)[0]
            clear_output(wait=False)
            corrected = load_sac(stafold /'corrected'/stanm/ '.'.join([stanm,ev.Name,tf]),rmresp=False)[0]
            clear_output(wait=False)
            raw[0].stats.location = 'Raw'
            corrected[0].stats.location = 'Corrected.'+label
            st = raw+corrected
            st.taper(0.001).filter('bandpass',freqmin=1/100,freqmax=1,zerophase=True,corners=4)
            st_hold = st_hold + st
    # For metrics
    elif type=='metrics':
        st_hold = Stream()
        Noise = load_pickle(list((evdir.parent/'AVG_STA'/stanm).glob('*sta.pkl'))[0])
        for evi,ev in enumerate(evmeta):
            raw = Stream([load_sac(stafold /'rmresp'/stanm/(ev.Name + '.'+c+'.SAC'),rmresp=False)[0][0] for c in ['H1','H2','HDH','HZ']])
            clear_output(wait=False)
            for i in range(len(raw)):raw[i].stats.location = 'Raw'
            corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
            corrected+=load_sac(stafold /'corrected'/stanm/ '.'.join([stanm,ev.Name,tf]),rmresp=False)[0][0]
            for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
            clear_output(wait=False)
            if not (len(corrected)==4) or not (len(raw)==4):print('Data missing');continue
            RtrZ=raw.select(channel='*Z')[0].copy()
            RtrZ.Metrics=OBSM.Metrics(tr1=raw.select(channel='*1')[0].copy(),tr2=raw.select(channel='*2')[0].copy(),trP=raw.select(channel='*DH')[0].copy(),trZ=raw.select(channel='*Z')[0].copy())
            CtrZ=corrected.select(channel='*Z')[0].copy()
            CtrZ.Metrics = OBSM.Metrics(tr1=corrected.select(channel='*1')[0].copy(),tr2=corrected.select(channel='*2')[0].copy(),trP=corrected.select(channel='*DH')[0].copy(),trZ=corrected.select(channel='*Z')[0].copy())
            CtrZ.Metrics = CtrZ.Metrics / RtrZ.Metrics
            RtrZ.Metrics = RtrZ.Metrics / CtrZ.Metrics
            CtrZ.Noise=Noise;RtrZ.Noise=Noise
            st_hold+=RtrZ;st_hold+=CtrZ
        st_hold.Noise=Noise
    return st_hold,evmeta
def get_event_traces(stanm,evdir,ev,tf='sta.ZP-21.HZ.SAC'):
    stafold = evdir / stanm;label=tf.replace('sta.','').replace('.SAC','')
    raw = Stream([load_sac(stafold / (ev.Name + '.'+c+'.SAC'),rmresp=True)[0][0] for c in ['H1','H2','HDH','HZ']])
    clear_output(wait=False)
    for i in range(len(raw)):raw[i].stats.location = 'Raw'
    corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
    corrected+=load_sac(stafold /  'CORRECTED' / '.'.join([stanm,ev.Name,tf]),rmresp=False)[0][0]
    for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
# # ======================================================================================================================================================
def write_pickle(file,var):
    import pickle
    with open(str(file), 'wb') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Saved to :' + str(file))
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        b = pickle.load(handle)
    return b
def distance(sta,ev,unit='deg'):
    origins=ev.origins[0]
    stalla,evlla=[sta.Latitude,sta.Longitude],[origins.latitude,origins.longitude]
    dist=locations2degrees(stalla[0],stalla[1],evlla[0],evlla[1])
    if unit.lower()=='km':dist=degrees2kilometers(dist)
    return dist
def unravel(lst):return list(itertools.chain.from_iterable(lst))
def get_noise(dirs,stanm):return load_pickle(list((dirs.SpectraAvg/stanm).glob('*sta.pkl'))[0])

def mirror_events(reports):
    nkeys = [n for n in list(reports[0].__dict__.keys()) if not n=='f']
    mirror = dict()
    for ni,n in enumerate(nkeys):
        skeys = list(reports[0][n].__dict__.keys())
        for si,s in enumerate(skeys):
            stanm = '.'.join([n,s]).replace('n','')
            ev0 = [k.replace('.','') for k in reports[0][n][s].events]
            ev1 = [k.replace('.','') for k in reports[1][n][s].events]
            mirror[stanm] = np.intersect1d(ev0,ev1)
    return mirror
def filter_event_catalog(cat,min_sta=15):
    EvMeta = cat.iloc[0].Events
    for sta in cat.iloc:EvMeta+=sta.Events
    SharedEvents =Catalog([EvMeta[ei] for ei in np.unique([e.Name for e in EvMeta],return_index=True)[1]])
    n_sta_events = Catalog([s for s in SharedEvents if len(s.Stations)>=min_sta])
    n_sta_events = Catalog([n_sta_events[ei] for ei in np.flip(np.argsort([s.magnitudes[0].mag for s in n_sta_events]))])
    stas=[];np.array([stas.extend(n.Stations) for n in n_sta_events]);stas=np.unique(stas)
    stations_not_used = cat[~np.any([s==cat.StaName for s in stas],axis=0)]
    return n_sta_events,stations_not_used

def get_event_P12(stas,evn,fold):
    st_p12 = Stream()
    for c in ['DH','1','2']:
        st_p12+=Stream([load_sac(fold/s/(evn + ('.H'+c+'.SAC')),rmresp=True)[0][0] for s in stas])
        clear_output(wait=False)
    return st_p12
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
    st=Stream([load_sac(datafolder/s/subfold/('*'+f),rmresp=rmresp)[0][0] for f,s in zip(files,stas)]);clear_output(wait=False)
    if P12_folder:
        st_P12=get_event_P12(stas,event.Name,P12_folder);st+=st_P12;clear_output(wait=False)
    st.taper(0.0001);st.trim(event.origins[0].time,event.origins[0].time+tlen,pad=True,fill_value=trim_gap)
    clear_output(wait=False)
    return st

# -----
def get_avg_noise(stanm,fold):return load_pickle(list((fold/stanm).glob('*sta.pkl'))[0])
# -----
def detect_outscale(raw,correct,vertical_scale=1.2,ylim=None,suppress=False):
    # Occasionally the ampltide changes (ie noise reduction) after correction is so 
    # significant that it makes it challenging to plot the raw and corrected on the 
    # same plot.This function detects when the scale differences exceed vertical_scale 
    # (120% by default) of the maximum value in the corrected trace. When it occurs, 
    # the distance in amplitudes the raw is from the corrected is shrunk to within this margin.
    if not ylim:ylim = np.array([np.max(np.abs(c.data))*vertical_scale for c in correct])
    out_scaled = np.array([np.max(np.abs(c.data),
    where=~(np.isinf(c.data)+np.isnan(c.data)),
    initial=0) for c in raw]) > ylim
    if np.any(out_scaled):
        if not suppress:print('Large amplitude scale differences detected in '+raw[0].id)
        for tr_ind,tr in enumerate(raw):
            if out_scaled[tr_ind]:tr.data = (tr.data/np.max(np.abs(tr.data)))*ylim[tr_ind]
    return raw

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
    st=Stream([load_sac(datafolder/s/subfold/('*'+f),rmresp=rmresp)[0][0] for f,s in zip(files,stas)]);clear_output(wait=False)

    # Metrics will approximately septuple the memory pressure for the code
    if P12_folder:st = stream_metrics(st,event.Name,P12_folder)

    st.taper(0.0001);st.trim(event.origins[0].time,event.origins[0].time+tlen,pad=True,fill_value=trim_gap)
    clear_output(wait=False)
    return st
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
def smooth(d,k=10):return np.convolve(d, np.ones(k) / k, mode='same')
def mirror(afold,bfold,events,comp='HZ'):
    mirrored=[ev for ev in events if 
    (len(list(afold.glob(f'*{ev.Name}*{comp}.SAC')))>0) 
    and (len(list(bfold.glob(f'*{ev.Name}*{comp}.SAC')))>0)]
    return Catalog(mirrored)
def hps_corrected_list(stanm,catalog):
    sta=catalog.loc[stanm]
    events=np.array([e.name.split(stanm+'.')[1].strip('.HZ.SAC') for e in list((dirs.Events_HPS/stanm/'CORRECTED').glob('*.HZ.SAC'))])
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    events =events[ib]
    events=events[np.where(np.array([np.sum([len(list((dirs.Events_HPS/stanm/'CORRECTED').glob(f'*{e}*{c}*.SAC'))) for c in ['Z','1','2','HDH']]) for e in events])==4)]
    c,ia,ib=np.intersect1d([e.Name for e in sta.Events],events,return_indices=True)
    events=Catalog([sta.Events[i] for i in ia])
    return events