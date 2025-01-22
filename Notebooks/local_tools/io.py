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
def save_tight(filename,fig=None,dpi=200,format=None):
        # Saves figure to PDF with no margins. Do not modify
        # plt.gca().set_axis_off()
        # plt.subplots_adjust(top = 1, bottom = 0.0, right = 1, left = 0,hspace = 0.07, wspace = 0.03)
        plt.margins(0.1,0.1)
        # plt.gca().xaxis.set_major_locator(plt.NullLocator())
        # plt.gca().yaxis.set_major_locator(plt.NullLocator())
        if fig is None:
                plt.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
                print('Complete')
        else:
                fig.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
                print('Complete')
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
def get_noise(dirs,stanm):return load_pickle(list((dirs.SpectraAvg/stanm).glob('*sta.pkl'))[0])
def get_event_traces(stanm,evdir,ev,tf='sta.ZP-21.HZ.SAC'):
    stafold = evdir / stanm;label=tf.replace('sta.','').replace('.SAC','')
    raw = Stream([load_sac(stafold / (ev.Name + '.'+c+'.SAC'),rmresp=True)[0][0] for c in ['H1','H2','HDH','HZ']])
    clear_output(wait=False)
    for i in range(len(raw)):raw[i].stats.location = 'Raw'
    corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
    corrected+=load_sac(stafold /  'CORRECTED' / '.'.join([stanm,ev.Name,tf]),rmresp=False)[0][0]
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
