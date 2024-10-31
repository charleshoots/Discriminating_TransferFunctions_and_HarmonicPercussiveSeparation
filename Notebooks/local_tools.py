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
    if not r=='ZP':r=''.join(sorted(r))
    AA=avg.power.__dict__['c'+r[0]+r[0]]
    BB=avg.power.__dict__['c'+r[1]+r[1]]
    AB=avg.cross.__dict__['c'+r[0]+r[1]]
    x=avg.f
    y=Meters[m](ab=AB,aa=AA,bb=BB)
    return x[x>=0],y[x>=0]
Meters={'Coherence':_calc_coherence,'Admittance':_calc_admittance,'Phase':_calc_phase}
def get_event_list(sta,evdir,tf = 'sta.ZP-21',mirror_fold=None):
    # sta,evdir,tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    stafold = evdir / sta.StaName
    evmeta = sta.EventMeta
    events = [f.name.replace(sta.StaName + '.','').replace('.'+tf,'') for f in list((stafold /  'CORRECTED').glob('*.'+tf))]
    if mirror_fold:
        mirror = np.unique([g[:g.find('.Z')] for g in 
        [f.name.replace(sta.StaName + '.','').replace('.sta','').replace('.day','').replace('.SA','').replace('.'+'HZ','') 
        for f in list((mirror_fold/sta.StaName/'CORRECTED').glob('*'+'HZ'+'.'+'SAC'))]])
        events = list(np.intersect1d(events,mirror))
    eind = [[np.intersect1d(events,[e.Name for e in evmeta],return_indices=True)[1]][0],[np.intersect1d(events,[e.Name for e in evmeta],return_indices=True)[2]][0]]
    evmeta = Catalog([evmeta[e] for e in eind[1]]);events = [events[e] for e in eind[0]]
    return evmeta
def get_station_events_hps(sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None):
    # sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    # --------------------------------------------- Event list
    evmeta = get_event_list(sta=sta,evdir=evdir[0],tf=tf,mirror_fold=mirror_fold)
    # --------------------------------------------- Load and store data
    label=tf.replace('sta.','').replace('.SAC','')
    if type=='stream':
        # For streams
        stafold = evdir[0] / sta.StaName
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            raw = load_sac(stafold / (ev.Name + '.HZ.SAC'),rmresp=True)[0]
            corrected = load_sac(stafold /  'CORRECTED' / '.'.join([sta.StaName,ev.Name,tf]),rmresp=False)[0]
            raw[0].stats.location = 'Raw'
            corrected[0].stats.location = 'Corrected.'+label
            st = raw+corrected
            st.taper(0.001).filter('bandpass',freqmin=1/100,freqmax=1,zerophase=True,corners=4)
            st_hold = st_hold + st
    elif type=='metrics':
        # For metrics
        st_hold = Stream()
        hpsev_fold = evdir[0] / sta.StaName;atacrevfold = evdir[1] / sta.StaName
        for evi,ev in enumerate(evmeta):
            raw = Stream([load_sac(hpsev_fold / (ev.Name + '.'+c+'.SAC'),rmresp=True)[0][0] for c in ['HZ']])
            raw2 = Stream([load_sac(atacrevfold / (ev.Name + '.'+c+'.SAC'),rmresp=True)[0][0] for c in ['H1','H2','HDH']])
            raw.taper(.001)
            raw.trim(raw2[0].stats.starttime,raw2[0].stats.endtime,pad=True,fill_value=0);raw+=raw2
            tlen = raw[0].stats.endtime-raw[0].stats.starttime
            for i in range(len(raw)):raw[i].stats.location = 'Raw'
            corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
            corrected+=load_sac(hpsev_fold /  'CORRECTED' / '.'.join([sta.StaName,ev.Name,tf]),rmresp=False)[0][0]
            corrected.taper(.001)
            corrected.trim(raw[0].stats.starttime,raw[0].stats.endtime,pad=True,fill_value=0)
            for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
            if not (len(corrected)==4) or not (len(raw)==4):print('Data missing');continue
            RtrZ=raw.select(channel='*Z')[0].copy()
            RtrZ.Metrics=OBSM.Metrics(tr1=raw.select(channel='*1')[0].copy(),tr2=raw.select(channel='*2')[0].copy(),trP=raw.select(channel='*H')[0].copy(),trZ=raw.select(channel='*Z')[0].copy())
            CtrZ=corrected.select(channel='*Z')[0].copy()
            CtrZ.Metrics = OBSM.Metrics(tr1=corrected.select(channel='*1')[0].copy(),tr2=corrected.select(channel='*2')[0].copy(),trP=corrected.select(channel='*H')[0].copy(),trZ=corrected.select(channel='*Z')[0].copy())
            CtrZ.Metrics = CtrZ.Metrics / RtrZ.Metrics
            RtrZ.Metrics = RtrZ.Metrics / CtrZ.Metrics
            st_hold+=RtrZ;st_hold+=CtrZ
        st_hold.Noise = load_pickle(list((atacrevfold.parent.parent/'AVG_STA'/sta.StaName).glob('*sta.pkl'))[0])
    return st_hold,evmeta
def get_station_events(sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None):
    # sta,evdir,type='stream',tf = 'sta.ZP-21.HZ.SAC',mirror_fold=None
    stafold = evdir / sta.StaName
    # --------------------------------------------- Event list
    evmeta = get_event_list(sta=sta,evdir=evdir,tf=tf,mirror_fold=mirror_fold)
    # --------------------------------------------- Load and store data
    label=tf.replace('sta.','').replace('.SAC','')
    # For streams
    if type=='stream':
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            raw = load_sac(stafold / (ev.Name + '.HZ.SAC'),rmresp=True)[0]
            corrected = load_sac(stafold /  'CORRECTED' / '.'.join([sta.StaName,ev.Name,tf]),rmresp=False)[0]
            raw[0].stats.location = 'Raw'
            corrected[0].stats.location = 'Corrected.'+label
            st = raw+corrected
            st.taper(0.001).filter('bandpass',freqmin=1/100,freqmax=1,zerophase=True,corners=4)
            st_hold = st_hold + st
    # For metrics
    elif type=='metrics':
        st_hold = Stream()
        for evi,ev in enumerate(evmeta):
            raw = Stream([load_sac(stafold / (ev.Name + '.'+c+'.SAC'),rmresp=True)[0][0] for c in ['H1','H2','HDH','HZ']])
            for i in range(len(raw)):raw[i].stats.location = 'Raw'
            corrected = Stream([raw.select(channel=c)[0].copy() for c in ['*1','*2','*H']]).copy()
            corrected+=load_sac(stafold /  'CORRECTED' / '.'.join([sta.StaName,ev.Name,tf]),rmresp=False)[0][0]
            for i in range(len(corrected)):corrected[i].stats.location = 'Corrected.'+label
            if not (len(corrected)==4) or not (len(raw)==4):print('Data missing');continue
            RtrZ=raw.select(channel='*Z')[0].copy()
            RtrZ.Metrics=OBSM.Metrics(tr1=raw.select(channel='*1')[0].copy(),tr2=raw.select(channel='*2')[0].copy(),trP=raw.select(channel='*DH')[0].copy(),trZ=raw.select(channel='*Z')[0].copy())
            CtrZ=corrected.select(channel='*Z')[0].copy()
            CtrZ.Metrics = OBSM.Metrics(tr1=corrected.select(channel='*1')[0].copy(),tr2=corrected.select(channel='*2')[0].copy(),trP=corrected.select(channel='*DH')[0].copy(),trZ=corrected.select(channel='*Z')[0].copy())
            CtrZ.Metrics = CtrZ.Metrics / RtrZ.Metrics
            RtrZ.Metrics = RtrZ.Metrics / CtrZ.Metrics
            st_hold+=RtrZ;st_hold+=CtrZ
        st_hold.Noise = load_pickle(list((evdir.parent/'AVG_STA'/sta.StaName).glob('*sta.pkl'))[0])
    return st_hold,evmeta

# # ======================================================================================================================================================
k = 1



# def get_noiseavg(stakey,noisedir,rebuild=False):
#     # -------Get good days
#     avg = load_pickle(list((noisedir.parent/'AVG_STA'/stakey).glob('*sta.pkl'))[0])
#     if rebuild:
#         gooddays = np.array(list((noisedir.parent/'Spectra'/stakey).glob('*.spectra.pkl')))[avg.gooddays]
#         days = np.array(list((noisedir / stakey).glob('*Z.SAC')))[avg.gooddays]
#         days = [d.name.replace('..','.').replace('.HZ.SAC','') for d in days]
#         overlap,window=load_pickle(gooddays[0]).overlap,load_pickle(gooddays[0]).window
#         goodwins = [load_pickle(g).goodwins for g in gooddays]
#         for w,d in zip(goodwins,days):
#             noise = Stream([load_sac(noisedir/stakey/''.join([d,'..',c,'.SAC']),rmresp=True)[0][0] for c in ['H1','H2','HZ','HDH']])
#             slides = [s for s in noise.slide(window_length=window,step=window*(1-overlap))]
#             slides = [s.split() for wi,s in enumerate(slides) if w[wi]]
#             noise = slides[0];st = Stream();chans=['*Z','*1','*2','*DH']
#             for s in slides:st+=s
#             avg=[np.mean(st.select(channel=c),axis=0) for c in chans]
#             for ci,c in enumerate(chans):noise.select(channel=c)[0].data=avg[ci]
#         return noise
#     else:return avg
# def get_noise_metrics(stakey,noisedir):
#     # -------Build noise
#     noise = get_noiseavg(stakey,noisedir)
#     # -------Build metrics
#     tr1=noise.select(channel='*1')[0];tr2=noise.select(channel='*2')[0]
#     trZ=noise.select(channel='*Z')[0];trP=noise.select(channel='*DH')[0]
#     NoiseM=OBSM.Metrics(tr1=tr1,tr2=tr2,trZ=trZ,trP=trP)
#     return NoiseM