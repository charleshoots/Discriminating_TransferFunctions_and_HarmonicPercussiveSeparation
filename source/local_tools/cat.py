# from imports import *
from modules import *
from local_tools.ObsQA.TOOLS import io
import local_tools
from local_tools import ObsQA
from local_tools.io import *
from local_tools.math import *
# from quick_class import *
# import warnings
# import fnmatch
# from obspy.core.inventory.inventory import read_inventory
# import operator
# from obspy import read
# from obspy.geodetics import locations2degrees
# # ========================================================================================================================================================
# from IPython.display import clear_output
# from matplotlib.collections import PatchCollection
# from matplotlib.patches import Rectangle
# from obspy import Stream
# from modules import *
IRIS_EVENT_LINK = lambda s: f'https://ds.iris.edu/ds/nodes/dmc/tools/event/{s.resource_id.id.split('=')[-1]}'
def cat_parser(cat,k,wins):
    sta_sets = AttribDict()
    k=k.replace('StaDepth','Water_Depth_m')
    if k=='event_magnitude':
        mags = [np.array([e.magnitudes[0].mag for e in s.Events]) for s in cat.iloc]
        for si,s in enumerate(cat.iloc):
            sta_sets[si]=AttribDict()
            for wi in range(len(wins)-1):
                sta_sets[si][wi]=AttribDict()
                w=[wins[wi],wins[wi+1]]
                sta_sets[si][wi]=(mags[si]>=w[0])&(mags[si]<w[1])
    else:
        sta_sets=AttribDict()
        for si,s in enumerate(cat.iloc):
            sta_sets[si]=AttribDict()
            if not np.isin(k,list(s.Deployment.keys())):sta_sets[si]=np.nan;continue
            sta_sets[si]=s.Deployment[k]
        if isinstance(wins[0], numbers.Number):
            for si in range(len(cat)):
                test=[]
                for wi in range(len(wins)-1):
                    w = [wins[wi],wins[wi+1]]
                    test.append((sta_sets[si]>=w[0])&(sta_sets[si]<w[1]))
                sta_sets[si]=test
        else:
            for si in range(len(cat)):
                test=[]
                for w in wins:test.append(sta_sets[si]==w)
                sta_sets[si]=test
    return sta_sets
def parsed_coh(cat,i,wins,notched=True,M='ATaCR',C='zp_21'):
    f=cat.iloc[0].Data.Coherence().f
    cohs=[]
    if isinstance(wins[0],str):nwins=len(wins)
    else:nwins=len(wins)-1
    for wi in range(nwins):
        cw=[]
        for si,s in enumerate(cat.iloc):
            if np.sum(i[si][wi])>0:
                if np.size(i[si][wi])==1:ind=np.repeat(i[si][wi],len(s.Events))
                else:ind=i[si][wi]
                if C=='H':
                    ci=np.mean([(np.atleast_2d((s.Data.Coherence()[M][comp].coh[ind,:].squeeze()))[:,(f<=fnotch(s.StaDepth)) if notched else (f>-1)]).mean(axis=0).mean(axis=-1) for comp in ['11','22']])
                else:
                    ci=(s.Data.Coherence()[M][C].coh[ind,:][:,(f<=fnotch(s.StaDepth)) if notched else (f>-1)]).mean(axis=0).mean(axis=-1)
            else: ci=np.nan
            cw.append(ci)
        cohs.append(cw)
    cohs = np.array(cohs).T
    return cohs
def evcoh(coh,ev):
    c=coh.copy()
    c.ATaCR.zp_21.coh=np.atleast_2d(c.ATaCR.zp_21.coh[c.ATaCR.zp_21.events==ev,:])
    c.ATaCR.zp_21.events=np.array(c.ATaCR.zp_21.events[c.ATaCR.zp_21.events==ev])
    c.HPS.zz.coh=np.atleast_2d(c.HPS.zz.coh[c.HPS.zz.events==ev,:])
    c.HPS.zz.events=np.array(c.HPS.zz.events[c.HPS.zz.events==ev])
    c.HPS['11'].coh=np.atleast_2d(c.HPS['11'].coh[c.HPS['11'].events==ev,:])
    c.HPS['11'].events=np.array(c.HPS['11'].events[c.HPS['11'].events==ev])
    c.HPS['22'].coh=np.atleast_2d(c.HPS['22'].coh[c.HPS['22'].events==ev,:])
    c.HPS['22'].events=np.array(c.HPS['22'].events[c.HPS['22'].events==ev])
    return c


def cat_to_evcat(cat):
    icat = cat.copy()
    keys=['Distance_from_Land_km',
    'Distance_to_Plate_Boundary_km',
    'Surface_Current_ms',
    'Crustal_Age_Myr',
    'Deployment_Length_days',
    'Seismometer']
    evskeys=['Name','StaName', 'Station', 'Network', 'Data','Origin', 'LaLo', 'Link', 'Distance', 'Magnitude',
    'Stations', 'Traces', 'Latitude', 'Longitude', 'Experiment', 'Environment',
    'Pressure_Gauge', 'StaDepth', 'Start', 'End','NoiseAverage',
    'Seismometer', 'Sediment_Thickness_m', 'Instrument_Design', 'Distance_from_Land_km', 'Distance_to_Plate_Boundary_km',
    'Surface_Current_ms', 'Crustal_Age_Myr', 'Deployment_Length_days']
    catevs = []
    for si,s in enumerate(icat.iloc):
        s = s.copy()
        for k in keys:s[k]=s.Deployment[k]
        events = s.Events.copy()
        for ei,ev in enumerate(events):
            snm=s.StaName
            evs = s.copy()
            IRIS_EVENT_LINK = lambda ev=ev: f'https://ds.iris.edu/ds/nodes/dmc/tools/event/{str(ev.resource_id).split('=')[-1]}'
            TF = evs.Data.TF;Noise = evs.Data.Noise
            Coherence = lambda snm=snm,evna=ev.Name:local_tools.io.get_reports(snm,evna=evna)
            Traces = lambda event=ev.Name,snm=snm:local_tools.io.get_traces(snm,event)
            evs['Name'] = ev.Name
            evs['Origin'] = ev.origins[0].time
            evs['LaLo'] = [ev.origins[0].latitude,ev.origins[0].longitude]
            evs['Link'] = IRIS_EVENT_LINK
            evs['Distance'] = distance(s,ev)
            evs['Magnitude'] = ev.magnitudes[0].mag
            evs['Stations'] = ev.Stations
            evs['Traces']=lambda:evs.Data.Traces(ev)
            evs['Data']=AttribDict({'TF':TF,
            'Noise':Noise,
            'Traces':Traces,
            'Coherence':Coherence})
            evs = evs[evskeys]
            catevs.append(evs.to_frame().T)
    EventsCat = pd.concat(catevs)
    EventsCat.set_index('Name', inplace=True,drop=False)
    return EventsCat


def RobustParityTests(cat):
    dirs = io.dir_libraries()
    D = dict()
    D['NumEv']=[];D['N_P_Z']=[];D['N_P_All']=[];D['N_P_Horiz']=[];D['Parity_Z']=[];D['Parity_All']=[];D['Parity_H_Horiz']=[];D['Events']=[]
    atacrfile = lambda c='*':f'{s.SN}.{evn}.sta.{tf}.{c}.SAC'
    hpsfile = lambda c='*':f'{s.SN}.{evn}.{c}.SAC'
    hpsdir = dirs.Events_HPS/'corrected'
    atacrdir = dirs.Events/'corrected'
    atacrchans = ['HZ']
    hpschans = ['H1','H2','HZ']
    tf = 'ZP-21'
    # tfs = ['ZP-21','ZP','Z2-1']
    s = cat.iloc[0]
    for s in cat.iloc:
        parity_Z={};parity_All={};parity_H_Horiz={}
        for ev in s.Events:
            evn = ev.Name
            a_test_all = np.array([(atacrdir/s.SN/atacrfile(c)).exists() for c in atacrchans]).sum()==1
            h_test_all = np.array([(hpsdir/s.SN/hpsfile(c)).exists() for c in ['H1','H2','HZ']]).sum()==3
            h_test_horiz = np.array([(hpsdir/s.SN/hpsfile(c)).exists() for c in ['H1','H2']]).sum()==2
            a_test_Z = (atacrdir/s.SN/atacrfile('HZ')).exists()
            h_test_Z = (hpsdir/s.SN/hpsfile('HZ')).exists()
            parity_Z[evn]=(h_test_Z & a_test_Z)
            parity_All[evn]=(h_test_all & a_test_all)
            parity_H_Horiz[evn]=(h_test_horiz)
        # D['SN'].append(s.SN)
        D['NumEv'].append(len(s.Events));D['Events'].append(s.Events)
        D['N_P_Z'].append(len(parity_Z))
        D['N_P_All'].append(len(parity_All))
        D['N_P_Horiz'].append(len(parity_H_Horiz))
        D['Parity_Z'].append((parity_Z))
        D['Parity_All'].append((parity_All))
        D['Parity_H_Horiz'].append(parity_H_Horiz)
    D = pd.DataFrame(D,index=cat.SN)
    PerfectMethodParity = (sum((D.N_P_Z/D.NumEv)==1)==len(cat))
    PerfectHorizParity = (sum((D.N_P_Horiz/D.NumEv)==1)==len(cat))
    PerfectParity =  PerfectMethodParity & PerfectHorizParity
    if PerfectParity:print(f'| Perfect parity currently exists | \nPerfect parity across both methods for Z and all of HPS for Horiz >> {str(PerfectParity).upper()}')
    else:
        print('Parity is not perfect. Return report')
        return D
    if len(D[D.NumEv<20])>0:
        print('\n| Warnings | \nThe following stations event lists are below 20 by this many:\n')
        print(20-D[D.NumEv<20].NumEv)

def update_noise_quarantine(cat):
    dirs=dir_libraries()
    for stanm in cat.StaName:
        fold = dirs.SpectraAvg/stanm
        f=load_pickle(list((fold).glob('*.avg_sta.pkl'))[0])
        gooddays = f.day_files[f.gooddays]
        baddays = f.day_files[~f.gooddays]
        if len(baddays)>0:
            fold=dirs.Spectra/stanm;(fold/'_quarantine').mkdir(exist_ok=True)
            for day in baddays:
                if (fold/day).exists():shutil.move(fold/day, fold/'_quarantine'/day)
            fold=dirs.Noise/'raw'/stanm;(fold/'_quarantine').mkdir(exist_ok=True)
            for day in baddays:
                d=day.split('.spectra')[0]
                bad=list(fold.glob(f'{d}.*.SAC'))
                [shutil.move(r, fold/'_quarantine'/r.name) for r in bad]
            fold=dirs.Noise/'rmresp'/stanm;(fold/'_quarantine').mkdir(exist_ok=True)
            for day in baddays:
                d=day.split('.spectra')[0]
                bad=list(fold.glob(f'{d}.*.SAC'))
                [shutil.move(r, fold/'_quarantine'/r.name) for r in bad]
    Ndays = {stanm:len(list((dirs.Spectra/stanm).glob('*spectra.pkl'))) for stanm in cat.StaName}
    return Ndays
def shared_events(events,dirs,chan='HZ',tf='sta.ZP-21',minsta=10,stanm=None):
    atacrfold = dirs.Events/'corrected'
    hpsfold = dirs.Events_HPS/'corrected'
    methods=['atacr','hps']
    if not chan=='HZ':atacrfold=hpsfold;methods=['hps','hps']
    good=[]
    for evm in events.copy():
        ev=evm.copy()
        check_both = lambda ev,sta:np.all([(fold/sta/f'{sta}.{ev.Name}.{f'{tf}.' if method=='atacr' else ''}{chan}.SAC').exists() for method,fold in zip(methods,[atacrfold,hpsfold])])
        ev.Stations = np.array(ev.Stations)[np.where([check_both(ev,sta) for sta in ev.Stations])[0]]
        if stanm is not None:
            if not np.isin(stanm,ev.Stations):continue
        if len(ev.Stations)>=minsta:
            good.append(ev)
    ev_cat_shared = Catalog(good)
    return ev_cat_shared
# A nice little script for efficiently quarantining bad data.
def quarantine_data(folder,reg_list,qfold='_quarantine'):
    qfold=(folder/qfold);qfold.mkdir(exist_ok=True)
    for r in reg_list:[shutil.move(f,qfold/f.name) for f in list(folder.glob(f'*{r}*'))]

def partition_extra_days(cat,dirs):
    qfold = 'extra.days'
    folders = [dirs.Noise/'raw',dirs.Noise/'rmresp',dirs.Spectra]
    for si,sta in cat.iloc:
        noise=sta.Data.Noise.Averaged()
        days=noise.day_files
        if sum(noise.gooddays)<10:print(f'{sta.StaName}::{sum(noise.gooddays)}');continue
        good=days[noise.gooddays]
        bad = days[~noise.gooddays]
        goodwins_n=np.array([sum(g) for g in noise.day_goodwins[noise.gooddays]])
        good = good[np.flip(np.argsort(goodwins_n))]
        goodwins_n=goodwins_n[np.flip(np.argsort(goodwins_n))]
        keep = good[:10]
        extra = ['.'.join(e.split('.')[:2]) for e in good[10:]]
        if len(extra)>0:[lt.cat.quarantine_data(folder/sta.StaName,extra,qfold=qfold) for folder in folders]

def banner(msg):
    w=30
    print('**'*w)
    print('__'*w)
    nmsg=f'< {msg} >'
    print(f'{'|'*int((w-len(msg))/1)}{nmsg}{'|'*int((w-len(msg))/1)}')
    print('__'*w)
    print('**'*w)

# q = [i.replace('QC','').replace('.png','')]
# for i in q:quarantine_data(dirs.Noise/'raw'/'.'.join(i.split('.')[:2]),['.'.join(i.split('.')[-2:])])
# print('...All raw data quarantined')
# for i in q:quarantine_data(dirs.Spectra/'.'.join(i.split('.')[:2]),['.'.join(i.split('.')[-2:])])
# print('...All day spectra quarantined')
# for i in q:quarantine_data(dirs.TransferFunctions/'.'.join(i.split('.')[:2]),['.'.join(i.split('.')[-2:])])
# print('...All transfer functions quarantined')

def inquarantine(stanm,day,noisefold):return np.all([len(list((noisefold/g/stanm/'_quarantine').glob(day+'*.SAC')))>0 for g in ['rmresp','raw']])

def update_noise_qurantine(dirs):
    badfolder = dirs.Spectra/'_QualityControls'/'AnalystBadDays'
    noisefolder = dirs.Noise
    files=list(badfolder.glob('*.png'))
    print(f'{len(files)} files to assert quarantine')
    sta=[];nq=[]
    for f in files:
        stanm='.'.join(f.name.split('.')[:2])
        sta.append(stanm)
        # .zfill(3)
        day='.'.join([s.zfill(3) for s in f.name.split('.')[2:4]])
        # day=('.'.join(f.name.split('.')[2:4])).zfill(3)
        for nfold in [noisefolder/'rmresp'/stanm,noisefolder/'raw'/stanm]:
            (nfold/'_quarantine').mkdir(exist_ok=True)
            [shutil.move(fn,fn.parent/'_quarantine'/fn.name) for fn in list(nfold.glob(day + '*'))]
        nq.append(len(list((nfold/'_quarantine').glob('*.SAC')))/4)
    print(f'...done')
    print('Stations with quarantine: ')
    _=[print(f'{str(pi+1)} | {p} | days quarantined : {n}') for pi,(n,p) in enumerate(zip(nq,sta))]
def update_rmresp_folder(datafold,reg='*SAC',ovr=False):
    rawfold=Path(datafold)/'raw'
    (Path(datafold)/'rmresp').mkdir(exist_ok=True)
    files=list(rawfold.rglob(reg))
    files=[f for f in files if (not f.parent.name=='_depreciated')]
    files=[f for f in files if (not f.parent.name=='_quarantine')]
    files=[f for f in files if (f.parent.parent.name=='raw')]
    for fi,f in enumerate(files):
        state=f'{np.round(100*((fi+1)/len(files)),3)}% | {str(fi+1).zfill(len(str(len(files))))}/{len(files)}'
        dest=Path(str(f).replace('/raw/','/rmresp/'))
        dest.parent.mkdir(exist_ok=True,parents=True)
        # clear_output()
        if dest.exists():
            if not ovr:print(f'{state} || Already exists');continue
        try:
            tr=io.load_sac(f,rmresp=True)
            if len(tr.data)>0:
                tr.write(str(dest))
                del tr
                print(f'{state} || Updating.')
            else:print(f'{state} || Failed IO')
        except:
            print(f'{state} || Failed IO')


def addstations(nets,stas,cat,HJan23):
    nets,stas=np.array(nets),np.array(stas)
    rows=[]
    for net,sta in zip(nets,stas):
        rowcopy=cat.iloc[0].copy()
        rowcopy
        for i in rowcopy.keys():rowcopy[i]=None
        candidate=HJan23[(HJan23.Network.isin([net])) & (HJan23.Station.isin([sta]))]
        assert candidate.shape[0]==1
        candidate=candidate.iloc[0]
        rowcopy.StaName=candidate.Network+'.'+candidate.Station
        rowcopy.Station=candidate.Station
        rowcopy.Network=candidate.Network
        rowcopy.Latitude=candidate.Latitude
        rowcopy.Longitude=candidate.Longitude
        rowcopy.Experiment=candidate.Experiment
        rowcopy.Environment=candidate.Environment
        rowcopy.Pressure_Gauge=candidate.Pressure_Gauge
        rowcopy.StaDepth=candidate.Water_Depth_m
        rowcopy.Start=candidate.Start
        rowcopy.End=candidate.End
        rowcopy.Good_Channels = candidate[-4:].sum()==4
        rowcopy.Deployment=AttribDict(candidate.to_dict())
        rows.append(rowcopy.to_frame().T)
    output=pd.concat(rows).reset_index(drop=True)
    output.set_index('StaName', inplace=True,drop=False)
    return output

def AuditEventFolder(icat,srcat,eventsfolder,parseby='*HZ.SAC',Minmag=6.0,Maxmag=8.0,nsta=10):
    client = Client()
    files=list(eventsfolder.rglob(parseby))
    evna=np.unique([f.name.split('.HZ.SAC')[0] for f in files])
    evna=np.unique(['.'.join(e.split('.')) for e in evna if e in srcat.Name.unique()])
    events=[]
    for ei,ev in enumerate(evna):
        print(f'{ei+1}/{len(evna)} Collecting catalog from IRIS')
        timedelta = 60
        start = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
        end = start + timedelta
        args=dict(starttime=start,endtime=end,
        minmagnitude=Minmag, maxmagnitude=Maxmag,
        orderby='magnitude')
        event=client.get_events(**args)
        assert len(event)>0, 'no events found'
        event=event[0]
        event.Name=ev
        events.append(event)
    print('done.')
    events = Catalog(events)
    for ei,e in enumerate(events):
        print(f'{ei+1}/{len(evna)} Associating events with stations')
        name=e.Name;e.Stations=[]
        for s in icat.StaName:
            nf=len(list((eventsfolder/s).glob(f'{name}*')))
            if nf>0:e.Stations.append(s)
    nsta_events = Catalog([e for e in events if len(e.Stations)>=nsta])
    nrem=len(events)-len(nsta_events)
    events=nsta_events
    nsr=sum([len(e.Stations) for e in events])
    print('done.')
    print(f'{ei+1}/{len(evna)} Associating stations with events')
    for ista,s in enumerate(icat.iloc):s.Events=Catalog([e for e in events if s.StaName in e.Stations])
    print('done.')
    print(f'{nrem} events removed for having less than {nsta} stations.')
    print(f'Number of events in new catalog: {len(events)}')
    print(f'Number of source-receiver pairs in new catalog: {nsr}')
    return icat


def update_inventory(cat,return_list=False):
    dirs=lt.io.dir_libraries()
    inv=[read_inventory(dirs.Events_HPS/'raw'/f'{s.StaName}'/f'{s.StaName}_inventory.xml') for si,s in enumerate(cat.iloc)]
    if return_list:return inv
    else:
        cat['Inventory']=inv
        return cat

# def mirror(afold,bfold,events,comp='HZ'):
#     mirrored=[ev for ev in events if 
#     (len(list(afold.glob(f'*{ev.Name}*{comp}.SAC')))>0) 
#     and (len(list(bfold.glob(f'*{ev.Name}*{comp}.SAC')))>0)]
#     return Catalog(mirrored)


def mirror(stanm,catalog,dirs):
    # print(stanm)
    sta=catalog.loc[stanm].copy()
    ref_events=[e.Name for e in sta.Events]
    for fi,fold in enumerate([dirs.Events_HPS,dirs.Events]):
        events=np.array([e.name.replace('.HZ.SAC','') for e in list((fold/'rmresp'/stanm).glob('*.HZ.SAC'))])
        c,ia,ib=np.intersect1d(ref_events,events,return_indices=True)
        events=np.array(['.'.join(e.name.replace('.HZ.SAC','').replace(f'{stanm}.','').split('.')[:4]) for e in list((fold/'corrected'/stanm).glob('*.HZ.SAC'))])
        c,ia,ib=np.intersect1d(c,events,return_indices=True)

        # c,ia,ib=np.intersect1d(events_pre,c,return_indices=True)

        c,ia,ib=np.intersect1d(ref_events,c,return_indices=True)
        if fi==1:
            c,ia,ib=np.intersect1d([e.Name for e in sta.Events],c,return_indices=True)
            events=Catalog([sta.Events[i] for i in ia])
        else:
            events=c
        ref_events=events.copy()
    return events

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
def unravel(lst):return list(itertools.chain.from_iterable(lst))
def unravel_cat(cat):
    events=unravel([e for e in cat.Events])
    events=Catalog([events[i] for i in np.unique([e.Name for e in events],return_index=True)[-1]])
    return events