# from imports import *
from modules import *
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
            tr=load_sac(f,rmresp=True)
            tr.write(str(dest))
            del tr
            print(f'{state} || Updating.')
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
def AuditEventFolder(cat,eventsfolder,parseby='*Z.SAC',Minmag=6.0,Maxmag=8.0,nsta=0,search=True):
    if search:
        client = Client()
        for ista,Station in enumerate(cat.iloc):
            stanm=Station.StaName
            stafolder = Path(eventsfolder) / stanm
            files = list(stafolder.glob(parseby))
            evna = [f.name for f in files]
            evna = [f.split(f'{stanm}.')[-1] for f in evna]
            evna = ['.'.join(f.split('.')[:4]) for f in evna]
            assert len(evna)==len(np.unique(evna))
            print(f'{ista+1}/{len(cat)} | {stanm} | {len(evna)} events found. Collecting metadata from IRIS..')
            events=[]
            for ev in evna:
                timedelta = 60
                start = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
                end = start + timedelta
                event = client.get_events(starttime=start, endtime=end,minmagnitude=Minmag, maxmagnitude=Maxmag,orderby='magnitude')[0]
                event.Name = ev
                events.append(event)
            events = Catalog(events)
            Station.Events = None
            Station.Events = events
    #--Event names parsed by station row
    evnames_persta=[[e.Name for e in Station.Events] for Station in cat.iloc]
    # --Associate with other stations that share event names.
    for ista,Station in enumerate(cat.iloc):
        evna = [e.Name for e in Station.Events]
        for Event in Station.Events:
            Event.Stations = list(np.array(cat.StaName[[i for i,sev in enumerate(evnames_persta) if np.isin(Event.Name,sev)]]))
    events=unravel_cat(cat)
    # --Cull events with low station density (>=nsta).
    events=[events[i] for i in np.where(np.array([len(e.Stations) for e in events])>=nsta)[0]]
    #--All remaining event names for the full catalog
    evnames_overall = [e.Name for e in events]
    #--Re-assert catalog for each station following the nsta density cull
    for ista,Station in enumerate(cat.iloc):Station.Events=Catalog([e for e in Station.Events if np.isin(e.Name,evnames_overall)])
    return cat
def update_inventory(cat):
    for s in cat.iloc:s.Inventory=read_inventory(dirs.Events/'raw'/f'{s.StaName}'/f'{s.StaName}_inventory.xml')
    for inv,s in zip(cat.Inventory,cat.iloc):s.Inventory=Inventory([inv[i] for i in np.unique([n.stations[0].channels[0]._code for n in inv],return_index=True)[1]])
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