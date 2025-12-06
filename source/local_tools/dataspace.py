from modules import *
octavg = lt.math.octave_average
os.system('cls' if os.name == 'nt' else 'clear')
# -----
import local_tools
from local_tools.io import *
from local_tools.math import *
from local_tools.ObsQA.TOOLS.io import event_stream_arrivals
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# ----------------------------------------------------------------------------------------
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
def get_inv(snm):
    dirs=lt.io.dir_libraries()
    f=dirs.ATaCR/f'EVENTS/raw/{snm}/{snm}_inventory.xml'
    return read_inventory(f)
def import_HJan23(hjanfile):
    # os.chdir(project_path)
    dirs=lt.io.dir_libraries()
    HJan23=pd.read_pickle(hjanfile)
    HJan23['Good_Channels']=HJan23.T[-4:].sum().T==4
    return HJan23
def import_cat(catfile,hjanfile,aggregate=False):
    # os.chdir(project_path)
    dirs=lt.io.dir_libraries()
    # ____________________________________________________________________________________
    # |||||||||||||||||||||||||||||||||||| Catalogs ||||||||||||||||||||||||||||||||||||||
    catalog = pd.read_pickle(catfile)
    hjancat=import_HJan23(hjanfile)
    catalog['Good_Channels']=np.array([hjancat[hjancat.StaName==sn].iloc[0].Good_Channels for sn in catalog.StaName])
    catalog.Experiment.replace('CASCADIA KECK','KECK',inplace=True)
    catalog.Experiment.replace('CASCADIA INITIATIVE','Cascadia',inplace=True)
    catalog.Experiment.replace('ALBACORE','Albacore',inplace=True)
    catalog.Experiment.replace('MARIANA','Mariana',inplace=True)
    catalog.Experiment.replace('PAPUA','Papua',inplace=True)
    catalog.Experiment.replace('LAU','Lau',inplace=True)
    catalog.Experiment.replace('NOMELT','No Melt',inplace=True)
    catalog.Experiment.replace('BLANCO','Blanco',inplace=True)
    catalog.set_index('StaName', inplace=True,drop=False)
    clear_output(wait=False)
    a=[]
    dirs=lt.io.dir_libraries()
    for sta in catalog.iloc:
        snm=sta.StaName
        SpecAVG=lambda file=_get_file(snm,dirs.SpectraAvg,'*.avg_sta.pkl')[0]:load_pickle(file)
        dayfiles=lambda file=_get_file(snm,dirs.Spectra,'*.spectra.pkl'): SpecAVG(file).day_files[SpecAVG(file).gooddays]
        DaySpec=lambda fo=dirs.Spectra/snm,dayfiles=dayfiles: [load_pickle(fo/fi) for fi in dayfiles()]
        TF=lambda file=_get_file(snm,dirs.TransferFunctions,'*-*transfunc.pkl')[0]:load_pickle(file)
        Coherence = lambda snm=snm:get_reports(snm)
        Traces = lambda event,snm=snm,channel='HZ',tf='sta.ZP-21',methods=['Original','NoiseCut','ATaCR'],preproc=True:get_traces(snm,event,channel=channel,tf=tf,methods=methods,preproc=preproc)
        IRIS_EVENT_LINK = lambda s: f'https://ds.iris.edu/ds/nodes/dmc/tools/event/{s.resource_id.id.split('=')[-1]}'
        Print_IRIS_Event_Page = lambda e: f'{e.magnitudes[0].magnitude_type} {e.magnitudes[0].mag} | {e.Name}: {IRIS_EVENT_LINK(e)}'
        IRISPages = lambda evs=Catalog([sta.Events[ii] for ii in np.flip(np.argsort([i.magnitudes[0].mag for i in sta.Events]))]):[Print_IRIS_Event_Page(e) for e in evs]
        inv = lambda snm=sta.StaName:get_inv(snm)
        Data=AttribDict()
        Noise=AttribDict()
        Noise.Day=DaySpec
        Noise.Averaged=SpecAVG
        Data.Noise=Noise
        Data.TF=TF
        Data.Coherence=Coherence
        Data.Traces=Traces
        Data.IRIS_Event_Pages = IRISPages
        Data.Inventory=inv
        
        a.append(Data)

    catalog['Inventory']=local_tools.cat.update_inventory(catalog,return_list=True) if aggregate else None

    catalog['Data']=a
    catalog=catalog[['StaName', 'Station', 'Network','Data', 'Latitude', 'Longitude', 'Experiment',
    'Environment', 'Pressure_Gauge', 'StaDepth', 'Start', 'End', 'Events','Deployment', 'Inventory']]

    catalog['NoiseAverage'] = None if not aggregate else [(10*np.log10(s.Data.Noise.Averaged().power.cZZ[s.Data.Noise.Averaged().f>0 & (s.Data.Noise.Averaged().f<=fnotch(s.StaDepth))])).mean() for s in catalog.iloc]
    catalog['Seismometer']=[i['Seismometer'] for i in catalog.Deployment]
    catalog['Sediment_Thickness_m']=[i['Sediment_Thickness_m'] if np.isin('Sediment_Thickness_m',list(i.keys())) else np.nan for i in catalog.Deployment]
    catalog['Instrument_Design']=[i['Instrument_Design'] for i in catalog.Deployment]
    for s in catalog.Events:
        for e in s:e.mag = e.magnitudes[0].mag;e.depth=e.origins[0].depth;e.time = e.origins[0].time;e.lla = [e.origins[0].latitude,e.origins[0].longitude]
    d=catalog.Sediment_Thickness_m.copy();d[np.isnan(d)]=np.mean(catalog[catalog.Network=='7D'].Sediment_Thickness_m);catalog['Sediment_Thickness_m']=d
    for i in range(len(catalog)):catalog.iloc[i].Deployment.Network = catalog.iloc[i].Network
    catalog['SN'] = catalog.StaName
    # ------------------
    return catalog
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# ----------------------------------------------------------------------------------------
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
def _get_file(snm,fold=lt.io.dir_libraries().SpectraAvg,search='*.avg_sta.pkl'):
    files=np.array(list((fold/snm).glob(search)))
    return files if len(files)>0 else [None]
def rcat_to_srcat(cat):
    icat = cat.copy()
    keys=['Distance_from_Land_km',
    'Distance_to_Plate_Boundary_km',
    'Surface_Current_ms',
    'Crustal_Age_Myr',
    'Deployment_Length_days',
    'Seismometer']
    evskeys=['Name','StaName', 'Station','Event','Network', 'Data','Origin', 'LaLo', 'Link','Distance','Phases','Magnitude',
    'Stations', 'Traces', 'Latitude', 'Longitude','Experiment', 'Environment',
    'Pressure_Gauge','StaDepth','EvDepth','Start','End','NoiseAverage',
    'Seismometer', 'Sediment_Thickness_m', 'Instrument_Design', 'Distance_from_Land_km', 'Distance_to_Plate_Boundary_km',
    'Surface_Current_ms', 'Crustal_Age_Myr', 'Deployment_Length_days','Inventory']
    catevs=[];sk=-1
    for si,s in enumerate(icat.iloc):
        s = s.copy()
        for k in keys:s[k]=s.Deployment[k]
        events = s.Events.copy()
        for ei,ev in enumerate(events):
            snm=s.StaName
            evs = s.copy()
            IRIS_EVENT_LINK=lambda ev=ev: f'https://ds.iris.edu/ds/nodes/dmc/tools/event/{str(ev.resource_id).split('=')[-1]}'
            TF=evs.Data.TF;Noise = evs.Data.Noise
            Phases = lambda a=[s.Latitude,s.Longitude,s.StaDepth],b=ev,phases=('ttall',):event_stream_arrivals(a,b,phases=phases)
            Coherence=lambda snm=snm,evna=ev.Name:local_tools.io.get_reports(snm,evna=evna)
            Traces=lambda event=ev.Name,snm=snm,channel='HZ',tf='sta.ZP-21',methods=['Original','NoiseCut','ATaCR'],preproc=True:local_tools.io.get_traces(snm,event,channel=channel,tf=tf,methods=methods,preproc=preproc)
            evs['Name'] = ev.Name
            evs['Origin'] = ev.origins[0].time
            evs['LaLo'] = [ev.origins[0].latitude,ev.origins[0].longitude]
            evs['Link'] = IRIS_EVENT_LINK
            evs['Distance'] = distance(s,ev)
            evs['Phases']=Phases
            evs['Magnitude'] = ev.magnitudes[0].mag
            evs['Type']=ev.magnitudes[0].magnitude_type
            evs['Stations'] = ev.Stations
            evs['EvDepth']=ev.origins[0].depth/1000
            evs['Data']=AttribDict({
            'TF':TF,
            'Noise':Noise,
            'Traces':Traces,
            'Coherence':Coherence})
            evs['Event']=ev
            Traces=evs.Data.Traces
            evs['Traces']=Traces
            evs = evs[evskeys]
            catevs.append(evs.to_frame().T)
            sk+=1
    EventsCat = pd.concat(catevs)

    # evs = s.copy()
    # IRIS_EVENT_LINK = lambda ev=ev: f'https://ds.iris.edu/ds/nodes/dmc/tools/event/{str(ev.resource_id).split('=')[-1]}'
    # TF = evs.Data.TF;Noise = evs.Data.Noise
    # Coherence = lambda snm=snm,evna=ev.Name:local_tools.io.get_reports(snm,evna=evna)
    # Traces = lambda event=ev.Name,snm=snm:local_tools.io.get_traces(snm,event)


    EventsCat.set_index('Name', inplace=True,drop=False)
    return EventsCat
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 978-1016
# ----------------------------------------------------------------------------------------
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
def rcat_to_scat(r):
        evs = local_tools.cat.unravel_cat(r)
        evs = pd.DataFrame.from_dict(evs)
        evs['Description']=[e[0].text for e in  evs.event_descriptions]
        evs['Type']=[e[0].magnitude_type for e in evs.magnitudes]
        evs['Origin']=[e[0].time for e in evs.origins]
        IRIS_EVENT_LINK = lambda ev: f'https://ds.iris.edu/ds/nodes/dmc/tools/event/{str(ev.resource_id).split('=')[-1]}'
        evs['EventPage']=[IRIS_EVENT_LINK(e) for e in evs.iloc]
        evs['LaLo']=evs.lla
        evs['Magnitude'] = evs['mag']
        usefulkeys = ['Name','Magnitude','Type', 'depth','Description','LaLo','Stations','EventPage','Origin']
        evs = evs[usefulkeys]
        return evs
# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# ----------------------------------------------------------------------------------------
# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
class dataspace:
    """Generic container class for signal-related data and metadata."""
    def __init__(self,sets=['all'],aggregate=True,demo=False,catfile=None,hjan=None,**kwargs):
        self.sets=sets
        self.sets=[i.lower().replace('hjan23','Janiszewski23'.lower()) for i in self.sets]
        self.aggregate=aggregate
        self.demo=demo
        self.catfile = lt.io.dir_libraries().Catalogs/'Catalog_041425.pkl' if catfile is None else Path(catfile)
        self.hjanfile = lt.io.dir_libraries().Catalogs/'Janiszewski_etal_2023_StationList.pkl' if hjan is None else Path(hjan)
        self.load()
        self.index()

        if self.demo:self._restrict_data()
        for k, v in kwargs.items():setattr(self, k, v)

        self._remove_depreciated_sr_pairs()

    def load(self):
        if 'all' in self.sets:
            self.sets=['sr','r','s','Janiszewski23'.lower()]
            self.__import_cat()
            self.__rcat_to_scat()
            self.__import_HJan23()
            self.__rcat_to_srcat()
        elif 'r' in self.sets:self.__import_cat()
        elif 's' in self.sets:self.__import_cat();self.__rcat_to_scat()
        elif 'Janiszewski23'.lower() in self.sets:self.__import_HJan23(self.hjanfile)
        elif 'sr' in self.sets:self.__import_cat();self.__rcat_to_srcat()
        if 'sr' in self.sets:
            self.sr['Index']=self.sr.Name;self.sr.set_index('Index', inplace=True, drop=True)
            self.sr.sort_values(by=['Magnitude','Network'],ascending=False,inplace=True)
        if 'r' in self.sets:
            self.r['Index']=self.r.StaName;self.r.set_index('Index', inplace=True, drop=True)
            self.r.sort_values(by=['StaDepth','Network'],ascending=True,inplace=True)
        if 's' in self.sets:
            self.s['Index']=self.s.Name;self.s.set_index('Index', inplace=True, drop=True)
            self.s.sort_values(by=['Magnitude','Name'],ascending=[False,False],inplace=True)
        if 'Janiszewski23'.lower() in self.sets:
            self.Janiszewski23['Index']=self.Janiszewski23.StaName
            self.Janiszewski23.set_index('Index', inplace=True, drop=True)
        if (self.aggregate)&('sr' in self.sets):
            dirs = lt.io.dir_libraries()
            df=load_pickle(dirs.Analysis/'BulkLoad.SR.Coherences_092625.pkl')
            df=df.iloc[np.where(np.isin(np.array([df.Name + '_' + df.StaName]),np.array([self.sr.Name + '_' + self.sr.StaName])))[1]]
            if self.aggregate:
                self.sr['Coherence']=[df.aloc[sr.Name].aloc[sr.StaName].iloc[0].Coherence for sr in self.sr.iloc]
                assert sum(np.array([d.Coherence.event==d.Name for d in self.sr.iloc])==False)==0
                assert sum(np.array([d.Coherence.StaName==d.StaName for d in self.sr.iloc])==False)==0
            load_snr=False
            if (self.aggregate)&load_snr:
                df=load_pickle(dirs.SNR/'SNR.Models'/'SNR_acausul.filter_V04_5s_bandwidth_100_bands.pkl')
                df=df.iloc[np.where(np.isin(np.array([df.Name + '_' + df.StaName]),np.array([self.sr.Name + '_' + self.sr.StaName])))[1]]
                if self.aggregate:
                    self.sr['SNR']=[df.aloc[sr.Name].aloc[sr.StaName].iloc[0].SNR for sr in self.sr.iloc]
                    pk=list(self.sr.iloc[0].SNR.keys())[-1]
                    assert sum(np.array([d.SNR[pk].stnm==d.StaName for d in self.sr.iloc])==False)==0

    def index(self):
        index_cols=['Name','StaName','Station','Network',
        'Experiment','Pressure_Gauge','Seismometer','Instrument_Design']
        # idx_dummy=[f"{c}_idx" for c in index_cols]
        idx_dummy=['_'*(ci+1) for ci,c in enumerate(index_cols)]
        for src, dst in zip(index_cols, idx_dummy):self.sr[dst]=self.sr[src]
        self.sr.set_index(idx_dummy,drop=True,inplace=True)
        index_cols=['StaName','Station', 'Network','Experiment',
        'Pressure_Gauge','Seismometer','Instrument_Design']
        # idx_dummy=[f"{c}_idx" for c in index_cols]
        idx_dummy=['_'*(ci+1) for ci,c in enumerate(index_cols)]
        for src, dst in zip(index_cols, idx_dummy):self.r[dst]=self.r[src]
        self.r.set_index(idx_dummy,drop=True,inplace=True)

    def _restrict_data(self,demo=None):
        if demo is None:demo=self.demo
        assert demo is not None, 'demo=None. No station constraint given.'
        demo=[self.demo] if isinstance(self.demo,str) else self.demo
        print(f'\nWarning. {demo} set as the demo station{'s' if len(demo)>1 else ''}.\nOnly data for this constraint will be shown.')
        self.r=self.r[self.r.StaName.isin(demo)]
        self.sr=self.sr[self.sr.StaName.isin(demo)]


    def copy(self):
        """Return a shallow copy of the Signal."""
        import copy
        return copy.copy(self)
    # --------------------------hidden supports--------------------------
    def __repr__(self):
        return '\n'.join(['Dataspace: A highly organized longitudinal data class.',
        '\nRun self.load() for default datasets.',
        'demo [string/iterable] | Constrains the catalog to only stations defined by the demo argument. Demo=False (default)',
        '\n',
        'aggregate [bool] | Loads all available data into a lambda pipeline for very fast use by the catalog.',
        'This includes all SNR, coherence, phase, admittance, tilt data, noise and event traces with quick plotting tools for each.',
        '\n',
        ])
    def __import_cat(self,**kwargs):self.r=import_cat(catfile=self.catfile,hjanfile=self.hjanfile,aggregate=False)
    def __rcat_to_srcat(self):self.sr=rcat_to_srcat(self.r)
    def __import_HJan23(self):self.Janiszewski23=import_HJan23(self.hjanfile)
    def __rcat_to_scat(self):self.s=rcat_to_scat(self.r)

    def _remove_depreciated_sr_pairs(self):
        #The following SR pairs were quarantined on April 2nd, 2025 for extreme outlier SNR reduction (negative)
        depreciated_ = np.array(['2014.109.13.28', '2010.222.05.23', '2011.175.03.09',
            '2014.101.07.07', '2011.097.14.32', '2013.242.16.25',
            '2015.142.21.45', '2015.208.21.41', '2011.232.18.19',
            '2011.236.17.46', '2011.246.22.55', '2014.144.09.25',
            '2010.216.22.01', '2013.033.14.17', '2015.208.04.49',
            '2013.134.00.32', '2013.289.10.31', '2014.215.00.22',
            '2011.130.08.55', '2014.320.22.33', '2013.274.03.38',
            '2013.144.14.56', '2015.116.07.09', '2015.191.04.12'])
        self.sr = self.sr[~self.sr.Name.isin(depreciated_)]
    
# D = []
# icat = cat.sr.copy()
# for sr in icat.iloc:
#     d = sr.Data.Coherence()
#     d.ATaCR.zp_21.coh
#     d.HPS.zz.coh
#     D.append(AttribDict({'TF':d.ATaCR.zp_21.coh,'HPS_Z':d.HPS.zz.coh,'HPS_H':np.dstack([d.HPS['11'].coh,d.HPS['22'].coh]).mean(axis=2),'event':d.ATaCR.zp_21.events[0],'StaName':sr.StaName}))
# icat['Coherence'] = D
# file = dirs.Analysis/'BulkLoad.SR.Coherences.pkl'
# icat.to_pickle(file)