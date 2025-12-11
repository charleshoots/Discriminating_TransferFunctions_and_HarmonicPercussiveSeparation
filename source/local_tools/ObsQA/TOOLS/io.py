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

from modules import *
from local_tools import ObsQA
import time  # Example usage for sleep
import multiprocessing
from obspy.core import UTCDateTime as UTCDateTime
import concurrent
import datetime
from concurrent.futures import wait
import time
import logging
from packages.NoiseCut import *
from IPython.display import clear_output
from scipy import signal
import local_tools as lt
from local_tools import math
octavg=lt.math.octave_average
from functools import reduce
# from modules import modules
import obstools
from obstools.scripts import comply_calculate, atacr_clean_spectra, atacr_correct_event, atacr_daily_spectra, atacr_download_data, atacr_download_event, atacr_transfer_functions
# unpacker=globals().update
# unpack=lambda Args,keys=None:[unpacker({k:Args[k]}) for k in [keys if keys is not None else list(Args.keys())][0]]
def fv(func):
        code=func.__code__;variables=code.co_varnames[:code.co_argcount]
        return list(variables)

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# from classes import OBSMetrics 
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ----------------------------------------------------------
#### ----------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ====================================================================================
def _valid_win_length_samples(win_length_samples, win_length, sampling_rate):
    _next_pow2 = lambda n:int(round(2**np.ceil(np.log2(n))))
    if win_length_samples is None and win_length is None:win_length_samples = _next_pow2(120*sampling_rate)
    elif win_length_samples is None and win_length is not None:win_length_samples = _next_pow2(win_length*sampling_rate)
    elif win_length is not None and win_length_samples is not None:raise ValueError('Parameters win_length and win_length_samples are both defined.')
    elif win_length_samples is not None and win_length is None:
        win_length_samples = int(win_length_samples)
        if win_length_samples != _next_pow2(win_length_samples):raise ValueError('Parameter win_length_samples must be a power of 2.')
    return win_length_samples
# ====================================================================================
def spectrogram(trace,fs=None,win_length=163.84,pow2db=False,mag=False):
    if fs==None:fs=trace.stats.sampling_rate
    if isinstance(trace,Trace):x=trace.data.astype(float)
    win_length_samples=None
    win_length_samples = _valid_win_length_samples(win_length_samples, win_length, fs)
    hop_length = win_length_samples // 4
    n_fft = win_length_samples
    S = librosa.stft(x,n_fft=n_fft,hop_length=hop_length,win_length=win_length_samples)
    if mag:S,phase=librosa.magphase(S)
    df = fs/win_length_samples
    f = np.arange(S.shape[0]) * df
    t = np.arange(S.shape[1]) * hop_length
    t = t/fs
    if pow2db:S=librosa.power_to_db(np.abs(S))


    return S,f,t
# ====================================================================================
def ax_pcolormesh(T,F,S,ax=None,cmap='magma',vlim=[-98,-18],ymax=1,figsize=(10,4)):
       ax.cla()
       if ax is None:fig,ax=plt.subplots(nrows=1,ncols=1,figsize=figsize)
       else:fig=ax.figure
       pcm=ax.pcolormesh(T,F,S,cmap=cmap,shading='auto',vmin=vlim[0],vmax=vlim[1]);ax.set_ylim(F[1],ymax);ax.set_yscale('log')
       cbar=fig.colorbar(pcm, ax=ax, pad= 0.01);cbar.ax.tick_params(labelsize=7)
def plot_spectrogram(S,F,T,ax=None,ymax=1,figsize=(10,4),cmap='magma',vlim=[-98,-18],cbar=True):
        # hps_res_cmap = mcolors.LinearSegmentedColormap.from_list("CustomDiverging", list(zip([0, 0.5, 1], [(1, 1, 1), (0.5, 0.5, 0.5),(0, 0, 0)])))
        units = ['seconds','minutes','hours']
        while T[-1]>7200:T = T/60;units.pop(0)
        T=T/3600
        if ax is None:fig,ax=plt.subplots(nrows=1,ncols=1,figsize=figsize)
        else:fig=ax.figure
        pcm=ax.pcolormesh(T,F,S,cmap=cmap,shading='auto',vmin=vlim[0],vmax=vlim[1]);ax.set_ylim(F[1],ymax);ax.set_yscale('log')
        plt.yticks(fontsize= 7)
        ax.set_xticks([])
        if cbar:cbar=fig.colorbar(pcm, ax=ax, pad= 0.01);cbar.ax.tick_params(labelsize=7);cbar.set_label('dB', fontsize=7)  # Customize font size and weight
        # plt.xlabel(units[0])
        return fig,pcm
# ====================================================================================
def basic_preproc(d,lowpass=1,percent=.03):
        d.detrend('linear');d.detrend('demean')
        if not percent==None:d.taper(percent)
        if not lowpass==None:d.filter('lowpass',freq=lowpass,zerophase=True)
        d.detrend('linear');d.detrend('demean')
        return d

def unravel_report(file,AVG=False):
    r = load_pickle(file)
    nets=[n for n in list(r.keys()) if n[0]=='n']
    d=AttribDict()
    for n in nets:
        stas=list(r[n].keys())
        for s in stas:
            stanm=f'{n[1:]}.{s}'
            d[stanm]=AttribDict()
            events=np.array([e.replace('.HDH.SAC','').replace('.H2.SAC','').replace('.HZ.SAC','').replace('.H1.SAC','') 
            for e in r[n][s].events])
            d[stanm].events=events
            d[stanm].coh=r[n][s].coh
            keep=np.isnan(d[stanm].coh).sum(axis=1)==0
            d[stanm].coh=d[stanm].coh[keep,:]
            d[stanm].events=d[stanm].events[keep]
            d.f=r.f
            d[stanm].f = d.f
            if AVG:d[stanm].coh=d[stanm].coh.mean(axis=0)
    return d
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        b = pickle.load(handle)
    return b

def reshapedict(c):
    n = np.array(list(c.keys()));n=n[n!='f']
    sn = unravel([[f'{ni[1:]}.{s}' for s in c[ni].keys()] for ni in n])
    r = AttribDict({f'{s}':c[f'n{s.split('.')[0]}'][s.split('.')[1]] for s in sn})
    r.f = c.f
    return r

def get_reports(sta=None,
    Archive=None,
    methods=['ATaCR','NoiseCut'],evna=None):
    if Archive is None:dirs=lt.io.dir_libraries();Archive=dirs.Data
    methods = np.array(methods)
    Report=AttribDict()
    AnalysisFolder = Archive/'Analysis'/'Coherence'
    lp = lambda f:reshapedict(load_pickle(f))
    atacrfile = lambda:AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{tf}.zz.coherence.pkl'
    hpsfile = lambda :AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{c}.coherence.pkl'
    unionevents = lambda a:reduce(np.intersect1d,a)
    methods = [m.replace('NoiseCut','HPS') for m in methods]
    for m in methods:
        if m.lower()=='atacr':
            tfs = ['zp_21','z2_1','zp']
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
        for s in stas:
            if s=='f':continue
            evs=unionevents([Report[m][list(Report[m].keys())[0]][s].events for m in methods])
            ind = np.isin(Report[m][c][s].events,evs)
            for i in Report[m].keys():
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
        Report  = lt.cat.evcoh(Report,evna)
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
        st.detrend('linear');st.detrend('demean')
        return st
def get_traces(stanm,event,channel='HZ',tf='sta.ZP-21',methods=['Original','NoiseCut','ATaCR'],preproc=True):
        if type(event)==obspy.core.event.event.Event:event=event.Name
        dirs=lt.io.dir_libraries()
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
                st=Stream([load_sac(f) for f in files])
                for s,m in zip(st,methods):s.stats.location=m
                if preproc:st=basic_preproc(st)
                return st
        else:print(f'Not all files exist : {stanm} | {event}');return None
# ====================================================================================

def audit_events(eventsfolder,Minmag=6.0,Maxmag=7.0):
    client = Client()
    timedelta = 60
    stations = ObsQA.TOOLS.io.get_event_catalog(eventsfolder)
    evlist = []
    [evlist.extend(stations.events[j]) for j in range(len(stations))]
    evlist = np.unique(evlist)
    evlist
    evaudit = dict()
    for ev in evlist:
        inds = [np.array([k for k in [np.where(np.array(stations.events[j])==ev)[0]]])[0] for j in range(len(stations))]
        for i in range(len(inds)):
            if len(inds[i])==0:
                inds[i] = None
            else:
                inds[i] = inds[i][0]
        idx = np.where([(j is not None) for j in inds])[0]
        stas = stations.Station[idx].tolist()
        evaudit[ev] = [stas]
        evaudit[ev] = [idx]
        ## evaudit['nsta'].append(len(stas))
    evaudit = pd.DataFrame.from_dict(evaudit,orient='index',columns=['Stations'])
    evaudit['Event'] = evaudit.index
    evaudit = evaudit.reset_index(drop=True)[['Event','Stations']]
    evaudit['nsta'] = [0 for _ in range(len(evaudit))]
    for j in range(len(evaudit)):
        evaudit.at[j,'Networks'] = np.array(list(stations.Network[evaudit.Stations[j]]))
        net = np.unique(np.array(list(stations.Network[evaudit.Stations[j]])))
        evaudit.at[j,'Stations'] = list(stations.Station[evaudit.Stations[j]])
        evaudit.at[j,'nsta'] = len(evaudit.iloc[j].Stations)
        
        print('[' + str(j) + '/' + str(len(evaudit)) + '] event: ' + evaudit.iloc[j].Event + ' || nsta: ' + str(evaudit.iloc[j].nsta))
        ev = evaudit.Event[j]
        starttime = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
        endtime = starttime + datetime.timedelta(seconds=timedelta)
        print('OBSPY Events: ' + str(starttime) + ' --> ' + str(endtime) + ' | Mags: m' + str(Minmag) +  '-' + str(Maxmag))
        cat = client.get_events(starttime=starttime, endtime=endtime,minmagnitude=Minmag, maxmagnitude=Maxmag)
        km = cat[0].origins[0].depth/1000
        mw = cat[0].magnitudes[0].mag
        lla = cat[0].origins[0].latitude,cat[0].origins[0].longitude
        evaudit.at[j,'MW'] = mw
        evaudit.at[j,'depth'] = km
        evaudit.at[j,'ev_lla'] = lla[0]
        evaudit.at[j,'ev_lon'] = lla[1]
    return evaudit
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def get_ATACR_HPS_Comp(d,atacr_TFs_used='ZP-21', win_length=163.84):
        RawP = d.Raw[0]['trP'].copy()
        RawZ = d.Raw[0]['trZ'].copy()
        PostATACR = d.Corrected[0][atacr_TFs_used].copy()
        PostATACR.stats.location = 'ATaCR (' + PostATACR.stats.location + ')'
        PostHPS, spectrograms = noisecut(RawZ.copy(), ret_spectrograms=True,win_length=win_length)
        PostBoth, spectrograms = noisecut(PostATACR.copy(), ret_spectrograms=True,win_length=win_length)
        return RawP,RawZ,PostATACR,PostHPS,PostBoth

#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
        times = [[a.name,a.time,a.incident_angle] for a in arrivals]
        return times
def event_stream_arrivals(R,S,phases=('ttall',),
        # shadow_phases=('PKIKP','SKS','SKIKSSKIKS'),
        ):
        if isinstance(R,Trace):
                if (R.stats.starttime>(S.origins[0].time+3600*5)) or (R.stats.endtime>(S.origins[0].time+3600*5)):
                       raise Exception('Event was not observed in this trace')
                stallaz=[R.stats.sac.stla,R.stats.sac.stlo,R.stats.sac.stel]
        else:stallaz=R
        if isinstance(S,obspy.core.event.event.Event):
                evllaz=[S.origins[0].latitude,S.origins[0].longitude,S.origins[0].depth/1000]
        else:evllaz=S
        # gcarc = locations2degrees(stallaz[0],stallaz[1],evllaz[0],evllaz[1])
        # if (gcarc>103) & (gcarc<141):
        #         if not (phases==('ttall',)):phases+=shadow_phases
        ar=get_arrivals(stallaz,evllaz,model='iasp91',phases=phases)
        ard = dict();[ard.update({a[0]:a[1:]}) for a in ar]
        return ard
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def obspyfft(d):
        signal = d.data
        fourier = np.abs(np.fft.rfft(signal))
        freq = np.abs(np.fft.rfftfreq(signal.size, d=1./d.stats.sampling_rate))
        return freq,fourier
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_event_catalog(eventsfolder,subfolder='CORRECTED',fmt = 'pkl'):
        evdict = dict()
        evdict['folder'] = [fld.split('/')[-2] for fld in g.glob(eventsfolder + '/*/')]
        evdict['Network'] = [fld.split('.')[0] for fld in evdict['folder']]
        evdict['Station'] = [fld.split('.')[1] for fld in evdict['folder']]
        evdict['StaName'] = [str(a) + '.' + str(b) for a,b in zip(evdict['Network'],evdict['Station'])]
        # evdict['folder'][0]
        evdict['n_events'] = list()
        evdict['events'] = list()
        for i in range(len(evdict['folder'])):
                staname =  evdict['StaName'][i]
                path = Path(eventsfolder)
                path = path / evdict['folder'][i] / subfolder
                files = list(path.glob('*.' + fmt))
                events = [str(f.name).replace(staname + '.','').replace('.sta','').replace('.day','').replace('.' + fmt,'') for f in files]
                # events = list(np.unique(['.'.join(files[g].split('.' + fmt)[0].split('.')[0:-1]) for g in range(len(files))]))
                evdict['events'].append(events)
                evdict['n_events'].append(len(events))
        catalog = pd.DataFrame(evdict)
        catalog = catalog.sort_values(by=['Network','Station'])
        catalog = catalog.reset_index(drop=True)
        return catalog
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def getstalist():
        # current_path = os.path.dirname(__file__)
        # excelfile = current_path + '/Janiszewski_etal_2023_StationList.xlsx'
        dirs=lt.io.dir_libraries()
        excelfile = dirs.Catalogs/'Janiszewski_etal_2023_StationList.xlsx'
        stations = pd.read_excel(excelfile)
        d = dict()
        d.update({a:b for a,b in zip(stations.columns,[c.replace('(','').replace(')','').replace(' ','_').replace('/','') for c in stations.columns])})
        stations = stations.rename(columns=d)
        stations.Station = [str(s) for s in stations.Station]
        stations.Network = [str(s) for s in stations.Network]
        stations['StaName'] = stations.Network.astype(str) + '.' + stations.Station.astype(str)
        allgood = np.in1d(stations[['Z_Is_Good','H1_Is_Good','H2_Is_Good','P_Is_Good']].sum(axis=1).tolist(),4)
        stations['Good_Channels'] = allgood
        stations = stations.assign(n_events=pd.Series())
        stations = stations.assign(Magnitude_mw=pd.Series())
        stations = stations.assign(Depth_KM=pd.Series())
        stations = stations.assign(Origin=pd.Series())
        stations = stations.assign(Metadata=pd.Series())
        stations = stations.assign(Averaging=pd.Series())
        stations = stations.assign(Events=pd.Series())
        stations = stations.assign(Files=pd.Series())
        stations = stations.sort_values(by=['Network','Station'])
        stations = stations.reset_index(drop=True)
        return stations

#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def mat2df(files):
        '''
        Converts all loaded matlab variables into a single dataframe
        '''
        if not isinstance(files,list):
                files = [files]
        files = sorted(files)
        out = []
        df = pd.DataFrame()
        FNUM = 0
        for f in files:
                df2 = pd.DataFrame()
                cfold = f[0:maxstrfind(f,'/')+1]
                cf = f[maxstrfind(f,'/')+1:len(f)]
                mat = spio.loadmat(f, simplify_cells=True)
                matvars = list(mat.keys())[3:len(list(mat.keys()))]
                oldkey = list(mat.keys())[-1] #Replaces Last Key with a single hardcoded name. Assumes a single variable (ie one struct) was saved to the matfile
                for k in range(3): #pop out the first three keys. they are artifacts from the mat import
                        mat.pop(list(mat.keys())[0])
                for k in matvars:
                        if isinstance(mat[k],list):
                                for dct in mat[k]:
                                        tmp = pd.DataFrame.from_dict(dct,orient='index').T
                                        tmp['FNUM'] = FNUM
                                        df2 = pd.concat([df2,tmp])
                                mat[k] = df2
                        else:
                                mat[k] = pd.DataFrame.from_dict(mat[k]).T
                                mat[k]['FNUM'] = FNUM
                if len(matvars)>1:
                        for i in range(len(matvars)-1):
                                mat[matvars[0]] = pd.merge(mat[matvars[0]],mat[matvars[i+1]],on='FNUM')
                fdf = mat[matvars[0]].set_index('FNUM')
                fdf['Folder'] = cfold
                fdf['File'] = cf
                out.append(fdf)
                FNUM += 1
        for d in out:
                df = pd.concat([df,d])
        df = df.groupby('FNUM',as_index=False,dropna=False).ffill().groupby('FNUM',as_index=False,dropna=False).bfill() #<---Fills data gaps using adjacent rows from only the same file.
        df = df.sort_values('File')
        return df
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def maxstrfind(s,p):
        '''
        Wild this doesnt exist in Python libraries. Finds LAST occurence of expression in a string.
        '''
        i = 0
        while i>-1:
                io = i
                i = s.find(p,i+1,len(s))
        return io
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def datenum_to_datetime64(dnum):
        '''
        Just some Matlab DateNum nonsense
        '''
        days = np.asarray(dnum) - 719529  # shift to unix epoch (1970-01-01)
        return np.round((days * 86400000)).astype("datetime64[ms]")
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_noise_metrics(atacr_parent,net,sta,avg='Day',update=True):
        atacr_parent = Path(atacr_parent)
        if avg.lower()=='sta':
                specfold = 'AVG_STA'
                filetag = 'avg_sta'
        elif avg.lower()=='day':
                specfold = 'SPECTRA'
                filetag = 'spectra'
        subfolder = atacr_parent / specfold / '.'.join([net,sta])
        if update==False:
                out = pd.read_pickle(subfolder / 'Metrics' /  ('.'.join([net,sta])+ '.' + avg + 'Metrics.pkl'))
        else:
                # print(subfolder)
                s = loadpickles(subfolder)
                out = []
                for o,f in zip(s.Output,s.File):
                        keys = list(o.__dict__.keys())
                        keys.remove('cross')
                        keys.remove('power')
                        if avg.lower()=='day':
                                keys.remove('ft1')
                                keys.remove('ft2')
                                keys.remove('ftZ')
                                keys.remove('ftP')
                        vals = [[o.__dict__[k]] for k in keys]
                        noise = dict()
                        noise['StaNoise'] = [o]
                        noise['File'] = f
                        for k,v in zip(keys,vals):
                                noise[k] = v
                        if avg.lower()=='sta':
                                year = [int(s.split('.')[0]) for s in f.split('.avg_sta.pkl')[0].split('-')]
                                jday = [int(s.split('.')[1]) for s in f.split('.avg_sta.pkl')[0].split('-')]
                                noise['year'] = [year]
                                noise['julday'] = [jday]
                        if avg.lower()=='day':
                                ft = dict()
                                ft['1'] = o.__dict__['ft1']
                                ft['2'] = o.__dict__['ft2']
                                ft['Z'] = o.__dict__['ftZ']
                                ft['P'] = o.__dict__['ftP']
                                tr = [np.mean(np.real(np.fft.ifft(ft[k][o.__dict__['goodwins']])),axis=0) for k in list(ft.keys())]
                                GoodNoise = dict()
                                for k,t in zip(list(ft.keys()),tr):
                                        GoodNoise[k] = t
                                noise['ft'] = [ft]
                                noise['GoodNoise'] = [GoodNoise]
                        noise_out = pd.DataFrame.from_dict(noise)
                        if avg.lower()=='sta':
                                spec = o.power.__dict__
                                spec.update(o.cross.__dict__)
                                spec['f'] = o.__dict__['f']
                                noise['CSD'] = [spec]
                                noise_trace = dict()
                                noise_trace_avg = dict()
                                gd = noise_out.gooddays[0]
                                gw = ObsQA.TOOLS.io.get_noise_metrics(atacr_parent,net,sta,avg='Day',update=True).goodwins[gd]
                                ft = ObsQA.TOOLS.io.get_noise_metrics(atacr_parent,net,sta,avg='Day',update=True).ft[gd]

                                for c in ['1','2','Z','P']:
                                        noise_trace[c] = [np.mean(spec[c][g],axis=0) for spec,g in zip(ft.iloc,gw.iloc)]
                                        noise_trace_avg[c] = np.real(np.fft.ifft(np.mean(np.array(noise_trace[c]),axis=0)))
                                noise['Noise_Average'] = [noise_trace_avg]
                                noise['Noise'] = [noise_trace]
                                noise_out = pd.DataFrame.from_dict(noise)
                        if len(out)==0:
                                out = noise_out
                        else:
                                out = pd.concat([out,noise_out])
                if avg.lower()=='day':
                        out = out[['StaNoise','GoodNoise','File','key','tkey','julday','year','QC','av',
                                        'window','overlap','dt','npts','fs','ncomp','tf_list',
                                        'goodwins','rotation','f','ft']]
                        out = out.sort_values(by=['year','julday'],ascending=False)
                else:
                        out = out[['key','StaNoise','Noise','Noise_Average','File', 'year', 'julday','f', 'nwins']]
        out.reset_index(drop=True,inplace=True)
        if avg.lower()=='sta':
                Metrics = dict()
                Metrics_Noise = ObsQA.OBSM.classes.OBSMetrics(csd=out.StaNoise[0],f=out.StaNoise[0].f)
                # Metrics['Noise'] = out
                Metrics['Noise'] = Metrics_Noise
                out = Metrics
        return out
def get_noise(stakey,noisedir,attempt_metrics=False):
        # _____________________________________________________________________
        # "attempt_metrics" is experimental. 
        # it recreates the raw traces and associated spectra used in the station average
        # that past the atacr QC/QA steps. 
        # this is an attempt to have the most raw form of the noise data instead of 
        # just these intermittent datatypes, namely channel csd's and psd's.
        # the benefit, aside from having direct access to the raw noise data, 
        # is it can now be pushed into an OBSMetrics object.
        # It works and is decimal accurate to the StaNoise objects saved mid-way in ATaCR.
        # Nonetheless, this it is still experimental. 
        # _____________________________________________________________________
        avg = load_pickle(list((noisedir.parent/'AVG_STA'/stakey).glob('*sta.pkl'))[0])
        if not attempt_metrics:
                return avg
        else:
                gooddays = np.array(list((noisedir.parent/'Spectra'/stakey).glob('*.spectra.pkl')))[avg.gooddays]
                days = np.array(list((noisedir / stakey).glob('*Z.SAC')))[avg.gooddays]
                days = [d.name.replace('..','.').replace('.HZ.SAC','') for d in days]
                overlap,window=load_pickle(gooddays[0]).overlap,load_pickle(gooddays[0]).window
                goodwins = [load_pickle(g).goodwins for g in gooddays]
                for d in days:
                        noise = Stream([load_sac(noisedir/stakey/''.join([d,'..',c,'.SAC']),rmresp=True)[0][0] for c in ['H1','H2','HZ','HDH']])
                        n = [s for s in noise.split(window_length=window,step=window*(1-overlap))]
                # -------Build noise
                noise = Stream([load_sac(noisedir/stakey/''.join([d,'..',c,'.SAC']),rmresp=True)[0][0] for c in ['H1','H2','HZ','HDH'] for d in days])
                NoiseTr=Stream([noise.select(channel='*Z')[0].copy(),noise.select(channel='*1')[0].copy(),noise.select(channel='*2')[0].copy(),noise.select(channel='*DH')[0].copy()]).copy()
                for c in ['H1','H2','HZ','HDH']:NoiseTr.select(channel='*'+c)[0].data=np.mean(noise.select(channel='*'+c),axis=0)
                # -------Build metrics
                NoiseM=OBSM.Metrics(tr1=NoiseTr.select(channel='*1')[0],tr2=NoiseTr.select(channel='*2')[0],trZ=NoiseTr.select(channel='*Z')[0],trP=NoiseTr.select(channel='*DH')[0])
                return NoiseM

#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def ClosestMLPreEventTF(tfdir,event_time):
        '''
        Replicates the ML ATaCR method for defining the 'nearest' TF to an event. It is a slightly different method than the Python version.
        '''
        files = g.glob(tfdir + '/*.mat')
        files = [ele for ele in files if '_AVERAGE_' not in ele] #<--Remove the station average from file list
        eventids = np.array(list(map(int,[f.split('/')[-1].split('_')[1] for f in files]))) #<Split file strings down to their event ids. Assumes format dir/stuff_EVENTID_stuff
        f = files[np.where(((np.array(event_time,dtype=int) - eventids)>0) & ((np.array(event_time,dtype=int) - eventids)==np.min((np.array(event_time,dtype=int) - eventids))))[0][0]]
        out = mat2df(f)
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def organize_evdata(evdata):
        '''
        Builds ObsPy trace objects inside the dataframes from the imported Matlab data
        '''
        itr1 = (evdata['channel'].squeeze().str.find('1')>0)
        itr2 = (evdata['channel'].squeeze().str.find('2')>0)
        itrZ = (evdata['channel'].squeeze().str.find('Z')>0)
        itrP = ~((evdata['channel'].squeeze().str.find('1')>0) + (evdata['channel'].squeeze().str.find('2')>0) + (evdata['channel'].squeeze().str.find('Z')>0))
        tmp = evdata[itr1]
        tr1 = Trace(data=tmp['data'].squeeze())
        tr1.stats.sampling_rate = tmp.sampleRate.squeeze()
        tr1.stats.network = tmp.network.squeeze()
        tr1.stats.station = tmp.station.squeeze()
        tr1.stats.channel = tmp.channel.squeeze()
        tr1.stats.location = 'ML-ATaCR-PreProc'
        tr1.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        tmp = evdata[itr2]
        tr2 = Trace(data=tmp['data'].squeeze())
        tr2.stats.sampling_rate = tmp.sampleRate.squeeze()
        tr2.stats.network = tmp.network.squeeze()
        tr2.stats.station = tmp.station.squeeze()
        tr2.stats.channel = tmp.channel.squeeze()
        tr2.stats.location = 'ML-ATaCR-PreProc'
        tr2.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        tmp = evdata[itrZ]
        trZ = Trace(data=tmp['data'].squeeze())
        trZ.stats.sampling_rate = tmp.sampleRate.squeeze()
        trZ.stats.network = tmp.network.squeeze()
        trZ.stats.station = tmp.station.squeeze()
        trZ.stats.channel = tmp.channel.squeeze()
        trZ.stats.location = 'ML-ATaCR-PreProc'
        trZ.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        tmp = evdata[itrP]
        trP = Trace(data=tmp['data'].squeeze())
        trP.stats.sampling_rate = tmp.sampleRate.squeeze()
        trP.stats.network = tmp.network.squeeze()
        trP.stats.station = tmp.station.squeeze()
        trP.stats.channel = tmp.channel.squeeze()
        trP.stats.location = 'ML-ATaCR-PreProc'
        trP.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        evdata = evdata.assign(Trace=0)
        evdata.iat[np.squeeze(np.where(itr1)).tolist(), -1] = tr1
        evdata.iat[np.squeeze(np.where(itr2)).tolist(), -1] = tr2
        evdata.iat[np.squeeze(np.where(itrZ)).tolist(), -1] = trZ
        evdata.iat[np.squeeze(np.where(itrP)).tolist(), -1] = trP
        evdata.tr1 = tr1
        evdata.tr2 = tr2
        evdata.trZ = trZ
        evdata.trP = trP
        return evdata
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def demean_detrend(stream):
        stream.detrend()
        for i in range(len(stream)):
            stream[i].data = stream[i].data - np.mean(stream[i].data)  
        stream.detrend()
        return stream
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_noisecut_output(net,sta,event,noisecut_datafolder):
        noisecut_datafolder = Path(noisecut_datafolder)
        file = '.'.join([net,sta,event,'pkl'])
        stafolder = '.'.join([net,sta])
        if len(list((noisecut_datafolder / stafolder).glob(file)))==0:
                print('{} not found'.format(file))
                return []
        else:
                st = pd.read_pickle(noisecut_datafolder / stafolder / file)
                for i,_ in enumerate(st.Raw[0]):
                        st.Raw[0][i].stats.location = 'Raw'
                return st
import warnings,fnmatch,operator,itertools
def unravel(lst):return list(itertools.chain.from_iterable(lst))
def get_noisecut_event(parentfolder,staname,event,channel=['*1','*2','*Z','*H'],len_hrs=2,pre_trim=False,post_trim=True,win_length=163.84,width=None,verbose=False):
        # datafold = 'Data'
        dateformat = '%Y.%j.%H.%M'
        event = UTCDateTime.strptime(str(event),dateformat)
        origin = event
        file_day = event
        channel = [cc.replace('**','*') for cc in ['*' + c for c in np.array([channel]).flatten().tolist()]]
        folder = Path(str(parentfolder)) / staname
        files = [list(folder.glob(f'{file_day.strftime(dateformat)}.{c}.SAC'))[0] for c in channel]
        if len(files)==0:print('||--24hr event data not found');return []
        raw = Stream([load_sac(f) for f in files])
        raw = Stream([raw.select(component=c)[0] for c in channel if len(raw.select(component=c))>0])
        if (event+2*3600-100)>raw[0].stats.endtime:
                print(staname+'| Data Gap')
                return []
        if pre_trim:raw.trim(origin,origin+3600*len_hrs)
        raw = demean_detrend(raw.copy())
        print(staname,'|',event.strftime(dateformat),'|','Begin NoiseCut')
        if sum([np.sum(np.isnan(tr.data)) for tr in raw])>raw[0].data.shape[0]:
                print('||--BAD EVENT TRACE--||'+f' {staname} - {event.strftime(dateformat)}')
                return []
        hps_out = [noisecut(s.copy(),ret_spectrograms=True,win_length=win_length,width=width,verbose=verbose) for s in raw]
        corrected = Stream([h[0].copy() for h in hps_out])
        if post_trim:corrected.trim(origin,origin+3600*len_hrs);raw.trim(origin,origin+3600*len_hrs)

        hps_spectrograms = [h[1] for h in hps_out]
        raw_collect = raw
        corrected_collect = corrected
        e = [event for _ in range(len(raw_collect))]
        n = [staname.split('.')[0] for _ in range(len(raw_collect))]
        s = [staname.split('.')[1] for _ in range(len(raw_collect))]
        c = [r.stats.channel for r in raw]
        out = pd.DataFrame({'Event':e,'Network':n,'Station':s,'Channel':c,'Raw':raw_collect,'Corrected':corrected_collect,'Spectrograms':hps_spectrograms})
        return out,hps_out

def run_noisecut(tr,win_length=163.84,width=None):
        hps, (S_full, S_background, S_hps, frequencies, times) = noisecut(tr, ret_spectrograms=True,win_length=win_length,width=width)
        return hps, (S_full, S_background, S_hps, frequencies, times)

def Get_ATaCR_CorrectedEvents(eventfolder,eventnames,net,sta,tfavg='sta',tf=''):
        if not isinstance(eventnames,list):
            eventnames = [eventnames]
        if not isinstance(net,list):
            net = [net]
        if not isinstance(sta,list):
            sta = [sta]
        dobspy,raw_collect,corrected = [],[],[]
        for i in range(len(eventnames)):
            neti = net[i]
            stai = sta[i]
            evi = eventnames[i]
            prefix = neti + '.' + stai
            folder = eventfolder / prefix / 'CORRECTED'
            f = list(folder.glob(prefix + '.' + evi + '*.pkl'))
            if len(f)==0:
                return []
            f = f[0]
            ri = pkl.load(open(f,'rb'))
            raw = {'tr1':ri.tr1.copy(),'tr2':ri.tr2.copy(),'trZ':ri.trZ.copy(),'trP':ri.trP.copy()}
            for k in list(raw.keys()):
                    raw[k].stats.location = 'Raw' #For better book-keeping. Label the raw traces.
            trcorr = {}
            raw_collect.append(raw)
            # rawz.append(ri.trZ.copy())
            dobspy.append(ri)
            for k in list(ri.correct.keys()): #Shape corrected traces into a list of ObsPy trace objects
                    tr = ri.trZ.copy()
                    tr.data = ri.correct[k]
                    tr.stats.location = k
                    trcorr[k] = tr
            corrected.append(trcorr)
        out = pd.DataFrame({'Event':eventnames,'Network':net,'Station':sta,'Raw':raw_collect,'Corrected':corrected})
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def GetML_EventData_and_TransferFunctions(event_time,network,sta,preprocfolder=None,correctedfolder=None,tffolder=None):
        evdata = None
        corrected_evdata = None
        ml_tf = None
        if preprocfolder is not None:
                folder = preprocfolder
                path = folder + '/' + event_time + '/*.mat'
                ml_files = g.glob(path)
                evdata = ObsQA.TOOLS.io.mat2df(ml_files)
                evdata = evdata[(evdata['Network']==network)&(evdata['Station']==sta)]
                evdata = ObsQA.TOOLS.io.organize_evdata(evdata)
        if correctedfolder is not None:
                folder = correctedfolder
                path = folder + '/' + network + '/' + sta + '/' + network + sta + '_' + str(event_time) + '_corrseis' + '.mat'
                ml_files = g.glob(path)
                corrected_evdata = ObsQA.TOOLS.io.mat2df(ml_files)
                corrected_evdata = corrected_evdata[(corrected_evdata['Network']==network)&(corrected_evdata['Station']==sta)]
                corrected_evdata = corrected_evdata.assign(Trace=0)
                for i in range(len(corrected_evdata)):
                        tmp = corrected_evdata.iloc[i]
                        tr = Trace(data=tmp['timeseries'].squeeze())
                        tr.stats.sampling_rate = 1/tmp['dt']
                        tr.stats.network = tmp['Network']
                        tr.stats.station = tmp['Station']
                        tr.stats.channel = tmp['label']
                        tr.stats.location = 'ML-ATaCR-Corrected'
                        if evdata is not None:
                                tr.stats.starttime = evdata.trZ.stats.starttime
                        corrected_evdata.iat[i, -1] = tr
        if tffolder is not None:
                ftf = tffolder + '/' + network + '/' + sta
                ml_tf = ObsQA.TOOLS.io.ClosestMLPreEventTF(ftf,event_time)
        ml_tf = ml_tf.reset_index()
        corrected_evdata = corrected_evdata.reset_index()
        evdata = evdata.reset_index()
        TF_list = {i : j for i, j in zip(corrected_evdata.label.to_list(), np.ones(len(corrected_evdata.label.to_list()),dtype=bool))}
        corrected_evdata['TF_list'] = TF_list
        for i in range(len(corrected_evdata)):
                corrected_evdata.at[i,'TF_list'] = TF_list
        corrected_evdata = corrected_evdata[['label','Network','Station','eventid','TFfilename', 'dt','NFFT', 'f', 'filop',
        'taxis', 'tf_op', 'spectrum', 'timeseries','isgood', 'Folder', 'File', 'TF_list','Trace']]
        return evdata,corrected_evdata,ml_tf
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def condensedata(input):
        collect = pd.DataFrame()
        for data in input:
                if isinstance(data,pd.DataFrame):
                        key = np.append(np.array(data.keys()[np.array(np.where(data.keys()=='channel')).flatten().tolist()]),np.array(data.keys()[np.array(np.where(data.keys()=='label')).flatten().tolist()]))[0]
                        labels = list(data[key])
                        if len(np.flatnonzero(np.core.defchararray.find(labels,'Z')!=-1))==1:
                                # isPreProc
                                data = data.rename(columns={'channel': 'label'})[['Trace','Network','Station','label','File']].copy().iloc[np.where(data[key]=='BHZ')]
                                data['State'] = 'PreProc'
                                data['CodeBase'] = 'Matlab'
                        else:
                                # isCorrected
                                data = data[['Trace','Network','Station','label','File']].copy()
                                data['State'] = 'Corrected'
                                data['CodeBase'] = 'Matlab'
                elif isinstance(data,obstools.atacr.classes.EventStream):
                        d1 = pd.DataFrame.from_dict({'Trace':data.trZ,'Network':data.trZ.stats.network,'Station':data.trZ.stats.station,'label':data.trZ.stats.channel,'File':data.prefix,'State':'PreProc','CodeBase':'Python'},orient='index').T
                        d2 = pd.DataFrame.from_dict({'Trace':list(evstream.correct.values()),'Network':data.trZ.stats.network,'Station':data.trZ.stats.station,'label':list(data.correct.keys()),'File':data.prefix,'State':'Corrected','CodeBase':'Python'})
                        data = pd.concat([d1,d2])
                collect = pd.concat([collect,data]).sort_values(by=['CodeBase'])
        return collect
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def MLtoDayNoise(folder,netlist,stalist):
        Data = pd.DataFrame({'Network':netlist,'Station':stalist,'DayNoise':[ [] for _ in range(len(stalist)) ],'NDays':np.zeros(len(stalist)),'Files':[ [] for _ in range(len(stalist)) ]})
        lls = []
        llf = []
        for ista in range(len(stalist)):
                csta=stalist[ista]
                cnet=netlist[ista]
                cfold = folder + '/' + cnet + '/' + csta
                df = mat2df(g.glob(cfold + '/*.mat'))
                files = list(df.File.unique())
                DayList = []
        for i in range(len(files)):
                cdf = organize_evdata(df.iloc[list(df.File==files[i])])
                window = cdf.iloc[0].Trace.stats.endtime - cdf.iloc[0].Trace.stats.starttime
                overlap = 0.3 #Default set in ML's ATaCR code (setup_parameter)
                key = cdf.iloc[0].Trace.stats.network + '.' + cdf.iloc[0].Trace.stats.station
                tr1 = list(cdf[list(cdf.channel=='BH1')].Trace)[0]
                tr2 = list(cdf[list(cdf.channel=='BH2')].Trace)[0]
                trZ = list(cdf[list(cdf.channel=='BHZ')].Trace)[0]
                trP = list(cdf[list(cdf.channel=='BDH')].Trace)[0]
                DN = DayNoise(tr1=tr1, tr2=tr2, trZ=trZ, trP=trP, window=window, key=key)
                DN.QC = True
                DN.av = True
                DayList.append(DN)
        Data.loc[ista,'NDays'] = len(DayList)
        lls.append(DayList)
        llf.append(files)
        Data['DayNoise'] = lls
        Data['Files'] = llf
        return Data
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def MLtoStaNoise(cnet,csta,preprocfolder,specfolder,level='b1'):
        DN = ObsQA.TOOLS.io.MLtoDayNoise(preprocfolder,[cnet],[csta])
        cfold = specfolder + '/' + cnet + '/' + csta + '/' + level.lower()
        files = g.glob(cfold + '/*.mat')
        df = ObsQA.TOOLS.io.mat2df(files)
        files = df.File.unique()
        udf = df.copy().iloc[0].to_frame().T.drop(0)
        for i in range(len(files)):
                udf = pd.concat([udf,df.iloc[list(df.File==files[i])].iloc[0].to_frame().T])
        df = udf.copy()
        for ifile in range(len(files)):
                cfile = files[ifile]
                cdf = df.iloc[list(df.File==cfile)]
                c11 = cdf.iloc[0].c11_stack.T
                c22 = cdf.iloc[0].c22_stack.T
                cZZ = cdf.iloc[0].czz_stack.T
                cPP = cdf.iloc[0].cpp_stack.T
                c12 = cdf.iloc[0].c12_stack.T
                c1Z = cdf.iloc[0].c1z_stack.T
                c1P = cdf.iloc[0].c1p_stack.T
                c2Z = cdf.iloc[0].c2z_stack.T
                c2P = cdf.iloc[0].c2p_stack.T
                cZP = cdf.iloc[0].cpz_stack.T
                cHH = cdf.iloc[0].chh_stack.T
                cHZ = cdf.iloc[0].chz_stack.T
                cHP = cdf.iloc[0].chp_stack.T
                tilt = cdf.iloc[0].rotor
                rotcoh = cdf.iloc[0].rotcoh
                f = cdf.iloc[0].f
                nwins = len(cdf.iloc[0].goodwins)
                direc = np.arange(0,360+10,10) #default in ML ATaCR code (b1)
                ncomp = sum(cdf.iloc[0].comp_exist)
                key = cnet + '.' + csta
                tf_list = {'ZP': True, 'Z1': True, 'Z2-1': True, 'ZP-21': True, 'ZH': True, 'ZP-H': True}
                power = obstools.atacr.classes.Power(c11 = c11, c22 = c22, cZZ = cZZ, cPP = cPP)
                cross = obstools.atacr.classes.Cross(c12=c12, c1Z=c1Z, c1P=c1P, c2Z=c2Z, c2P=c2P, cZP=cZP)
                rotation = obstools.atacr.classes.Rotation(cHH=cHH, cHZ=cHZ, cHP=cHP, coh=None, ph=None, tilt=tilt, coh_value=rotcoh, phase_value=None, direc=direc)
                goodwins = np.array(cdf.iloc[0].goodwins,dtype=bool)
                ind = np.where(np.char.replace(DN.iloc[0].Files,'.mat','_' + level + '_spectra.mat')==cfile)[0][0]
                DN.iloc[0].DayNoise[ind].power = power
                DN.iloc[0].DayNoise[ind].cross = cross
                DN.iloc[0].DayNoise[ind].rotation = rotation
                DN.iloc[0].DayNoise[ind].f = cdf.iloc[0].f
                DN.iloc[0].DayNoise[ind].goodwins = goodwins
        SN = StaNoise(daylist=DN.iloc[0].DayNoise)
        SN.init()
        SN.key = 'ML-ATaCR: ' + SN.key
        return SN
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def loadpickles(path):
        if not isinstance(path,Path):
                path = Path(path)
        py_files = list(path.glob('*.pkl'))
        if len(py_files)==0:
                raise Exception('Folder contains no .pkl files')
        out = pd.DataFrame({'Output':[ [] for _ in range(len(py_files)) ],'File':[ [] for _ in range(len(py_files)) ]})
        for i,f in enumerate(py_files):
                file = open(f, 'rb')
                pydata = pkl.load(file)
                out.iloc[i]['Output'] = pydata
                out.iloc[i]['File'] = f.name
                file.close()
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def randomdays(TStart,TEnd,seed='MESSI_22FIFA_WORLD_CUP',days=10,dateformat = '%Y-%m-%d %H:%M:%S'):
        TStart = datetime.datetime(year=TStart.year,month=TStart.month,day=TStart.day)
        TEnd = datetime.datetime(year=TEnd.year,month=TEnd.month,day=TEnd.day)
        Starts = [TStart + datetime.timedelta(days=j) for j in range((TEnd - TStart).days)]
        Ends = [TStart + datetime.timedelta(days=j+1) for j in range((TEnd - TStart).days)]
        if isinstance(seed,str):
                seed = int(''.join([str(ord(a)) for a in seed]))
        else:
                seed = int(str(seed).replace('.',''))
        rng = np.random.default_rng(seed=seed)
        idx = rng.choice([i for i in range(len(Starts))], size=days, replace=False)
        Starts = np.array(Starts)[idx].tolist()
        Ends = np.array(Ends)[idx].tolist()
        Starts = [str(UTCDateTime.strptime(str(a),dateformat)) for a in Starts]
        Ends = [str(UTCDateTime.strptime(str(a),dateformat)) for a in Ends]
        return Starts,Ends
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def fv(func):
        code=func.__code__
        variables=code.co_varnames[:code.co_argcount]
        return list(variables)

def build_staquery(d=None,chan='H',ATaCR_Parent=None,staquery_output='./sta_query.pkl'):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        out = dict()
        if not isinstance(staquery_output,Path):
                if (ATaCR_Parent is None) & (staquery_output is not None):
                        staquery_output = os.getcwd() + staquery_output.replace('./','/')
                else:
                        if staquery_output is not None:
                                os.chdir(ATaCR_Parent)
                                staquery_output = str(ATaCR_Parent) + '/' + staquery_output.replace('./','')
        else:
                if staquery_output is not None:
                        # staquery_output.parent.mkdir(parents=True,exist_ok=True)
                        os.chdir(str(staquery_output.parent))

        for csta in d.iloc:
                net = csta.Network
                sta = csta.Station
                key = net + '.' + sta
                csta_dict = {'station':sta, 'network':net, 'altnet':[],'channel':chan, 'location':['--'], 'latitude':csta['Latitude'], 'longitude':csta['Longitude'], 'elevation':-csta['StaDepth']/1000, 'startdate':UTCDateTime(csta.Start), 'enddate':UTCDateTime(csta.End), 'polarity':1.0, 'azcorr':0.0, 'status':'open'}
                out[key] = csta_dict
        output = pd.DataFrame.from_dict(out)
        if staquery_output is not None:
                output.to_pickle(staquery_output)
                print('Station Query File Written to: ' + str(staquery_output))
        else:
                return output
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DownloadEvents(catalog=None,ATaCR_Parent=None,netsta_names=None,day_mode=False,Minmag=6.3,Maxmag=6.7,limit=1000,pre_event_min_aperture=1,logoutput_subfolder='',log_prefix = '',staquery_output = './sta_query.pkl',chan='H',event_window=7200,channels='Z,P,12',ovr=False):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        # sys.stdout.flush()
        dirs=lt.io.dir_libraries()
        logfilename = '_Step_2_7_EventDownload_logfile.log'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix
        else:
                logoutput = '_Step_2_7_EventDownload_logfile.log'
        datafolder = './EVENTS/raw'
        dateformat = '%Y.%j.%H.%M'
        print('----Begin Event Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                if os.path.isdir(datafolder + staname)==False:
                        os.system('mkdir ' + datafolder + '/' + staname)
                ObsQA.TOOLS.io.build_staquery(d=catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                log_fout = logoutput  + logfilename
                # log_fout
                logoutput_subfolder = Path(logoutput_subfolder)
                logoutput_subfolder.mkdir(parents=True,exist_ok=True)
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        ev = csta.Events[j]
                        print('---'*20)
                        print(staname + '| S:[' +str(i+1) + '/' + str(len(catalog)) + '] | E:[' + str(j+1) + '/' + str(len(csta.Events)) + '] | ' + str(ev.Name),flush=True)
                        print('---'*20)
                        EventStart = UTCDateTime.strptime(ev.Name,dateformat)
                        EventEnd = UTCDateTime.strptime(ev.Name,dateformat) + 60
                        print('Event ' + str(j+1) + '/' + str(len(csta.Events)) + ' | Start -> End : ' + str(EventStart) + ' -> ' + str(EventEnd))
                        if ATaCR_Parent is not None:
                                os.chdir(ATaCR_Parent)
                        # args = [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]
                        args = [str(staquery_output),'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--window={}'.format(event_window),'--channels={}'.format(channels)]
                        if ovr:args.append('--overwrite')
                        if day_mode:args.append(f'--pre_window={22*3600}');args.append(f'--eventpath={dirs.Events_HPS}')
                        # [print(ii) for ii in [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]]
                        # with open(log_fout, 'w') as sys.stdout:
                        atacr_download_event.main(atacr_download_event.get_event_arguments(args))
        print(' ')
        print('----Event Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DayNoiseWhileLoop(catalog=None,NoiseFolder=None,ATaCR_Parent=None,days=10,attempts=50,seed='MESSI_22FIFA_WORLD_CUP!',staquery_output = './sta_query.pkl',ovr=False,channels='Z,P,12'):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        ndays = lambda fold: min([len(list(fold.glob(f'*{c}.SAC'))) for c in ['Z','1','2','H']])
        NoiseFolder = Path(NoiseFolder)
        cat = catalog.copy()
        from modules import load_pickle
        stafolders = [NoiseFolder / s for s in cat.StaName]
        load_sta_atacrnoise = lambda stanm: load_pickle(list((ATaCR_Parent/'AVG_STA'/Station.StaName).glob('*.pkl'))[0])
        # gooddays = load_sta_atacrnoise(stanm).gooddays
        if isinstance(seed,str):seed = int(''.join([str(ord(a)) for a in seed]))
        else:seed = int(str(seed).replace('.',''))
        for ista,cstaf in enumerate(stafolders):
                print('--<>--'*20)
                nfiles = len(list(cstaf.glob('*Z.SAC')))
                nstart = nfiles
                print('Folder: [' + cstaf.name + '] (' + str(ista) + '/' + str(len(stafolders)) + ')'  + '  ||  Current day count: ' + str(nfiles))
                queries = 0
                if nfiles<days:
                        while (queries<attempts) and (nfiles<days):
                                queries+=1
                                seed*=2
                                query_day_len = days - nfiles
                                icatalog = cat[cat.StaName==cstaf.name]
                                Station=icatalog.iloc[0]
                                print('__'*50)
                                print(' ')
                                print(f'Attempt:{str(queries)} , Days needed:{str(query_day_len)}')
                                print(' ')
                                print('__'*50)

                                # Random days
                                baddays = sum([not i for i in load_sta_atacrnoise(Station.StaName).gooddays])
                                nfiles=ndays(cstaf)-baddays
                                nget=days-nfiles
                                if nfiles>0:
                                        Starts,Ends = ObsQA.TOOLS.io.randomdays(Station.Start,Station.End,seed=seed,days=nget)
                                        Addendum = [Station.Start + datetime.timedelta(days=90) + datetime.timedelta(days=int((Station.End - Station.Start + datetime.timedelta(days=90)).days/5))*i for i in range(4)]
                                        Starts.extend([d.strftime('%Y-%m-%dT00:00:00.000000Z') for d in Addendum])
                                        Starts = list(np.unique(Starts))
                                        Ends = [(UTCDateTime(d) + 3600*24).strftime('%Y-%m-%dT00:00:00.000000Z') for d in Starts]
                                        run_atacr_daynoise(Station,Starts,Ends,staquery_output=staquery_output,event_mode=False,ovr=False)
                                        # ObsQA.TOOLS.io.DownloadDayNoise(icatalog, randomloop=False, ATaCR_Parent = NoiseFolder.parent,log_prefix=icatalog.StaName.iloc[0],days=query_day_len,seed=seed)
                                nfiles=ndays(cstaf)-baddays
                                if nfiles<days:print(' || Day requirements not satisfied. Attempting new seed. ||')
                                else:print('[' + str(ista) + '/' + str(len(stafolders)) + '][DAY REQUIREMENTS SATISFIED] | Days gained: ' + str(nfiles-nstart))
                else:print('[' + str(ista) + '/' + str(len(stafolders)) + '][DAY REQUIREMENTS ALREADY SATISFIED] | Skipping.')
                if (nfiles<days) & (queries>=attempts):print('Attempt ' + str(attempts) + '/' + str(attempts) + '. Maximum number of attempts exceeded. Skipping station.')
        failedattempts = [(cstaf.parts[-1], len([f for f in cstaf.glob('*Z.SAC')])) for cstaf in stafolders if len([f for f in cstaf.glob('*Z.SAC')])<days]        
        print('Complete')
        if len(failedattempts)>0:
                print('Stations that still do not meet requirements:')
                print(failedattempts)
        else:print('All station folders now meet day requirements')
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DownloadDayNoise(catalog=None,days=[],message=None,randomloop=True,end_delta=24,event_mode=False,seed='MESSI_22FIFA_WORLD_CUP!',ATaCR_Parent=None,netsta_names=None,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',channels='Z,P,12',ovr=False):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        logfilename = log_prefix + '_Step_3_7_NoiseDownload_logfile.log'
        dateformat = '%Y.%j.%H.%M'
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        datafolder = Path('Data')
        if logoutput_subfolder is not None:
                logfolder = Path(logoutput_subfolder)
                logfolder.mkdir(exist_ok=True,parents=True)
        else:
                logfolder = datafolder
        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i,Station in enumerate(catalog.iloc):
                if event_mode:
                        Ends = [UTCDateTime.strptime(e.Name,dateformat) + 3600*2 for e in Station.Events]
                        Starts = [e-(3600*24) for e in Ends]
                elif isinstance(days,int):
                        # Random days
                        Starts,Ends = ObsQA.TOOLS.io.randomdays(Station.Start,Station.End,seed=seed,days=days)
                        Addendum = [Station.Start + datetime.timedelta(days=90) + datetime.timedelta(days=int((Station.End - Station.Start + datetime.timedelta(days=90)).days/5))*i for i in range(4)]
                        Starts.extend([d.strftime('%Y-%m-%dT00:00:00.000000Z') for d in Addendum])
                        Starts = list(np.unique(Starts))
                        Ends = [(UTCDateTime(d) + 3600*24).strftime('%Y-%m-%dT00:00:00.000000Z') for d in Starts]
                else:
                        # Specific days
                        Starts = days
                        Ends = [datetime.datetime(year=s.year,month=s.month,day=s.day) + datetime.timedelta(hours=end_delta) for s in Starts]
                # Log file stdout
                staname = Station.StaName
                sta_folder = Path('Data') /'raw' / staname
                
                # sta_folder.mkdir(exist_ok=True,parents=True)
                if logoutput_subfolder is not None:
                        log_fout = logfolder / logfilename
                else:
                        log_fout = sta_folder / logfilename
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)

                ObsQA.TOOLS.io.build_staquery(d=Station.to_frame().T,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)

                if randomloop & isinstance(days,int) & (not event_mode):
                        NoiseFolder = Path(ATaCR_Parent) / 'Data' / 'raw'
                        ObsQA.TOOLS.io.DayNoiseWhileLoop(Station.to_frame().T,NoiseFolder,ATaCR_Parent,days=days,attempts=100,ovr=ovr,channels=channels)
                else:
                        run_atacr_daynoise(Station,Starts,Ends,staquery_output=staquery_output,event_mode=event_mode,ovr=ovr)
        print(' ')
        print('----Noise Download Complete----')
        # sys.stdout = original

def run_atacr_daynoise(Station,Starts,Ends,channels='Z,P,12',staquery_output='./sta_query.pkl',event_mode=False,ovr=False,message=None):
        for j,(NoiseStart,NoiseEnd) in enumerate(zip(Starts,Ends)):
                print('<||>'*30)
                # print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Day ' + str(j+1) + '/' + str(len(Starts)),flush=True)
                if event_mode:
                        print('____'*10);print('>>'*15 +'|| 24-hr EVENT MODE ||'+'<<'*15)
                        args = [str(staquery_output),'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd),'--channels={}'.format(channels),'--evn']
                else:
                        args = [str(staquery_output),'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd),'--channels={}'.format(channels)]
                if ovr:args.append('--overwrite')
                status = f'|| [STATUS] | {Station.StaName} | Day:{j+1}/{len(Starts)} | '
                if message:status+=message
                print(status)
                print('____'*10)
                # dateformat = '%Y.%j.%H.%M'
                # jkljk = hjkhkj
                if not ovr:
                       fls=len(list((lt.io.dir_libraries().Noise/'raw'/Station.StaName).rglob(f'{UTCDateTime(NoiseStart).strftime('%Y.%j')}*.SAC')))
                       if fls>=len(''.join(channels.split(','))):continue
                fls=len(list((lt.io.dir_libraries().Noise/'raw'/'_quarantine'/Station.StaName).rglob(f'{UTCDateTime(NoiseStart).strftime('%Y.%j')}*.SAC')))
                if fls>0:continue
                atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
                clear_output(wait=False);os.system('cls' if os.name == 'nt' else 'clear')

#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DownloadEVNoise(catalog=None,ATaCR_Parent=None,netsta_names=None,pre_event_day_aperture=10,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H'):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        # sys.stdout.flush()
        logfilename = log_prefix + '_Step_3_7_NoiseDownload_logfile.log'
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        datafolder = Path('./Data/')
        if logoutput_subfolder is not None:
                logfolder = Path(logoutput_subfolder)
                logfolder.mkdir(exist_ok=True,parents=True)
        else:
                logfolder = datafolder
        dateformat = '%Y.%j.%H.%M'

        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                sta_folder = Path(datafolder / staname)
                sta_folder.mkdir(exist_ok=True,parents=True)
                ObsQA.TOOLS.io.build_staquery(d=catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                if logoutput_subfolder is not None:
                        log_fout = logfolder / logfilename
                else:
                        log_fout = sta_folder / logfilename
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Event ' + str(j+1) + '/' + str(len(csta.Events)),flush=True)
                        ev = csta.Events[j]
                        NoiseStart = UTCDateTime.strptime(ev,dateformat) - datetime.timedelta(days=pre_event_day_aperture)
                        NoiseStart = NoiseStart - datetime.timedelta(hours = NoiseStart.hour, minutes = NoiseStart.minute, seconds=NoiseStart.second) #rounds down to the nearest day
                        NoiseEnd = UTCDateTime.strptime(ev,dateformat)
                        NoiseEnd = NoiseEnd - datetime.timedelta(hours = NoiseEnd.hour, minutes = NoiseEnd.minute, seconds=NoiseEnd.second) #rounds down to the nearest day
                        args = [str(staquery_output),'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd)]
                        atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
        print(' ')
        print('----Noise Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DailySpectra(catalog=None,ATaCR_Parent=None,netsta_names=None,extra_flags = '--figQC --figAverage --figCoh --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1,ovr=False):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        # sys.stdout.flush()
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_4_7_QCSpectra_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_4_7_QCSpectra_logfile.log'
                log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [str(staquery_output)]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        if ovr:args.append('--overwrite')
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        # print('----Begin Daily Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        args = atacr_daily_spectra.get_dailyspec_arguments(args)
        if not fork:
                if ATaCR_Parent is not None:os.chdir(ATaCR_Parent)
                atacr_daily_spectra.main(args)
        else:
                # if __name__ == '__main__':
                
                TaskMaster()
                # if __name__ == '__main__':
                #         # Specify the log file and setup logging in the main process
                #         setup_logging(str(log_fout))  
                #         # Create and start the multiprocessing process
                #         p = multiprocessing.Process(target=worker, args=(atacr_daily_spectra.main, args, str(log_fout)))
                #         print("Starting child process.")
                #         p.start()
                #         # Wait for the child process to finish
                #         p.join()
                #         if p.is_alive():print("post-Child process is still running unexpectedly.")
                #         else:print("post-Child process has terminated as expected.")
                #         print("Child process has completed.")
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def CleanSpectra(catalog=None,ATaCR_Parent=None,netsta_names=None,extra_flags = '--figQC --figAverage --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1,ovr=False):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        # sys.stdout.flush()
        #  --figQC --figAverage --save-fig
        # datafolder = './Data/'
        # if logoutput_subfolder is not None:
                # logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_5_7_CleanSpectra_logfile.log'
                # log_fout = logoutput
        # else:
                # logoutput = '_Step_5_7_CleanSpectra_logfile.log'
                # log_fout = datafolder + logoutput
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # log_fout = Path(log_fout)
        # log_fout.parent.mkdir(parents=True,exist_ok=True)
        # sys.stdout = open(str(log_fout),'w+')
        # logging.basicConfig(stream=log_fout)
        # print('----Begin Clean Spectra----')
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [str(staquery_output)]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        if ovr:args.append('--overwrite')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_clean_spectra.get_cleanspec_arguments(args)
        # atacr_clean_spectra.main(args)
        if not fork:
                atacr_clean_spectra.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        # Start the load operations and mark each future with its URL
                        print('----Begin Clean Spectra----')
                        future = process_executor.submit(atacr_clean_spectra.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Clean Spectra Complete----')
                # sys.stdout = original
        # print('----Clean Spectra Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def TransferFunctions(catalog=None,ATaCR_Parent=None,netsta_names=None,extra_flags = '--figTF --save-fig'
        ,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1,ovr=False):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        # sys.stdout.flush()
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_6_7_CalcTFs_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_6_7_CalcTFs_logfile.log'
                log_fout = datafolder + logoutput
        args = [str(staquery_output)]
        extra_flags = extra_flags # + ' --taper ' + str(taper_mode)
        [args.append(flg) for flg in extra_flags.split()]
        if ovr:args.append('--overwrite')
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_transfer_functions.get_transfer_arguments(args)
        # atacr_transfer_functions.main(args)
        if not fork:
                atacr_transfer_functions.main(args)
        else:
                print('----Begin Transfer Functions----')
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        future = process_executor.submit(atacr_transfer_functions.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Transfer Functions Complete----')
                # sys.stdout = original
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def CorrectEvents(catalog=None,ATaCR_Parent=None,netsta_names=None,extra_flags = '--figRaw --figClean --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1,ovr=False):
        # [exec(f'{k}=args.{k}') for k in list(args.keys())]
        datafolder = './Data/'
        # staquery_output = str(staquery_output)
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_7_7_CorrectEvents_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_7_7_CorrectEvents_logfile.log'
                log_fout = datafolder + logoutput
        args = [str(staquery_output)]
        [args.append(flg) for flg in extra_flags.split()]
        if ovr:args.append('--overwrite')
        # with open(log_fout, 'w') as sys.stdout:
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_correct_event.get_correct_arguments(args)
        if not fork:
                atacr_correct_event.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        future = process_executor.submit(atacr_correct_event.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Correct Events Complete----')
                # sys.stdout = original
        # atacr_correct_event.main(args)
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def Run_ATaCR(args):
        # message=None,ATaCR_Parent = None, STEPS=[1,2,3,4,5,6,7], netsta_names=None, chan='H', Minmag=6.3, Maxmag=6.7,
        # limit=1000, pre_event_min_aperture=1, pre_event_day_aperture=30,
        # dailyspectra_flags='--figQC --figAverage --figCoh --save-fig',
        # cleanspectra_flags='--figQC --figAverage --figCoh --figCross --save-fig',
        # tf_flags='--figTF --save-fig', 
        # correctevents_flags='--figRaw --figClean --save-fig',
        # logoutput_subfolder=None,log_prefix = '',staquery_output = './sta_query.pkl',fork=True,days=10,
        # seed='MESSI_22FIFA_WORLD_CUP!',max_workers=1,event_mode=False,event_dt=None,
        # event_window=7200,channels='Z,P,12',ovr=False):
        # _=unpack(args,fv(CorrectEvents)[1:])
        dirs = lt.io.dir_libraries()

        if 1 in args.STEP:
                print('Step 1/7 - BEGIN: Station Metadata')
                curargs = dict({k:args[k] for k in fv(build_staquery) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.build_staquery(**curargs)
                print('Step 1/7 - COMPLETE: Station Metadata')
        if 2 in args.STEP:
                print('Step 2/7 - BEGIN: Download Event Data')
                args.logoutput_subfolder = str(dirs.Logs) + '/2_7'
                curargs = dict({k:args[k] for k in fv(DownloadEvents) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.DownloadEvents(**curargs)
                print('Step 2/7 - COMPLETE: Download Event Data')
        if 3 in args.STEP:
                print('Step 3/7 - BEGIN: Download Day Data')
                args.logoutput_subfolder = str(dirs.Logs) + '/3_7'
                curargs = dict({k:args[k] for k in fv(DownloadDayNoise) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.DownloadDayNoise(**curargs)
                print('Step 3/7 - COMPLETE: Download Day Data')
        if 4 in args.STEP:
                print('Step 4/7 - BEGIN: Quality Control Noise Data')
                args.logoutput_subfolder = str(dirs.Logs) + '/4_7'
                curargs = dict({k:args[k] for k in fv(DailySpectra) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.DailySpectra(**curargs)
                print('Step 4/7 - COMPLETE: Quality Control Noise Data')
        if 5 in args.STEP:
                print('Step 5/7 - BEGIN: Spectral Average of Noise Data')
                args.logoutput_subfolder = str(dirs.Logs) + '/5_7'
                curargs = dict({k:args[k] for k in fv(CleanSpectra) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.CleanSpectra(**curargs)
                print('Step 5/7 - COMPLETE: Spectral Average of Noise Data')
        if 6 in args.STEP:
                print('Step 6/7 - BEGIN: Calculate Transfer Functions')
                args.logoutput_subfolder = str(dirs.Logs) + '/6_7'
                curargs = dict({k:args[k] for k in fv(TransferFunctions) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.TransferFunctions(**curargs)
                print('Step 6/7 - COMPLETE: Calculate Transfer Functions')
        if 7 in args.STEP:
                print('Step 7/7 - BEGIN: Correct Event Data')
                args.logoutput_subfolder = str(dirs.Logs) + '/7_7'
                curargs = dict({k:args[k] for k in fv(CorrectEvents) if np.isin(k,list(args.keys()))})
                ObsQA.TOOLS.io.CorrectEvents(**curargs)
                print('Step 7/7 - COMPLETE: Correct Event Data')
        plt.close('all')


def get_data(datapath, tstart, tend,seismic_pre_filt=[0.001, 0.002, 45.0, 50.0], pressure_pre_filt=[0.001, 0.002, 45.0, 50.0],seismic_units="DISP",pressure_units="DEF",pressure_water_level=None,seismic_water_level=60):
    """
    Function to grab all available noise data given a path and data time range

    Parameters
    ----------
    datapath : str
        Path to noise data folder
    tstart : :class:`~obspy.class.UTCDateTime`
        Start time for query
    tend : :class:`~obspy.class.UTCDateTime`
        End time for query

    Returns
    -------
    tr1, tr2, trZ, trP : :class:`~obspy.core.Trace` object
        Corresponding trace objects for components H1, H2, HZ and HP. Returns
        empty traces for missing components.

    """

    # Define empty streams
    trN1 = Stream()
    trN2 = Stream()
    trNZ = Stream()
    trNP = Stream()

    # Time iterator
    t1 = tstart
    
    # Cycle through each day within time range
    while t1 < tend:

        # Time stamp used in file name
        tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

        # Cycle through directory and load files
        p = datapath.glob('*.*')
        files = [x for x in p if x.is_file()]
        with warnings.catch_warnings():
                warnings.filterwarnings("ignore");
                for file in files:
                        inv = read_inventory(Path(str(file)).parent / '*_inventory.xml')
                        if fnmatch.fnmatch(str(file), '*' + tstamp + '*1.SAC'):
                                tr = read(str(file))
                                remargs=dict(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
                                tr.remove_response(**remargs)
                                trN1.append(tr[0])
                        elif fnmatch.fnmatch(str(file), '*' + tstamp + '*2.SAC'):
                                tr = read(str(file))
                                remargs=dict(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
                                tr.remove_response(**remargs)
                                trN2.append(tr[0])
                        elif fnmatch.fnmatch(str(file), '*' + tstamp + '*Z.SAC'):
                                tr = read(str(file))
                                remargs=dict(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
                                tr.remove_response(**remargs)
                                trNZ.append(tr[0])
                        elif fnmatch.fnmatch(str(file), '*' + tstamp + '*H.SAC'):
                                tr = read(str(file))
                                remargs=dict(inventory=inv,pre_filt=pressure_pre_filt,output=pressure_units,water_level=pressure_water_level,hide_sensitivity_mismatch_warning=True)
                                tr.remove_response(**remargs)
                                trNP.append(tr[0])

        # Increase increment
        t1 += 3600.*24.

    # Fill with empty traces if components are not found
    ntr = len(trNZ)
    if not trN1 and not trN2:
        for i in range(ntr):
            trN1.append(Trace())
            trN2.append(Trace())
    if not trNP:
        for i in range(ntr):
            trNP.append(Trace())

    if ntr > 0:
        # Check that all sampling rates are equal - otherwise resample
        if trNZ[0].stats.sampling_rate != trNP[0].stats.sampling_rate:

            # These checks assume that all seismic data have the same sampling
            if trNZ[0].stats.sampling_rate < trNP[0].stats.sampling_rate:
                trNP.resample(trNZ[0].stats.sampling_rate, no_filter=False)
            else:
                trNZ.resample(trNP[0].stats.sampling_rate, no_filter=False)
                if trN1:
                    trN1.resample(trNP[0].stats.sampling_rate, no_filter=False)
                if trN2:
                    trN2.resample(trNP[0].stats.sampling_rate, no_filter=False)
    return trN1, trN2, trNZ, trNP

def load_sac(file,rmresp=False,inv=[],
        pressure_pre_filt=[0.001, 0.002, 45.0, 50.0],
        seismic_pre_filt=[0.001, 0.002, 45.0, 50.0],
        seismic_units="DISP",pressure_units="DEF",
        pressure_water_level=None,seismic_water_level=60):
        # ----------------------
    if rmresp:inv=read_inventory(Path(str(file)).parent / '*_inventory.xml')
    try:
        with warnings.catch_warnings():
                warnings.filterwarnings("ignore");
                if fnmatch.fnmatch(str(file),'*1.SAC') or fnmatch.fnmatch(str(file),'*2.SAC') or fnmatch.fnmatch(str(file),'*Z.SAC') :
                        tr = read(str(file))
                        if rmresp:tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
                elif fnmatch.fnmatch(str(file),'*H.SAC'):
                        tr = read(str(file))
                        if rmresp:tr.remove_response(inventory=inv,pre_filt=pressure_pre_filt,output=pressure_units,water_level=pressure_water_level,hide_sensitivity_mismatch_warning=True)
                return tr[0]
    except:
        print('WARNING:Load Error')
        return Trace() #,inv

def get_event(eventpath, tstart=None, tend=None,evlist=None,seismic_pre_filt=[0.001, 0.002, 45.0, 50.0], pressure_pre_filt=[0.001, 0.002, 45.0, 50.0],seismic_units="DISP",pressure_units="DEF",pressure_water_level=None,seismic_water_level=60):
    """
    Function to grab all available earthquake data given a path and data time
    range
    Parameters
    ----------
    eventpath : str
        Path to earthquake data folder
    tstart : :class:`~obspy.class.UTCDateTime`
        Start time for query
    tend : :class:`~obspy.class.UTCDateTime`
        End time for query
    Returns
    -------
    tr1, tr2, trZ, trP : :class:`~obspy.core.Trace` object
        Corresponding trace objects for components H1, H2, HZ and HP. Returns
        empty traces for missing components.
    """
    # Define empty streams
    tr1,tr2,trP,trZ = Stream(),Stream(),Stream(),Stream()
    # Cycle over all available files in time range
    if evlist is None:
        # Find out how many events from Z.SAC files
        eventfiles = list(eventpath.glob('*Z.SAC'))
        if not eventfiles:
            raise(Exception("No event found in folder "+str(eventpath)))
        # Extract events from time stamps
        prefix = [file.name.split('.') for file in eventfiles]
        evstamp = [p[0]+'.'+p[1]+'.'+p[2]+'.'+p[3]+'.' for p in prefix]
        evDateTime = [UTCDateTime(p[0]+'-'+p[1]+'T'+p[2]+":"+p[3]) for p in prefix]
        for event, tstamp in zip(evDateTime, evstamp):
            if event >= tstart and event <= tend:
                warnings.filterwarnings("ignore")
                p = list(eventpath.glob('*Z.SAC'))
                files = [x for x in p if x.is_file()]
                for file in files:
                    tr1=load_sac(Path(str(file).replace('Z.SAC','*1.SAC')))
                    tr2=load_sac(Path(str(file).replace('Z.SAC','*2.SAC')))
                    trZ=load_sac(Path(str(file).replace('Z.SAC','*Z.SAC')))
                    trP=load_sac(Path(str(file).replace('Z.SAC','*H.SAC')))
                    tr1, tr2, trZ, trP = check_fs(tr1,tr2,trZ,trP)
                    yield tr1[0], tr2[0], trZ[0], trP[0]
    else:
        for ev in evlist:
            tr1=load_sac(Path(str(eventpath / ev)+'*1.SAC'))
            tr2=load_sac(Path(str(eventpath / ev)+'*2.SAC'))
            trZ=load_sac(Path(str(eventpath / ev)+'*Z.SAC'))
            trP=load_sac(Path(str(eventpath / ev)+'*H.SAC'))
            tr1, tr2, trZ, trP = check_fs(tr1,tr2,trZ,trP)
            yield tr1[0], tr2[0], trZ[0], trP[0]

def check_fs(tr1,tr2,trZ,trP):
    if trZ[0].stats.sampling_rate != trP[0].stats.sampling_rate:
        if trZ[0].stats.sampling_rate < trP[0].stats.sampling_rate:
            trP.resample(trZ[0].stats.sampling_rate, no_filter=False)
    else:
        trZ.resample(trP[0].stats.sampling_rate, no_filter=False)
        if tr1:
            tr1.resample(trP[0].stats.sampling_rate, no_filter=False)
        if tr2:
            tr2.resample(trP[0].stats.sampling_rate, no_filter=False)
    return tr1,tr2,trZ,trP

def distance(sta,ev,unit='deg'):
    origins=ev.origins[0]
    stalla,evlla=[sta.Latitude,sta.Longitude],[origins.latitude,origins.longitude]
    dist=locations2degrees(stalla[0],stalla[1],evlla[0],evlla[1])
    if unit.lower()=='km':dist=degrees2kilometers(dist)
    return dist

def fnotch(d):
        '''The frequency knotch root function described in Crawford et al., 1998.
        depth (d) is in meters. Returned (f) is in Hz.'''
        g = 9.80665
        f = (g/(2*np.pi*d))**0.5
        return f

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
def meta_hist_plot(catalog):
    full_event_catalog=Catalog(unravel([s.Events for s in catalog.iloc]))
    full_event_catalog=Catalog([full_event_catalog[i] for i in np.unique([e.Name for e in full_event_catalog],return_index=True)[1]])
    ev=full_event_catalog[0]
    MetaHold = dict()
    MetaHold['Event_depths_km'] = [ev.origins[0].depth/1000 for ev in full_event_catalog]
    MetaHold['Event_Magnitude_M'] = [ev.magnitudes[0].mag for ev in full_event_catalog]
    MetaHold['Event_Distances'] = unravel([[distance(catalog[catalog.StaName==s].iloc[0],ev) for s in ev.Stations if np.isin(s,catalog.StaName)] for ev in full_event_catalog])
    MetaHold['Stations_per_Event'] = [len(ev.Stations) for ev in full_event_catalog]
    MetaHold['Deployment_Depths'] = catalog.StaDepth
    MetaHold['Pressure_Gauges'] = catalog.Pressure_Gauge
    MetaHold['Experiment'] = '('+catalog.Network+')'+catalog.Experiment
    MetaHold['Noise_Notch_period'] = 1/fnotch(catalog.StaDepth)
    MetaHold['Regions'] = catalog.Environment
    deployment_metas = [
    'Distance_from_Land_km','Distance_to_Plate_Boundary_km',
    'Sediment_Thickness_m','Surface_Current_ms','Crustal_Age_Myr',
    'Deployment_Length_days','Instrument_Design','Seismometer']
    _ = [MetaHold.update({m:unravel([[s[m] for s in catalog.Deployment]])}) for m in deployment_metas]
    # labels=list(MetaHold.keys())
    labels=['Event_depths_km',
    'Event_Magnitude_M',
    'Event_Distances',
    'Stations_per_Event',
    'Deployment_Depths',
    'Pressure_Gauges',
    'Experiment',
    'Noise_Notch_period',
    'Distance_from_Land_km',
    'Distance_to_Plate_Boundary_km',
    'Sediment_Thickness_m',
    'Surface_Current_ms',
    'Crustal_Age_Myr',
    'Instrument_Design',
    'Seismometer']
    n_ax=len(list(MetaHold.keys()))
    fig,axes=plt.subplots(ncols=3,nrows=int(np.ceil(n_ax/3)),figsize=(20,8),layout='constrained',squeeze=True)
    axes=axes.reshape(-1)
    bins={'Event_depths_km':[0,100,200,300,400,500,600,700],
    'Stations_per_Event':np.arange(1,40,2),
    'Event_Magnitude_M':[6,6.5,7.0,7.5,7.9],
    'Event_Distances':np.arange(0,180,30),
    #  'Deployment_Depths':np.arange(0,6000,1000),
    'Deployment_Depths':[0,1500,3000,4000,6000],
    'Pressure_Gauges':2,
    'Noise_Notch_period':11,
    'Experiment':11,
    'Regions':10,
    'Distance_from_Land_km':None,
    'Distance_to_Plate_Boundary_km':None,
    'Sediment_Thickness_m':None,
    'Surface_Current_ms':None,
    'Crustal_Age_Myr':np.arange(0,160,30),
    'Deployment_Length_days':None,
    'Instrument_Design':8,
    'Seismometer':3}

    txi=0
    for axi,(ax,label) in enumerate(zip(axes,labels)):
        if 'Event' in label:ylabel='Events'
        if 'Station' in label:ylabel='Stations'
        else:ylabel='Stations'
        xlabel = label.replace('_period',', period').replace('_km',', km').replace('_m',', m').replace('_ms',', m/s').replace('_Myr',', Myr').replace('_days',', days').replace('_',' ')
        xlabel= xlabel[0].upper() +xlabel[1:].lower()
        xlabel = xlabel.replace('Deployment depths','Deployment depths, m')
        if xlabel=='Event magnitude m':xlabel=('Event magnitude, M')
        if xlabel=='Event distances':xlabel=f'{xlabel}, °'
        isnumber= not isinstance(MetaHold[label][0],str)
        if not isnumber:bn=len(np.unique(MetaHold[label]))
        else:bn=bins[label]
        n,b,p=ax.hist(MetaHold[label],bins=bn,edgecolor='k',rwidth=0.9,facecolor='darkgrey')
        ax.set_xticks(b[:-1] + np.diff(b)[0]/2)
        if isnumber:
            if label=='Surface_Current_ms':ax.set_xticklabels([str(np.round(float(f),2)) for f in [tx._text for tx in ax.get_xticklabels()]])
            elif not label=='Event_Magnitude_M':ax.set_xticklabels([str(int(float(f))) for f in [tx._text for tx in ax.get_xticklabels()]])
        else:ax.set_xticklabels(ax.get_xticklabels(),rotation=[0,32,15,0,0,0][txi]);txi+=1
        if label=='Deployment_Depths':aa=[n,b,p,ax.get_xticks()];ax.set_xticks([750,2250,3500,5000]);ax.set_xticklabels([750,2250,3500,5000])
        ax.set_xlabel(xlabel,fontweight='bold')
        if np.any(np.array(['Event_Magnitude_M','Event_Distances','Event_depths_km'])==label):ax.set_ylabel('Events')
        else:ax.set_ylabel(ylabel)
    [ax.axis('off') for axi,ax in enumerate(axes) if (axi+1)>len(labels)]
    clear_output(wait=False);os.system('cls' if os.name == 'nt' else 'clear')
    return fig


def get_reports(sta=None,
    Archive=None,
    methods=['ATaCR','NoiseCut'],evna=None):
    methods = np.array(methods)
    Report=AttribDict()
    if Archive is None:dirs=lt.io.dir_libraries();Archive=dirs.Data
    AnalysisFolder = Archive/'Analysis'/'Coherence'
    lp = lambda f:reshapedict(load_pickle(f))
    atacrfile = lambda:AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{tf}.zz.coherence.pkl'
    hpsfile = lambda :AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{c}.coherence.pkl'
    unionevents = lambda a:reduce(np.intersect1d,a)
    methods = [m.replace('NoiseCut','HPS') for m in methods]
    for m in methods:
        if m.lower()=='atacr':
            tfs = ['zp_21','z2_1','zp']
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
        for s in stas:
            if s=='f':continue
            evs=unionevents([Report[m][list(Report[m].keys())[0]][s].events for m in methods])
            ind = np.isin(Report[m][c][s].events,evs)
            for i in Report[m].keys():
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
        Report  = lt.cat.evcoh(Report,evna)
    return Report
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# eof