# from obspy import read
# import numpy as np
# from scipy.signal import stft, detrend
# import warnings
# import fnmatch
# import obspy
# from obspy import read_inventory
# from imports import *
from modules import *
from scipy.signal import coherence,welch
from modules import *
# ----
def get_sac(file,seismic_pre_filt=[0.001, 0.002, 45.0, 50.0], pressure_pre_filt=[0.001, 0.002, 45.0, 50.0],seismic_units="DISP",pressure_units="DEF",pressure_water_level=None,seismic_water_level=60):
    warnings.filterwarnings("ignore")
    pressure = np.any([fnmatch.fnmatch(str(file), '*' + s + '.SAC') for s in ['H']])
    seis = np.any([fnmatch.fnmatch(str(file), '*' + s + '.SAC') for s in ['1','2','Z']])
    inv = read_inventory(list(file.parent.glob('*_inventory.xml'))[0])
    if seis:
        tr = read(str(file))
        warnings.filterwarnings("ignore")
        tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level)
    elif pressure:
        tr = read(str(file))
        warnings.filterwarnings("ignore")
        tr.remove_response(inventory=inv,pre_filt=pressure_pre_filt,output=pressure_units,water_level=pressure_water_level)
    return tr[0]

def pull_cohphadm(stanm,cat,UncorrectedFold,CorrectedFold,
        tf='ZP-21',
        corrected_comp='HZ',orig_comp='HDH',reverse=False):
        # from modules import load_sac
        # from modules import basic_preproc
        correvpath = CorrectedFold
        origevpath = UncorrectedFold
        if 'HPS_Data' in str(CorrectedFold):corrmethod='NoiseCut' 
        elif 'rmresp' in str(CorrectedFold):corrmethod='Original'
        else:corrmethod='ATaCR'
        if 'HPS_Data' in str(UncorrectedFold):origmethod='NoiseCut' 
        elif 'rmresp' in str(UncorrectedFold):origmethod='Original'
        else:origmethod='ATaCR'
        if len(tf)>0:fclip = '.sta.' + tf
        else:fclip = ''
        correvs = [f.name for f in list(correvpath.glob('*' + fclip + f'.{corrected_comp}.SAC'))]
        origevs = [f.name for f in list(origevpath.glob('*' + f'.{orig_comp}.SAC'))]
        evs=np.intersect1d([r.split(f'.{corrected_comp}.SAC')[0].replace(stanm+'.','').replace('.sta.'+tf,'') for r in correvs],[r.split(f'.{orig_comp}.SAC')[0] for r in origevs])
        c=evs.copy();correvs=[f'{stanm}.{i}.{corrected_comp}.SAC' for i in c]
        c=evs.copy();origevs=[f'{i}.{orig_comp}.SAC' for i in c]
                #     events = [c.replace('.sta','').replace('.HZ.SAC','').split(stanm + '.')[-1].split('.' + tf)[0] for c in correvs]
        events = [c.replace('.sta','').replace('.HZ.SAC','').replace('.H1.SAC','').replace('.H2.SAC','').replace(stanm+'.','') for c in correvs]
        if len(tf)>0:
                events=[c.replace('.'+tf,'') for c in events]
                correvs=[c.replace('.HZ.','.sta.'+tf+'.HZ.') for c in correvs]
        cpa_list,ev_list = [],[]

        #data_test: Tests if each event in the catalog was found.
        data_test = [(np.isin(e.Name,events)).tolist() for e in cat.loc[stanm].Events.iloc[0]]
        #data_test_conv: Tests if each event found is also in the catalog.
        data_test_conv = [(np.isin(e,[o.Name for o in cat.loc[stanm].iloc[0].Events])).tolist() for e in events]

        log_warnings=[]
        #I don't want to track whether ZP coherence post-correction is collected
        if (not np.all(data_test)) & (not np.isin('HDH',[orig_comp,corrected_comp])):
                [print(f'{e.Name}: {t}') for t,e in zip(data_test,cat.loc[stanm].iloc[0].Events) if t==False]
                log_warnings=[f'{stanm}:{e.Name}:Catalog event not found' for t,e in zip(data_test,cat.loc[stanm].iloc[0].Events) if t==False]
        #        assert np.all(data_test)
        corr,orig=[],[]
        for conv_test,ev,r,c in zip(data_test_conv,events,origevs,correvs):
                #Only collect coherences for events currently listed in the data volume.
                origst=cat.loc[stanm].iloc[0].Data.Traces(ev,channel=orig_comp,methods=[origmethod])
                corrst=cat.loc[stanm].iloc[0].Data.Traces(ev,channel=corrected_comp,methods=[corrmethod])
                if np.any([origst==None,corrst==None]):msg=f'{stanm}|{ev}| Data missing or corrupt';log_warnings.append(msg);print(msg);continue
                if np.any([len(origst)<1,len(corrst)<1]):msg=f'{stanm}|{ev}| Data missing or corrupt';log_warnings.append(msg);print(msg);continue
                origst=origst.select(location=origmethod)[0]
                corrst=corrst.select(location=corrmethod)[0]
                corr.append(corrst);orig.append(origst)
        corr=Stream(corr);orig=Stream(orig)
        for ev,origst,corrst in zip(events,orig,corr):
                if reverse:cpa = cohphadm(corrst,origst)
                else:cpa = cohphadm(origst,corrst)
                cpa.Event = ev
                cpa_list.append(cpa)
                ev_list.append(ev)
        _=[print(s) for s in log_warnings]
        return ev_list,cpa_list,log_warnings

class DataHold(object):
    def __init__(self):self

class cohphadm(object):
        def __init__(self,A, B=None,overlap=0.2,csd=None,f=None,fs=None,preproc=True):
                self.preproc=preproc
                if B is None:B=A.copy()
                self.overlap = overlap
                self.A = A.copy()
                self.B = B.copy()

                # --Settings--
                self.dt = np.unique([s.stats.delta for s in [self.A,self.B] if isinstance(s,obspy.core.trace.Trace)])[0]
                self.fs = 1/self.dt
                assert self.fs>1, 'Bad samplerate'
                # self.tlen = np.size(self.A.data)/self.fs
                self.depth = abs(self.A.stats.sac.stel)*1000
                self.net = self.A.stats.network
                self.sta = self.A.stats.station
                self.chan = self.A.stats.channel
                self.window = 600 #window length in seconds. Default = 600 (10min)
                self.lowpass=None
                self.percent=.02 #Hanning taper percent. Default = 0.02 (2%)
                self.tlen = 7200 #Trace length. Default 7200s (2hrs)
                self.npts = int(self.window*self.fs)
                self.csd = csd
                self.f = f
                self.nperseg=self.npts
                self.noverlap=self.npts//int(1/self.overlap)
                if self.preproc:self.basic_preproc()

        def basic_preproc(self):
                A,B = self.A,self.B
                seconds = lambda A:A.hour*3600+A.minute*60+A.second
                start = max([B.stats.starttime,A.stats.starttime])
                end = min([B.stats.endtime,A.stats.endtime])
                A.trim(start,end);B.trim(start,end)
                tolerance=1 # seconds
                tend = np.min([A.stats.endtime,B.stats.endtime])
                A.trim(tend-self.tlen,tend)
                B.trim(tend-self.tlen,tend)
                trim_n = min([A.data.shape[0],B.data.shape[0]])
                B.data=B.data[:trim_n];A.data=A.data[:trim_n]
                assert np.allclose(seconds(B.stats.endtime),seconds(A.stats.endtime),tolerance),'End times do not match'
                assert np.allclose(seconds(B.stats.starttime),seconds(A.stats.starttime),tolerance),'Start times do not match'
                for d in [A,B]:
                        # d.detrend('constant');d.detrend('demean')
                        if not self.percent==None:d.taper(self.percent)
                        if not self.lowpass==None:d.filter('lowpass',freq=self.lowpass,zerophase=True);d.detrend('constant');d.detrend('demean')
                self.A,self.B = A,B

        def fnotch(self,d=None,n=1):
                '''The frequency notch root function described in Crawford et al., 1998.
                depth (d) is in meters. Returned (f) is in Hz.'''
                if d==None:d = self.depth
                g = 9.80665
                f = (g/(2*np.pi*d*n))**0.5
                return f
        def COH(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                coh = self._calc_coherence(ab,aa,bb) #(MSC from Bell, 2015)
                # f, coh = coherence(self.A.data,self.B.data, fs=self.fs, nperseg=self.nperseg,noverlap=self.noverlap)
                return f[f>0],coh[f>0]
        def PH(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                ph = self._calc_phase(ab)
                return f[f>0],ph[f>0]
        def ADM(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                adm = self._calc_admittance(ab,bb)
                return f[f>0],adm[f>0]

        def _calc_coherence(self,ab,aa,bb):
                #Magnitude-Squared Coherence (MSC) function (Bell,2015):
                coh = ( np.abs(ab)**2 )  /  ( np.abs(aa) * np.abs(bb) )
                # ------------------------------
                #Complex-Valued Coherence (CVC) (Crawford&Webb,2000):
                # coh = (  (ab)  /  ((aa*bb)**0.5) )
                return np.abs(coh)

        def _calc_phase(self,ab):
                ph = np.angle(ab,deg=True)
                return ph
        def _calc_admittance(self,ab,bb):
                ad = np.abs(np.abs(ab)/bb)
                return ad
        def _window(self,window=None,overlap=None):
                if overlap is None:overlap=self.overlap
                if window is None:window=self.window
                # Points in window
                ws = int(window/self.dt)
                # Number of points to overlap
                ss = int(window*overlap/self.dt)
                # hanning window
                hanning = np.hanning(2*ss)
                wind = np.ones(ws)
                wind[0:ss] = hanning[0:ss]
                wind[-ss:ws] = hanning[ss:ws]
                return wind,ws,ss
        def _powercross(self,a,b):
                f,ab = self._csd_helper(a,b,return_onesided=False) #csd
                f,aa = self._csd_helper(a,a,return_onesided=False) #psd
                f,bb = self._csd_helper(b,b,return_onesided=False) #psd
                return f,ab,aa,bb
        def _csd_helper(self,a,b,window=None,overlap=None,return_onesided=False):
                f,_t,a_ft = self._stft(a,window=window,overlap=overlap,return_onesided=return_onesided) # ( time-windows x frequency-bins )
                f,_t,b_ft = self._stft(b,window=window,overlap=overlap,return_onesided=return_onesided) # ( time-windows x frequency-bins )
                cab = np.mean(self._calc_csd(a_ft,b_ft),axis=0) #Welch method
                return f,cab
        def _calc_csd(self,a_ft,b_ft):
                cab = np.conj(a_ft)*b_ft
                return cab
        def _stft(self,tr,scaling='spectrum',window=None,overlap=None,return_onesided=False):
                wind,ws,ss = self._window(window=window,overlap=overlap)
                _f, _t, ft = stft(tr,self.fs, return_onesided=return_onesided, boundary=None,padded=True, 
                window=wind, nperseg=ws, noverlap=ss,scaling=scaling)
                _f = _f.reshape(-1)
                ft = np.atleast_2d(ft)
                return _f.T,_t.T,ft.T
# ------------------------------------------------      ------------------------------------------------                ------------------------------------------------
# ------------------------------------------------------------------------------------------------              ------------------------------------------------
# ------------------------------------------------      ------------------------------------------------                ------------------------------------------------
class Signal(object):
        def __init__(self,A, B=None,overlap=0.2,csd=None,f=None,fs=None,preproc=True):
                self.preproc=preproc
                if B is None:B=A.copy()
                self.overlap = overlap
                self.A = A.copy()
                self.B = B.copy()
                # --Settings--
                self.dt = np.unique([s.stats.delta for s in [self.A,self.B] if isinstance(s,obspy.core.trace.Trace)])[0]
                self.fs = 1/self.dt
                assert self.fs>1, 'Bad samplerate'
                # self.tlen = np.size(self.A.data)/self.fs
                self.depth = abs(self.A.stats.sac.stel)*1000
                self.net = self.A.stats.network
                self.sta = self.A.stats.station
                self.chan = self.A.stats.channel
                self.window = 600 #window length in seconds. Default = 600 (10min)
                self.lowpass=None
                self.percent=.02 #Hanning taper percent. Default = 0.02 (2%)
                self.tlen = 7200 #Trace length. Default 7200s (2hrs)
                self.npts = int(self.window*self.fs)
                self.csd = csd
                self.f = f
                self.nperseg=self.npts
                self.noverlap=self.npts//int(1/self.overlap)
                if self.preproc:self.basic_preproc()
        def basic_preproc(self):
                A,B = self.A,self.B
                seconds = lambda A:A.hour*3600+A.minute*60+A.second
                start = max([B.stats.starttime,A.stats.starttime])
                end = min([B.stats.endtime,A.stats.endtime])
                A.trim(start,end);B.trim(start,end)
                tolerance=1 # seconds
                tend = np.min([A.stats.endtime,B.stats.endtime])
                A.trim(tend-self.tlen,tend)
                B.trim(tend-self.tlen,tend)
                trim_n = min([A.data.shape[0],B.data.shape[0]])
                B.data=B.data[:trim_n];A.data=A.data[:trim_n]
                assert np.allclose(seconds(B.stats.endtime),seconds(A.stats.endtime),tolerance),'End times do not match'
                assert np.allclose(seconds(B.stats.starttime),seconds(A.stats.starttime),tolerance),'Start times do not match'
                for d in [A,B]:
                        # d.detrend('constant');d.detrend('demean')
                        if not self.percent==None:d.taper(self.percent)
                        if not self.lowpass==None:d.filter('lowpass',freq=self.lowpass,zerophase=True);d.detrend('constant');d.detrend('demean')
                self.A,self.B = A,B
        def fnotch(self,d=None,n=1):
                '''The frequency notch root function described in Crawford et al., 1998.
                depth (d) is in meters. Returned (f) is in Hz.'''
                if d==None:d = self.depth
                g = 9.80665
                f = (g/(2*np.pi*d*n))**0.5
                return f
        def coherence(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                coh = self._calc_coherence(ab,aa,bb) #(MSC from Bell, 2015)
                # f, coh = coherence(self.A.data,self.B.data, fs=self.fs, nperseg=self.nperseg,noverlap=self.noverlap)
                return f[f>0],coh[f>0]
        
        def phase(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                ph = self._calc_phase(ab)
                return f[f>0],ph[f>0]
        def admittance(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                adm = self._calc_admittance(ab,bb)
                return f[f>0],adm[f>0]
        def _calc_coherence(self,ab,aa,bb):
                #Magnitude-Squared Coherence (MSC) function (Bell,2015):
                coh = ( np.abs(ab)**2 )  /  ( np.abs(aa) * np.abs(bb) )
                # ------------------------------
                #Complex-Valued Coherence (CVC) (Crawford&Webb,2000):
                # coh = (  (ab)  /  ((aa*bb)**0.5) )
                return np.abs(coh)
        def _calc_phase(self,ab):
                ph = np.angle(ab,deg=True)
                return ph
        def _calc_admittance(self,ab,bb):
                ad = np.abs(np.abs(ab)/bb)
                return ad
        def _window(self,window=None,overlap=None):
                if overlap is None:overlap=self.overlap
                if window is None:window=self.window
                # Points in window
                ws = int(window/self.dt)
                # Number of points to overlap
                ss = int(window*overlap/self.dt)
                # hanning window
                hanning = np.hanning(2*ss)
                wind = np.ones(ws)
                wind[0:ss] = hanning[0:ss]
                wind[-ss:ws] = hanning[ss:ws]
                return wind,ws,ss
        def _powercross(self,a,b):
                f,ab = self._csd_helper(a,b,return_onesided=False) #csd
                f,aa = self._csd_helper(a,a,return_onesided=False) #psd
                f,bb = self._csd_helper(b,b,return_onesided=False) #psd
                return f,ab,aa,bb
        def _csd_helper(self,a,b,window=None,overlap=None,return_onesided=False):
                f,_t,a_ft = self._stft(a,window=window,overlap=overlap,return_onesided=return_onesided) # ( time-windows x frequency-bins )
                f,_t,b_ft = self._stft(b,window=window,overlap=overlap,return_onesided=return_onesided) # ( time-windows x frequency-bins )
                cab = np.mean(self._calc_csd(a_ft,b_ft),axis=0) #Welch method
                return f,cab
        def _calc_csd(self,a_ft,b_ft):
                cab = np.conj(a_ft)*b_ft
                return cab
        def _stft(self,tr,scaling='spectrum',window=None,overlap=None,return_onesided=False):
                wind,ws,ss = self._window(window=window,overlap=overlap)
                _f, _t, ft = stft(tr,self.fs, return_onesided=return_onesided, boundary=None,padded=True, 
                window=wind, nperseg=ws, noverlap=ss,scaling=scaling)
                _f = _f.reshape(-1)
                ft = np.atleast_2d(ft)
                return _f.T,_t.T,ft.T