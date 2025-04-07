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
        gs=True,
        corrected_comp='HZ',raw_comp='HDH'):
        correvpath = CorrectedFold
        rawevpath = UncorrectedFold
        # raw_stanm_subfolder=None,corrected_stanm_subfolder=None,
        #     correvpath = CorrectedFold / stanm /'CORRECTED'
        #     rawevpath = UncorrectedFold / stanm
        #     if raw_stanm_subfolder is not None:rawevpath=rawevpath/raw_stanm_subfolder
        #     if corrected_stanm_subfolder is not None:correvpath=correvpath/corrected_stanm_subfolder
        if len(tf)>0:fclip = '.sta.' + tf
        else:fclip = ''
        correvs = [f.name for f in list(correvpath.glob('*' + fclip + f'.{corrected_comp}.SAC'))]
        rawevs = [f.name for f in list(rawevpath.glob('*' + f'.{raw_comp}.SAC'))]
        evs=np.intersect1d([r.split(f'.{corrected_comp}.SAC')[0].replace(stanm+'.','').replace('.sta.'+tf,'') for r in correvs],[r.split(f'.{raw_comp}.SAC')[0] for r in rawevs])
        c=evs.copy();correvs=[f'{stanm}.{i}.{corrected_comp}.SAC' for i in c]
        c=evs.copy();rawevs=[f'{i}.{raw_comp}.SAC' for i in c]
                #     events = [c.replace('.sta','').replace('.HZ.SAC','').split(stanm + '.')[-1].split('.' + tf)[0] for c in correvs]
        events = [c.replace('.sta','').replace('.HZ.SAC','').replace('.H1.SAC','').replace('.H2.SAC','').replace(stanm+'.','') for c in correvs]
        if len(tf)>0:
                events=[c.replace('.'+tf,'') for c in events]
                correvs=[c.replace('.HZ.','.sta.'+tf+'.HZ.') for c in correvs]
        cpa_list,ev_list = [],[]

        #data_test: Tests if each event in the catalog was found.
        data_test = [(np.isin(e.Name,events)).tolist() for e in cat.loc[stanm].Events]
        #data_test_conv: Tests if each event found is also in the catalog.
        data_test_conv = [(np.isin(e,[o.Name for o in cat.loc[stanm].Events])).tolist() for e in events]

        log_warnings=[]
        #I don't want to track whether ZP coherence post-correction is collected
        if (not np.all(data_test)) & (not np.isin('HDH',[raw_comp,corrected_comp])):
                [print(f'{e.Name}: {t}') for t,e in zip(data_test,cat.loc[stanm].Events) if t==False]
                log_warnings=[f'{stanm}:{e.Name}:Catalog event not found' for t,e in zip(data_test,cat.loc[stanm].Events) if t==False]
        #        assert np.all(data_test)
        corr,raw=[],[]
        for conv_test,ev,r,c in zip(data_test_conv,events,rawevs,correvs):
                #Only collect coherences for events currently listed in the data volume.
                if not conv_test:
                        log_warnings.append(f'{stanm}:{ev}:Event not found in catalog')
                        continue #These two lines do the exact same thing.
                rawst = load_sac(rawevpath/r,rmresp=False)
                if len(rawst)==0:print(f'{stanm}|{ev}| Data missing or corrupt');continue
                corrst = load_sac(correvpath/c,rmresp=False)
                if len(corrst)==0:print(f'{stanm}|{ev}| Data missing or corrupt');continue
                corr.append(corrst);raw.append(rawst)
        corr=Stream(corr);raw=Stream(raw)
        for ev,rawst,corrst in zip(events,raw,corr):
                tend = np.min([rawst.stats.endtime,corrst.stats.endtime])
                rawst.trim(tend-7200,tend)
                corrst.trim(tend-7200,tend)
                trim_n = min([corrst.data.shape[0],rawst.data.shape[0]])
                rawst.data=rawst.data[:trim_n]
                corrst.data=corrst.data[:trim_n]
                if len(np.unique([len(rawst.data),len(corrst.data)]))>1:raise Exception('Trace lenghts not equal')
        corr=basic_preproc(corr);raw=basic_preproc(raw)
        for ev,rawst,corrst in zip(events,raw,corr):
                cpa = cohphadm(rawst,corrst)
                cpa.Event = ev
                cpa_list.append(cpa)
                ev_list.append(ev)
        _=[print(s) for s in log_warnings]
        return ev_list,cpa_list,log_warnings

class DataHold(object):
    def __init__(self):self

# corrected_trace = corrected_traces[ei]
# original_trace = original_traces[ei]
# npts = original_trace.stats.npts
# npts = 1024
# # npts=10*int(round(2**np.ceil(np.log2(120*original_trace.stats.sampling_rate))))
# # npts=2**int(np.ceil(np.log2(original_trace.stats.npts)))
# npts=int(np.ceil(np.log2(original_trace.stats.npts)))*1024
# fs = original_trace.stats.sampling_rate
# fvals, coh = coherence(original_trace.data, corrected_trace.data, fs=fs, nperseg=npts)
class cohphadm(object):
        def __init__(self,A, B=None,overlap=0.3,csd=None,f=None,fs=None):
                if B is None:A=B.copy()
                self.overlap = overlap
                self.A = A.copy()
                self.B = B.copy()
                self._meta()
                self.csd = csd
                self.f = f
        def fnotch(self,d=None):
                '''The frequency knotch root function described in Crawford et al., 1998.
                depth (d) is in meters. Returned (f) is in Hz.'''
                if d==None:
                        d = self.depth
                g = 9.80665
                f = (g/(2*np.pi*d))**0.5
                return f
        # def COH(self):
        #         f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
        #         coh = self._calc_coherence(ab,aa,bb)
        #         return f[f>0],coh[f>0]
        def COH(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                coh = self._calc_coherence(ab,aa,bb)
                fs=self.B.stats.sampling_rate
                npts = int(np.ceil(np.log2(self.B.stats.npts)))*1024
                f, coh = coherence(self.A.data,self.B.data, fs=fs, nperseg=npts)
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
                coh = ((abs(ab)**2)/(np.abs(aa)*np.abs(bb)))
                # ------------------------------
                #Complex-Valued Coherence (CVC) (Crawford&Webb,2000):
                #Only useful if I want to retain phase or slightly bias output towards 1 (unity).
                # coh = (  (ab)  /  ((aa*bb)**0.5) )
                return coh

        def _calc_phase(self,ab):
                ph = np.angle(ab,deg=True)
                return ph
        def _calc_admittance(self,ab,bb):
                ad = np.abs(ab)/bb
                return ad
        def _meta(self):
                self.dt = np.unique([s.stats.delta for s in [self.A,self.B] if isinstance(s,obspy.core.trace.Trace)])[0]
                self.fs = 1/self.dt
                self.tlen = np.size(self.A.data)/self.fs
                self.depth = abs(self.A.stats.sac.stel)*1000
                self.net = self.A.stats.network
                self.sta = self.A.stats.station
                self.chan = self.A.stats.channel
                # if self.tlen>7200:self.window = 7200
                self.window = 500
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
                f,ab = self._csd_helper(a,b,return_onesided=True)
                f,aa = self._csd_helper(a,a,return_onesided=True)
                f,bb = self._csd_helper(b,b,return_onesided=True)
                return f,ab,aa,bb
        def _csd_helper(self,a,b,window=None,overlap=None,return_onesided=False):
                f,_t,a_ft = self._stft(a,window=window,overlap=overlap,return_onesided=return_onesided)
                f,_t,b_ft = self._stft(b,window=window,overlap=overlap,return_onesided=return_onesided)
                cab = np.mean(self._calc_csd(a_ft,b_ft),axis=0)
                return f,cab
        def _calc_csd(self,a_ft,b_ft):
                cab = a_ft*np.conj(b_ft)
                return cab
        def _stft(self,tr,scaling='spectrum',window=None,overlap=None,return_onesided=False):
                wind,ws,ss = self._window(window=window,overlap=overlap)
                _f, _t, ft = stft(tr,self.fs, return_onesided=return_onesided, boundary=None,padded=False, window=wind, nperseg=ws, noverlap=ss,detrend='constant',scaling=scaling)
                _f = _f.reshape(-1)
                ft = np.atleast_2d(ft)
                ft = ft*ws
                return _f.T,_t.T,ft.T
# ------------------------------------------------      ------------------------------------------------                ------------------------------------------------
# ------------------------------------------------------------------------------------------------              ------------------------------------------------
# ------------------------------------------------      ------------------------------------------------                ------------------------------------------------