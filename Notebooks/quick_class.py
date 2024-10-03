from obspy import read
import numpy as np
from scipy.signal import stft, detrend
import warnings
import fnmatch
import obspy
from obspy import read_inventory
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
# def pull_cohphadm(stanm,EvFolder,tf='ZP-21'):
#     correvpath = EvFolder / stanm / 'CORRECTED'
#     rawevpath = EvFolder / stanm
#     fclip = '.sta.' + tf
#     correvs = [f.name for f in list(correvpath.glob('*' + fclip + '.HZ.SAC'))]
#     rawevs = [f.replace(fclip,'').replace(stanm + '.','') for f in correvs]
#     events = [c.replace('.sta','').replace('.HZ.SAC','').split(stanm + '.')[-1].split('.' + tf)[0] for c in correvs]
#     cpa_list = []
#     for r,c in zip(rawevs,correvs):
#         # r = rawevs[0]
#         # c = correvs[0]
#         rawst = read(rawevpath / r)[0]
#         corrst = read(correvpath / c)[0]
#         cpa = cohphadm(rawst,corrst)
#         cpa_list.append(cpa)
#         # f,coh = cpa.COH()
#         # plt.scatter(f,coh,s=0.05)
#         # plt.axvline(cpa.fnotch(),linewidth=0.2,c='k')
#         # plt.xscale('log')
#         # plt.title('backup')
#     return events,cpa_list
def pull_cohphadm(stanm,EvFolder,tf='ZP-21',g=True):
    correvpath = EvFolder / stanm / 'CORRECTED'
    rawevpath = EvFolder / stanm
    fclip = '.sta.' + tf
    correvs = [f.name for f in list(correvpath.glob('*' + fclip + '.HZ.SAC'))]
    rawevs = [f.replace(fclip,'').replace(stanm + '.','') for f in correvs]
    events = [c.replace('.sta','').replace('.HZ.SAC','').split(stanm + '.')[-1].split('.' + tf)[0] for c in correvs]
    cpa_list = []
    for r,c in zip(rawevs,correvs):
        # r = rawevs[0]
        # c = correvs[0]
        if g:
            rawst = get_sac(rawevpath / r)
        else:
            rawst = read(rawevpath / r)[0]
        corrst = read(correvpath / c)[0]

        cpa = cohphadm(rawst,corrst)
        cpa_list.append(cpa)
        # f,coh = cpa.COH()
        # plt.scatter(f,coh,s=0.05)
        # plt.axvline(cpa.fnotch(),linewidth=0.2,c='k')
        # plt.xscale('log')
        # plt.title('backup')
    return events,cpa_list
class cohphadm(object):
        def __init__(self,A=None, B=None,overlap=0.3,csd=None,f=None,fs=None):
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
        def COH(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                coh = self._calc_coherence(ab,aa,bb)
                return f,coh
        def PH(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                ph = self._calc_phase(ab)
                return f,ph
        def ADM(self):
                f,ab,aa,bb = self._powercross(self.A.data,self.B.data)
                adm = self._calc_admittance(ab,bb)
                return f,adm
        def _calc_coherence(self,ab,aa,bb):
                coh = ((abs(ab)**2)/abs(aa*bb))
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
                if self.tlen>7200:self.window = 7200
                else:self.window = 500
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
# --------