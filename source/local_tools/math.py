from modules import *
from scipy.stats import iqr

def fnotch(d,n=1):
        '''The frequency knotch root function described in Crawford et al., 1998.
        depth (d) is in meters. Returned (f) is in Hz.'''
        assert (0.5<=n) & (n<=2), f'n={n} violates 0.5<=n<=2'
        g = 9.80665
        f = (g/(2*np.pi*n*d))**0.5
        return f
def smooth(d,k=10):return np.convolve(d, np.ones(k) / k, mode='same')
def distance(sta,ev,unit='deg'):
    origins=ev.origins[0]
    stalla,evlla=[sta.Latitude,sta.Longitude],[origins.latitude,origins.longitude]
    dist=locations2degrees(stalla[0],stalla[1],evlla[0],evlla[1])
    if unit.lower()=='km':dist=degrees2kilometers(dist)
    return dist
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

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
def spectra(trace,fs=None,
    win_length=163.84,win_length_samples=None, #length of fft window
    db=True, #Converts either amplitude or power to dB, depending on options set for density or power.
    density=True,power=True, #Power and density yield the same output. just a matter of what the user prefers to call it
    taper=None, #Fraction of trace length to taper before the transform
    n_fft=2, #Padding for each window. Measured in units of integer lenghts of the fft window.
    center=True, #Centers (True) window on the frequency bins or to the left of each (False).
    complex=False, #Returns as a complex value. not affected by other options
    spectrogram=False #Returns a spectrogram. if False, this just averages over the time-axis to yield a single trace spectra.
    ):
    def next_pow2(win_length_samples, win_length, sampling_rate): 
        #Overrides the window length to be the next multiple of two.
        #This ensures more explicit control over things like padding aswell as an overall faster computation.
        _next_pow2 = lambda n:int(round(2**np.ceil(np.log2(n))))
        if win_length_samples is None:
            if win_length is None:win_length_samples=_next_pow2(120*sampling_rate)
            elif win_length is not None:win_length_samples=_next_pow2(win_length*sampling_rate)
        elif win_length_samples is not None:
            win_length_samples = int(win_length_samples)
            if win_length_samples!=_next_pow2(win_length_samples):raise ValueError('Parameter win_length_samples must be a power of 2.')
        return win_length_samples
    density=not ((not power) or (not density)) #This bool ensures correct use of both options no matter which is used.
    if fs==None:fs=trace.stats.sampling_rate
    trace.detrend('simple');trace.detrend('demean')
    if isinstance(trace,Trace):x=trace.data.copy() #demean & detrend
    x=x.astype(float)
    if taper:trace.taper(taper)
    win_length_samples = next_pow2(win_length_samples, win_length, fs)
    hop_length = win_length_samples // 4

    n_fft = int(n_fft)*win_length_samples #Pads the time window that has a length of win_length_samples to equal a length of n_fft
    #Amplitude spectra
    S=librosa.stft(x,n_fft=n_fft,hop_length=hop_length,win_length=win_length_samples,center=center)

    S,phase=librosa.magphase(S) #Separates amplitude (real) from phase (complex)

    f=(np.arange(S.shape[0]) * (fs/win_length_samples))
    t=((np.arange(S.shape[1]) * hop_length)/fs)
    n=2*(1/trace.stats.sampling_rate)/trace.data.size

    #Power spectra
    if density:S=(abs(S)**2)*n #Power Density
    else:S=S=(abs(S)**1)*n #Power

    d=S
    # if db:d=10*np.log10(d)
    if db:d=PowDisp_to_AcceldB(f,d)
    if complex:d=d*phase
    if not spectrogram:d=d.mean(axis=1)

    return f,d

def PowDisp_to_AcceldB(f,y,smooth=False):
    y=20*np.log10(y) #decibels
    if smooth:y=utils.smooth(y,50)
    disp_to_accel = 40*np.log10(2*np.pi*f,where=f>0.)
    y=y+disp_to_accel #displacement to acceleration
    y=y-np.mean(abs(y),axis=0)+40 #demean
    return y
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

def octave_average(d, f, fmin=1/100, fmax=1, fraction=8, domain="geo", N=291):
    f = np.asarray(f); d = np.atleast_2d(d)
    if fmin is None: fmin = f[f>0].min()
    if fmax is None: fmax = f.max()
    if domain == "geo":
        N = int(N)
        q = (fmax/fmin)**(1.0/(N-1))
        edges = np.geomspace(fmin, fmax, N)# centers (N,)
        fc=edges[1:]
        # edges = fmin * q**(np.arange(N+1) - 0.5)               # N+1 edges, half-step in log-space
    else:
        k = np.arange(np.log2(fmin), np.log2(fmax) + 1e-12, 1.0/fraction)
        fc = 2.0**k                                            # centers (B,)
        r = 2.0**(1.0/(2.0*fraction))
        edges = np.c_[fc/r, fc*r]                              # temp edges per band → flatten below
        edges = np.r_[edges[:,0], edges[-1,1]]                 # (B+1,)
    # build band masks from edges
    lo, hi = edges[:-1], edges[1:]
    m = (f[None,:] >= lo[:,None]) & (f[None,:] <= hi[:,None])  # (B,F)
    keep = m.any(1); m = m[keep]; fc = fc[keep]
    # averaging domain
    if domain == "db":      x = 10**(d/10.0)
    elif domain == "log":   x = np.exp(d)
    else:                   x = d
    X = np.where(m[None,:,:], x[:,None,:], np.nan)             # (Nsig,B,F)
    avg = np.nanmean(X, axis=2)
    if domain == "db":      avg = 10*np.log10(avg)
    elif domain == "log":   avg = np.log(avg)
    return fc, avg

def cohstats(coh,margin=1,axis=0):
    coh = coh.squeeze()
    coh = coh[~np.any(np.isnan(coh),axis=1)]
    q1q3iqr = lambda coh,axis=0:np.array((np.percentile(coh,25,axis=axis),np.percentile(coh,75,axis=axis),iqr(coh,rng=(0.25,0.75),axis=axis)))
    S=q1q3iqr(coh)
    lower_whisker=S[0,:]-1.5*S[2,:]
    upper_whisker=S[1,:]+1.5*S[2,:]

    inliers = (coh<=(np.array(upper_whisker)*margin))&(coh>=(np.array(lower_whisker)/margin))
    outliers = ~inliers
    verify=np.sum([np.any((c[o]<=np.max(c[i]))&(c[o]>=np.min(c[i]))) for c,i,o in zip(coh.T,inliers.T,outliers.T)])
    # median=np.array([np.median(c[i]) for c,i in zip(coh.T,inliers.T)])
    mean_inliers = np.array([np.mean(c[i]) for c, i in zip(coh.T, inliers.T)])
    # inliers=[c[(c>=(lw/margin))&(c<=(uw*margin))] for lw,uw,c in zip(lower_whisker,upper_whisker,coh.T)]
    # outliers=[c[(c<(lw/margin))+(c>(uw*margin))] for lw,uw,c in zip(lower_whisker,upper_whisker,coh.T)]
    # upper_whisker=[np.max(a) if len(a)>0 else w for a,w in zip(inliers,upper_whisker)]
    # lower_whisker=[np.min(a) if len(a)>0 else w for a,w in zip(inliers,lower_whisker)]
    lower=lower_whisker
    upper=upper_whisker
    # upper=np.array([np.max(c[i]) for c,i in zip(coh.T,inliers.T)])
    # lower=np.array([np.min(c[i]) for c,i in zip(coh.T,inliers.T)])
    # return {'upper':upper_whisker,'lower':lower_whisker,'median':median,'outliers':outliers}
    return upper, lower, mean_inliers, outliers, inliers

def smooth(data, nd, axis=0):
    if np.any(data):
        if data.ndim > 1:
            filt = np.zeros(data.shape)
            for i in range(data.shape[::-1][axis]):
                if axis == 0:filt[:, i] = np.convolve(data[:, i], np.ones((nd,))/nd, mode='same')
                elif axis == 1:filt[i, :] = np.convolve(data[i, :], np.ones((nd,))/nd, mode='same')
        else:filt = np.convolve(data, np.ones((nd,))/nd, mode='same')
        return filt
    else:return None
def Stream_LogPSD(st,smoothed=True,window=7200,overlap=0.3):
    fs=st[0].stats.sampling_rate
    dt=1/fs
    days = [[tr.stats.starttime.strftime('%Y.%j') for tr in sst] for sst in [st.select(channel='*1'),st.select(channel='*2'),st.select(channel='*Z'),st.select(channel='*DH')]]
    ws = int(window/dt)
    ss = int(window*overlap/dt)
    hanning = np.hanning(2*ss)
    wind = np.ones(ws)
    wind[0:ss] = hanning[0:ss]
    wind[-ss:ws] = hanning[ss:ws]
    _stft = lambda tr:stft(tr.data, fs, return_onesided=False, boundary=None,
    padded=False, window=wind, nperseg=ws, noverlap=ss,detrend='constant')
    f, t, _ = _stft(st[0])
    psd = np.array([(_stft(tr.data)[-1]*ws) for tr in st])
    psd = abs(psd)**2*2/dt
    logpsd = 10*np.log10(psd,where=(psd>0.))
    # smoothed = True if logpsd.shape[-1]>50 else False
    smoothed = True
    if smoothed:logpsd=np.array([smooth(sl,50,axis=0) for sl in logpsd])
    # logpsd = logpsd.mean(axis=2).T
    faxis = int(len(f)/2)
    f = f[0:faxis]
    logpsd = logpsd[:,0:faxis,:]
    # logpsd = 10*np.log10(psd,where=(psd>0.))
    disp_to_accel=40*np.log10(2*np.pi*f,where=f>0.).reshape(-1,1);disp_to_accel[f==0]=0
    faxis=f>0 #Extract positive frequencies only
    disp_to_accel=disp_to_accel[faxis]
    f=f[faxis]
    sls = logpsd[:,faxis,:]
    logpsd=sls
    comps = [tr.stats.channel.replace('HDH','HP')[-1] for tr in st]
    sls=logpsd
    if ~np.isin('P',comps):
        sls = sls+disp_to_accel
    sls = np.array([sl-abs(sl[f<1,:]).mean(axis=0) for sl in sls])
    sls = sls[:,:,~np.isinf(sls.sum(axis=0).sum(axis=0))]
    sls=sls.mean(axis=-1)
    return comps,days,f,t,sls