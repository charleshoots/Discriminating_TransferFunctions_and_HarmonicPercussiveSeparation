from imports import *
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
def get_noise(dirs,stanm):return load_pickle(list((dirs.SpectraAvg/stanm).glob('*sta.pkl'))[0])
reportfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/Analysis/NetworkCoherences')
plotfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_FigureArchive/_GEN6')
bands = ['1-10','10-30','30-100']
def smooth(d,k=10):return np.convolve(d, np.ones(k) / k, mode='same')
NoiseColors = [mcolors.to_hex(m) for m in [cm.__dict__[e].resampled(30).resampled(6).colors for e in ['devon_categorical']][0]]
plotfold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_FigureArchive/_GEN6/StationEventPages/low_rez')
cat = catalog.copy()
def unravel(lst):return list(itertools.chain.from_iterable(lst))


def get_stft(tr,window=7200,overlap=0.3):
    fs=tr.stats.sampling_rate;dt=1/fs
    ws = int(window/dt)
    ss = int(window*overlap/dt)
    hanning = np.hanning(2*ss)
    wind = np.ones(ws)
    wind[0:ss] = hanning[0:ss]
    wind[-ss:ws] = hanning[ss:ws]
    _f, _t, ft = stft(tr.data,fs, 
    return_onesided=False, 
    boundary=None,padded=False, 
    window=wind, nperseg=ws, noverlap=ss,
    detrend='constant')
    return _f,ft

# ================================================================================================================

def get_noise_trace(dirs,stanm):
    ### Aprox. 40seconds to run per station requested
    fchans = ['Z','DH','1','2']
    noise = get_noise(dirs,stanm)
    files = [Path(f) for f in noise.day_files[noise.gooddays]]
    goodwins = noise.day_goodwins
    del noise
    days = [f.name.split('.spectra.pkl')[0] for f in files]
    #41-seconds for 1 station of noise
    trchans = ['HZ', 'HDH', 'H1', 'H2'];d = AttribDict()
    [d.update({c:[]}) for c in trchans]
    day_noise = d.copy()
    for gw,day in zip(goodwins,days):
        dayst = Stream([load_sac(dirs.Noise/stanm/f'{day}..H{c}.SAC')[0][0] for c in fchans])
        clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
        # [day_noise[tr.stats.channel].append(np.mean(get_stft(tr)[1][:,gw],axis=1)) for tr in dayst]
        [day_noise[tr.stats.channel].append(get_stft(tr)[1][:,gw]) for tr in dayst]
        if day==days[0]:
            save_st = dayst.copy()
        del dayst
    noise=d.copy()
    [noise.update({c:np.mean(np.hstack(day_noise[c]),axis=1)}) for c in trchans]
    for c in trchans:save_st.select(channel=c)[0].data = noise[c]
    del day_noise,noise
    noise = NoiseWrapper(save_st)
    return noise

stanm = cat.iloc[0].StaName
noise = get_noise_trace(dirs,stanm)

s = catalog.iloc[0]
stanm = s.StaName
evm = Catalog([s.Events[0]])
metrics = get_station_events(stanm,dirs.Events,type='metrics',evmeta=evm)