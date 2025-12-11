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

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *


# get wins
def get_wins(Ph,p_lead,s_lead,u,noisewin,bn,s,wlen):
    # start value
    start = lambda k: Ph[k]-(p_lead if k in ['P','Pdiff'] else 0)  - (s_lead if k in ['S','Sdiff'] else 0)
    # stop value
    stop = lambda k: Ph[k]+wlen[k][bn] - (p_lead if k in ['P','Pdiff'] else 0) - (s_lead if k in ['S','Sdiff'] else 0)
    # wins value
    wins = {k:[start(k),stop(k)] for k in Ph.keys()}
    wins.update({'Rg':[s.Origin+i for i in RgWin(x,u)],
    'Noise':noisewin})
    if 'P' in Ph.keys():wins['Noise'] = [min(wins['P'])-list(p_lead+abs(np.diff(noisewin)))[0],min(wins['P'])-p_lead]
    if 'Pdiff' in Ph.keys():wins['Noise'] = [min(wins['Pdiff'])-list(p_lead+abs(np.diff(noisewin)))[0],min(wins['Pdiff'])-p_lead]
    return AttribDict(wins)
# get noise win
def get_noise_win(iSR,p_lead,s_lead,wlen,prenoise_buffer):
    # minph value
    minph = iSR[iSR.Distance==iSR.Distance.min()].iloc[0].Phases(phases=('P','S'))
    # minph value
    minph={k:[minph[k][0]-(p_lead if k in ['P','Pdiff'] else 0)  - (s_lead if k in ['S','Sdiff'] else 0),minph[k][0]+wlen[k]['30_100'] - (p_lead if k in ['P','Pdiff'] else 0) - (s_lead if k in ['S','Sdiff'] else 0)] for k in minph.keys()}
    # if 'P' in minph.keys():noisewin=[iSR.Origin[0] + (prenoise_buffer*7200) , iSR.Origin[0]+(min(minph['P'])-bleed_buffer)]
    # else:
    noisewin=[iSR.Origin[0]  , iSR.Origin[0] + noisewlen]
    return noisewin
# baz value
baz=lambda:obspy.geodetics.base.gps2dist_azimuth(s.Latitude,s.Longitude,s.LaLo[0],s.LaLo[1])[1]
# function zne2lqt
def zne2lqt(st,baz,inc,rgtime):
    z=st.select(channel='*Z')[0].data
    # n value
    n=st.select(channel='*1')[0].data
    # e value
    e=st.select(channel='*2')[0].data
    # ns value
    ns=min([len(z),len(n),len(e)]) #sometimes its one sample off after trimming
    z,n,e=z[:ns],n[:ns],e[:ns]
    l,q,t=obspy.signal.rotate.rotate_zne_lqt(z,n,e,baz,inc)
    for i,chan,d in zip(st,['L','Q','T'],[l,q,t]):
        i.data=d;i.stats.channel=chan
        i.stats.location=i.stats.location.split('.')[0]+f'.{chan}'
        i.stats.baz=baz;i.stats.inc=inc
    return st
# collect traces
def collect_traces(s):
    # tr value
    tr=s.Traces()

    # h value
    h=s.Traces(channel='H1',methods=['Original','NoiseCut'])
    # h2 value
    h2=s.Traces(channel='H2',methods=['Original','NoiseCut'])
    [h.append(i) for i in h2]

    for m in ['Original','ATaCR']:
        # tmp value
        tmp=h.select(location='Original').copy()
        # loop over tmp
        for i in tmp:i.stats.location=m
        [tr.append(i) for i in tmp]
    [tr.append(i) for i in h.select(location='NoiseCut')]
    # loop over tr
    for i in tr:i.stats.location=f'{i.stats.location}.{i.stats.channel}'
    # loop over tr
    for i in tr:i.stats.location=i.stats.location.replace('ATaCR','TF').replace('NoiseCut','HPS')
    return tr
# standardbands value
standardbands = ['1_10','10_30','30_100']
# bandmap value
bandmap=lambda bn:np.array(standardbands)[np.min(np.where(np.min([float(i) for i in bn.split('_')])<=np.array([[float(i) for i in bi.split('_')] for bi in standardbands])[:,1]))]
# --------------------------------------------------------------------------------
P_snr_wlen=AttribDict({'1_10':150,'10_30':150,'30_100':150,})
# S snr wlen value
S_snr_wlen=AttribDict({'1_10':150,'10_30':150,'30_100':150,})
# Pdiff snr wlen value
Pdiff_snr_wlen=AttribDict({'1_10':150,'10_30':150,'30_100':150,})
# Sdiff snr wlen value
Sdiff_snr_wlen=AttribDict({'1_10':150,'10_30':150,'30_100':150,})
# Rg snr wlen value
Rg_snr_wlen=AttribDict({'1_10':[2.0,4.2],'10_30':[2.0,4.2],'30_100':[2.0,4.2],})

# note value
note = 'V04'
fold = dirs.SNR/'SNR.Models'
cat = catalog.copy()

# cat.sr.reset_index()
SR=cat.sr.sort_values(by=['Name','Distance']).copy()
SR.sort_values(by='Magnitude',inplace=True,ascending=False)
wlen = {'P':P_snr_wlen,'S':S_snr_wlen,
'Sdiff':Sdiff_snr_wlen,'Pdiff':Pdiff_snr_wlen,
'Rg':Rg_snr_wlen}

def make_bands(N, width=None, lo=1.0, hi=100.0):
    if width is None:width=(hi-lo)/N
    if width<=0 or width>(hi-lo): raise ValueError("width must be in (0, hi-lo]")
    s=np.linspace(lo, hi-width, N)
    return np.c_[s, s+width]


# bands = ['1_10','10_30','30_100'] #standard
# bands = ['1_10','10_20','20_30','30_100'] #4 bands, a bit wider
# bands = ['1_10','10_15','15_30','30_100'] #ideal for test but risky

# bands=np.arange(0,102,2)

width=5;N=100;bands=make_bands(N,width);note=f'V04_5s_bandwidth'
bands=[f'{b[0]:.2f}_{b[1]:.2f}' for b in bands]

noisewlen=300;bleed_buffer = 24
prenoise_buffer=.01
p_lead=pdiff_lead=20;s_lead=sdiff_lead=20

snrlabel='SNR'
SR[snrlabel]=[[] for _ in SR.iloc];SR['wins']=[[] for _ in SR.iloc]
evn=SR.Name.unique()
RgWin=lambda x,u:sorted([x/i for i in u])
bandkeys=np.array(list(wlen['P'].keys()))
bandkey_bands=np.array([[float(t) for t in i.split('_')] for i in wlen['P'].keys()])
methods=['Original','NoiseCut','ATaCR']

# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------



filter_type = 'acausul' #Filter type used.
LQT = False #Experimental. Not finished. Don't use. Designed to calculate SNR after rotation into ray coordinates.


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

for ei,ev in enumerate(evn):
    iSR = SR[SR.Name==ev].copy()
    for si,s in enumerate(iSR.iloc):
        siPh=s.Phases(phases=('P','S','Sdiff','Pdiff'))
        retphases = list(siPh.keys())
        if ('Pdiff' in retphases)&('P' in retphases):siPh.remove('Pdiff')
        if ('Sdiff' in retphases)&('S' in retphases):siPh.remove('Sdiff')
        retphases = list(siPh.keys())
        Ph={k:s.Origin+siPh[k][0] for k in siPh.keys()}
        inclination=np.mean([siPh[k][1] for k in siPh.keys()]) if len(siPh.keys())>0 else 0


        tr=collect_traces(s)
        mkeys=[i.stats.location for i in tr]

        if tr==None:continue
        tr=tr.copy().trim(s.Origin,s.Origin+7200)


        x=degrees2kilometers(s.Distance)
        print(f' {ei+1}/{len(evn)} | {ev} | {si+1}/{len(iSR)} | Phases: {', '.join(retphases)}')
        bsnr=AttribDict({})
        bwins = AttribDict({})
        noisewin = get_noise_win(iSR,p_lead,s_lead,wlen,prenoise_buffer)
        for b in bands:
            band=[float(i) for i in b.split('_')]
            bn=b
            bn=np.array(bandmap(bn))
            assert bn.size==1
            bn=str(bn)
            itr=tr.copy()
            nproc=1
            for _ in range(nproc):itr.detrend('demean').detrend('linear')
            itr.taper(prenoise_buffer)
            itr.filter('bandpass',freqmax=1/min(band),freqmin=1/max(band),zerophase={'causul':False,'acausul':True}[filter_type])
            for _ in range(nproc):itr.detrend('demean').detrend('linear')
            itr.taper(prenoise_buffer)

            u=wlen['Rg'][bn]
            wins=get_wins(Ph,p_lead,s_lead,u,noisewin,bn,s,wlen)

            if LQT:
                rot_lqt=Stream(unravel([zne2lqt(itr.select(location=f'{m}*').copy(),baz(),inclination,wins.Rg) for m in ['Original','HPS','TF']]))
                itr=rot_lqt
            # test = rot_lqt.copy()
            # test=Stream([test[j] for j in np.where(np.isin([i.stats.location for i in test],['Original.L','Original.Q']))[0]])

            # I'm using a post-filter taper for the first 1% to protect from any risk of gibbs in my noise window. I've seen it in a few traces on specific stations (XN). 

            # ----Window tests
            # assert (tr[0].stats.starttime-min(wins['Noise']))<=0,'Noise window preceeds trace'
            if (('S' in wins.keys())&('P' in wins.keys())):assert (max(wins['P']) - min(wins['S']))<0,'P and S windows overlap'
            if (('S' in wins.keys())&('Rg' in wins.keys())):assert (max(wins['S']) - min(wins['Rg']))<0,'S and Rg windows overlap'
            if (('Sdiff' in wins.keys())&('Pdiff' in wins.keys())):assert (max(wins['Pdiff']) - min(wins['Sdiff']))<0,'P and S windows overlap'
            if (('Sdiff' in wins.keys())&('Rg' in wins.keys())):assert (max(wins['Sdiff']) - min(wins['Rg']))<0,'S and Rg windows overlap'
            # if (('P' in wins.keys())&('Noise' in wins.keys())):assert ( ((max(wins['Noise']))+5)-min(wins['P']))<=0,'Noise window overlaps with P'
            # Aggregate noise
            # rms = lambda y:(y**2).mean()**.5 #<----I was originally doing it as this, the rms of amplitude.
            # rmse = lambda y, y_hat: (((y - y_hat)**2).mean())**.5 #<---General definition of RMSE
            rmse=lambda y:( (  ( abs(y)-abs(y).mean() )**2  ).mean())**.5 #<----Standard (mean) definition of (ABS) RMSE
            # rmse = lambda y:abs(y).std(ddof=0) #<----the same thing when ddof is 0
            noise = lambda k: itr.select(location=f'*{k}').copy().trim(min(wins['Noise']),max(wins['Noise']))[0].data

            signalpeak = lambda i:abs(i.data).max() #the original

            LT={i.stats.location:
                rmse(i.data)
                for i in itr.copy().trim(min(wins['Noise']),max(wins['Noise']))}
            # Aggregate signal
            ST= {k: #<---Phase
                {i.stats.location: #<---Method
                signalpeak(i)
                for i in itr.copy().trim(min(wins[k]),max(wins[k]))} for k in wins.keys() if k!='Noise'}
            snr={k: #<---Phase
                {i.stats.location: #<---Method
                ST[k][i.stats.location]/LT[i.stats.location] #<----SNR
                for i in itr} for k in wins.keys() if k!='Noise'}
            snr.update({'ST':ST,'LT':LT,'stnm':s.StaName}) #Keep the num/denom for later.
            wins.update({'stnm':s.StaName})
            bsnr.update({b:snr})
            bwins.update({b:wins})
        idx=np.where((SR.StaName==s.StaName)&(SR.Name==ev))[0][0]
        SR.iat[idx, SR.columns.get_loc(snrlabel)]=bsnr
        SR.iat[idx, SR.columns.get_loc('wins')]=bwins

safekeys=['Name', 'StaName', 'Station', 'Event', 'Network', 'Origin',
'LaLo', 'Distance', 'Magnitude', 'Stations',
'Latitude', 'Longitude', 'Experiment', 'Environment', 'Pressure_Gauge',
'StaDepth', 'Start', 'End', 'NoiseAverage', 'Seismometer',
'Sediment_Thickness_m', 'Instrument_Design', 'Distance_from_Land_km',
'Distance_to_Plate_Boundary_km', 'Surface_Current_ms',
'Crustal_Age_Myr', 'Deployment_Length_days', 'Inventory','wins', snrlabel]


fold = dirs.SNR/'SNR.Models'
file=f'{snrlabel}_{filter_type}.filter_{note}.pkl'

if LQT:file=file.replace('.pkl','.LQT.pkl')

if len(bands)>3:file=file.replace('.pkl',f'_{len(bands)}_bands.pkl')
fold.mkdir(exist_ok=True)
SR[safekeys].to_pickle(fold/file)

