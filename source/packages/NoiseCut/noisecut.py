# -----------------------------------------------------------------------------
# NoiseCut
#
# This file is part of the NoiseCut library. For licensing information see the
# accompanying file `LICENSE`.
# -----------------------------------------------------------------------------

import math
import numpy as np
import librosa
from obspy import Trace
from obspy.core.util.attribdict import AttribDict
import librosa.display
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from pathlib import Path
import os

def _next_pow2(n):return int(round(2**np.ceil(np.log2(n))))
def _valid_win_length_samples(win_length_samples, win_length, sampling_rate):
    if win_length_samples is None and win_length is None:
        # fully automatic window length
        win_length_samples = _next_pow2(120*sampling_rate)
    elif win_length_samples is None and win_length is not None:
        win_length_samples = _next_pow2(win_length*sampling_rate) #<----DEFAULT.
    elif win_length_samples is not None and win_length is not None:
        raise ValueError(
            'Parameters win_length and win_length_samples are mutually '
            'exclusive.')
    elif win_length_samples is not None and win_length is None:
        # check win_length_samples is a power of 2
        win_length_samples = int(win_length_samples)
        if win_length_samples != _next_pow2(win_length_samples):
            raise ValueError('Parameter win_length_samples must be a power of 2.')
    return win_length_samples
def plot_noisecut_spectrograms(S_full, S_background, S_hps, frequencies, times, fig=None,ymax=1,figsize=(15,9)):
    show = False
    units = ['seconds','minutes','hours']
    while times[-1]>7200:
        times = times/60
        units.pop(0)
    if fig is None:
        fig= plt.figure(figsize=figsize)
        # PLOT-1::----
        axs=fig.add_subplot(311)
        pcm=axs.pcolormesh(times, frequencies, librosa.power_to_db(np.abs(S_full)), cmap = 'magma', shading= 'auto')
        plt.title('Full spectrogram', fontsize=14)
        plt.ylabel('Frequency (Hz)', fontsize=14)
        plt.yticks (fontsize= 14)
        axs.set_xticks([])
        cbar=fig.colorbar(pcm, ax=axs, pad= 0.01)
        cbar.ax.tick_params(labelsize=14)
        plt.yscale('log')
        plt.ylim(0,ymax)
        # PLOT-2::----
        axs=fig.add_subplot(312)
        pcm=axs.pcolormesh(times, frequencies, librosa.power_to_db(np.abs(S_background)), cmap = 'magma', shading= 'auto')
        plt.title('Noise spectrogram', fontsize=14)
        plt.ylabel('Frequency (Hz)', fontsize=14)
        plt.yticks (fontsize= 14)
        axs.set_xticks([])
        cbar=fig.colorbar(pcm, ax=axs, pad= 0.01)
        cbar.ax.tick_params(labelsize=14)
        plt.yscale('log')
        plt.ylim(0,ymax)
        # PLOT-3::----
        axs=fig.add_subplot(313)
        pcm=axs.pcolormesh(times, frequencies, librosa.power_to_db(np.abs(S_hps)), cmap = 'magma', shading= 'auto')
        plt.title('Denoised spectrogram', fontsize=14)
        plt.ylabel('Frequency (Hz)', fontsize=14)
        plt.yticks (fontsize= 14)
        cbar=fig.colorbar(pcm, ax=axs, pad= 0.01)
        cbar.ax.tick_params(labelsize=14)
        plt.yscale('log')
        plt.ylim(0,ymax)
        # labels at the end
        plt.xlabel(units[0])
        # plt.tight_layout()
        fig.savefig ('NoiseCut spectrograms.png', dpi=100)
def QAPlot(S,f,t,clim=None,ttl='No title',stats=None,log=True,cmap=None,pabs=False,usedbS=True):
    plt.close()
    if cmap:
        cmap = mcolors.LinearSegmentedColormap.from_list("CustomDiverging", 
        list(zip([0, 0.5, 1], [(1, 1, 1), (0.5, 0.5, 0.5),(0, 0, 0)])))
        cmap = [cmap,[-100,-40]]
    else:vlim=[None,None]
    fig,axes=plt.subplots(2,1,figsize=(5.5,9))
    fig.suptitle(ttl)
    if clim is not None:clim = librosa.power_to_db(clim)
    if cmap is None:cmap='magma';vlim=clim
    else:cmap,vlim=cmap;clim=vlim
    if usedbS:dbS1=librosa.power_to_db(np.abs(S))
    else:dbS1=10*np.log10(S);vlim=[-100,-30];clim=vlim
    flims=f[(~(S.sum(axis=1)==0))&(f>0)]
    tlim=[86400-7200,86400]
    yt=[1,1/5,1/10,1/30,1/60,1/100]
    for axi,ax in zip(range(2),axes):
        s = ax.pcolormesh(t, f, dbS1 if (pabs&(axi==1)) else dbS1, cmap = cmap if axi==1 else 'magma', shading= 'auto',clim=vlim if axi==1 else clim)
        ax.set_xlabel('sec')
        ax.set_ylabel('Hz')
        ax.set_xlim(tlim)
        ax.set_yscale('linear')
        if axi==1:ax.set_yscale('log');ax.set_yticks(yt);ax.set_ylim(flims.min(),flims.max())
        if axi==0:ax.set_yscale('log');ax.set_ylim([f[1],1.0])
        ax.set_xlim(tlim)
        plt.colorbar(s)
    ttl = ttl[:ttl.find('\n')]
    file=ttl.replace(',','').replace('(','').replace(')','').replace('<','_').replace('=','').replace(' ','_').replace('+','')
    stanm = f'{stats.network}.{stats.station}'
    evn = f'{stats.starttime.strftime('%Y.%j.%H.%M')}'
    staevname = f'{stanm}.{evn}.{stats.channel}'
    file=staevname + '_' + file
    plotfold = Path(os.getcwd()) / 'QAOutput' / stanm;plotfold.mkdir(exist_ok=True,parents=True)
    fig.savefig(plotfold/(file+'.png'), dpi=300)
    return fig
def noisecut(
        trace,
        ret_spectrograms=False,
        win_length_samples=None,
        win_length=163.84,resample_factor=1.0,width=None,kernel_size=80,overlap=0.75,verbose=False,margin=5):
    '''
    Reduce noise from all the components of the OBS data using HPS noise
    reduction algorithms.
    :param win_length_samples:
        Window length in samples. Must be a power of 2. Alternatively it can be
        set with `win_length`.
    :type win_length_samples:
        int
    :param win_length:
        Window length [s]. Alternatively it can be set with
        `win_length_samples`.
    :type win_length:
        float
    :returns:
        The HPS trace and the spectrograms of the original, noise, and hps
        trace as well as an array with the frequencies.
    :return_type:
        tuple ``(hps_trace, (s_original, s_noise, s_hps, frequencies))``
    '''
    verbose=np.array([verbose]).ravel()
    if True in verbose:verbose=np.array([0,1,2,3,4,5,6,7,8,9,10])
    if resample_factor!=1.0:
        print('Resampling')
        trace.resample(np.ceil(trace.stats.sampling_rate*resample_factor))
    x = trace.data.astype(float)

    win_length_samples = _valid_win_length_samples(
        win_length_samples, win_length, trace.stats.sampling_rate)

    windows = int(1/(1-overlap))
    hop_length = win_length_samples // windows
    n_fft = win_length_samples
    if 0 in verbose:
        print('Building raw spectrogram from STFT | hop_length={hop_length}, n_fft={n_fft}, win_length_samples={win_length_samples}'.format(hop_length=hop_length,n_fft=n_fft,win_length_samples=win_length_samples))
    # Compute the spectrogram amplitude and phase
    S_full, phase = librosa.magphase(librosa.stft(
        x,n_fft=n_fft,
        hop_length=hop_length,win_length=win_length_samples))
    fq=(np.arange(S_full.shape[0]) * (trace.stats.sampling_rate/win_length_samples))
    t=((np.arange(S_full.shape[1]) * hop_length)/trace.stats.sampling_rate)
    clim = [S_full.min(),S_full.max()]

    # Concerning win_length:
    # Smaller values improve the temporal resolution of the STFT
    # (i.e. the ability to discriminate impulses that are closely spaced in time) at the expense
    # of frequency resolution (i.e. the ability to discriminate pure tones that are closely spaced
    # in frequency). This effect is known as the time-frequency localization trade-off and needs to
    # be adjusted according to the properties of the input signal

    if 1 in verbose:fig1=QAPlot(S_full,fq,t,clim,'S01 S_full (Original)',trace.stats)

    l1 = math.floor((0.1 * win_length_samples) / trace.stats.sampling_rate)
    l2 = math.ceil((1 * win_length_samples) / trace.stats.sampling_rate)

    # We consider the frequency range of [0.1-1] Hz for the second step
    S_full2 = np.zeros((S_full.shape[0], S_full.shape[1]))
    S_full2[l1:l2, :] = S_full[l1:l2, :]
    if 2 in verbose:fig2=QAPlot(S_full2,fq,t,clim,'S02 S_full2, (0.1<=f<1.0)\n X'' used for MED in Zali''23',trace.stats)
    # We consider the frequency range outside of [0.1-1] Hz for the first step
    S_full1 = np.zeros((S_full.shape[0], S_full.shape[1]))
    S_full1[:l1, :] = S_full[:l1, :]
    S_full1[l2:, :] = S_full[l2:, :]

    if 3 in verbose:fig3=QAPlot(S_full1,fq,t,clim,'S03 S_full1, not (0.1<=f<1.0)\n V used for SIM in Zali''23 (eq.1)',trace.stats)
    # We'll compare frames using cosine similarity, and aggregate similar
    # frames by taking their (per-frequency) median value.

    if width is None: # Waiting Factor 
        wait_factor = 7200 #User defined wait factor, in seconds.
        # width: 
        #   The waiting factor (width) is the minimum number of samples between filtered data in the SIM stage of NoiseCut.
        #   Zali '23 default setting: width=200. Set to produce a 2-hr wait factor for a 100hz trace but it's closer to 2.23 hours.

        # width = ((((S_full1.shape[-1] - 1) // 2) - 1) // 5) - 10 #Hoots: This reproduces Zali's original width=200 for a 100Hz trace and tries to retain the same width for lower rate data.
        width = np.where(t>=wait_factor)[0].min() #This is a more direct way of giving the wait factor suggested to me by Janiszewski.

    if 0 in verbose:print(f'(SIM) Match-Filter | width={width} samples ({int((t[width]/3600))} hours)')

    S_filter = librosa.decompose.nn_filter(S_full1,aggregate=np.median,metric='cosine', width=width)
    
    S_filter = np.minimum(np.abs(S_full1), np.abs(S_filter)) # The output of the filter shouldn't be greater than the input
    if 4 in verbose:fig4=QAPlot(S_filter,fq,t,clim,'S04 S_filter, not (0.1<=f<1.0) \n W-hat (from W) in Zali''23',trace.stats)
    margin_i = 1;power = 2

    # Once we have the masks, simply multiply them with the input spectrogram
    # to separate the components.
    mask_i = librosa.util.softmask(S_filter,margin_i * (S_full1 - S_filter),power=power)

    S_background = mask_i * S_full1

    if 5 in verbose:fig5=QAPlot(S_background,fq,t,clim,'S05 S_background, not (0.1<=f<1.0) \n [R] Repeating spectrogram (R, eq.5 in Zali''23)',trace.stats)
    if 0 in verbose:print('HPSS Median-Filter | kernel_size='+str(kernel_size) + ', margin=' + str(margin))
    # In the second step we apply a median filter
    D_harmonic, D_percussive = librosa.decompose.hpss(
        S_full2,
        kernel_size=kernel_size,
        margin=margin)
    if 6 in verbose:fig6=QAPlot(D_harmonic,fq,t,clim,'S06 D_harmonic, (0.1<=f<1.0) \n [H] Harmonic from MED in Zali''23',trace.stats)
    if 7 in verbose:fig7=QAPlot(D_percussive,fq,t,clim,'S07 D_percussive, (0.1<=f<1.0) \n Percussive component of MED never used in Zali''23',trace.stats)
    S_background = S_background + D_harmonic
    if 8 in verbose:fig8=QAPlot(S_background,fq,t,clim=None,ttl='S08 S_background + D_harmonic \n N = R + H in Zali''23',stats=trace.stats,usedbS=False)
    f = S_background * phase
    L = x.shape[0]
    new = librosa.istft(f,hop_length=hop_length,win_length=win_length_samples,window='hann',length=L) #Noise model returned to time-domain
    z = x - new #Noise, in the time-domain, subtracted from the original input.

    stats = trace.stats
    if len(stats.location)>0:stats.location = stats.location + '->HPS'
    else:stats.location = 'HPS' #Adds a note to the metadata that this trace was processed in NoiseCut.

    hps_trace = Trace(data=z, header=stats)

    # hps_trace.write( 'NoiseCut.mseed', format='MSEED', encoding=5, reclen=4096)
    if 9 in verbose:fig10=QAPlot(S_full - S_background,fq,t,clim,'S10 Output \n T = V - N = V - (R + H) in Zali''23',trace.stats)
    if ret_spectrograms:
        S_hps = S_full - S_background
        df = trace.stats.sampling_rate/win_length_samples
        frequencies = np.arange(S_hps.shape[0]) * df
        times = np.arange(S_hps.shape[1]) * hop_length
        times = times/trace.stats.sampling_rate #Time axis in samples doesn't translate intuitively in plot_noisecut_spectrograms. Keep in seconds.
        specs = AttribDict()
        specs.Original = trace
        specs.Full = S_full
        specs.Phase = phase
        specs.HPS = S_hps
        specs.Background = S_background
        specs.Removed = new
        specs.f = frequencies
        specs.t = times
        return hps_trace,specs
    else:
        return hps_trace