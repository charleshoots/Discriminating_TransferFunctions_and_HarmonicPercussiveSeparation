# Copyright 2019 Pascal Audet & Helen Janiszewski
#
# This file is part of OBStools.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
:mod:`~obstools.plot` contains several functions for plotting the results
of the analysis at various final and intermediate steps.

"""

import numpy as np
from matplotlib import pyplot as plt
from obstools.atacr import utils
from obspy import Trace
from obspy import Stream,UTCDateTime
from obspy.clients.syngine import Client as SynClient

def get_synthetic(stanm,event,components='Z',units='velocity',dt=0.2,evformat='%Y.%j.%H.%M'):
    synclient = SynClient()
    net,sta = stanm.split('.')
    evid=''.join(['GCMT:','C',UTCDateTime.strptime(event,evformat).strftime('%Y%m%d%H%M'),'A'])
    print(''.join(['Downloading synthetic event: ',evid,',  Components: ',components, ', Units: ', units]))
    st = synclient.get_waveforms(model='prem_a_2s', network=net, station=sta,eventid=evid,label=''.join([stanm,'_',event]),components=components,units=units,dt=dt)
    return st

def fig_QC(f, power, gooddays, ncomp, key='',diff2accel=True,mode='DayNoise'):
    """
    Function to plot the Quality-Control step of the analysis. This function
    is used in both the `atacr_daily_spectra` or `atacr_clean_spectra` scripts.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    power : :class:`~obstools.classes.Power`
        Container for the Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not
        (False)
    ncomp : int
        Number of components used in analysis (can be 2, 3 or 4)
    key : str
        String corresponding to the station key under analysis

    """

    sl_c11 = power.c11
    sl_c22 = power.c22
    sl_cZZ = power.cZZ
    sl_cPP = power.cPP

    if ncomp == 2:
        comps = ['Z','P']
        sls = [sl_cZZ, sl_cPP]
        title = ['HZ component, Station: '+key,
                 'HP component, Station: '+key]
    elif ncomp == 3:
        comps = ['1','2','Z']
        sls = [sl_c11, sl_c22, sl_cZZ]
        title = ['H1 component, Station: '+key,
                 'H2 component, Station: '+key,
                 'HZ component, Station: '+key]
    else:
        comps = ['1','2','Z','P',]
        sls = [sl_c11, sl_c22, sl_cZZ, sl_cPP]
        title = ['H1 component, Station: '+key,
                 'H2 component, Station: '+key,
                 'HZ component, Station: '+key,
                 'HP component, Station: '+key]

    # Extract only positive frequencies
    faxis = f > 0
    f = f[faxis]
    sls = [cc[faxis] if cc is not None else cc for cc in sls]
    if mode=='StaNoise':
        decrement = 2*10*np.log10(10) #aka...20
        sls = [cc + decrement if cc is not None else cc for cc in sls]
    # disp_to_accel = 40*np.log10(2*np.pi*f).reshape(-1,1)
    # sls = [g-np.mean(abs(g),axis=0) if g is not None else g for i,g in enumerate(sls)]
    if diff2accel:
        disp_to_accel = 40*np.log10(2*np.pi*f,where=f>0.).reshape(-1,1)
        sls = [g+disp_to_accel if ((not np.any(g==None)) and (comps[i]!='P')) else g for i,g in enumerate(sls)]
    sls = [g-np.mean(abs(g),axis=0) if (not np.any(g==None)) else g for i,g in enumerate(sls)]

    fig = plt.figure(6)
    for i, sl in enumerate(sls):
        good=sl[:, gooddays];good=good[:,~np.isinf(good.sum(axis=0))]
        ax = fig.add_subplot(ncomp, 1, i+1)
        ax.semilogx(f, good, 'k', lw=0.5)
        ax.set_xlim(1/150,5)
        if comps[i]=='P':
            ax.set_ylim([-150,150])
            ax.set_ylim(-max(np.abs([sl[~np.isinf(sl)].min(),sl[~np.isinf(sl)].max()])),max(np.abs([sl[~np.isinf(sl)].min(),sl[~np.isinf(sl)].max()])))
        elif np.isin(comps[i],['1','2']):ax.set_ylim([-200,5]);#ax.set_ylim([-150,0])
        else:ax.set_ylim([-200,5])
        # ax.set_ylim(-max(np.abs([sl.min(),sl.max()])),max(np.abs([sl.min(),sl.max()])))
        if np.sum(~gooddays) > 0:
            bad=sl[:, ~gooddays];bad=bad[:,~np.isinf(bad.sum(axis=0))]
            ax.semilogx(f,bad, 'r', lw=0.5)
        ax.set_title(title[i], fontdict={'fontsize': 8})
        if i == len(sls)-1:
            ax.set_xlabel('Frequency (Hz)', fontdict={'fontsize': 8})
        # ax.set_yticks(np.arange(ax.get_ylim()[0],ax.get_ylim()[1],tick_increment))
    plt.tight_layout()
    fig.set_size_inches([4,10])
    plt.tight_layout()

    return plt


def fig_average(f, power, bad, gooddays, ncomp, key='',diff2accel=True):
    """
    Function to plot the averaged spectra (those qualified as 'good' in the
    QC step). This function is used
    in both the `atacr_daily_spectra` or `atacr_clean_spectra` scripts.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    power : :class:`~obstools.classes.Power`
        Container for the Power spectra
    bad : :class:`~obstools.classes.Power`
        Container for the *bad* Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not
        (False)
    ncomp : int
        Number of components used in analysis (can be 2, 3 or 4)
    key : str
        String corresponding to the station key under analysis

    """

    c11 = power.c11
    c22 = power.c22
    cZZ = power.cZZ
    cPP = power.cPP
    bc11 = bad.c11
    bc22 = bad.c22
    bcZZ = bad.cZZ
    bcPP = bad.cPP


    if ncomp == 2:
        comps = ['Z','P']
        ccs = [cZZ, cPP]
        bcs = [bcZZ, bcPP]
        title = ['Average HZ, Station: '+key,
                 'Average HP, Station: '+key]
    elif ncomp == 3:
        comps = ['1','2','Z']
        ccs = [c11, c22, cZZ]
        bcs = [bc11, bc22, bcZZ]
        title = ['Average H1, Station: '+key,
                 'Average H2, Station: '+key,
                 'Average HZ, Station: '+key]
    else:
        comps = ['P','1','2','Z']
        ccs = [cPP, c11, c22, cZZ]
        bcs = [bcPP, bc11, bc22, bcZZ]
        title = ['Average HP, Station: '+key,
                 'Average H1, Station: '+key,
                 'Average H2, Station: '+key,
                 'Average HZ, Station: '+key]

    # Extract only positive frequencies
    faxis = f > 0
    f = f[faxis]
    decrement = 2*10*np.log10(10) #aka...20
    disp_to_accel = 40*np.log10(2*np.pi*f,where=f>0.)
    ccs = [cc[faxis] if cc is not None else cc for cc in ccs]
    bcs = [cc[faxis] if cc is not None else cc for cc in bcs]
    # Smooth prior to plotting
    ccs = [10*np.log10(cc,where=(cc > 0.)) if cc is not None else cc for cc in ccs]
    bcs = [10*np.log10(cc,where=(cc > 0.)) if cc is not None else cc for cc in bcs]
    ccs = [utils.smooth(cc, 50) if cc is not None else cc for cc in ccs]
    bcs = [utils.smooth(cc, 50) if cc is not None else cc for cc in bcs]
    if diff2accel:
        ccs = [g+disp_to_accel if g is not None and comps[i]!='P' else g for i,g in enumerate(ccs)]
        bcs = [g+disp_to_accel if g is not None and comps[i]!='P' else g for i,g in enumerate(bcs)]
    ccs = [c - np.mean(abs(c),axis=0) if c is not None else c for i,c in enumerate(ccs)]
    bcs = [c - np.mean(abs(c),axis=0) if c is not None else c for i,c in enumerate(bcs)]
    ccs = [g+decrement if g is not None and comps[i]!='P' else g for i,g in enumerate(ccs)]
    bcs = [g+decrement if g is not None and comps[i]!='P' else g for i,g in enumerate(bcs)]
    fig = plt.figure()
    for i, (cc, bc) in enumerate(zip(ccs, bcs)):
        ax = fig.add_subplot(ncomp, 1, i+1)
        ax.semilogx(f,cc, 'k', lw=0.5)
        ax.set_xlim(1/150,5)
        if np.sum(~gooddays) > 0:
            ax.semilogx(f, bc, 'r', lw=0.5)
        ax.set_title(title[i], fontdict={'fontsize': 8})
        if comps[i]=='P':
            ax.set_ylim(top=150)
            ax.set_ylim(bottom=-100)
            tick_increment = 50
        else:
            ax.set_ylim(top=-40)
            ax.set_ylim(bottom=-200)
            tick_increment = 20
        if i == len(ccs)-1:
            plt.xlabel('Frequency (Hz)', fontdict={'fontsize': 8})
        ax.set_yticks(np.arange(ax.get_ylim()[0],ax.get_ylim()[1],tick_increment))
    plt.tight_layout()
    fig.set_size_inches([4,10])
    plt.tight_layout()

    return plt


def fig_av_cross(f, field, gooddays, ftype, ncomp, key='',
                save=False, fname='', form='png', **kwargs):
    """
    Function to plot the averaged cross-spectra (those qualified as 'good' in
    the QC step). This function is used in the `atacr_daily_spectra` script.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    field : :class:`~obstools.classes.Rotation`
        Container for the Power spectra
    gooddays : List
        List of booleans representing whether a window is good (True) or not
        (False)
    ftype : str
        Type of plot to be displayed. If ftype is Admittance, plot is loglog.
        Otherwise semilogx
    key : str
        String corresponding to the station key under analysis
    **kwargs : None
        Keyword arguments to modify plot

    """

    # Extract only positive frequencies
    faxis = f > 0

    if ncomp == 2:
        fieldZP = field.cZP.T
        fields = [fieldZP]
        title = [': ZP']
        fig = plt.figure(figsize=(6, 2.667))
    elif ncomp == 3:
        field12 = field.c12.T
        field1Z = field.c1Z.T
        field2Z = field.c2Z.T
        fields = [field12, field1Z, field2Z]
        title = [': 12', ': 1Z', ': 2Z']
        fig = plt.figure(figsize=(6, 4))
    else:
        fieldZP = field.cZP.T
        field12 = field.c12.T
        field1Z = field.c1Z.T
        field2Z = field.c2Z.T
        field1P = field.c1P.T
        field2P = field.c2P.T
        fields = [field12, field1Z, field1P, field2Z, field2P, fieldZP]
        title = [': 12', ': 1Z', ': 1P', ': 2Z', ': 2P', ': ZP']
        fig = plt.figure(figsize=(6, 8))
    fields = [c[faxis] for c in fields]
    for i, field in enumerate(fields):
        ax = fig.add_subplot(len(fields), 1, i+1)
        # Extact field
        if ftype == 'Admittance':
            ax.loglog(
                f[faxis], field[:, gooddays], color='gray', **kwargs)
            if np.sum(~gooddays) > 0:
                ax.loglog(
                    f[faxis], field[:, ~gooddays], color='r', **kwargs)
        else:
            ax.semilogx(
                f[faxis], field[:, gooddays], color='gray', **kwargs)
            if np.sum(~gooddays) > 0:
                ax.semilogx(
                    f[faxis], field[:, ~gooddays], color='r', **kwargs)
        plt.ylabel(ftype, fontdict={'fontsize': 8})
        plt.title(key+' '+ftype+title[i], fontdict={'fontsize': 8})
        if i == len(fields)-1:
            plt.xlabel('Frequency (Hz)', fontdict={'fontsize': 8})

    plt.tight_layout()

    return plt


def fig_coh_ph(coh, ph, direc):
    """
    Function to plot the coherence and phase between the rotated H and Z
    components, used to characterize the tilt direction.

    Parameters
    ----------
    coh : :mod:`~numpy.ndarray`
        Coherence between rotated H and Z components
    ph : :mod:`~numpy.ndarray`
        Phase between rotated H and Z components
    direc : :mod:`~numpy.ndarray`
        Directions considered in maximizing coherence between H and Z

    """

    colors = plt.cm.cividis(np.linspace(0, 1, coh.shape[0]))

    if coh.ndim > 1:
        f, (ax1, ax2) = plt.subplots(1, 2)
        for i, (co, p) in enumerate(zip(coh, ph)):
            ax1.plot(direc, co, c=colors[i])
            ax2.plot(direc, p*180./np.pi, c=colors[i])
        ax1.set_ylabel('Coherence')
        ax1.set_ylim((0, 1.))
        ax2.set_ylabel('Phase')
        ax1.set_xlabel('Angle from H1')
        ax2.set_xlabel('Angle from H1')
        plt.tight_layout()

    else:
        plt.figure()
        plt.subplot(121)
        plt.plot(direc, coh, c=colors[0])
        plt.ylim((0, 1.))
        plt.subplot(122)
        plt.plot(direc, ph*180./np.pi, c=colors[0])
        plt.tight_layout()

    return plt


def fig_TF(f, day_trfs, day_list, sta_trfs, sta_list, skey='',append_title=''):
    """
    Function to plot the transfer functions available.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    day_trfs : Dict
        Dictionary containing the transfer functions for the daily averages
    day_list : Dict
        Dictionary containing the list of daily transfer functions
    sta_trfs : Dict
        Dictionary containing the transfer functions for the station averages
    sta_list : Dict
        Dictionary containing the list of average transfer functions
    skey : str
        String corresponding to the station key under analysis

    """

    import matplotlib.ticker as mtick
    # Extract only positive frequencies
    faxis = f > 0

    # Get max number of TFs to plot
    ntf = max(sum(day_list.values()), sum(sta_list.values()))

    # Define all possible compbinations
    tf_list = {'ZP': True, 'Z1': True, 'Z2-1': True,
               'ZP-21': True, 'ZH': True, 'ZP-H': True}

    if ntf == 1:
        fig = plt.figure(figsize=(6, 1.75))
    else:
        fig = plt.figure(figsize=(6, 1.33333333*ntf))

    j = 1
    for key in tf_list:

        if not day_list[key] and not sta_list[key]:
            continue

        ax = fig.add_subplot(ntf, 1, j)

        if day_list[key]:
            for i in range(len(day_trfs)):
                ax.loglog(f[faxis],np.abs(day_trfs[i][key]['TF_'+key][faxis]),'gray', lw=0.5)
        if sta_list[key]:
            ax.loglog(f[faxis],np.abs(sta_trfs[key]['TF_'+key][faxis]),'k', lw=0.5)
        if key == 'ZP':
            ax.set_ylim(1.e-10, 1.e-5)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(append_title + skey+' Transfer Function: ZP',fontdict={'fontsize': 8})
        elif key == 'Z1':
            ax.set_ylim(1.e-5, 1.e0)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(append_title + skey+' Transfer Function: Z1',fontdict={'fontsize': 8})
        elif key == 'Z2-1':
            ax.set_ylim(1.e-5, 1.e0)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(append_title + skey+' Transfer Function: Z2-1',fontdict={'fontsize': 8})
        elif key == 'ZP-21':
            ax.set_ylim(1.e-10, 1.e-5)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(append_title + skey+' Transfer Function: ZP-21',fontdict={'fontsize': 8})
        elif key == 'ZH':
            ax.set_ylim(1.e-5, 1.e0)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(append_title + skey+' Transfer Function: ZH',fontdict={'fontsize': 8})
        elif key == 'ZP-H':
            ax.set_ylim(1.e-10, 1.e-5)
            ax.set_xlim(1.e-4, 2.5)
            ax.set_title(append_title + skey+' Transfer Function: ZP-H',fontdict={'fontsize': 8})

        j += 1

    ax.set_xlabel('Frequency (Hz)')
    plt.tight_layout()

    return plt


def fig_comply(f, day_comps, day_list, sta_comps, sta_list, skey=None,
    elev=-1000., f_0=None):
    """
    Function to plot the transfer functions available.

    Parameters
    ----------
    f : :mod:`~numpy.ndarray`
        Frequency axis (in Hz)
    day_comps : Dict
        Dictionary containing the compliance functions for the daily averages
    day_list : Dict
        Dictionary containing the list of daily transfer functions
    sta_comps : Dict
        Dictionary containing the compliance functions for the station averages
    sta_list : Dict
        Dictionary containing the list of average transfer functions
    skey : str
        String corresponding to the station key under analysis
    elev : float
        Station elevation in meters (OBS stations have negative elevations)
    f_0 : float
        Lowest frequency to consider in plot (Hz)

    """

    import matplotlib.ticker as mtick
    import matplotlib.pyplot as plt

    # Extract only positive frequencies
    faxis = f > 0

    # Positive station elevation for frequency limit calc
    elev = -1.*elev

    # Calculate theoretical frequency limit for infra-gravity waves
    f_c = np.sqrt(9.81/np.pi/elev)/2.

    # Define all possible combinations
    comp_list = {'ZP': True, 'ZP-21': True, 'ZP-H': True}

    # Get max number of subplot
    nkeys_day = sum(day_list[key] for key in comp_list)
    nkeys_sta = sum(sta_list[key] for key in comp_list)
    ncomps = max(nkeys_day, nkeys_sta)

    if ncomps == 1:
        fig = plt.figure(figsize=(6, 1.75))
    else:
        fig = plt.figure(figsize=(6, 1.33333333*ncomps))

    for j, key in enumerate(comp_list):

        if not day_list[key] and not sta_list[key]:
            continue

        ax = fig.add_subplot(ncomps, 2, j*2+1)
        ax.tick_params(labelsize=8)
        ax.yaxis.get_offset_text().set_fontsize(8)

        if day_list[key]:
            for i in range(len(day_comps)):
                compliance = np.abs(day_comps[i][key][0])
                ax.plot(
                    f[faxis],
                    compliance[faxis],
                    'gray', alpha=0.3, lw=0.5)
                ax.set_xlim(f_0, f_c)
                ytop = np.max(compliance[(f > f_0) & (f < f_c)])
                ybot = np.min(compliance[(f > f_0) & (f < f_c)])
                ax.set_ylim(ybot, ytop)

        if sta_list[key]:
            ax.plot(
                f[faxis],
                np.abs(sta_comps[key][0][faxis]),
                'k', lw=0.5)

        if key == 'ZP':
            ax.set_title(skey+' Compliance: ZP',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-21':
            ax.set_title(skey+' Compliance: ZP-21',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-H':
            ax.set_title(skey+' Compliance: ZP-H',
                         fontdict={'fontsize': 8})

        if f_0:
            ax.axvline(f_0, ls='--', c='k', lw=0.75)
        ax.axvline(f_c, ls='--', c='k', lw=0.75)

        ax = fig.add_subplot(ncomps, 2, j*2+2)
        ax.tick_params(labelsize=8)

        if day_list[key]:
            for i in range(len(day_comps)):
                ax.semilogx(
                    f[faxis],
                    np.abs(day_comps[i][key][1][faxis]),
                    'gray', alpha=0.3, lw=0.5)
        if sta_list[key]:
            ax.semilogx(
                f[faxis], 
                np.abs(sta_comps[key][1][faxis]),
                'k', lw=0.5)

        if key == 'ZP':
            ax.set_title(skey+' Coherence: ZP',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-21':
            ax.set_title(skey+' Coherence: ZP-21',
                         fontdict={'fontsize': 8})
        elif key == 'ZP-H':
            ax.set_title(skey+' Coherence: ZP-H',
                         fontdict={'fontsize': 8})

        if f_0:
            ax.axvline(f_0, ls='--', c='k', lw=0.75)
        ax.axvline(f_c, ls='--', c='k', lw=0.75)

    axes = plt.gcf().get_axes()
    axes[-2].set_xlabel('Frequency (Hz)', fontsize=8)
    axes[-1].set_xlabel('Frequency (Hz)', fontsize=8)

    plt.tight_layout()

    return plt


def fig_event_raw(evstream, fmin=1./150., fmax=2.):
    """
    Function to plot the raw (although bandpassed) seismograms.

    Parameters
    ----------
    evstream : :class:`~obtsools.classes.EventStream`
        Container for the event stream data
    fmin : float
        Low frequency corner (in Hz)
    fmax : float
        High frequency corner (in Hz)

    """

    from obspy import Stream

    # Unpack traces
    tr1 = evstream.tr1.copy()
    tr2 = evstream.tr2.copy()
    trZ = evstream.trZ.copy()
    trP = evstream.trP.copy()
    st = Stream(traces=[tr for tr in [tr1, tr2, trZ, trP] if np.any(tr.data)])

    st.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    stats = trZ.stats
    sr = stats.sampling_rate

    taxis = np.arange(0., stats.npts/sr, 1./sr)

    fig = plt.figure(figsize=(6, 6))

    ax = fig.add_subplot(4, 1, 1)
    ax.plot(taxis, trZ.data, 'k', lw=0.5)
    ax.set_title(evstream.key + ' ' + evstream.tstamp +
                 ': Z', fontdict={'fontsize': 8})
    ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                        scilimits=(-3, 3))
    ax.set_xlim((0., stats.npts/sr))

    if len(tr1.data > 0):
        ax = fig.add_subplot(4, 1, 2)
        ax.plot(taxis, tr1.data, 'k', lw=0.5)
        ax.set_xlim((0., 7200.))
        ax.set_title(evstream.tstamp + ': 1', fontdict={'fontsize': 8})
        ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                            scilimits=(-3, 3))

        ax = fig.add_subplot(4, 1, 3)
        ax.plot(taxis, tr2.data, 'k', lw=0.5)
        ax.set_xlim((0., stats.npts/sr))
        ax.set_title(evstream.tstamp + ': 2', fontdict={'fontsize': 8})
        ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                            scilimits=(-3, 3))

    if len(trP.data > 0):
        if len(tr1.data > 0):
            ax = fig.add_subplot(4, 1, 4)
        else:
            ax = fig.add_subplot(4, 1, 2)
        ax.plot(taxis, trP.data, 'k', lw=0.5)
        ax.ticklabel_format(axis='y', style='sci', useOffset=True,
                            scilimits=(-3, 3))
        ax.set_xlim((0., stats.npts/sr))
        ax.set_title(evstream.tstamp + ': P', fontdict={'fontsize': 8})

    plt.xlabel('Time since earthquake (sec)')
    plt.tight_layout()

    return plt


# def fig_event_corrected(evstream, TF_list, fmin=1./250., fmax=2.,taper=0.1):
#     """
#     Function to plot the corrected vertical component seismograms.
#     Parameters
#     ----------
#     evstream : :class:`~obtsools.classes.EventStream`
#         Container for the event stream data
#     Tf_list : list
#         List of Dictionary elements of transfer functions used
#         for plotting the corrected vertical component.
#     """
#     # trZ = evstream.trZ.copy()
#     # trZ.taper(taper,side='both')
#     # trZ.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
#     # corrected = dict()
#     # [corrected.update({key:Trace(data=evstream.correct[key],header=stats).copy()}) for key in list(evstream.correct.keys())]
#     # [corrected[c].taper(taper,side='both') for c in list(corrected.keys())]
#     # [corrected[c].filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True) for c in list(corrected.keys())]
#     # Unpack vertical trace and filter
#     keys = ['Z1','Z2-1','ZP-21','ZH','ZP-H','ZP']
#     # Aggregate
#     st = Stream()
#     st.append(evstream.trZ.copy())
#     [st.append(Trace(data=evstream.correct[key].copy(),header=st[0].stats).copy()) for key in list(evstream.correct.keys())]

#     # Preproc
#     # st.taper(taper,side='both')
#     st.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)

#     # Disagregate
#     raw,st = st[0],st[1:]
#     corrected = dict()
#     dummyindex = raw.copy()
#     dummyindex.data = dummyindex.data*np.nan
#     [corrected.update({key:dummyindex}) for key in keys]
#     [corrected.update({key:tr}) for key,tr in zip(list(evstream.correct.keys()),st)]
#     taxis = raw.times()

#     # fig, axes = plt.subplots(nrows=6, ncols=1,figsize=(8, 8),layout='constrained',squeeze=True,sharey='col',sharex='col')
#     fig, axes = plt.subplots(nrows=6, ncols=1,figsize=(8, 8),squeeze=True,sharey='col',sharex='col')
#     # ylim = [raw.data.min(),raw.data.max()]
#     # ylim = [np.min(np.array([corrected[c].data for c in list(evstream.correct.keys())]).min(axis=0)),max(np.array([corrected[c].data for c in list(evstream.correct.keys())]).max(axis=0))]
#     # [a.set_ylim(ylim) for a in axes]
#     [a.set_xlim((0., taxis[-1])) for a in axes]
#     [a.plot(taxis, raw.data, 'lightgray', lw=0.5) for a in axes]
#     [a.plot(taxis, corrected[c].data, 'k', lw=0.5) for c,a in zip(keys,axes)]
#     [a.set_title(evstream.key + ' ' + evstream.tstamp +':' + key, fontdict={'fontsize': 8}) for key,a in zip(keys,axes)]
#     [a.ticklabel_format(axis='y', style='sci', useOffset=True,scilimits=(-3, 3)) for a in axes]
#     plt.xlabel('Time since earthquake (sec)')
#     plt.tight_layout()
#     return fig


def fig_event_corrected(evstream, TF_list, fmin=1./250., fmax=2,synthetics=False,taper = 0.05):
    """
    Function to plot the corrected vertical component seismograms.

    Parameters
    ----------
    evstream : :class:`~obtsools.classes.EventStream`
        Container for the event stream data
    Tf_list : list
        List of Dictionary elements of transfer functions used
        for plotting the corrected vertical component.

    """
    # Unpack vertical trace and filter
    # fig,axes = plt.subplots(nrows=len(st), ncols=1,figsize=(8,13),layout='constrained',squeeze=False,sharey='all',sharex='all')
    # [ax.plot(s.times(),s.data) for ax,s in zip(axes,st)]

    trZ = evstream.trZ.copy()
    stats = trZ.stats
    depth = round(abs(trZ.stats.sac.stel)*1000)
    stanm = '.'.join([stats.network,stats.station])
    trimS,trimE = stats.starttime,stats.endtime
    if synthetics:
        synth = get_synthetic(stanm,stats.starttime.strftime('%Y.%j.%H.%M'))[0]
        fmin,fmax = 1/100,0.5
        synth.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
        trimS,trimE = stats.starttime,stats.starttime+3600
        synth.trim(trimS,trimE)
    # ----------------------------------
    keys = ['Z1','Z2-1','ZP-21','ZH','ZP-H','ZP']
    # Aggregate
    st = Stream()
    st.append(evstream.trZ.copy())
    [st.append(Trace(data=evstream.correct[key].copy(),header=st[0].stats).copy()) for key in list(evstream.correct.keys())]
    # Preproc
    if taper>0:
        st.taper(taper,side='both')
    st.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    st.trim(trimS,trimE)
    # Disagregate
    trZ,st = st[0],st[1:]
    corrected = dict()
    dummyindex = trZ.copy()
    dummyindex.data = dummyindex.data*np.nan
    [corrected.update({key:dummyindex}) for key in keys]
    [corrected.update({key:tr}) for key,tr in zip(list(evstream.correct.keys()),st)]

    # ylim = [j([j(corrected[k].data) for k in list(evstream.correct.keys())]) for j in [np.min,np.max]]
    ylim = [j(trZ.data) for j in [np.min,np.max]]
    # if synthetics:
    #     ylim = [j(synth.data) for j in [np.min,np.max]]
    # ----------------------------------

    taxis = st[0].times()
    tstamp = evstream.tstamp

    plt.figure(figsize=(6, 8))

    for row,key in enumerate(keys):
        plt.subplot(611 + row)
        plt.plot(taxis, trZ.data, 'lightgray', lw=0.5)
        plt.plot(taxis,corrected[key].data,'k',lw=0.5)
        if synthetics:
            plt.plot(synth.times(),synth.data,'r:',lw=0.5)
        plt.title(stanm + ' ,' + str(depth) + 'm |' + tstamp + ': ' + key, fontdict={'fontsize': 8})
        plt.gca().ticklabel_format(axis='y', style='sci', useOffset=True,scilimits=(-3, 3))
        plt.xlim(taxis[0],taxis[-1])
        plt.ylim(ylim)
    plt.xlabel('Time since earthquake (sec)')
    plt.tight_layout()

    return plt