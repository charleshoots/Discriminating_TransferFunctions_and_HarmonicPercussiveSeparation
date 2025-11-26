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
:mod:`~obstools.atacr.utils` contains several functions that are used in the
class methods of `~obstools.atacr.classes`.

"""


import os
import math
import numpy as np
import fnmatch
from matplotlib import pyplot as plt
from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
from obspy.core.inventory.inventory import read_inventory
from pathlib import Path
import warnings,itertools,re
warnings.filterwarnings('ignore')

def unravel(lst):return list(itertools.chain.from_iterable(lst))
def atacr_status(step,station,stage=None,note=None):
    s = f'[{str(step)}] || [{str(station)}]'
    if stage:s=f'{s} | {str(stage)}'
    if note:s=f'{s} || {str(note)}'
    print(s)

def traceshift(trace, tt):
    """
    Function to shift traces in time given travel time


    Parameters
    ----------

    trace : :class:`~obspy.core.Trace` object
        Trace object to update
    tt : float
        Time shift in seconds

    Returns
    -------

    rtrace : :class:`~obspy.core.Trace` object
        Updated trace object


    """

    # Define frequencies
    nt = trace.stats.npts
    dt = trace.stats.delta
    freq = np.fft.fftfreq(nt, d=dt)

    # Fourier transform
    ftrace = np.fft.fft(trace.data)

    # Shift
    for i in range(len(freq)):
        ftrace[i] = ftrace[i]*np.exp(-2.*np.pi*1j*freq[i]*tt)

    # Back Fourier transform and return as trace
    rtrace = trace.copy()
    rtrace.data = np.real(np.fft.ifft(ftrace))

    # Update start time
    rtrace.stats.starttime -= tt

    return rtrace


def QC_streams(start, end, st):
    """
    Function for quality control of traces, which compares the
    start and end times that were requested, as well as the total n
    length of the traces.


    Parameters
    ----------

    start : :class:`~obspy.core.UTCDateTime` object
        Start time of requested stream
    end : :class:`~obspy.core.UTCDateTime` object
        End time of requested stream
    st : :class:`~obspy.core.Stream` object
        Stream object with all trace data

    Returns
    -------

    (pass): bool
        Whether the QC test has passed
    st : :class:`~obspy.core.Stream` object
        Updated stream object


    """

    # Check start times
    if not np.all([tr.stats.starttime == start for tr in st]):
        print("* Start times are not all close to true start: ")
        [print("*   "+tr.stats.channel+" " +
               str(tr.stats.starttime)+" " +
               str(tr.stats.endtime)) for tr in st]
        print("*   True start: "+str(start))
        print("* -> Shifting traces to true start")
        delay = [tr.stats.starttime - start for tr in st]
        st_shifted = Stream(
            traces=[traceshift(tr, dt) for tr, dt in zip(st, delay)])
        st = st_shifted.copy()

    # Try trimming
    dt = st[0].stats.delta
    try:
        st.trim(start, end-dt, fill_value=0., pad=True)
    except Exception:
        print("* Unable to trim")
        print("* -> Skipping")
        print("**************************************************")
        return False, None

    # Check final lengths - they should all be equal if start times
    # and sampling rates are all equal and traces have been trimmed
    sr = st[0].stats.sampling_rate
    if not np.allclose([tr.stats.npts for tr in st[1:]], st[0].stats.npts):
        print("* Lengths are incompatible: ")
        [print("*     "+str(tr.stats.npts)) for tr in st]
        print("* -> Skipping")
        print("**************************************************")

        return False, None

    elif not np.allclose([st[0].stats.npts], int((end - start)*sr), atol=1):
        print("* Length is too short: ")
        print("*    "+str(st[0].stats.npts) +
              " ~= "+str(int((end - start)*sr)))
        print("* -> Skipping")
        print("**************************************************")

        return False, None
    elif np.any([np.any(np.isnan(tr.data)) for tr in st]):
        print("* Dead traces(NaNs reported): -> Skipping")
        return False, None

    else:
        return True, st


def update_stats(tr, stla, stlo, stel, cha, evla=None, evlo=None):
    """
    Function to include SAC metadata to :class:`~obspy.core.Trace` objects
    Parameters
    ----------
    tr : :class:`~obspy.core.Trace` object
        Trace object to update
    stla : float
        Latitude of station
    stlo : float
        Longitude of station
    stel : float
        Station elevation (m)
    cha : str
        Channel for component
    evla : float, optional
        Latitude of event
    evlo : float, optional
        Longitute of event
    Returns
    -------
    tr : :class:`~obspy.core.Trace` object
        Updated trace object
    """

    tr.stats.sac = AttribDict()
    tr.stats.sac.stla = stla
    tr.stats.sac.stlo = stlo
    tr.stats.sac.stel = stel
    tr.stats.sac.kcmpnm = cha
    tr.stats.channel = cha
    if evla is not None and evlo is not None:
        tr.stats.sac.evla = evla
        tr.stats.sac.evlo = evlo
    return tr

def get_event(eventpath,ovr=None):
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
    corrpath=Path('EVENTS')/'corrected'/eventpath.name
    # if np.isin('rmresp',[e.name for e in list(eventpath.parents)]):rmresp=False
    # else:rmresp=True
    # if evlist is None:
    # Find out how many events from Z.SAC files
    # Extract events from time stamps
    stanm=eventpath.name
    files = [x for x in list(eventpath.glob('*Z.SAC')) if x.is_file()]
    if not files:raise(Exception("No event found in folder "+str(eventpath)))
    for evi,file in enumerate(files):
        event='.'.join(file.name.split('.')[:-2])
        print(f'CorrectEvent (A7) |{stanm}{'='*20} EVENT [ {str(evi)}/{str(len(files))}: {event}{'='*20}')
        if not ovr:
            if np.any(len(list(corrpath.glob(f'*{'.*'.join(file.name.split('.'))}')))>0):
                print(f'-- Output exists already. Skipping event.')
                continue
        repstr='.'.join((str(file).split('.')[-2:]))
        trfiles = [Path(str(file).replace(repstr,f'*{c}.SAC')) for c in ['1','2','Z','H']]
        tr1,tr2,trZ,trP = load_data(trfiles)
        yield tr1,tr2,trZ,trP

def check_fs(st):
    rates=[tr.stats.sampling_rate for tr in st]
    inequal_rates=len(np.unique([tr.stats.sampling_rate for tr in st]))>1
    if inequal_rates:
        new_fs=min(rates) #choose the lowest sample rate from all given traces
        [tr.resample(new_fs, no_filter=False) for tr in st]
    return st

def calculate_tilt(ft1, ft2, ftZ, ftP, f, goodwins, tiltfreq=[0.005, 0.035]):
    """
    Determines tilt direction from maximum coherence between rotated H1 and Z.

    Parameters
    ----------
    ft1, ft2, ftZ, ftP : :class:`~numpy.ndarray`
        Fourier transform of corresponding H1, H2, HZ and HP components
    f : :class:`~numpy.ndarray`
        Frequency axis in Hz
    goodwins : list
        List of booleans representing whether a window is good (True) or not
        (False). This attribute is returned from the method
        :func:`~obstools.atacr.classes.DayNoise.QC_daily_spectra`
    tiltfreq : list, optional
        Two floats representing the frequency band at which the tilt is
        calculated

    Returns
    -------
    cHH, cHZ, cHP : :class:`~numpy.ndarray`
        Arrays of power and cross-spectral density functions of components HH
        (rotated H1 in direction of maximum tilt), HZ, and HP
    coh : :class:`~numpy.ndarray`
        Coherence value between rotated H and Z components, as a function of
        directions (azimuths)
    ph : :class:`~numpy.ndarray`
        Phase value between rotated H and Z components, as a function of
        directions (azimuths)
    direc : :class:`~numpy.ndarray`
        Array of directions (azimuths) considered
    tilt : float
        Direction (azimuth) of maximum coherence between rotated H1 and Z
    coh_value : float
        Coherence value at tilt direction
    phase_value : float
        Phase value at tilt direction

    """

    direc = np.arange(0., 360., 10.)
    coh = np.zeros(len(direc))
    ph = np.zeros(len(direc))
    cZZ = np.abs(np.mean(ftZ[goodwins, :] *
                         np.conj(ftZ[goodwins, :]), axis=0))[0:len(f)]

    for i, d in enumerate(direc):

        # Rotate horizontals
        ftH = rotate_dir(ft1, ft2, d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins, :] *
                             np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins, :] *
                      np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ph = phase(cHZ)

        # Calculate coherence over frequency band
        coh[i] = np.mean(Co[(f > tiltfreq[0]) & (f < tiltfreq[1])])
        ph[i] = np.pi/2. - np.mean(Ph[(f > tiltfreq[0]) & (f < tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(coh == coh.max())

    # Phase and direction at maximum coherence
    phase_value = ph[ind[0]][0]
    coh_value = coh[ind[0]][0]
    tilt = direc[ind[0]][0]

    # Refine search
    rdirec = np.arange(direc[ind[0]][0]-10., direc[ind[0]][0]+10., 1.)
    rcoh = np.zeros(len(direc))
    rph = np.zeros(len(direc))

    for i, d in enumerate(rdirec):

        # Rotate horizontals
        ftH = rotate_dir(ft1, ft2, d)

        # Get transfer functions
        cHH = np.abs(np.mean(ftH[goodwins, :] *
                             np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
        cHZ = np.mean(ftH[goodwins, :] *
                      np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]

        Co = coherence(cHZ, cHH, cZZ)
        Ph = phase(cHZ)

        # Calculate coherence over frequency band
        rcoh[i] = np.mean(Co[(f > tiltfreq[0]) & (f < tiltfreq[1])])
        rph[i] = np.pi/2. - np.mean(Ph[(f > tiltfreq[0]) & (f < tiltfreq[1])])

    # Index where coherence is max
    ind = np.argwhere(rcoh == rcoh.max())

    # Phase and direction at maximum coherence
    phase_value = rph[ind[0]][0]
    coh_value = rcoh[ind[0]][0]
    tilt = rdirec[ind[0]][0]

    # Phase has to be close to zero - otherwise add pi
    if phase_value > 0.5*np.pi:
        tilt += 180.
    if tilt > 360.:
        tilt -= 360.

    # print('Maximum coherence for tilt = ', tilt)

    # Now calculate spectra at tilt direction
    ftH = rotate_dir(ft1, ft2, tilt)

    # Get transfer functions
    cHH = np.abs(np.mean(ftH[goodwins, :] *
                         np.conj(ftH[goodwins, :]), axis=0))[0:len(f)]
    cHZ = np.mean(ftH[goodwins, :]*np.conj(ftZ[goodwins, :]), axis=0)[0:len(f)]
    if np.any(ftP):
        cHP = np.mean(ftH[goodwins, :] *
                      np.conj(ftP[goodwins, :]), axis=0)[0:len(f)]
    else:
        cHP = None

    return cHH, cHZ, cHP, coh, ph, direc, tilt, coh_value, phase_value


def smooth(data, nd, axis=0):
    """
    Function to smooth power spectral density functions from the convolution
    of a boxcar function with the PSD

    Parameters
    ----------
    data : :class:`~numpy.ndarray`
        Real-valued array to smooth (PSD)
    nd : int
        Number of samples over which to smooth
    axis : int, optional
        axis over which to perform the smoothing

    Returns
    -------
    filt : :class:`~numpy.ndarray`, optional
        Filtered data

    """
    if np.any(data):
        if data.ndim > 1:
            filt = np.zeros(data.shape)
            for i in range(data.shape[::-1][axis]):
                if axis == 0:
                    filt[:, i] = np.convolve(
                        data[:, i], np.ones((nd,))/nd, mode='same')
                elif axis == 1:
                    filt[i, :] = np.convolve(
                        data[i, :], np.ones((nd,))/nd, mode='same')
        else:
            filt = np.convolve(data, np.ones((nd,))/nd, mode='same')
        return filt
    else:
        return None


def admittance(Gxy, Gxx):
    """
    Calculates admittance between two components
    Parameters
    ---------
    Gxy : :class:`~numpy.ndarray`
        Cross spectral density function of `x` and `y`
    Gxx : :class:`~numpy.ndarray`
        Power spectral density function of `x`
    Returns
    -------
    : :class:`~numpy.ndarray`, optional
        Admittance between `x` and `y`
    """
    if np.any(Gxy) and np.any(Gxx):return np.abs(Gxy)/Gxx
    else:return None

def coherence(Gxy, Gxx, Gyy):
    """
    Calculates coherence between two components
    Parameters
    ---------
    Gxy : :class:`~numpy.ndarray`
        Cross spectral density function of `x` and `y`
    Gxx : :class:`~numpy.ndarray`
        Power spectral density function of `x`
    Gyy : :class:`~numpy.ndarray`
        Power spectral density function of `y`
    Returns
    -------
    : :class:`~numpy.ndarray`, optional
        Coherence between `x` and `y`
    """
    if np.any(Gxy) and np.any(Gxx) and np.any(Gxx):return np.abs(Gxy)**2/(Gxx*Gyy)
    else:return None

def phase(Gxy):
    """
    Calculates phase angle between two components
    Parameters
    ---------
    Gxy : :class:`~numpy.ndarray`
        Cross spectral density function of `x` and `y`
    Returns
    -------
    : :class:`~numpy.ndarray`, optional
        Phase angle between `x` and `y`
    """
    if np.any(Gxy):return np.angle(Gxy)
    else:return None

def rotate_dir(tr1, tr2, direc):

    d = -direc*np.pi/180.+np.pi/2.
    rot_mat = np.array([[np.cos(d), -np.sin(d)],
                        [np.sin(d), np.cos(d)]])

    v12 = np.array([tr2, tr1])
    vxy = np.tensordot(rot_mat, v12, axes=1)
    tr_2 = vxy[0, :]
    tr_1 = vxy[1, :]

    return tr_1

def ftest(res1, pars1, res2, pars2):

    from scipy.stats import f as f_dist

    N1 = len(res1)
    N2 = len(res2)

    dof1 = N1 - pars1
    dof2 = N2 - pars2

    Ea_1 = np.sum(res1**2)
    Ea_2 = np.sum(res2**2)

    Fobs = (Ea_1/dof1)/(Ea_2/dof2)

    P = 1. - (f_dist.cdf(Fobs, dof1, dof2) - f_dist.cdf(1./Fobs, dof1, dof2))

    return P

def date_from_filestr(file,stanm=None):
    warnings.filterwarnings('ignore')
    file=file.replace('.SAC','')#Remove format
    if stanm:file=file.replace(stanm,'') #Remove station names
    file=[file.replace(r,'') for r in ['.HZ','.HDH','.H1','.H2'] if file.find(r)>-1][0] #Remove channel names
    file=re.sub(".[^0-9+.]\D+", '',file) #Parse out the date
    return file

def search_files(datapath, tstart, tend):
    """
    Works on unique DAY-LONG files.
    Returns a list of 4-length tuples (1 for each comp) containing the files that fit within tstart,tend
    -CHoots, 2023
    """
    # Define empty streams
    trN1_files = [];trN2_files = [];trNZ_files = [];trNP_files = []
    t1 = tstart # Time iterator
    files = [x for x in list(datapath.glob('*.*')) if x.is_file()] #Get a general list of all files
    while t1 < tend: # Cycle through each day within time range
        # Time stamp used in file name
        tstamp=t1.strftime('%Y.%j.')
        # Cycle through directory and load files
        f=[i for i in [f for f in files if ((str(f).find('1.SAC')>=0) & (str(f).find(tstamp)>=0))]]
        if len(f)==1:trN1_files.append(f)
        f=[i for i in [f for f in files if ((str(f).find('2.SAC')>=0) & (str(f).find(tstamp)>=0))]]
        if len(f)==1:trN2_files.append(f)
        f=[i for i in [f for f in files if ((str(f).find('Z.SAC')>=0) & (str(f).find(tstamp)>=0))]]
        if len(f)==1:trNZ_files.append(f)
        f=[i for i in [f for f in files if ((str(f).find('HDH.SAC')>=0) & (str(f).find(tstamp)>=0))]]
        if len(f)==1:trNP_files.append(f)
        # Increase increment
        t1 += 3600.*24.
    if not np.all([len([date_from_filestr(str(f),datapath.name) for f in c]) for c in [trN1_files,trN2_files,trNZ_files,trNP_files]]):raise Exception('File mismatch')
    files = [[tr1[0],tr2[0],trZ[0],trP[0]] for tr1,tr2,trZ,trP in zip(trN1_files,trN2_files,trNZ_files,trNP_files)]
    return files

def run_rmresp(tr,inv):
    #Runs OBSPY's rmresp function following the deignated defaults below for either pressure or seismic data.
    #input: OBSPY Trace (tr) and Inventory (inv) objects.
    #output: Trace object.
    RespConfig=AttribDict();RespConfig.Seismic=AttribDict();RespConfig.Pressure=AttribDict()
    RespConfig.Seismic.pre_filt=[0.001, 0.002, 45.0, 50.0]
    RespConfig.Seismic.units="DISP"
    RespConfig.Seismic.water_level=60
    RespConfig.Seismic.zero_mean=True   
    RespConfig.Pressure.pre_filt=[0.001, 0.002, 45.0, 50.0]
    RespConfig.Pressure.units="DEF"
    RespConfig.Pressure.water_level=None
    RespConfig.Pressure.zero_mean=False
    if tr.stats.channel=='HDH':Config=RespConfig.Pressure
    else:Config=RespConfig.Seismic
    tr.remove_response(inventory=inv,pre_filt=Config.pre_filt,output=Config.units,water_level=Config.water_level,zero_mean=Config.zero_mean)
    return tr

def load_data(files,rmresp=None):
    if not isinstance(files,list):files=[files]
    """
    Loads four files defined by the tuple, files.
    -CHoots, 2023
    """
    # Implicitly define (was already done if a parent folder is named 'rmresp') whether remove response is needed.
    if rmresp is None:
        if np.any(np.isin('rmresp',unravel([[e.name for e in list(f.parents)] for f in files]))):
            rmresp=False
        else:rmresp=True;inv=read_inventory(Path(str(files[0])).parent / '*_inventory.xml')
    #Load data
    st = [read(str(f))[0] for f in files]
    #Remove response
    if rmresp:st = [run_rmresp(tr,inv) for tr in st]
    #Confirm equal sample rates.
    st=check_fs(st)
    return st

def get_data_generator(datapath, tstart, tend,skipfiles=None,rmresp=None):
    """
    The same as get_data but now factored as a generator object such that the data is never loaded until its index is called.
    Generators can only be iterated, and activated by use in a for loop
    Significantly more memory efficient.
    -CHoots, 2023
    """
    file_list = search_files(datapath, tstart, tend)
    for files in file_list:
        if (skipfiles is not None) & np.any(np.isin(files,skipfiles)):continue
        yield load_data(files,rmresp=rmresp)

def make_inventory(stats,response):
    from obspy.core.inventory import Inventory, Network, Station, Channel, Site
    if 'sac' in list(stats.__dict__.keys()):instrument_key = 'sac'
    else: instrument_key = 'mseed'
    # We'll first create all the various objects. These strongly follow the
    # hierarchy of StationXML files.
    inv = Inventory(
        # We'll add networks later.
        # resource_id='station_resource_id',
        networks=[],
        # The source should be the id whoever create the file.
        source="")
    net = Network(
        # This is the network code according to the SEED standard.
        code=stats.network,
        # A list of stations. We'll add one later.
        stations=[],
        description="",
        # Start-and end dates are optional.
        start_date=stats.starttime)
    sta = Station(
        # This is the station code according to the SEED standard.
        code=stats.station,
        latitude=stats[instrument_key].stla,
        longitude=stats[instrument_key].stlo,
        elevation=stats[instrument_key].stel*1000,
        creation_date=stats.starttime,
        site=Site(name=""))
    cha = Channel(
        # This is the channel code according to the SEED standard.
        code=stats.channel,
        # This is the location code according to the SEED standard.
        location_code=stats.location,
        # Note that these coordinates can differ from the station coordinates.
        latitude=stats[instrument_key].stla,
        longitude=stats[instrument_key].stlo,
        elevation=stats[instrument_key].stel*1000,
        depth=abs(stats[instrument_key].stel*1000),
        azimuth=0.0,
        dip=-90.0,
        sample_rate=1/stats.delta)
    # Now tie it all together.
    cha.response = response
    sta.channels.append(cha)
    net.stations.append(sta)
    inv.networks.append(net)
    return inv

def save_inventory(filename,st):
    tr = st.select(channel='*Z')
    if len(tr)==0:tr = st[0]
    else:tr = tr[0]
    inventory = make_inventory(tr.copy().stats,tr.copy().stats.response)
    for tr in st:
        inventory+= make_inventory(tr.copy().stats,tr.copy().stats.response)
    print('...saving station inventory')
    inventory.write(filename,format="STATIONXML")






# depreciated codes

# def load_sac(file,rmresp=False,inv=[],
#         pressure_pre_filt=[0.001, 0.002, 45.0, 50.0],
#         seismic_pre_filt=[0.001, 0.002, 45.0, 50.0],
#         seismic_units="DISP",pressure_units="DEF",
#         pressure_water_level=None,seismic_water_level=60):
#         # ----------------------
#     if rmresp:inv=read_inventory(Path(str(file)).parent / '*_inventory.xml')
#     try:
#         if fnmatch.fnmatch(str(file),'*1.SAC'):
#                 tr = read(str(file))
#                 if rmresp:
#                        for t in tr:t.stats.location=''
#                        tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
#         elif fnmatch.fnmatch(str(file),'*2.SAC'):
#                 tr = read(str(file))
#                 if rmresp:
#                        for t in tr:t.stats.location=''
#                        tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
#         elif fnmatch.fnmatch(str(file),'*Z.SAC'):
#                 tr = read(str(file))
#                 if rmresp:
#                        for t in tr:t.stats.location=''
#                        tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level,hide_sensitivity_mismatch_warning=True)
#         elif fnmatch.fnmatch(str(file),'*H.SAC'):
#                 tr = read(str(file))
#                 if rmresp:
#                        for t in tr:t.stats.location=''
#                        tr.remove_response(inventory=inv,pre_filt=pressure_pre_filt,output=pressure_units,water_level=pressure_water_level,hide_sensitivity_mismatch_warning=True)
#         return tr,inv
#     except:
#         print('WARNING:Load Error')
#         return Stream(),inv

# def get_data(datapath, tstart, tend,seismic_pre_filt=[0.001, 0.002, 45.0, 50.0], pressure_pre_filt=[0.001, 0.002, 45.0, 50.0],seismic_units="DISP",pressure_units="DEF",pressure_water_level=None,seismic_water_level=60):
#     """
#     Function to grab all available noise data given a path and data time range

#     Parameters
#     ----------
#     datapath : str
#         Path to noise data folder
#     tstart : :class:`~obspy.class.UTCDateTime`
#         Start time for query
#     tend : :class:`~obspy.class.UTCDateTime`
#         End time for query

#     Returns
#     -------
#     tr1, tr2, trZ, trP : :class:`~obspy.core.Trace` object
#         Corresponding trace objects for components H1, H2, HZ and HP. Returns
#         empty traces for missing components.

#     """

#     # Define empty streams
#     trN1 = Stream()
#     trN2 = Stream()
#     trNZ = Stream()
#     trNP = Stream()

#     # Time iterator
#     t1 = tstart
    
#     # Cycle through each day within time range
#     while t1 < tend:

#         # Time stamp used in file name
#         tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'

#         # Cycle through directory and load files
#         p = datapath.glob('*.*')
#         files = [x for x in p if x.is_file()]
#         for file in files:
#             inv = read_inventory(Path(str(file)).parent / '*_inventory.xml')
#             if fnmatch.fnmatch(str(file), '*' + tstamp + '*1.SAC'):
#                 tr = read(str(file))
#                 tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level)
#                 trN1.append(tr[0])
#             elif fnmatch.fnmatch(str(file), '*' + tstamp + '*2.SAC'):
#                 tr = read(str(file))
#                 tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level)
#                 trN2.append(tr[0])
#             elif fnmatch.fnmatch(str(file), '*' + tstamp + '*Z.SAC'):
#                 tr = read(str(file))
#                 tr.remove_response(inventory=inv,pre_filt=seismic_pre_filt, output=seismic_units,water_level=seismic_water_level)
#                 trNZ.append(tr[0])
#             elif fnmatch.fnmatch(str(file), '*' + tstamp + '*H.SAC'):
#                 tr = read(str(file))
#                 tr.remove_response(inventory=inv,pre_filt=pressure_pre_filt,output=pressure_units,water_level=pressure_water_level)
#                 trNP.append(tr[0])

#         # Increase increment
#         t1 += 3600.*24.

#     # Fill with empty traces if components are not found
#     ntr = len(trNZ)
#     if not trN1 and not trN2:
#         for i in range(ntr):
#             trN1.append(Trace())
#             trN2.append(Trace())
#     if not trNP:
#         for i in range(ntr):
#             trNP.append(Trace())

#     if ntr > 0:
#         # Check that all sampling rates are equal - otherwise resample
#         if trNZ[0].stats.sampling_rate != trNP[0].stats.sampling_rate:

#             # These checks assume that all seismic data have the same sampling
#             if trNZ[0].stats.sampling_rate < trNP[0].stats.sampling_rate:
#                 trNP.resample(trNZ[0].stats.sampling_rate, no_filter=False)
#             else:
#                 trNZ.resample(trNP[0].stats.sampling_rate, no_filter=False)
#                 if trN1:
#                     trN1.resample(trNP[0].stats.sampling_rate, no_filter=False)
#                 if trN2:
#                     trN2.resample(trNP[0].stats.sampling_rate, no_filter=False)
#     return trN1, trN2, trNZ, trNP