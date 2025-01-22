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
        [print("*   "+tr.stats.channel+" " +str(tr.stats.starttime)+" " +str(tr.stats.endtime)) for tr in st]
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
        print("*    "+str(st[0].stats.npts)+" ~= "+str(int((end - start)*sr)))
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

def get_event(eventpath, tstart=None, tend=None,evlist=None,skipfiles=None):
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
    if np.isin('rmresp',[e.name for e in list(eventpath.parents)]):rmresp=False
    else:rmresp=True
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
                    if (skipfiles is not None) & np.any(np.isin(file,skipfiles)):continue
                    tr1=load_data(Path(str(file).replace('Z.SAC','*1.SAC')))[0]
                    tr2=load_data(Path(str(file).replace('Z.SAC','*2.SAC')))[0]
                    trZ=load_data(Path(str(file).replace('Z.SAC','*Z.SAC')))[0]
                    trP=load_data(Path(str(file).replace('Z.SAC','*H.SAC')))[0]
                    tr1,tr2,trZ,trP = check_fs([tr1,tr2,trZ,trP])
                    yield tr1,tr2,trZ,trP
    else:
        for ev in evlist:
            files=[Path(str(eventpath / ev)+f'*{c}.SAC') for c in ['1','2','Z','H']]
            if (skipfiles is not None) & np.any(np.isin(files,skipfiles)):continue
            tr1,tr2,trZ,trP=[load_data(f)[0] for f in files]
            tr1,tr2,trZ,trP = check_fs([tr1,tr2,trZ,trP])
            yield tr1,tr2,trZ,trP

def check_fs(st):
    rates=[tr.stats.sampling_rate for tr in st]
    inequal_rates=len(np.unique([tr.stats.sampling_rate for tr in st]))>1
    if inequal_rates:
        new_fs=min(rates) #choose the lowest sample rate from all given traces
        [tr.resample(new_fs, no_filter=False) for tr in st]
    return st


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

def load_data(files):
    if not isinstance(files,list):files=[files]
    """
    Loads four files defined by the tuple, files.
    -CHoots, 2023
    """
    # Implicitly define (was already done if a parent folder is named 'rmresp') whether remove response is needed.
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

def get_data_generator(datapath, tstart, tend,skipfiles=None):
    """
    The same as get_data but now factored as a generator object such that the data is never loaded until its index is called.
    Generators can only be iterated, and activated by use in a for loop
    Significantly more memory efficient.
    -CHoots, 2023
    """
    file_list = search_files(datapath, tstart, tend)
    for files in file_list:
        if (skipfiles is not None) & np.any(np.isin(files,skipfiles)):continue
        yield load_data(files)

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
