#!/usr/bin/env python

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
# Import modules and functions
import numpy as np
import os.path
import pickle
import stdb
from obspy.clients.fdsn import Client
from obspy import Stream, UTCDateTime
from obstools.atacr import utils
atacr_status = utils.atacr_status
from pathlib import Path

from argparse import ArgumentParser
from os.path import exists as exist
from numpy import nan


def get_daylong_arguments(argv=None):
    """
    Get Options from :class:`~optparse.OptionParser` objects.

    Calling options for the script `obs_download_data.py` that accompany this
    package.

    """

    parser = ArgumentParser(
        usage="%(prog)s [options] <indb>",
        description="Script used " +
        "to download and pre-process up to four-component " +
        "(H1, H2, Z and P), day-long seismograms to use in " +
        "noise corrections of vertical component of OBS data. " +
        "Data are requested from the internet using the client " +
        "services framework for a given date range. The stations " +
        "are processed one by one and the data are stored to disk.")
    parser.add_argument(
        "indb",
        help="Station Database to process from.",
        type=str)

    # General Settings
    parser.add_argument(
        "--keys",
        action="store",
        type=str,
        dest="stkeys",
        default="",
        help="Specify a comma-separated list of station keys " +
        "for which to perform the analysis. These must be " +
        "contained within the station database. Partial keys " +
        "will be used to match against those in the dictionary. " +
        "For instance, providing IU will match with all stations " +
        "in the IU network. " +
        "[Default processes all stations in the database]")
    parser.add_argument(
        "-C", "--channels",
        action="store",
        type=str,
        dest="channels",
        default="Z,P,1,2",
        help="Specify a comma-separated list of channels for " +
        "which to perform the transfer function analysis. " +
        "Possible options are '12' (for horizontal channels 1 and 2) " +
        "and/or 'P' (for pressure channel). Specifying '12' allows " +
        "for tilt correction. Specifying 'P' allows for compliance " +
        "correction. [Default '12,P' looks for both horizontal and " +
        "pressure and allows for both tilt AND compliance corrections]")
    parser.add_argument(
        "-O", "--overwrite",
        action="store_true",
        dest="ovr",
        default=False,
        help="Force the overwriting of pre-existing data. " +
        "[Default False]")
    parser.add_argument(
        "-e", "--evn",
        action="store_true",
        dest="evn",
        default=False,
        help="Forces the download manager to write downloaded files to event specific filenames." +
        "[Default False]")
    # Server Settings
    ServerGroup = parser.add_argument_group(
        title="Server Settings",
        description="Settings associated with which "
        "datacenter to log into.")
    ServerGroup.add_argument(
        "-S", "--Server",
        action="store",
        type=str,
        dest="Server",
        default="IRIS",
        help="Specify the server to connect to. Options include: " +
        "BGR, ETH, GEONET, GFZ, INGV, IPGP, IRIS, KOERI, LMU, NCEDC, " +
        "NEIP, NERIES, ODC, ORFEUS, RESIF, SCEDC, USGS, USP. " +
        "[Default IRIS]")
    ServerGroup.add_argument(
        "-U", "--User-Auth",
        action="store",
        type=str,
        dest="UserAuth",
        default="",
        help="Enter your IRIS Authentification Username and Password " +
        "(--User-Auth='username:authpassword') to access and download " +
        "restricted data. [Default no user and password]")
    """
    # Database Settings
    DataGroup = parser.add_argument_group(
        title="Local Data Settings",
        description="Settings associated with defining " +
        "and using a local data base of pre-downloaded day-long SAC files.")
    DataGroup.add_argument(
        "--local-data",
        action="store",
        type=str,
        dest="localdata",
        default=None,
        help="Specify a comma separated list of paths containing " +
        "day-long sac files of data already downloaded. " +
        "If data exists for a seismogram is already present on disk, " +
        "it is selected preferentially over downloading " +
        "the data using the Client interface")
    DataGroup.add_argument(
        "--no-data-zero",
        action="store_true",
        dest="ndval",
        default=False,
        help="Specify to force missing data to be set as zero, " +
        "rather than default behaviour which sets to nan.")
"""

    # Constants Settings
    FreqGroup = parser.add_argument_group(
        title='Frequency Settings',
        description="Miscellaneous frequency settings")
    FreqGroup.add_argument(
        "--sampling-rate",
        action="store",
        type=float,
        dest="new_sampling_rate",
        default=10.,
        help="Specify new sampling rate (float, in Hz). [Default 5.]")
    FreqGroup.add_argument(
        "--units",
        action="store",
        type=str,
        dest="units",
        default="DISP",
        help="Choose the output seismogram units. Options are: " +
        "'DISP', 'VEL', 'ACC'. [Default 'DISP']")
    FreqGroup.add_argument(
        "--pre-filt",
        action="store",
        type=str,
        dest="pre_filt",
        default=None,
        help="Specify four comma-separated corner frequencies " +
        "(float, in Hz) for deconvolution pre-filter. " +
        "[Default 0.001,0.005,45.,50.]")

    # Event Selection Criteria
    DaysGroup = parser.add_argument_group(
        title="Time Search Settings",
        description="Time settings associated with searching " +
        "for day-long seismograms")
    DaysGroup.add_argument(
        "--start",
        action="store",
        type=str,
        dest="startT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start day for the data search. This will override any " +
        "station start times. " +
        "[Default start date for each station in database]")
    DaysGroup.add_argument(
        "--end",
        action="store",
        type=str,
        dest="endT",
        default="",
        help="Specify a UTCDateTime compatible string representing " +
        "the start time for the event search. This will override any " +
        "station end times [Default end date for each station in database]")

    args = parser.parse_args(argv)

    # Check inputs
    if not exist(args.indb):
        parser.error("Input file " + args.indb + " does not exist")

    # create station key list
    if len(args.stkeys) > 0:
        args.stkeys = args.stkeys.split(',')
    # create channel list
    if len(args.channels) > 0:
        args.channels = args.channels.split(',')
    else:
        args.channels = ["1","2", "P"]
    # for cha in args.channels:
    #     if cha not in ["12", "P","Z"]:
    #         parser.error("Error: Channel not recognized " + str(cha))

    # construct start time
    if len(args.startT) > 0:
        try:
            args.startT = UTCDateTime(args.startT)
        except Exception:
            parser.error(
                "Error: Cannot construct UTCDateTime from start time: " +
                args.startT)
    else:
        args.startT = None

    # construct end time
    if len(args.endT) > 0:
        try:
            args.endT = UTCDateTime(args.endT)
        except Exception:
            parser.error(
                "Error: Cannot construct UTCDateTime from end time: " +
                args.endT)
    else:
        args.endT = None

    # Parse User Authentification
    if not len(args.UserAuth) == 0:
        tt = args.UserAuth.split(':')
        if not len(tt) == 2:
            parser.error(
                "Error: Incorrect Username and Password Strings for " +
                "User Authentification")
        else:
            args.UserAuth = tt
    else:
        args.UserAuth = []

    # # Parse Local Data directories
    # if args.localdata is not None:
    #     args.localdata = args.localdata.split(',')
    # else:
    #     args.localdata = []

    # # Check NoData Value
    # if args.ndval:
    #     args.ndval = 0.0
    # else:
    #     args.ndval = nan

    if args.units not in ['DISP', 'VEL', 'ACC']:
        raise(Exception(
            "Error: invalid --units argument. Choose among " +
            "'DISP', 'VEL', or 'ACC'"))
    if args.pre_filt is None:
        args.pre_filt = [0.001, 0.005, 45., 50.]
        # args.pre_filt = None
    else:
        args.pre_filt = [float(val) for val in args.pre_filt.split(',')]
        args.pre_filt = sorted(args.pre_filt)
        if (len(args.pre_filt)) != 4:
            raise(Exception(
                "Error: --pre-filt should contain 4 comma-separated floats"))

    return args


def main(args=None):
    status=lambda stage=None,note=None:atacr_status(step,stanm,stage,note)
    step='DownloadNoise'
    # stage=None;note=None;status()
    if args is None:
        # Run Input Parser
        args = get_daylong_arguments()
    # Load Database
    # stdb>0.1.3
    try:
        db, stkeys = stdb.io.load_db(fname=args.indb, keys=args.stkeys)
    # stdb=0.1.3
    except Exception:
        db = stdb.io.load_db(fname=args.indb)
        # Construct station key loop
        allkeys = db.keys()
        sorted(allkeys)
        # Extract key subset
        if len(args.stkeys) > 0:
            stkeys = []
            for skey in args.stkeys:
                stkeys.extend([s for s in allkeys if skey in s])
        else:
            stkeys = db.keys()
            sorted(stkeys)
    # Loop over station keys
    for stkey in list(stkeys):
        # Extract station information from dictionary
        sta = db[stkey]
        stanm = '['+'.'.join([sta.network, sta.station])+']'
        status()
        # Define path to see if it exists
        datapath = Path('Data')/'raw'/stkey
        if not datapath.is_dir():
            stage=None;note=None;status(note=f'Path to {str(datapath)} doesn`t exist - creating it')
            datapath.mkdir(parents=True)
        # Establish client
        if len(args.UserAuth) == 0:
            client = Client(args.Server)
        else:
            client = Client(
                args.Server, user=args.UserAuth[0], password=args.UserAuth[1])
        # Get catalogue search start time
        if args.startT is None:
            tstart = sta.startdate
        else:
            tstart = args.startT
        # Get catalogue search end time
        if args.endT is None:
            tend = sta.enddate
        else:
            tend = args.endT
        if tstart > sta.enddate or tend < sta.startdate:
            continue
        # Temporary print locations
        tlocs = sta.location
        if len(tlocs) == 0:
            tlocs = ['']
        for il in range(0, len(tlocs)):
            if len(tlocs[il]) == 0:
                tlocs[il] = "--"
        sta.location = tlocs
        # Update Display
        status(note=f'Channels:{sta.channel},S:{tstart.strftime("%Y-%m-%d")},E:{tend.strftime("%Y-%m-%d")}')
        # Split into 24-hour long segments
        dt = 3600.*24.
        t1 = tstart
        t2 = tstart + dt
        while t2 <= tend:
            # Time stamp
            tstamp = str(t1.year).zfill(4)+'.'+str(t1.julday).zfill(3)+'.'
            status(note=f'Day {str(t1.year)}.{str(t1.julday)}')
            # Define file names (to check if files already exist)
            chans = [];[chans.extend(k) for k in [a for a in args.channels]];args.channels=chans
            files={};_=[files.update({c:datapath / (tstamp+'.'+sta.channel+c.replace('P','DH')+'.SAC')}) for c in args.channels]
            if args.evn:
                files={};_=[files.update({c:datapath / ((t2-2*3600).strftime('%Y.%j.%H.%M')+'.'+sta.channel+c.replace('P','DH')+'.SAC')}) for c in args.channels]
            else:
                files={};_=[files.update({c:datapath / (tstamp+'.'+sta.channel+c.replace('P','DH')+'.SAC')}) for c in args.channels]
            all_exist = np.all([files[k].exists() for k in files.keys()])
            chan_fstr = lambda chan,comp: f'E{chan.upper()}{comp.upper()},H{chan.upper()}{comp.upper()},B{chan.upper()}{comp.upper()}'
            channels=','.join([chan_fstr(sta.channel,c.replace('P','H')) for c in [c for c in ''.join(args.channels)]])
            channels=channels.replace('P','H').replace('HHH','HDH').replace('BHH','BDH').replace('EHH','EDH')
            pressure_channels = ','.join([c for c in channels.split(',') if ('HDH' in c) or ('BDH' in c) or ('EDH' in c)])
            seismic_channels = ','.join([c for c in channels.split(',') if ('HDH' not in c) and ('BDH' not in c) and ('EDH' not in c)])
            status(note=f'Pressure:{pressure_channels}  Seismic:{seismic_channels}')
            stp=Stream();sth=Stream()
            if all_exist & (not args.ovr):
                status(note=f'{tstamp} File exists, continuing')
                t1 += dt;t2 += dt
                continue
            if len(seismic_channels):
                try:
                    stage='[SEISMIC]'
                    status(stage,note=f'{tstamp} Downloading')
                    sth = client.get_waveforms(network=sta.network, station=sta.station,location=sta.location[0],
                    channel=seismic_channels,
                    starttime=t1, endtime=t2,
                    minimumlength=(24*60*60)-1,
                    attach_response=True) #minimumlength=24*3600-1 <---Add to client.get_waveforms
                    status(stage='[SEISMIC]',note="*      ...done")
                except Exception:
                    status(stage,note=" Error: Unable to download ?H? components - : continuing")
                    t1 += dt
                    t2 += dt
                    continue
            if len(pressure_channels):
                try:
                    stage='[PRESSURE]'
                    status(stage,note=f'{tstamp} Downloading')
                    stp = client.get_waveforms(network=sta.network, station=sta.station,location=sta.location[0],
                    channel=pressure_channels,
                    starttime=t1, endtime=t2,
                    minimumlength=(24*60*60)-1,
                    attach_response=True) #minimumlength=24*3600-1 <---Add to client.get_waveforms
                    status(stage='[PRESSURE]',note="*      ...done")
                except Exception:
                    status(stage,note=" Error: Unable to download ?H? components - : continuing")
                    t1 += dt
                    t2 += dt
                    continue
                if len(stp) > 1:
                    status(stage,note=f'WARNING: There are more than one ?DH trace, Keeping the highest sampling rate,Renaming channel to {sta.channel[0]}DH')
                    if stp[0].stats.sampling_rate >stp[1].stats.sampling_rate:stp = Stream(traces=stp[0])
                    else:stp = Stream(traces=stp[1])
            if len(sth)>0:
                sth.detrend('demean');sth.detrend('linear')
                if not sth[0].stats.sampling_rate==args.new_sampling_rate:
                    sth.filter('lowpass', freq=0.5*args.new_sampling_rate,corners=2, zerophase=True)
                    sth.resample(args.new_sampling_rate)
            if len(stp)>0:
                stp.detrend('demean');stp.detrend('linear')
                if not stp[0].stats.sampling_rate==args.new_sampling_rate:
                    stp.filter('lowpass', freq=0.5*args.new_sampling_rate,corners=2, zerophase=True)
                    stp.resample(args.new_sampling_rate)
            st = sth.merge() + stp.merge()
            all_channels_received = len(st)==len(''.join(args.channels))
            # Check streams
            is_ok, st = utils.QC_streams(t1, t2, st)
            if not is_ok:status(note='QC_Streams says "not ok"...continuing');t1 += dt;t2 += dt;continue
            if not all_channels_received:t1 += dt;t2 += dt;status("Missing components. Skipping day.");continue
            # Extract traces - P
            st_sac = Stream()
            if len(st.select(component='H')):
                fileP=files['P']
                trP = st.select(component='H')[0]
                trP = utils.update_stats(trP, sta.latitude, sta.longitude, sta.elevation,sta.channel[0]+'DH')
                status(note='Saving: ' + str(fileP))
                trP.write(str(fileP), format='SAC')
                st_sac+=trP
            if len(st.select(component='Z')):
                fileZ=files['Z']
                # Extract traces - Z
                trZ = st.select(component='Z')[0]
                trZ = utils.update_stats(trZ, sta.latitude, sta.longitude, sta.elevation,sta.channel+'Z')
                status(note='Saving: ' + str(fileZ))
                trZ.write(str(fileZ), format='SAC')
                
                st_sac+=trZ
            # Extract traces - H
            if len(st.select(component='1')):
                file1=files['1']
                tr1 = st.select(component='1')[0]
                tr1 = utils.update_stats(tr1, sta.latitude, sta.longitude, sta.elevation,sta.channel+'1')
                status(note='Saving: ' + str(file1))
                tr1.write(str(file1), format='SAC')
                st_sac+=tr1
            if len(st.select(component='2')):    
                file2=files['2']
                tr2 = st.select(component='2')[0]
                tr2 = utils.update_stats(tr2, sta.latitude, sta.longitude, sta.elevation,sta.channel+'2')
                status(note='Saving: ' + str(file2))
                tr2.write(str(file2), format='SAC')
                st_sac+=tr2
            inventory_file = datapath / (datapath.name + '_inventory.xml')
            # Writing response files to station inventory
            if not inventory_file.exists():utils.save_inventory(str(inventory_file),st_sac)
            t1 += dt
            t2 += dt


if __name__ == "__main__":
    # Run main program
    main()
