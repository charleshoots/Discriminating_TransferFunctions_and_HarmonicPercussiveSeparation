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

from obspy.clients.fdsn import Client
from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
client = Client('IRIS')

def update_stats(trace, stla=None, stlo=None, stel=None, cha=None, evla=None, evlo=None):
    """
    Function to include SAC metadata to :class:`~obspy.core.Trace` objects
    Parameters
    ----------
    trace : :class:`~obspy.core.Trace` object
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
    trace : :class:`~obspy.core.Trace` object
        Updated trace object
    """
    trace.stats.sac = AttribDict()
    if stla is None:
        inventory = client.get_stations(network=trace.stats.network, station=trace.stats.station,starttime=trace.stats.starttime,endtime=trace.stats.endtime)
        stla,stlo,stel = [inventory[0][0].__dict__[k] for k in ['_latitude','_longitude','_elevation']]

    trace.stats.sac.stla = stla
    trace.stats.sac.stlo = stlo
    trace.stats.sac.stel = stel

    if cha is not None:
        trace.stats.sac.kcmpnm = cha
        trace.stats.channel = cha
    if evla is not None and evlo is not None:
        trace.stats.sac.evla = evla
        trace.stats.sac.evlo = evlo
    return trace