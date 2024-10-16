from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
# ========================================================================================================================================================
from IPython.display import clear_output
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
def write_pickle(file,var):
    import pickle
    with open(str(file), 'wb') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Saved to :' + str(file))
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        b = pickle.load(handle)
    return b
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
reportfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/Analysis/NetworkCoherences')
bands = ['1-10','10-30','30-100']

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# s1 Build event-station metadata catalogue
# Using Obspy Station Objects and obspy Event Objects inside a Catalog Object
# ------------------------------------------------------
# # class Meta(object)
# __init__(self, sta, event, gacmin=5., gacmax=30., depmax=1000.)
# via..:
# ------------------------------------------------------
# bngdata = BNG(sta)
# # Add event to object
# accept = bngdata.add_event(
# ev, gacmin=args.mindist, gacmax=args.maxdist,
# depmax=args.maxdep, returned=True)
# 
# self.data = Stream(traces=[trZ, tr1, tr2])
# 
# ------------------------------------------------------


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# dldata = DL(sta)
# # Add event to object
# accept = dldata.add_event(ev, gacmin=args.mindist, gacmax=args.maxdist,depmax=args.maxdep, returned=True)
# has_data = dldata.download_data(client=wf_client, stdata=args.localdata,ndval=args.ndval, new_sr=2., t1=t1, t2=t2,returned=True, verbose=args.verb)
# dldata.calc(showplot=False)
# dldata = DL(sta)
# # Add event to object
# accept = dldata.add_event(ev, gacmin=args.mindist, gacmax=args.maxdist,depmax=args.maxdep, returned=True)
# has_data = dldata.download_data(client=wf_client, stdata=args.localdata,ndval=args.ndval, new_sr=2., t1=t1, t2=t2,returned=True, verbose=args.verb)
# dldata.calc(showplot=False)

# ev_intputfolder = Path(dirs['Py_CorrectedTraces'])
ev_intputfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/DLOPY_Data/EVENTS')
stakey = '7D.J44C'
dateformat = '%Y.%j.%H.%M'

cat = catalog[catalog.StaName==stakey]
staquery_output = build_staquery(cat,staquery_output=None)
# Load catalog
events = [ev[0] if isinstance(ev,obspy.core.event.catalog.Catalog) else ev for ev in cat.iloc[0].Metadata]
evcat = Catalog();evcat.extend(events)
evtimes = [e.origins[0].time.strftime(dateformat) for e in evcat]

dl_list = []
for evi,(ev,evstr) in enumerate(zip(events,evtimes)):
    # Load data
    st,inv = Stream(),[]
    files = [evstr+'.*'+c+'.SAC' for c in ['Z','1','2']]
    data = [load_sac(ev_intputfolder / stakey / f) for f in files]
    if any([d is None for d in data]):print(str(evi+1),'/',str(len(events)),' None returned');continue
    inv = data[0][1]
    _ = [st.extend(s[0]) for s in data]
    if len(st)<3:print(str(evi+1),'/',str(len(events)),' Missing data');continue
    # [(st.append(tr[0]),inv.append(invi)) if len(tr)>0 for tr,invi in [load_sac(ev_intputfolder / stakey / f) for f in files]];inv = inv[0]
    # fig = inv.plot(show=False,color='r',block=True,draw=False)
    # evcat.plot(fig=fig)
    # plt.show()
    dldata = DL(staquery_output[staquery_output.keys()[0]])
    dldata.add_event(ev, returned=True,gacmax=175)
    dl_list.append(dldata)
    if dldata.meta.accept:
        print(str(evi+1),'/',str(len(events)),' Event accepted')
        dldata.data = st
        dldata.calc(showplot=True)
    else:print(str(evi+1),'/',str(len(events)),' Event not accepted')
    plt.close('all')
    k=1
k = 1

# Plot event list
# evplot = evcat.plot()
# plt.show()
# plt.close('all')
# clear_output(wait=False)