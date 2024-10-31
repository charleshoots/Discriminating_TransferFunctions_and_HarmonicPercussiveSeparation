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
# ------------------------------------------------------

# ev_intputfolder = Path(dirs['Py_CorrectedTraces'])
plotfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/DLOPY_Data/Figures')
ev_intputfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/DLOPY_Data/EVENTS')
stakey = '7D.J44C'
dateformat = '%Y.%j.%H.%M'
tf = 'ZP-21'

# cat = catalog[catalog.StaName==stakey]
# cat = catalog[catalog.Station.isin(['J44C','G03A','M08A'])]
cat = catalog[catalog.Station.isin(['M07A'])]

dl_err = dict()
dl_err['7D.M07A'] = (93.16,5.38)
dl_err['7D.M08A'] = (4.52,14.88)
dl_err['7D.J44C'] = (295.31,3.78)
dl_err['7D.G03A'] = (155.77,11.94)


goodevents = dict()
goodevents['7D.M08A'] = ['2012.181.21.07','2011.301.18.54',
 '2012.024.00.52','2012.033.13.34',
 '2012.074.09.08','2012.074.21.13',
 '2012.080.17.56','2012.080.18.02']


goodevents['7D.M07A'] = ['2012.112.01.16',
'2012.033.13.34',
'2012.074.09.08',
'2012.112.01.16',
'2012.033.13.34',
'2012.074.09.08',
'2012.080.17.56',
'2012.080.18.02',
'2011.294.17.57']


# for sta in cat.iloc:
#     staframe = sta.to_frame().T
#     events = sta.EventMeta
#     evtimes = [e.origins[0].time.strftime(dateformat) for e in events]
#     for evi,(ev,evstr) in enumerate(zip(events,evtimes)):
#         datahold = []
#         for isource in ['raw','atacr']:
#             if isource=='raw':
#                 subfolder = ''
#                 rmresp=True
#             else:
#                 rmresp = False
#                 subfolder = 'CORRECTED'
#             stakey = sta.StaName
#             # ------------------------------------------------------------------------------
#             staquery_output = build_staquery(staframe,staquery_output=None)
#             # Load catalog
#             args = AttribDict()
#             args.showplot=True;args.saveplot=True
#             args.cc = 0.6
#             dl_list = []
#             R1phi = []; R1cc = []; R2phi = []; R2cc = []
#             # Load data
#             st,inv = Stream(),[]
#             if isource=='raw':files_z = [evstr+'.*'+c+'.SAC' for c in ['Z']]
#             else:files_z = [stakey +'.'+evstr+'.sta.'+tf + '.*'+c+'.SAC' for c in ['Z']]
#             if len(list((ev_intputfolder / stakey / subfolder).glob(files_z[0])))<1:
#                 print(stakey + ' | ' + evstr + ' | Event correction not found. Skipping.')
#                 continue
#             # ----
#             data = [load_sac(ev_intputfolder / stakey / subfolder / f,rmresp=rmresp) for f in files_z]
#             inv = data[0][1];data = [data[0][0][0]]
#             files = [evstr+'.*'+c+'.SAC' for c in ['1','2']]
#             data2 = [load_sac(ev_intputfolder / stakey / f,rmresp=rmresp)[0][0] for f in files]
#             data.extend(data2)
#             for di,d in enumerate(data):data[di].stats.loc = isource #''.join([data[di].id,':',isource])
#             datahold.extend(data)
#             st = Stream(datahold)
#             # ----
#             if any([d is None for d in data]):print(str(evi+1),'/',str(len(events)),' None returned');continue
#             if len(st)<3:print(str(evi+1),'/',str(len(events)),' Missing data');continue
#         if len(list((ev_intputfolder / stakey / subfolder).glob(files_z[0])))<1:continue
#         z = st.select(component='*Z').copy()
#         for di,d in enumerate(z):z[di].id = z[di].id+':'+z[di].stats.loc
#         z.taper(0.05)
#         z.filter('bandpass',freqmin=1/200,freqmax=1,zerophase=True)
#         plot = z.plot(show=False,draw=False)
#         (plotfolder/'traces'/stakey).mkdir(parents=True,exist_ok=True)
#         save_tight(plotfolder /'traces'/stakey/(evstr + '.png'),plot)
#         k = 1

# LPF = [0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.005]
# HPF = [0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015]
# for L,H in zip(LPF,HPF):
#     m = z.copy();m.filter('bandpass',freqmin=L,freqmax=H,zerophase=True);m.plot()


# _______________________________________________________________________|
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx||
# ----------# ----------# ----------# ----------# ----------# ----------||
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx||
# ----------# ----------# ----------# ----------# ----------# ----------||
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx||
# _______________________________________________________________________|

ev_source = 'raw'
# ev_source = 'atacr'
# -------------------
# ev_sourcei = ['atacr','atacr','atacr','raw','raw','raw','raw','raw','raw']
for sta in cat.iloc:
    stakey = sta.StaName
    # ------------------------------------------------------------------------------
    staquery_output = build_staquery(cat,staquery_output=None)
    # Load catalog
    events = cat.iloc[0].EventMeta
    evtimes = [e.origins[0].time.strftime(dateformat) for e in events]
    args = AttribDict()
    args.showplot=True;args.saveplot=True
    args.cc = 0.75
    dl_list = []
    R1phi = []; R1cc = []; R2phi = []; R2cc = []
    ni = -1
    evtimes = evtimes[:10]
    events = events[:10]
    for evi,(ev,evstr) in enumerate(zip(events,evtimes)):
        # if not np.isin(evstr,goodevents[stakey]):
        #     continue
        # else:
        #     print('--'*30 + 'Event: ' + evstr + '--'*30)
        #     print('=='*70)
        # ni+=1
        # ev_source = ev_sourcei[ni]
        print('--'*30 + 'Event: ' + evstr + '--'*30);print('=='*70)
        if ev_source=='raw':
            subfolder = '';rmresp=True
        else:
            subfolder = 'CORRECTED';rmresp = False
        # Load data
        st,inv = Stream(),[]
        if ev_source=='raw':files = [evstr+'.*'+c+'.SAC' for c in ['Z']]
        else:files = [stakey +'.'+evstr+'.sta.'+tf + '.*'+c+'.SAC' for c in ['Z']]
        if len(list((ev_intputfolder / stakey / subfolder).glob(files[0])))<1:
            print(stakey + ' | ' + evstr + ' | Event correction not found. Skipping.')
            continue
        # ----
        data = [load_sac(ev_intputfolder / stakey / subfolder / f,rmresp=rmresp) for f in files]
        inv = data[0][1];data = [data[0][0][0]]
        files = [evstr+'.*'+c+'.SAC' for c in ['1','2']]
        data2 = [load_sac(ev_intputfolder / stakey / f,rmresp=rmresp)[0][0] for f in files]
        data.extend(data2)
        st = Stream(data)
        # ----
        if any([d is None for d in data]):print(str(evi+1),'/',str(len(events)),' None returned');continue
        if len(st)<3:print(str(evi+1),'/',str(len(events)),' Missing data');continue
        # [(st.append(tr[0]),inv.append(invi)) if len(tr)>0 for tr,invi in [load_sac(ev_intputfolder / stakey / f) for f in files]];inv = inv[0]
        # fig = inv.plot(show=False,color='r',block=True,draw=False)
        # evcat.plot(fig=fig)
        # plt.show()
        dldata = DL(staquery_output[staquery_output.keys()[0]])
        dldata.add_event(ev, returned=True,gacmax=175)
        dldata.meta.stakey = stakey
        if dldata.meta.accept:
            print(str(evi+1),'/',str(len(events)),' Event accepted')
            dldata.data = st
            dldata.calc(showplot=True)
        else:print(str(evi+1),'/',str(len(events)),' Event not accepted');continue
        dl_list.append(dldata)
        plt.close('all')
        k=1;print('Max cc:', str(np.max(np.concatenate((dldata.meta.R1cc, dldata.meta.R2cc), axis=None))))

        # for folder in os.listdir(indir):
        #     # Load meta data
        #     filename = indir / folder / "Meta_data.pkl"
        # if not filename.is_file():continue
        # meta = pickle.load(open(filename, 'rb'))np.max(np.concatenate((dldata.meta.R1cc, dldata.meta.R2cc), axis=None))
        R1phi.append(dldata.meta.R1phi);R2phi.append(dldata.meta.R2phi)
        R1cc.append(dldata.meta.R1cc);R2cc.append(dldata.meta.R2cc)
        del dldata
    R1phi = np.array(R1phi).flatten()
    R1cc = np.array(R1cc).flatten()
    R2phi = np.array(R2phi).flatten()
    R2cc = np.array(R2cc).flatten()
    phi = np.concatenate((R1phi, R2phi), axis=None)
    cc = np.concatenate((R1cc, R2cc), axis=None)
    ind = cc > args.cc
    val, err = orientpy.utils.estimate(phi, ind)
    print("|    D-L mean, error, data included: " + "{0:.2f}, {1:.2f}, {2}".format(val, err, np.sum(ind)))
    print("|    D-L CC level: {0:.1f}".format(args.cc))
    print()
    # if np.sum(np.isnan(np.array([val, err])))>0:
    #     continue
    if args.showplot or args.saveplot:
        plot = orientpy.plotting.plot_dl_results(stakey, R1phi, R1cc, R2phi, R2cc, ind,
        val, err, phi, cc, args.cc)
        plot.get_axes()[0].axhline(dl_err[stakey][0],c='gray',linestyle='-',linewidth=0.2)
        plot.get_axes()[0].text(0.0,dl_err[stakey][0],'Doran & Laske (2017)='+str(np.round(dl_err[stakey][0]*10)/10),horizontalalignment='left',verticalalignment='top',bbox=dict(facecolor='w',alpha=0.2,edgecolor='none'))
        plot.get_axes()[0].text(1.0,np.mean(val),'$\phi$='+str(np.round(np.mean(val)*10)/10),horizontalalignment='right',verticalalignment='bottom',bbox=dict(facecolor='w',alpha=0.2,edgecolor='none'))
    file = '.'.join([stakey,ev_intputfolder.parent.name,ev_source,'png'])
    save_tight(plotfolder / file,plot)

k = 1

# Plot event list
# evplot = evcat.plot()
# plt.show()
# plt.close('all')
# clear_output(wait=False)