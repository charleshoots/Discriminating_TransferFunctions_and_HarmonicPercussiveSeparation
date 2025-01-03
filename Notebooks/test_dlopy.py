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

# ev_input_folder = Path(dirs['Py_CorrectedTraces'])
# plotfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/DLOPY_Data/Figures')
plotfolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/DLOPY_Data/Figures')
ev_input_folder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/DLOPY_Data/EVENTS')
stakey = '7D.J44C'
dateformat = '%Y.%j.%H.%M'
tf = 'ZP-21'

# cat = catalog[catalog.StaName==stakey]
# cat = catalog[catalog.Station.isin(['J44C','G03A','M08A'])]
# cat = catalog[catalog.Station.isin(['M07A'])]
cat=catalog[catalog.Network.isin(['7D'])].copy()


doran_laske_calcs=Path('/Users/charlesh/Desktop/711_Proposal/figures/2016165_esupp_Table_S1.csv')
janiszewski23_meta=Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/Catalogs/Janiszewski_etal_2023_StationList.xlsx')
dl_data = pd.read_csv(doran_laske_calcs)
haj_dat= pd.read_excel(janiszewski23_meta)
haj_dat=haj_dat[haj_dat.Station.isin(dl_data.Station)]
dl_data=dl_data[dl_data.Station.isin(haj_dat.Station)]
dl_data=dl_data.sort_values(by='Station')
haj_dat=haj_dat.sort_values(by='Station')
haj_dat=haj_dat.rename(columns={
'Latitude (deg)':'Latitude','Longitude (deg)':'Longitude',
'Pressure Gauge':'Pressure_Gauge','Water Depth (m)':'StaDepth','Instrument Design':'Instrument_Design'})[['Station', 'Network',
'Latitude','Longitude','Experiment',
'Environment','Pressure_Gauge','StaDepth','Start','End','Instrument_Design','Seismometer']]
dl_data=dl_data[['Station','Orientation','Uncertainty','Total Phases','Unique Events','Include High Fqs','CC used']]
dl_data.reset_index(drop=True)
haj_dat.reset_index(drop=True)
dl_data=pd.merge(dl_data,haj_dat)
dl_data.set_index('Station', inplace=True,drop=False)


dl_err = dict()
for stakey in cat.StaName:dl_err[stakey]=dl_data.loc[stakey.split('.')[1]].Orientation,dl_data.loc[stakey.split('.')[1]].Uncertainty


# _______________________________________________________________________|
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx||
# ----------# ----------# ----------# ----------# ----------# ----------||
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx||
# ----------# ----------# ----------# ----------# ----------# ----------||
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx||
# _______________________________________________________________________|
dists=[]
ev_source = 'raw'
ev_source = 'atacr'
# -------------------
# ev_sourcei = ['atacr','atacr','atacr','raw','raw','raw','raw','raw','raw']
ev_input_folder = dirs.Events/'rmresp'
# ev_input_folder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/DLOPY_Data/EVENTS')

for ev_source in ['atacr','raw']:
    DL_OUTPUTS = AttribDict()
    for sta in cat.iloc:
        stakey = sta.StaName
        # if not stakey=='7D.M07A':
        #     continue
        # else:
        #     goodevents = np.array(['2012.112.01.16',
        #     '2012.033.13.34',
        #     '2012.074.09.08',
        #     '2012.112.01.16',
        #     '2012.033.13.34',
        #     '2012.074.09.08',
        #     '2012.080.17.56',
        #     '2012.080.18.02',
        #     '2011.294.17.57'])
        #     events=Catalog([e for e in sta.Events if np.any(e.Name==goodevents)])
        # ------------------------------------------------------------------------------
        staquery_output = build_staquery(cat,staquery_output=None)
        # Load catalog
        events = sta.Events
        evtimes = [e.Name for e in events]
        args = AttribDict()
        args.showplot=True;args.saveplot=True
        args.cc = 0.75
        dl_list = []
        R1phi = []; R1cc = []; R2phi = []; R2cc = [];LPF=[];HPF=[]
        for evi,(ev,evstr) in enumerate(zip(events,evtimes)):
            print('--'*30 + 'Event: ' + evstr + '--'*30);print('=='*70)
            # Load data
            st,inv = Stream(),[]
            if ev_source=='raw':
                files = [evstr+'.*'+c+'.SAC' for c in ['HZ','H1','H2']]
                files=[g[0] for g in [list((ev_input_folder /stakey ).glob(f)) for f in files] if len(g)==1]
            else:
                files = [stakey +'.'+evstr+'.*'+c+'.SAC' for c in ['sta.ZP-21.HZ']]
                files.extend([evstr+'.*'+c+'.SAC' for c in ['H1','H2']])
                files=[g[0] for g in [list((fo /stakey ).glob(f)) for fo,f in zip([dirs.Events/'corrected',dirs.Events/'rmresp',dirs.Events/'rmresp'],files)] if len(g)==1]
                # files=[g[0] for g in [list((fo).glob(f)) for fo,f in zip([ev_input_folder/stakey/'CORRECTED',ev_input_folder/stakey,ev_input_folder/stakey],files)] if len(g)==1]
            if len(files)<3:
                print(stakey + ' | ' + evstr + ' | Event data not found. Skipping.')
                continue
            # ----
            if ev_source=='raw':
                data = Stream([load_sac(ev_input_folder / stakey / f,rmresp=False)[0][0] for f in files])
                # data = Stream([load_sac(ev_input_folder / stakey / f,rmresp=True)[0][0] for f in files])
            else:
                data = Stream([load_sac(f,rmresp=False)[0][0] for f in files])
                # data = Stream([load_sac(f,rmresp=rmresp)[0][0] for f,rmresp in zip(files,[False,True,True])])
                # data.resample(40)
            # ----
            if any([d is None for d in data]):print(str(evi+1),'/',str(len(events)),' None returned');continue
            if len(data)<3:print(str(evi+1),'/',str(len(events)),' Missing data');continue

            dldata = DL(staquery_output[staquery_output.keys()[0]])
            dldata.add_event(ev, returned=True,gacmax=175)
            dldata.meta.stakey = stakey
            if dldata.meta.accept:

                Rearth = 6371.25
                circE = 2.*np.pi*Rearth
                dist2 = circE - dldata.meta.epi_dist
                dists.append(dist2)

                print(str(evi+1),'/',str(len(events)),' Event accepted')
                dldata.data = data
                dldata.calc(showplot=True)
            else:print(str(evi+1),'/',str(len(events)),' Event not accepted');continue
            dl_list.append(dldata)
            plt.close('all')
            k=1;print('Max cc:', str(np.max(np.concatenate((dldata.meta.R1cc, dldata.meta.R2cc), axis=None))))


            R1phi.append(dldata.meta.R1phi);R2phi.append(dldata.meta.R2phi)
            R1cc.append(dldata.meta.R1cc);R2cc.append(dldata.meta.R2cc)
            LPF.append(np.array(dldata.meta.LPF))
            HPF.append(np.array(dldata.meta.HPF))
            del dldata
        if len(R1phi)==0:
            print('No station data. Skipping.')
        else:
            R1phi = np.array(R1phi).flatten()
            R1cc = np.array(R1cc).flatten()
            R2phi = np.array(R2phi).flatten()
            R2cc = np.array(R2cc).flatten()

            LPF = np.array(LPF).flatten()
            HPF = np.array(HPF).flatten()
            R2_LPF=LPF[~np.isnan(R2cc)]
            R2_HPF=HPF[~np.isnan(R2cc)]
            R1_LPF=LPF[~np.isnan(R1cc)]
            R1_HPF=HPF[~np.isnan(R1cc)]

            R1phi=R1phi[~np.isnan(R1phi)]
            R1cc=R1cc[~np.isnan(R1cc)]
            R2phi=R2phi[~np.isnan(R2phi)]
            R2cc=R2cc[~np.isnan(R2cc)]

            dct=AttribDict()
            dct.R1=AttribDict()
            dct.R2=AttribDict()
            dct.R1.LPF=R1_LPF
            dct.R1.HPF=R1_HPF
            dct.R2.LPF=R2_LPF
            dct.R2.HPF=R2_HPF
            dct.R1.phi=R1phi
            dct.R1.cc=R1cc
            dct.R2.phi=R2phi
            dct.R2.cc=R2cc
            DL_OUTPUTS[stakey]=dct
            fname=Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/DLOPY_Data')
            write_pickle(fname/f'DL_OUTPUTS.{ev_source}.pkl',DL_OUTPUTS)

            phi = np.concatenate((R1phi, R2phi), axis=None)
            cc = np.concatenate((R1cc, R2cc), axis=None)
            # args.cc=0.01
            ind = cc > args.cc
            val, err = orientpy.utils.estimate(phi, ind)
            print("|    D-L mean, error, data included: " + "{0:.2f}, {1:.2f}, {2}".format(val, err, np.sum(ind)))
            print("|    D-L CC level: {0:.1f}".format(args.cc))
            print()

            if args.showplot or args.saveplot:
                plot = orientpy.plotting.plot_dl_results(stakey, R1phi, R1cc, R2phi, R2cc, ind,
                val, err, phi, cc, args.cc)
                plot.set_figwidth(11)
                plot.set_figheight(3)
                if np.any(np.array(list(dl_err.keys()))==stakey):
                    plot.get_axes()[0].axhline(dl_err[stakey][0],c='gray',linestyle='-',linewidth=0.2)
                    plot.get_axes()[0].text(0.0,dl_err[stakey][0],'Doran & Laske (2017)='+str(np.round(dl_err[stakey][0]*10)/10),horizontalalignment='left',verticalalignment='top',bbox=dict(facecolor='w',alpha=0.2,edgecolor='none'))
                plot.get_axes()[0].text(1.0,np.mean(val),'$\phi$='+str(np.round(np.mean(val)*10)/10),horizontalalignment='right',verticalalignment='bottom',bbox=dict(facecolor='w',alpha=0.2,edgecolor='none'))
            file = '.'.join([stakey,ev_input_folder.parent.name,ev_source,'png'])
            save_tight(plotfolder / file,plot)
            k=1

k = 1
