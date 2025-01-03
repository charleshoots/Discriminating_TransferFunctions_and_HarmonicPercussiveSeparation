from imports import *
from modules import *
from quick_class import *
# import warnings
# import fnmatch
# from obspy.core.inventory.inventory import read_inventory # from IPython.display import clear_output
# import os
from obspy.core.util.attribdict import AttribDict as attr
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cat = catalog.copy()

# EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
tf,method = '','hps'

# EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# tf,method = 'ZP-21','atacr'

nets = cat.Network.unique()
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences')
ovr = True
coh_report = attr()

# corrected_comp='HZ';raw_comp = 'HDH' #ZP
# corrected_comp='HZ';raw_comp = 'HZ' #ZZ
# corrected_comp='H1';raw_comp = 'H1' #H1H1
# corrected_comp='H2';raw_comp = 'H2' #H2H2
corrected_comp='HDH';raw_comp = 'HDH' #H2H2

for ni,n in enumerate(nets):
    N = 'n'+n
    s_corrected_comp = corrected_comp.replace('HZ','Z').replace('H1','1').replace('H2','2').replace('HDH','P')
    s_raw_comp=raw_comp.replace('HZ','Z').replace('H1','1').replace('H2','2').replace('HDH','P')
    s_comps=''.join([s_corrected_comp,s_raw_comp])
    file = str(savefolder / method / f'{n}_{method.lower()}.{s_comps}_coh.report.pkl')
    if not ovr:
        if Path(file).exists():
            print(n + ' Network report complete')
            continue
    coh_report[N] = attr()
    icat = cat[cat.Network==n].copy()
    for si,sta in enumerate(icat.iloc):
        S = sta.Station
        coh_report[N][S] = attr()
        stanm = '.'.join([sta.Network,sta.Station])
        print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')
        corr_fold = CorrFolder/stanm/'CORRECTED'
        raw_fold = EvFolder/stanm/'rmresp'
        events,cpa_list = pull_cohphadm(stanm,
        EvFolder=raw_fold,CorrFolder=corr_fold,tf=tf,
        corrected_comp=corrected_comp,raw_comp=raw_comp)
        os.system('clear')
        print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')

        if len(cpa_list)==0:raise Exception('No data')
        if si==0:f,coh = cpa_list[0].COH();coh_report.f = f

        coh_report[N][S].coh = np.array([c.COH()[-1] for c in cpa_list])
        # coh_report[N][S].cophadm = cpa_list
        coh_report[N][S].events = events
        stafolder = savefolder / method / N[1:] / 'Stations' / s_comps ;stafolder.mkdir(parents=True,exist_ok=True)
        stafile = f'{N[1:]}.{S}.{s_comps}.pkl'
        write_pickle(stafolder/stafile,cpa_list)
    netfolder = stafolder.parent.parent
    write_pickle(netfolder / Path(file).name,coh_report[N].copy())
methodfolder = netfolder.parent
write_pickle(str(methodfolder /(f'complete_{method.lower()}.{s_comps}_coh.report.pkl')),coh_report)