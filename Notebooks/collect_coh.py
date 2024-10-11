from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory # from IPython.display import clear_output
import os
from obspy.core.util.attribdict import AttribDict as attr
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cat = catalog.copy()

EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data/Data')
CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data/Output')
tf,method = '','hps'

# EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# tf,method = 'ZP-21','atacr'

nets = cat.Network.unique()
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/Analysis/NetworkCoherences')
ovr = True
coh_report = attr()
for ni,n in enumerate(nets):
    N = 'n'+n
    file = str(savefolder / method / (n + '_' + method.lower() + '_coh.report.pkl'))
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
        events,cpa_list = pull_cohphadm(stanm,EvFolder=EvFolder,CorrFolder=CorrFolder,tf=tf)
        os.system('clear')
        print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')

        if len(cpa_list)==0:raise Exception('No data')
        if si==0:f,coh = cpa_list[0].COH();coh_report.f = f

        coh_report[N][S].zzcoh = np.array([c.COH()[-1] for c in cpa_list])
        coh_report[N][S].events = events
    write_pickle(file,coh_report[N].copy())
write_pickle(str(savefolder / method /('complete_'+method.lower()+'_coh.report.pkl')),coh_report)