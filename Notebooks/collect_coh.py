from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory # from IPython.display import clear_output
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
cat = catalog.copy()
# EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data/Data')
CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data/Output')
CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data/_depreciated_Output')
# tf = 'ZP-21'
tf = ''
nets = cat.Network.unique()
method = 'hps'
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/Analysis/NetworkCoherences')
g = True
# nets = ['X9']
for ni,n in enumerate(nets):
    coh_report = dict()
    icat = cat[cat.Network==n].copy()
    for si,sta in enumerate(icat.iloc):
        # if len(list((EvFolder / sta.StaName / 'CORRECTED').glob('*' + '.sta.' + tf + '.HZ.SAC')))>=10:
        #     g = True
        stanm = '.'.join([sta.Network,sta.Station])
        print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')
        events,cpa_list = pull_cohphadm(stanm,EvFolder=EvFolder,CorrFolder=CorrFolder,tf=tf,g=g)
        # clear_output(wait=False)
        if len(cpa_list)==0:
            events,cpa_list = pull_cohphadm(stanm,EvFolder,tf=tf,g=g)
            raise Exception('No data')
            print('No Data')
            continue
        print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')
        # lkjj = jkl
        f,coh = cpa_list[0].COH()
        coh_report['f'] = f
        coh_report[sta.StaName] = np.array([c.COH()[-1] for c in cpa_list])
    file = str(savefolder / (n + '_' + method.lower() + '_coh.report.pkl'))
    write_pickle(file,coh_report)
    nsave = 1
    del coh_report