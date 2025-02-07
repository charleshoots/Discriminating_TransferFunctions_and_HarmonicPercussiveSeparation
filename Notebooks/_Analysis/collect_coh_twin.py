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
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences')


EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
tf,method = '','hps';coh_stage='post'

# EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# tf,method = 'ZP-21','aNotebooks/collect_coh.pytacr';coh_stage='post'

# EvFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# CorrFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/EVENTS')
# tf,method = '','atacr';coh_stage='pre'


nets = cat.Network.unique()
ovr = True;coh_report = attr()

comps=['HDH','HZ','H1','H2']
for compi,comp in enumerate(comps):
    for ni,n in enumerate(nets):
        N = 'n'+n
        s_comp = comp.replace('HZ','Z').replace('H1','1').replace('H2','2').replace('HDH','P')
        if coh_stage=='post':file=str(savefolder / method / (n + '_' + method.lower()+'.Z' + s_comp+'_coh.report.pkl'))
        else:file=str(savefolder / 'uncorrected' / (n + '_' + 'Uncorrected'+'.Z' + s_comp+'_coh.report.pkl'))
        if not ovr:
            if Path(file).exists():
                print(n + ' Network report complete')
                continue
        coh_report[N] = attr()
        icat = cat[cat.Network==n].copy()
        for si,sta in enumerate(icat.iloc):
            status = f'|Components {comp} {compi+1}/{len(comps)}| Network {n} {ni+1}/{len(nets)}| Station {sta.Station} {si+1}/{len(icat)}|'
            clear_output()
            print(status)
            S = sta.Station
            coh_report[N][S] = attr()
            stanm = '.'.join([sta.Network,sta.Station])
            print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')

            if coh_stage=='post':
                correvpath = CorrFolder / stanm /'CORRECTED'
                rawevpath = EvFolder / stanm
                events,cpa_list = pull_cohphadm(stanm,EvFolder=rawevpath,CorrFolder=correvpath,tf=tf,comp=comp)
            else:
                corr_rmresp=True
                correvpath = EvFolder / stanm
                rawevpath = EvFolder / stanm
                events,cpa_list = pull_cohphadm(stanm,EvFolder=rawevpath,CorrFolder=correvpath,tf='',comp=comp,corr_rmresp=corr_rmresp)

            os.system('clear')
            print(sta.StaName+' [Net:'+str(ni+1)+'/'+str(len(nets))+']'+' [Sta:'+str(si+1)+'/'+str(len(icat))+']')
            if len(cpa_list)==0:raise Exception('No data')
            if si==0:f,coh = cpa_list[0].COH();coh_report.f = f
            coh_report[N][S].coh = np.array([c.COH()[-1] for c in cpa_list])
            # coh_report[N][S].cophadm = cpa_list
            coh_report[N][S].events = events

            if coh_stage=='post':stafolder = savefolder / method / N[1:] / 'Stations' / ('Z'+s_comp) ;stafolder.mkdir(parents=True,exist_ok=True)
            else:stafolder = savefolder / 'uncorrected' / N[1:] / 'Stations' / ('Z'+s_comp) ;stafolder.mkdir(parents=True,exist_ok=True)

            if coh_stage=='post':stafile = f'{N[1:]}.{S}.Z{s_comp}.pkl'
            else:stafile = f'{N[1:]}.{S}.Z{s_comp}.Uncorrected.pkl'
            write_pickle(stafolder/stafile,cpa_list)
        netfolder = stafolder.parent.parent
        write_pickle(netfolder / Path(file).name,coh_report[N].copy())
    methodfolder = netfolder.parent
    if coh_stage=='post':write_pickle(str(methodfolder /('complete_'+method.lower()+'.Z'+s_comp+'_coh.report.pkl')),coh_report)
    else:write_pickle(str(methodfolder /('complete_'+'Uncorrected'+'.Z'+s_comp+'_coh.report.pkl')),coh_report)