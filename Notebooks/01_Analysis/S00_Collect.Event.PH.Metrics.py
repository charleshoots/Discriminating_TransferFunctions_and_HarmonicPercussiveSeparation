import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
import os,sys;os.chdir('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/Notebooks')
sys.path.append('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/Notebooks')
from imports import *
from modules import *
from quick_class import *
from obspy.core.util.attribdict import AttribDict as attr


cat = catalog.copy()
nets = cat.Network.unique()
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkPhase')


method_ind=0
ovr = True

# ------------------------------------------------------------------------------------------------
# base_pairs=[['HZ','HDH'],['HDH','HZ'],['HZ','HZ'],['H1','H1'],['H2','H2']]#Every auto-phase, only useful for NoiseCut
base_pairs=[['HZ','HZ']]#Defaults: ZZ,ZP,Z1,Z2


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

log_list = []
# 'hps'
for method_ind,method in enumerate(['hps','atacr']):
    ph_sets=base_pairs.copy()
    if method=='hps':tf,method = '','hps';CorrectedFold=dirs.Events_HPS/'corrected';UncorrectedFold=dirs.Events_HPS/'rmresp'
    else:tf,method = 'ZP-21','atacr';CorrectedFold = dirs.Events/'corrected';UncorrectedFold = dirs.Events/'rmresp'
    if method.lower()=='hps':ph_sets.extend([['H1','H1'],['H2','H2']])
    for set_i,(corrected_comp,raw_comp) in enumerate(ph_sets):
        ph_report = attr()
        if method.lower()=='atacr':
            if (corrected_comp=='H2') or (corrected_comp=='H1'):continue
        os.system('clear')
        for ni,n in enumerate(nets):

            N = 'n'+n
            s_corrected_comp = corrected_comp.replace('HZ','Z').replace('H1','1').replace('H2','2').replace('HDH','P')
            s_raw_comp=raw_comp.replace('HZ','Z').replace('H1','1').replace('H2','2').replace('HDH','P')
            s_comps=''.join([s_corrected_comp,s_raw_comp])
            file = str(savefolder / method / f'{n}_{method.lower()}.{s_comps}_ph.report.pkl')

            if not ovr:
                if Path(file).exists():print(n + ' Network report complete');continue

            ph_report[N] = attr()
            icat = cat[cat.Network==n].copy()
            for si,sta in enumerate(icat.iloc):
                state = lambda:f'PHASE|{method.upper()}-{s_comps} ({set_i+1}/{len(ph_sets)}) | {sta.StaName} | N {str(ni+1)}/{str(len(nets))} | S {str(si+1)}/{str(len(icat))}'
                S = sta.Station
                ph_report[N][S] = attr()
                stanm = '.'.join([sta.Network,sta.Station])
                print(state())
                corrected_sta_fold = CorrectedFold/stanm
                uncorrected_sta_fold = UncorrectedFold/stanm
                if method.lower()=='atacr':
                    if corrected_comp=='HZ':tf='ZP-21';corrected_sta_fold=CorrectedFold/stanm
                    else:corrected_sta_fold=uncorrected_sta_fold;tf=''

                events,cpa_list,log_warnings = pull_cohphadm(stanm,icat,
                UncorrectedFold=uncorrected_sta_fold,CorrectedFold=corrected_sta_fold,tf=tf,
                corrected_comp=corrected_comp,raw_comp=raw_comp)
                log_list.append([f'{method}:{en}' for en in log_warnings])

                if len(cpa_list)==0:
                    print(state()+'|| No data')
                    continue
                    # raise Exception('No data')
                if si==0:f,ph = cpa_list[0].PH();ph_report.f = f

                ph_report[N][S].ph = np.array([np.abs(c.PH()[-1]) for c in cpa_list])
                # ph_report[N][S].cophadm = cpa_list
                ph_report[N][S].events = events
                stafile = f'{N[1:]}.{S}.{s_comps}.pkl'
            # netfolder = netfolder = savefolder / method / 'Networks'
            # netfolder = netfolder / N[1:]
            # netfolder.mkdir(exist_ok=True,parents=True)
            # write_pickle(netfolder / Path(file).name,ph_report[N].copy())
        completefold = savefolder / method / 'Complete'
        completefold.mkdir(exist_ok=True,parents=True)
        write_pickle(str(completefold /(f'complete_{method.lower()}.{s_comps}_ph.report.pkl')),ph_report)
        nr=1
    m=1
pd.to_pickle(lt.cat.unravel(log_list),dirs.Analysis/'Latest_DataAvailabilityWarnings.pkl')
print('--Complete--')