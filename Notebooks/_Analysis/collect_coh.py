import os,sys;os.chdir('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/Notebooks')
sys.path.append('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/Notebooks')
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
nets = cat.Network.unique()
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences')



ovr = True

# tf,method = '','hps';CorrectedFold=dirs.Events_HPS/'corrected';UncorrectedFold=dirs.Events_HPS/'rmresp'
tf,method = 'ZP-21','atacr';CorrectedFold = dirs.Events/'corrected';UncorrectedFold = dirs.Events/'rmresp'
# ------------------------------------------------------------------------------------------------
# coh_sets=[['HZ','HDH'],['HDH','HZ'],['HZ','HZ'],['H1','H1'],['H2','H2']]#Every auto-coherence, only useful for NoiseCut
coh_sets=[['HZ','HZ'],['HZ','HDH'],['HZ','H1'],['HZ','H2']]#Defaults: ZZ,ZP,Z1,Z2

if method.lower()=='hps':coh_sets.extend([['HDH','HDH'],['H1','H1'],['H2','H2']])
# ------------------------------------------------------------------------------------------------


for corrected_comp,raw_comp in coh_sets:
    coh_report = attr()
    if method.lower()=='atacr':
        if (corrected_comp=='H2') or (corrected_comp=='H1'):continue
    os.system('clear')
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
            state = lambda:f'{method.upper()} - {s_comps} || {sta.StaName} | NET: {str(ni+1)} / {str(len(nets))} | STA: {str(si+1)} / {str(len(icat))} '
            S = sta.Station
            coh_report[N][S] = attr()
            stanm = '.'.join([sta.Network,sta.Station])
            print(state())
            corrected_sta_fold = CorrectedFold/stanm
            uncorrected_sta_fold = UncorrectedFold/stanm
            if method.lower()=='atacr':
                if corrected_comp=='HZ':tf='ZP-21';corrected_sta_fold=CorrectedFold/stanm
                else:corrected_sta_fold=uncorrected_sta_fold;tf=''

            events,cpa_list = pull_cohphadm(stanm,
            UncorrectedFold=uncorrected_sta_fold,CorrectedFold=corrected_sta_fold,tf=tf,
            corrected_comp=corrected_comp,raw_comp=raw_comp)

            if len(cpa_list)==0:
                print(state()+'|| No data')
                continue
                # raise Exception('No data')
            if si==0:f,coh = cpa_list[0].COH();coh_report.f = f

            coh_report[N][S].coh = np.array([c.COH()[-1] for c in cpa_list])
            # coh_report[N][S].cophadm = cpa_list
            coh_report[N][S].events = events
            stafile = f'{N[1:]}.{S}.{s_comps}.pkl'
        netfolder = netfolder = savefolder / method / 'Networks'
        netfolder = netfolder / N[1:]
        netfolder.mkdir(exist_ok=True,parents=True)
        write_pickle(netfolder / Path(file).name,coh_report[N].copy())
    completefold = savefolder / method / 'Complete'
    completefold.mkdir(exist_ok=True,parents=True)
    write_pickle(str(completefold /(f'complete_{method.lower()}.{s_comps}_coh.report.pkl')),coh_report)
k=1