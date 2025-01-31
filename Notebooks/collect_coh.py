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


tf,method = '','hps';CorrFolder=dirs.Events_HPS/'corrected';EvFolder=dirs.Events_HPS/'rmresp'
# tf,method = 'ZP-21','atacr';EvFolder = dirs.Events/'rmresp';CorrFolder = dirs.Events/'corrected'

nets = cat.Network.unique()
savefolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences')
ovr = True
coh_report = attr()



coh_sets=[['HZ','HDH'],['HDH','HZ'],['HZ','HZ'],['H1','H1'],['H2','H2']]#[2];coh_sets=[coh_sets]

for corrected_comp,raw_comp in coh_sets:
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
            corr_fold = CorrFolder/stanm
            raw_fold = EvFolder/stanm
            if method.lower()=='atacr':
                if corrected_comp=='HZ':tf='ZP-21';corr_fold=CorrFolder/stanm
                else:corr_fold=raw_fold;tf=''

            events,cpa_list = pull_cohphadm(stanm,
            EvFolder=raw_fold,CorrFolder=corr_fold,tf=tf,
            corrected_comp=corrected_comp,raw_comp=raw_comp)

            if len(cpa_list)==0:
                print(state()+'|| No data')
                continue
                # raise Exception('No data')
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
k=1