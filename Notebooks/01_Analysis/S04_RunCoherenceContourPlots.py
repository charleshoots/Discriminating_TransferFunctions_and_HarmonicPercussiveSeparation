import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
from imports import *
from modules import *
hpsfold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
# hpsfold = Path('/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data')
atacrfold = dirs.Events

cat = catalog.copy()
# cat = cat.iloc[:3]

# Comps = ['ZZ','11','22','PP']
Comps = ['ZZ']

Methods = ['ATaCR','NoiseCut']

plotfold = dirs.P01.S04


cat.sort_values(by='StaDepth',inplace=True)

# for ci,Comp in enumerate(Comps):
for mi,method in enumerate(Methods):
    # All_AvgCorrectedCoh,All_AvgRawCoh = AttribDict(),AttribDict()
    # [All_AvgCorrectedCoh.update({c:[]}) for c in Comps];[All_AvgRawCoh.update({c:[]}) for c in Comps]
    # for si,sta in enumerate(cat.iloc):
    #     if method=='ATaCR':tf='sta.ZP-21.HZ.SAC';evdir=atacrfold;mirror_fold=hpsfold
    #     elif method=='Noisecut':tf='HZ.SAC';evdir=[hpsfold,dirs.Events];mirror_fold=atacrfold
    #         # ======
    #         # ------
    #     # clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
    #     # print(f'[{si+1}/{len(cat)}]{sta.StaName}|[{mi+1}/{len(Methods)}]]{method}')
    #     if method.lower()=='noisecut':
    #         st_hold,evmeta = lt.io.get_station_events_hps(sta.StaName,evdir,tf=tf,type='metrics',evmeta=sta.Events)
    #         raw=st_hold.select(location='*Raw*').copy()
    #         corrected=st_hold.select(location='*Corrected*').copy()
    #     else:
    #         st_hold,evmeta = lt.io.get_station_events(sta.StaName,evdir,tf=tf,type='metrics',evmeta=sta.Events)
    #         raw=st_hold.select(location='*Raw*').copy()
    #         corrected=st_hold.select(location='*Corrected*').copy()
    #     # clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
    #     for ci,Comp in enumerate(Comps):
    #         # clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
    #         print(f'Loading ... [{si+1}/{len(cat)}]{sta.StaName} | [{mi+1}/{len(Methods)}]] [{method}] | :: [{ci+1}/{len(Comps)}]]{Comp}')
    #         f=raw[0].Metrics.Coherence(Comp)[0]
    #         RawCoh=np.array([tr.Metrics.Coherence(Comp)[1] for tr in raw])
    #         CorrectedCoh=np.array([tr.Metrics.Coherence(Comp)[1] for tr in corrected])
    #         AvgRawCoh=RawCoh.mean(axis=0)
    #         AvgCorrectedCoh=CorrectedCoh.mean(axis=0)
    #         del RawCoh,CorrectedCoh
    #         All_AvgRawCoh[Comp].append(AvgRawCoh)
    #         All_AvgCorrectedCoh[Comp].append(AvgCorrectedCoh)
    #         # ======
    #         # ------
    #         # fig,ax = plt.subplots(figsize=(13,6))
    #         # [ax.scatter(f,c,c='k',s=0.5) for c in CorrectedCoh]
    #         # ax.set_xscale('log');ax.set_ylim(0,1);ax.set_xlim(1/500,1)

    # --------------------------------------------------------
    # for c in Comps:All_AvgCorrectedCoh[c]=np.array(All_AvgCorrectedCoh[c])
    # for c in Comps:All_AvgRawCoh[c]=np.array(All_AvgRawCoh[c])
    vlim = [0,1];levels=np.linspace(0,1,40)
    for ci,Comp in enumerate(Comps):
        if (method=='ATaCR')&(not Comp=='ZZ'):continue
        report = get_reports(Comp,catalog,dirs.Archive,dirs,AVG=True,methods=[method])
        c_report = report[method]
        f=report.f
        # c_report[f'n{n}'][s].coh.shape[np.intersect1d(c_report[f'n{n}'][s].events,[e.Name for e in lt.cat.mirror('7D.FS42D',cat,dirs)],return_indices=True)[1],:]

        # np.intersect1d(c_report[f'n{'7D'}'][s].events,[e.Name for e in lt.cat.mirror('7D.FS42D',cat,dirs)],return_indices=True)[1]
        n,s=cat.Network[0],cat.Station[0]
        coh=np.array([np.nanmean(
        c_report[stnm].coh[np.intersect1d(c_report[stnm].events,[e.Name for e in lt.cat.mirror(stnm,cat,dirs)],return_indices=True)[1],:],
        axis=0) for stnm in cat.StaName])

        z=[z for z in cat.StaDepth]
        # clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
        print(f'Plotting ... [{mi+1}/{len(Methods)}]] [{method}] | :: [{ci+1}/{len(Comps)}]]{Comp}')
        # --------------------------------------------------------
        # if mi==0:
        #     file=plotfold/f'Stations.Contour.{Comp}.Uncorrected.StationEventAverage.png'
        #     z=cat.StaDepth;coh=All_AvgRawCoh[Comp]
        #     fig,ax,cnt = lt.plots.dataset_averaged_coherence_plot(f,z,coh,
        #     title=f'Station Event Average Uncorrected {Comp} Coherences',
        #     vlim=vlim,levels=levels)
        #     save_tight(file,fig)
        # --------------------------------------------------------
        file=plotfold/f'Contour.{method}.{Comp}_Coherence.StationEventAverage.png'
        # _AltCode
        # z=cat.StaDepth;coh=All_AvgCorrectedCoh[Comp]
        fig = lt.plots.dataset_averaged_coherence_plot(f,z,coh,
        title=f'Station Event Average {method} {Comp} Coherences',
        vlim=vlim,levels=levels)
        save_tight(file,fig)
        # --------------------------------------------------------
        plt.close('all')