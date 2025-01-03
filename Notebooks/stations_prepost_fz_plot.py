from imports import *
from modules import *
hpsfold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/HPS_Data/Data')
# hpsfold = Path('/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data')
atacrfold = dirs.Events

cat = catalog.copy()
# cat = cat.iloc[:3]

Comps = ['ZZ','ZP','Z1','Z2']
Methods = ['ATaCR','Noisecut']

for mi,method in enumerate(Methods):
    All_AvgCorrectedCoh,All_AvgRawCoh = AttribDict(),AttribDict()
    [All_AvgCorrectedCoh.update({c:[]}) for c in Comps];[All_AvgRawCoh.update({c:[]}) for c in Comps]
    for si,sta in enumerate(cat.iloc):
        if method=='ATaCR':tf='sta.ZP-21.HZ.SAC';evdir=atacrfold;mirror_fold=hpsfold
        elif method=='Noisecut':tf='HZ.SAC';evdir=[hpsfold,dirs.Events];mirror_fold=atacrfold
            # ======
            # ------
        clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
        print(f'[{si+1}/{len(cat)}]{sta.StaName}|[{mi+1}/{len(Methods)}]]{method}')
        if method.lower()=='noisecut':
            st_hold,evmeta = get_station_events_hps(sta.StaName,evdir,tf=tf,mirror_fold=mirror_fold,type='metrics',evmeta=sta.Events)
            raw=st_hold.select(location='*Raw*').copy()
            corrected=st_hold.select(location='*Corrected*').copy()
        else:
            st_hold,evmeta = get_station_events(sta.StaName,evdir,tf=tf,mirror_fold=mirror_fold,type='metrics',evmeta=sta.Events)
            raw=st_hold.select(location='*Raw*').copy()
            corrected=st_hold.select(location='*Corrected*').copy()
        clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
        for ci,Comp in enumerate(Comps):
            clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
            print(f'[{si+1}/{len(cat)}]{sta.StaName}|[{mi+1}/{len(Methods)}]]{method}|[{ci+1}/{len(Comps)}]]{Comp}')
            f=raw[0].Metrics.Coherence(Comp)[0]
            RawCoh=np.array([tr.Metrics.Coherence(Comp)[1] for tr in raw])
            CorrectedCoh=np.array([tr.Metrics.Coherence(Comp)[1] for tr in corrected])
            AvgRawCoh=RawCoh.mean(axis=0)
            AvgCorrectedCoh=CorrectedCoh.mean(axis=0)
            del RawCoh,CorrectedCoh
            All_AvgRawCoh[Comp].append(AvgRawCoh)
            All_AvgCorrectedCoh[Comp].append(AvgCorrectedCoh)
            # ======
            # ------
            # fig,ax = plt.subplots(figsize=(13,6))
            # [ax.scatter(f,c,c='k',s=0.5) for c in CorrectedCoh]
            # ax.set_xscale('log');ax.set_ylim(0,1);ax.set_xlim(1/500,1)

    # --------------------------------------------------------
    for c in Comps:All_AvgCorrectedCoh[c]=np.array(All_AvgCorrectedCoh[c])
    for c in Comps:All_AvgRawCoh[c]=np.array(All_AvgRawCoh[c])
    vlim = [0,1];levels=np.linspace(0,1,40)
    for ci,Comp in enumerate(Comps):
        clear_output(wait=False);os.system('cls' if os.name=='nt' else 'clear')
        print(f'[{si+1}/{len(cat)}]{sta.StaName}|[{mi+1}/{len(Methods)}]]{method}|[{ci+1}/{len(Comps)}]]{Comp}')
        # --------------------------------------------------------
        if mi==0:
            file=dirs.Plots/'_Plots'/f'Stations.Contour.{Comp}.Uncorrected.StationEventAverage_AltCode.png'
            z=cat.StaDepth;coh=All_AvgRawCoh[Comp]
            fig,ax,cnt = dataset_averaged_coherence_plot(f,z,coh,
            title=f'Station Event Average Uncorrected {Comp} Coherences',
            vlim=vlim,levels=levels)
            save_tight(file,fig)
        # --------------------------------------------------------
        file=dirs.Plots/'_Plots'/f'Stations.Contour.{Comp}.{method}.StationEventAverage_AltCode.png'
        z=cat.StaDepth;coh=All_AvgCorrectedCoh[Comp]
        fig,ax,cnt = dataset_averaged_coherence_plot(f,z,coh,
        title=f'Station Event Average {method} {Comp} Coherences',
        vlim=vlim,levels=levels)
        save_tight(file,fig)
        # --------------------------------------------------------
    del All_AvgRawCoh,All_AvgCorrectedCoh