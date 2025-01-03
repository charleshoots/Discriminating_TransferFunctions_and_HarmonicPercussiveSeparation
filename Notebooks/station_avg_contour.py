Comps = ['ZZ','ZP','Z1','Z2']
for comp in Comps:
    Report = get_reports(comp,Archive)
    f = Report.f
    AVG = True
    for stanm in catalog.StaName:
        if AVG:
            Report.Uncorrected[stanm].coh=Report.Uncorrected[stanm].coh.mean(axis=0)
            Report.ATaCR[stanm].coh=Report.ATaCR[stanm].coh.mean(axis=0)
            Report.Noisecut[stanm].coh=Report.Noisecut[stanm].coh.mean(axis=0)
    for method in ['Uncorrected','ATaCR','Noisecut']:
        file=dirs.Plots/'_Plots'/f'Stations.Contour.{comp}.{method}.StationEventAverage.png'
        coh = np.array([Report[method][stanm].coh for stanm in catalog.StaName])
        z=np.array(catalog.StaDepth.to_list())
        vlim = [0,1];levels=np.linspace(0,1,300)
        fig,ax,cnt = dataset_averaged_coherence_plot(f,z,coh,
        title=f'Station Event Average {method} {comp} Coherences',
        vlim=vlim,levels=levels)
        save_tight(file,fig)