from imports import *
# |||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||| Nx2 . COHERENCE RECORDS . PLOT CODE ||||||||
# |||||||||||||||||||||||||||||||||||||||||| =|||||||||
nsta_bar = 10
# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(catfolder / 'event_record_audit.pkl')
evaudit = evaudit[evaudit.nsta>=nsta_bar]
# evaudit = evaudit[evaudit.MW>=7.1
tapers = [0]
methods = ['PostATACR','PostHPS']
correction_method = 'PostATACR'
ysep_scl = 1.3
figsize = (10,10)
for correction_method in methods:
    coh_comp = correction_method.replace('PostHPS','HPS').replace('PostATACR','ATaCR')
    if correction_method=='PostHPS':
      return_hps = True
    else:
      return_hps = False

    OutFolder = Path(plotfolder)
    SubFolders = Path('EventRecords') / correction_method / 'coherence'
    OutFolder = OutFolder / SubFolders
    OutFolder.mkdir(parents=True,exist_ok=True)
    for evi,ev in enumerate(evaudit.iloc):
        event = ev.Event
        File = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + event + '_' + correction_method.replace('Post','') + '_COH'
        title = File.replace('_',' | ').replace('z','z: ').replace('m','mag: m')
        print('[' + str(evi) + '/' + str(len(evaudit)) + '] ' + File)
        stations = ev.Stations
        networks = ev.Networks.tolist()
        evdepth = ev.depth
        post_record = Stream()
        pre_record = Stream()
        Metrics = []
        for i,(net,sta) in enumerate(zip(networks,stations)):
          try:
            M,Comp = get_metrics_comp(net,sta,archive,event,return_hps=return_hps,events_folder='EVENTS')
            M['Noise'] = get_Noise(archive,net,sta,'sta')['Noise']
            Metrics.append(M.copy())
          except:
            _ = stations.pop(i)
            _ = networks.pop(i)
            continue
          post_record += Comp[correction_method].copy()
          pre_record += Comp['RawZ']
        if len(Metrics)==0:
          continue
        evstream = post_record.copy()
        evstream_back = pre_record.copy()
        title = event
        sortindex = list(np.argsort([np.abs(st.stats.sac.stel*1000) for st in evstream]))
        sortindex.reverse()
        Metrics = [Metrics[s] for s in sortindex]
        fig, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
        for ci,cmp in enumerate(['ZP','ZZ']):
          ax = axes[0,ci]
    # -------------
          if cmp=='ZP':
            [ax.scatter(M['Noise'].Coherence(cmp)[0],M['Noise'].Coherence(cmp)[1] + ysep*ysep_scl,c='gray',s=1) for ysep,M in enumerate(Metrics)]
            [ax.plot(M['Raw'].Coherence(cmp)[0],M['Raw'].Coherence(cmp)[1] + ysep*ysep_scl,c='b',linewidth=0.7) for ysep,M in enumerate(Metrics)]
          if cmp=='ZP':
            [ax.plot(M[coh_comp].Coherence(cmp)[0],M[coh_comp].Coherence(cmp)[1] + ysep*ysep_scl,c='r',linewidth=0.1) for ysep,M in enumerate(Metrics)]
          else:
            [ax.plot(M[coh_comp].Coherence(cmp)[0],M[coh_comp].Coherence(cmp)[1] + ysep*ysep_scl,c='r',linewidth=0.8) for ysep,M in enumerate(Metrics)]
    # -------------
          [ax.plot(M[coh_comp].Coherence(cmp)[0],M[coh_comp].Coherence(cmp)[1]*0 + ysep*ysep_scl,c='k',linewidth=0.1) for ysep,M in enumerate(Metrics)]
          fn = [fnotch(np.round(1000*abs(M['Raw'].traces['Z'][0].stats.sac.stel))) for M in Metrics]
          yticks = [ysep*ysep_scl for ysep,_ in enumerate(Metrics)]
          [ax.plot([f,f],[y,y+1],c='k',linewidth=0.4) for f,y in zip(fn,yticks)]
          [ax.text(f,y+1,'fn',horizontalalignment='center',fontsize=7,fontweight='bold') for f,y in zip(fn,yticks)]
          yticklabels = [str(round(1000*abs(M['Raw'].traces['Z'][0].stats.sac.stel))) + 'm [' + str( M['Raw'].traces['Z'][0].stats.network) + '] ' + str( M['Raw'].traces['Z'][0].stats.station) for M in Metrics]
          ax.set_xscale('log')
          f = Metrics[0][coh_comp].Coherence(cmp)[0]
          ax.set_xlim(f[1],f[-1])
          ax_title = cmp + ' Coherence'
          ax.set_title(ax_title,fontweight='bold')
          ax.set_yticks(yticks)
          ax.set_yticklabels(yticklabels)
          ax.set_ylim(yticks[0],yticks[-1]+1.3)
        fig.suptitle(File.replace('_',' | '),fontweight='bold',fontsize=15)
        save_tight(OutFolder / (File + '.eps' ),dpi=700)
        plt.close('all')
        # || 26min to complete