# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||| 1x3 . SEPARATED BANDS . EVENT RECORDS . PLOT CODE |||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
nsta_bar = 10
# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(catfolder / 'event_record_audit.pkl')
evaudit = evaudit[evaudit.nsta>=nsta_bar]
evaudit = evaudit.iloc[90:]
# evaudit = evaudit[evaudit.MW>=7.1]
display(evaudit)
# folder = 'grouped_bands'
folder = 'separated_bands'

ysep_scl = 1.3
figsize = (14,13)
bands = [(1/10,1),(1/30,1/10),(1/100,1/30)]
trim = (10,7200)
tapers = [3]
# methods = ['PostATACR','PostHPS']
methods = ['PostATACR']
# methods = ['PostHPS']
for correction_method in methods:
  coh_comp = correction_method.replace('PostHPS','HPS').replace('PostATACR','ATaCR')
  if correction_method=='PostHPS':
    return_hps = True
  else:
    return_hps = False
  for taper_mode in tapers:
    OutFolder = Path(plotfolder)
    SubFolders = Path('EventRecords') / ('Taper_' + str(taper_mode)) / correction_method / folder
    OutFolder = OutFolder / SubFolders
    OutFolder.mkdir(parents=True,exist_ok=True)
    for evi,ev in enumerate(evaudit.iloc):
        event = ev.Event
        stations = ev.Stations
        networks = ev.Networks
        event = ev.Event
        stations = ev.Stations
        networks = ev.Networks.tolist()
        evdepth = ev.depth
        post_record = Stream()
        pre_record = Stream()
        print('=='*25)
        print('|| [' + str(evi) + '/' + str(len(evaudit)) + '] ' + event + ' | ')
        print('||---Begin load')
        for i,(net,sta) in enumerate(zip(networks,stations)):
          # try:
          try:
            Metrics,Comp = get_metrics_comp(net,sta,datafolder,event,return_hps=return_hps,events_folder='EVENTS')
          except:
            _ = stations.pop(i)
            _ = networks.pop(i)
            continue
          post_record += Comp[correction_method].copy()
          pre_record += Comp['RawZ'].copy()
          del Metrics
          del Comp
          if len(post_record)==0:
            continue
        print('||---Load complete')
        phases = ('P','S','SKS','PKiKP','SKiKS','SKSSKS',)
        # phases = ('P','S',)
        # phases=('ttall',)
        evstream = post_record.copy()
        evstream_back = pre_record.copy()
        facecolor=('b','r')
        title = event
        sortindex = None
        normscale = 0.7
        residual_fraction = 0.5
        for bandi,(band,s) in enumerate(zip(bands,[[1,1],[1,1],[1,1]])):
          band_sec = np.sort([1/b for b in band])
          print('|| Plotting band: ' + str(band_sec[0]) + ' to ' + str(band_sec[-1]) + 's')
          fig, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
          # fig, axes = plt.subplots(nrows=1, ncols=3,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
          # norms='postset'
          norms='trace'
          # norms = 'col'
          # norms = np.array([[abs(c.data).max() for c in s] for s in [pre_record,post_record]]).T.max(axis=1).tolist() #global max norm
          # norms = np.array([[abs(c.data).max() for c in s] for s in [pre_record]]).T.max(axis=1).tolist() #pre max norm
          # norms = np.array([[abs(c.data).max() for c in s] for s in [post_record]]).T.max(axis=1).tolist() #post max norm. >This works well at >m6.5
          for record_i,evstream in enumerate([pre_record,post_record]):
            # File = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + event + '_' + str(band_sec[0]) + 'to' + str(band_sec[-1]) + 's_' + correction_method.replace('Post','') + '_T'
            File = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + event + '_' + str(band_sec[0]) + 'to' + str(band_sec[-1]) + 's_' + correction_method.replace('Post','')
            fig_title = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + '-'.join(event.split('.')[:2]) + ' ' + ':'.join(event.split('.')[2:]) + '_' + str(band_sec[0]) + 'to' + str(band_sec[-1]) + 's_Corrected using' + correction_method.replace('Post','')
            ax = axes[0,record_i]
            # linewidth = [0.05,0.05,0.05][bandi]
            linewidth = [0.2,0.2,0.2][bandi]
            ax_title = ['Pre-Correction \n ','Post-Correction \n '][record_i]
            title = File.replace('_',' | ').replace('z','z: ').replace('m','mag: m')
            # ---------------
            # More reasonable black and white plots:
            ax = event_record_plot(evstream=evstream,evstream_back=None,norm=norms,scales=s,linewidth=linewidth,figsize=figsize,band=band,trim=trim,facecolor=facecolor,evdepth=evdepth,phases=phases,title=title,sortindex = sortindex,ax=ax,normscale=normscale,residual_fraction=residual_fraction)
            # Use this for those ugly fill mode plots:
            # ax = event_record_plot(evstream=evstream,evstream_back=None,norm=norm,scales=s,linewidth=linewidth,figsize=figsize,band=band,trim=trim,facecolor=facecolor,evdepth=evdepth,phases=phases,title=title,sortindex = sortindex,ax=ax,normscale=normscale,residual_fraction=residual_fraction)
            # ---------------

            ax.set_title(ax_title,fontweight='bold')
            
# --------------------------------------------- Switch between grouped-band and separated bands plots
          # ##### Separated bands
          save_format = 'svg'
          fig.suptitle(fig_title.replace('_',' | ').replace('to',' to '),fontweight='bold',fontsize=15)
          File = File + '.' + save_format
          outfile = OutFolder / File
          save_tight(outfile,format=save_format,dpi=500)
          # fig.savefig(outfile,format=save_format,dpi=400)
        # ##### Grouped bands
        # fig.suptitle(File.replace('_',' | ').replace('to',' to '),fontweight='bold',fontsize=15)
        # save_tight(OutFolder / (File + '.eps' ),dpi=600)
# --------------------------------------------
          plt.close('all')
  print('000'*30)
  print(correction_method + ' ||---EV RECORDS COMPLETE---||')
  print('000'*30)