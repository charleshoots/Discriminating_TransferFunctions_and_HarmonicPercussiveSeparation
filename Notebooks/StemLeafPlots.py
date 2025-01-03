# # ||||||||||||||||||||||||||||||||||||||
# # ||||||||||||||| STEM/LEAF PLOTS ||||||
# # ||||||||||||||||||||||||||||||||||||||
# cat = catalog.copy()
# taper_mode = 0
# return_noise = True
# return_hps = False
# events_folder = 'EVENTS_Taper_' + str(taper_mode)
# figsize = (20,14)
# coh_comps = ['ATaCR_ZZ','Noise_ZP','Raw_ZP','ATaCR_ZP','Raw_1Z','ATaCR_1Z','Raw_2Z','ATaCR_2Z']
# nrows = len(coh_comps) + 1
# fig, axes = plt.subplots(nrows=nrows, ncols=1,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
# for stai,Station in enumerate(cat.iloc):
#   Metrics = []
#   for evi,event in enumerate(Station.Events):
#     notify = Station.StaName + ' ' + str(stai+1) + '/' + str(len(cat)) + ' || ' + str(evi+1) + '/' + str(len(Station.Events))
#     try:
#       M,Comp = get_metrics_comp(Station.Network,Station.Station,datafolder,event,return_hps=return_hps,return_noise=return_noise,events_folder=events_folder)
#       Metrics.append(M.copy())
#       print(notify)
#     except:
#       print(notify + ' | 404 Error')
#       continue
#   X = [stai]
#   stats = metric_stats(Metrics)
#   pre_stalta = stats['log_rms_stalta'][:,0]
#   post_stalta = stats['log_rms_stalta'][:,1]
#   ax = axes[0,0]
#   metric_title = 'STA/LTA'
#   for pi,p in enumerate([pre_stalta,post_stalta]):
#     if stai==0:
#       label = ['Pre','Post'][pi]
#     else:
#       label = None
#     Y = [p.mean()]
#     c = ['b','r'][pi]
#     ls = ['-','--'][pi]
#     minmax = np.atleast_2d(np.array([p.min(),p.max()])).T
#     ax.scatter(X,Y,c=c,s=800,marker='_',label=label)
#     ax.scatter(X,Y-p.std(),c=c,s=100,marker='_')
#     ax.scatter(X,Y+p.std(),c=c,s=100,marker='_')
#     # ax.errorbar(X,Y,yerr=minmax,ecolor=c,ls=ls)
#     ax.vlines(X,minmax[0],minmax[1],colors = c,linestyles=ls,linewidth=3)
#     # ax.errorbar(X,Y,yerr=p.std(), ecolor='r',capthick=0)
#   ax.set_title(metric_title + '\n Min/Max : Vertical Lines ' + ' '*10 + 'σ : Short Horiz. Lines' + ' '*10 + 'μ : Long Horiz. Lines',fontweight='bold')
#   ax.set_ylim(0)
#   ax.legend(ncol=2,loc='upper left')
#   band_colors = '#901f62','#588a7d','#284ea5'
#   bands = [[1/100,1/30],[1/30,1/10],[1/10,1]]
#   bands = 1/np.sort(bands)
#   for axi,coh in enumerate(coh_comps):
#     ax = axes[axi+1,0]
#     [ax.stem(X,np.mean((stats['coh_mu_rms_bands'][coh][:,bi])**2)**0.5,markerfmt=band_colors[bi]) for bi in range(len(bands))]
#     ax.set_title(coh.replace('_','-') + ' Coherence RMS')
#     if stai==0:
#       if axi==0:
#         [ax.stem(X,np.mean((stats['coh_mu_rms_bands'][coh][:,bi])**2)**0.5,markerfmt=band_colors[bi],label=str(bands[bi][1]) + '-' + str(bands[bi][0]) + 's',linefmt=band_colors[bi],basefmt='w') for bi in range(len(bands))]
#     ax.legend(ncol=3,loc='upper left')
#     ax.set_ylim(0)
# ax = axes[-1,0]
# xticks = [xsep for xsep,_ in enumerate(cat.iloc)]
# xticklabels = [c.StaName for c in cat.iloc]
# ax.set_xticks(xticks)
# ax.set_xticklabels(xticklabels,rotation = 90)
# folder = Path(plotfolder) / 'EventRecords' / 'Stats'
# folder.mkdir(exist_ok=True)
# files = 'signal_metrics_taper_' + str(taper_mode)
# save_tight(folder / (files + '.eps' ),dpi=500)
# # || 77min to complete