from imports import *
# ========================================================================================================================================================
# ================================================================CODE SNIPPETS===========================================================================
# ========================================================================================================================================================
# for a in cm._cmap_names_categorical:
#     display(cm.__dict__[a].resampled(4))
# # ======================================================================================================================================================
# # ======================================================================================================================================================
# for (Event,Station,Metrics,Comp) in OBS_Generator(catalog,dirs['Py_DataParentFolder']):
#     print(Event)
# for i,(Event,Station,Metrics,Comp) in zip(range(1),OBS_Generator(catalog,dirs['Py_DataParentFolder'])):
#     print(Event)
# # ======================================================================================================================================================
# # ======================================================================================================================================================
def smooth(d,k=10):return np.convolve(d, np.ones(k) / k, mode='same')
# Station,evi = catalog.iloc[22],3
# Event = Station.Events[evi]
# display(catalog)
# ========================================================================================================================================================
# ^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^^_^
# ========================================================================================================================================================

NoiseColors = [mcolors.to_hex(m) for m in [cm.__dict__[e].resampled(30).resampled(6).colors for e in ['devon_categorical']][0]]
# [display(c) for c in [cm.__dict__[e].resampled(70).resampled(5) for e in ['nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]]
# np.array([[mcolors.to_hex(c) for c in r] for r in [cm.__dict__[e].resampled(70).resampled(5).colors for e in ['nuuk_categorical','devon_categorical','hawaii_categorical','imola_categorical','lapaz_categorical']]])

# import Nx2_CoherenceRecords

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++ CONSTRUCTOR AREA ++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++ CONSTRUCTOR AREA ++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
methods = ['PostATACR']
atacrdatafolder = archive / 'ATaCR_Data' / 'ATaCR_Python'
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
    for station in catalog.iloc:
       stations = [station.Station]
       networks = [station.Network]
       events = station.Events
    #    events = ['2012.069.07.09','2012.181.21.07']
    for i,(net,sta) in enumerate(zip(networks,stations)):
        Metrics = []
        for evi,event in enumerate(events):
            depth = round(station.Metadata[evi].origins[0].depth/1000)
            mag = station.Metadata[evi].magnitudes[0].mag
            File = 'm' + str(mag) + '_z' + str(depth) + 'km' + '_' + event + '_' + correction_method.replace('Post','') + '_COH'
            title = File.replace('_',' | ').replace('z','z: ').replace('m','mag: m')
            print('[' + str(evi) + '/' + str(len(events)) + '] ' + File)
            post_record = Stream()
            pre_record = Stream()

            M,Comp = get_metrics_comp(net,sta,atacrdatafolder,event,return_hps=return_hps,events_folder='EVENTS')
            # M['Noise'] = get_Noise(atacrdatafolder,net,sta,'sta')['Noise']
            Metrics.append(M.copy())


pairs = ['ZP']
meters = ['psd','Coherence','Phase','Admittance']
fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(8,12),layout='constrained',squeeze=False,sharey='all',sharex='all')
axes = axes.reshape(-1)


Pre = Metrics[0]['Raw']
Post = Metrics[0]['ATaCR']
Noise = Metrics[0]['Noise']
for pi,(ax,m) in enumerate(zip(axes,meters)):
    evf,prey = Pre.__getattribute__(m)(pairs[0])
    evf,posty = Post.__getattribute__(m)(pairs[0])
    if m=='psd':
        noisef,noisey = Noise.f,Noise.StaNoise.power.__dict__['cZZ'].shape
    ax.set_xlabel('Frequency')
    ax.set_ylabel('Coherence')
    ax.set_title(pair + '-Coherence',fontweight='bold')
    if not pair=='ZZ':
        ax.scatter(noisef,noisecoh,s=2,c='gray')
    ax.scatter(evf,precoh,s=0.5,c='b')
    ax.scatter(evf,postcoh,s=0.5,c='r')
    ax.set_xscale('log')
# plt.tight_layout()
