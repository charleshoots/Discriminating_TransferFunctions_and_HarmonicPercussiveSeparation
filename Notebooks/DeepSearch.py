from imports import *
from scipy.stats import iqr
from local_tools.math import cohstats
from scipy.ndimage import gaussian_filter
from local_tools.plots import ax_sta_metrics
cat = catalog.copy()
octavg=lt.math.octave_average
# TF-data
# _______________________________________________________________________________________________________________________
report=get_reports('ZZ',catalog,dirs.Archive,dirs,AVG=True)
get_tf = lambda stanm:load_pickle(list((dirs.TransferFunctions/stanm).glob('*-*.pkl'))[0])
get_day_tfs = lambda stanm: [load_pickle(i) for i in [f for f in list((dirs.TransferFunctions/stanm).glob('*.pkl')) if str(f).find('-')<0]]
# get_day_tfs = lambda stanm,days: [load_pickle(dirs.TransferFunctions/stanm/i) for i in days]
# _______________________________________________________________________________________________________________________

minmag,maxmag=0,0
flim = [1/230,1]
figsize=(6.5,3)
octave_av=False
nsta=15


# fold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_FigureArchive/_Papers/Ch1/temp')
# d = np.array([[i.split('.')[-2],i.split('.Mw')[0]] for i in [f.name.replace('_psd.and.psd_residual.png','') for f in list(fold.glob('*.png'))]])
# stas = d[:,0]
# evs = d[:,1]
# evs = ['2010.251.11.37','2012.181.21.07','2011.042.20.05','2013.039.11.12','2015.140.22.48','2012.312.16.35','2015.127.07.10']

fold = Path('/Users/charlesh/Downloads/untitled')
files = [f.name.replace('_psd.and.psd_residual.png','').replace('_psd.and.psd_residual_octave.png','') for f in list(fold.rglob('*.png'))]
evs = np.unique([f.split('.Mw')[0] for f in files])
stas = np.unique(['.'.join(f.split('.Mw')[1].split('_')[-1].split('.')[:-1]) for f in files])

# _______________________________________________________________________________________________________________________
magfold = f'{int(minmag*10)}_{int(maxmag*10)}'
cat=catalog.copy()

deepstations = cat.sort_values(by='StaDepth')[-nsta:].copy()
shallowstations=cat.sort_values(by='StaDepth')[:nsta].copy()
deepevents = [e.Name for e in lt.cat.unravel_cat(deepstations)]
shallowevents = [e.Name for e in lt.cat.unravel_cat(shallowstations)]
cat = pd.concat([deepstations,shallowstations])



# _______________________________________________________________________________________________________________________
# outlierprops={'color':'r','s':10,'alpha':0.09}
# inlierprops={'color':'dodgerblue','s':13,'alpha':0.2} #'#7370cb' royalblue dodgerblue
# whiskerprops={'linewidth':0.5,'color':'k'}
# whiskerwidth=0.002;margin=1.02
# midlineprops={'linewidth':1.5,'color':'k','alpha':0.8}
# loglogprops={'linewidth':1.5,'color':'k','alpha':0.8,'nonpositive':'mask'}
# midscatterprops={'s':10,'color':'k'} #args.csd_pairs[pair]
notchprops = {'alpha':0.4,'linewidth':1,'color':'w','linestyle':'-.'}
labelprops = {'fontsize':5}
# methodlabelprops = {'fontweight':'bold','fontsize':5,'verticalalignment':'bottom','horizontalalignment':'right'}
# _______________________________________________________________________________________________________________________
status = lambda:print(f'M: {mi+1}/{4} | S: {si+1}/{len(cat)} | E: {evi+1}/{len(sta.Events)} :: {method}')
colors = [(1, 1, 1), (0.5, 0.5, 0.5),(0, 0, 0)]  # White -> Black -> Gray
positions = [0, 0.5, 1]  # Positions for black, white, and gray
custom_cmap_abs = mcolors.LinearSegmentedColormap.from_list("CustomDiverging", list(zip(positions, colors)))
colors = [(0.5, 0.5, 0.5),(1, 1, 1),'r']  # White -> Gray -> Red
positions = [0, 0.5, 1]  # Positions for black, white, and gray
custom_cmap = mcolors.LinearSegmentedColormap.from_list("CustomDiverging", list(zip(positions, colors)))
# _______________________________________________________________________________________________________________________
cbar_dy = 0.13
titlefontszie = 7
win_length=163.84
xlim=[0+(2*win_length/7200),2-(2*win_length/7200)]

# _______________________________________________________________________________________________________________________

# for octave_av in [False,True]:
deep_shallow_parity_test = lambda event:(np.isin(event.Name,shallowevents)) & (np.isin(event.Name,deepevents))
for mi,(minmag,maxmag) in enumerate([[6.0,6.3],[6.3,6.8],[6.8,7.4],[7.4,8.0]]):
    # ------------------------------------
    for si,(sta,ns) in enumerate(zip(cat.iloc,range(len(cat)))):



        magfold = f'{int(minmag*10)}_{int(maxmag*10)}'
        stanm=sta.StaName
        events = sta.Events
        events = [e for e in events if e.magnitudes[0].mag>minmag]
        events = [e for e in events if e.magnitudes[0].mag<=maxmag]
        events = [events[i] for i in np.flip(np.argsort([e.magnitudes[0].mag for e in events]))] # Sort by magnitude
        events = [e for e in events if deep_shallow_parity_test(e)] # Subset by only events that have both a deep and shallow observing station
        nev=50 #Reduce this number to cap the number of events explored.
        for evi,(event,ne) in enumerate(zip(events,range(nev))):
            # if not np.isin(stanm,stas):continue
            # if not np.isin(event.Name,evs):continue
            if not deep_shallow_parity_test(event):continue


            for method_ind,method in enumerate(['NoiseCut','ATaCR']):
                status()
                try:
                    traces = get_traces(stanm,event.Name)
                    specs = {tr.stats.location:spectrogram(tr.copy(),win_length=win_length,pow2db=False) for tr in traces}
                    if octave_av:
                        for k in specs.keys():
                            S,f,t=specs[k]
                            foct,Soct=octavg(S.T,f)
                            specs[k] = Soct,foct,t

                    for k in ['ATaCR','NoiseCut']:
                        temp= list(specs['Raw']);temp[0]=specs[k][0] - temp[0] # RESIDUAL = CORRECTED-RAW
                        # temp= list(specs['Raw']);temp[0]=temp[0]-specs[k][0]
                        specs[k+'_Residual'] = temp
                        del temp
                    # vlim=[np.min([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()]),np.max([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()])]
                    vlim=np.array([-100,-40])
                except:
                    print(f'{event.Name} : {stanm} : Data?')
                    continue


                # ------- FIGURE SETUP
                fig = plt.figure(figsize=figsize)
                gapwidth=.02 #Gap between subplots
                cohwidth=0.20 #Relative Width of coherence subplot
                spectwidth=0.40 #Relative Width of spectrograms
                gs = gridspec.GridSpec(1, 5, width_ratios=[cohwidth, gapwidth,spectwidth, gapwidth, spectwidth], wspace=0.09)
                ax1 = fig.add_subplot(gs[0, 0])  # First subplot
                gap1 = fig.add_subplot(gs[0, 1]);gap1.set_visible(False)
                ax2 = fig.add_subplot(gs[0, 2])  # Second subplot
                gap2 = fig.add_subplot(gs[0, 3]);gap2.set_visible(False)
                ax3 = fig.add_subplot(gs[0, 4])  # Third subplot
                axes = [ax1, ax2, ax3]

                SUBPLOT_A = True
                SUBPLOT_A_Add_Residual_Avergae = True
                if SUBPLOT_A:
                    # ------- a) COHERENCE STATION AVERAGE (RED LINE + BLACK STEMS)
                    ax=axes[0]
                    ev_name = UTCDateTime.strptime(event.Name,'%Y.%j.%H.%M').strftime('%Y-%m-%d %H:%M:00')
                    ax_sta_metrics(ax,report,sta,method,event_name=event.Name,flim=flim,octave_av=octave_av)
                    ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
                    ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
                    ax.set_title(f'a)',loc='left', pad=0,y=1.1,x=-.1,fontsize=titlefontszie)
                    ax.set_xlabel('Coherence',fontsize=titlefontszie,labelpad=11)
                    if SUBPLOT_A_Add_Residual_Avergae:
                        # ------- a) RESIDUAL AVERAGE (BLUE LINE)
                        # Residual average (BLUE LINE), or perhaps a Compensated-Coherence? (e.g. Coherence * exp(-((STFT(a)-STFT(b))**2)))
                        # Create a second x-axis
                        ax02 = ax.twiny()
                        S_Raw,spec_f,spec_t=specs['Raw']
                        S_Post,spec_f,spec_t=specs[method]
                        S=np.abs(S_Post-S_Raw)**2
                        S=librosa.power_to_db(S)
                        # y=np.median(S,axis=1)
                        y=S.mean(axis=1)
                        # y=S.max(axis=1)
                        ax02.plot(y,spec_f,**{'linewidth':0.75,'color':'blue','alpha':0.8})
                        ax02.scatter(y,spec_f,c='blue',s=0.5)
                        if method_ind==0:res_avg_vlim=[int(S.mean(axis=1).min()), int(0.7*(S.mean(axis=1).max()))]
                        ax02.set_xlim(res_avg_vlim)
                        # ax02.set_yscale('log')
                        ax02.set_xticklabels(ax02.get_xticklabels(),**labelprops)
                        # ax02.set_yticklabels(ax02.get_yticklabels(),**labelprops)
                        ax02.set_xlabel('Residual Average',fontsize=titlefontszie,labelpad=11)
                        ax02.xaxis.set_label_position('top')  # Move label to the top
                        ax02.xaxis.set_ticks_position('top')  # Move ticks to the top
                        ax02.invert_xaxis()
                
                SUBPLOT_B = True
                if SUBPLOT_B:
                    # ------- b) CORRECTED SPECTROGRAM
                    if method=='ATaCR':ax=axes[1];k='ATaCR'
                    else:ax=axes[1];k='NoiseCut'
                    S,spec_f,spec_t=specs[k]
                    S = np.abs(S)**2
                    S=librosa.power_to_db(np.abs(S))
                    # ---------------------------------PLOT
                    specfig,pcm=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=vlim,cbar=False)
                    # cbar_ax1 = fig.add_axes([ax2.get_position().x0,ax2.get_position().y1+0.02,ax2.get_position().width,0.02])
                    cbar_ax1 = fig.add_axes([ax2.get_position().x0,ax2.get_position().y0-0.03,ax2.get_position().width,0.02])
                    cbar1 = fig.colorbar(pcm,cax=cbar_ax1,pad= 0.01,orientation='horizontal',location='bottom')
                    cbar1.ax.tick_params(labelsize=labelprops['fontsize'])
                    cbar1.set_label(f'Power Spectral Density, dB\nPost-{method}', fontsize=7)  # Cust
                    cbarticks=[vlim[0],0 if np.all([min(vlim)<0,max(vlim>0)]) else int(np.mean(vlim)) ,vlim[-1]]
                    cbar1.set_ticks(cbarticks)
                    cbar1.set_ticklabels(cbarticks,**labelprops)
                    ax.set_ylim(flim[0],flim[1])
                    ax.axhline(fnotch(sta.StaDepth),**notchprops)
                    ax.set_title(f'b)',loc='left', pad=0,y=1.1,x=-.1,fontsize=titlefontszie)
                    ax.set_yticklabels([])
                    ax.set_xlabel('Time after origin, hours',**labelprops)
                    ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800)/3600)
                    ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
                    ax.set_xlim(xlim)
                    # cbar1.ax.xaxis.set_ticks_position('top')
                    # cbar1.ax.xaxis.set_label_position('top')
                    ax.xaxis.set_ticks_position('top')
                    ax.xaxis.set_label_position('top')

                SUBPLOT_C = True
                if SUBPLOT_C: # c) RESIDUAL SPECTROGRAM
                    # | DATA | ---------------------------------
                    # ---------------------------------Calculate Residual Density
                    S_Raw,spec_f,spec_t=specs['Raw']
                    S_Post,spec_f,spec_t=specs[k]
                    S = S_Post-S_Raw
                    S = librosa.power_to_db(np.abs(S)**2)
                    note='psd.and.psd_residual'
                    cbartitle = f'$|\Delta|$ Residual Spectral Density, dB\n|Post-{method} - Original|'
                    cmap=custom_cmap_abs
                    if method_ind==0: #Forces NoiseCut and ATaCR to share the exact same clim for comparability
                        # residual_vlim = [-100,-40];subfold = 'Original'
                        residual_vlim=[int(S.mean(axis=1).min()), int(0.7*(S.mean(axis=1).max()))];subfold = 'Adjusted_clim'
                    vlim=residual_vlim
                    # | PLOT | ---------------------------------
                    vlim=np.array(vlim)
                    ax=axes[2]
                    specfig,pcm=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=vlim,cbar=False,cmap=cmap)
                    # cbar_ax2 = fig.add_axes([ax3.get_position().x0,ax3.get_position().y1+0.02,ax3.get_position().width,0.02])
                    cbar_ax2 = fig.add_axes([ax3.get_position().x0,ax3.get_position().y0-0.03,ax3.get_position().width,0.02])
                    cbar2 = fig.colorbar(pcm,cax=cbar_ax2,pad= 0.01,orientation='horizontal',location='bottom') 
                    cbar2.ax.tick_params(labelsize=labelprops['fontsize'])
                    cbar2.set_label(cbartitle, fontsize=7)  # Customize font size and weight
                    if vlim[0] is not None and vlim[1] is not None:cbarticks=[vlim[0],0 if np.all([min(vlim)<0,max(vlim>0)]) else int(np.mean(vlim)) ,vlim[-1]]
                    else:cbarticks=np.linspace(cbar2.get_ticks().min(),cbar2.get_ticks().max(),4)
                    cbar2.set_ticks(cbarticks)
                    cbar2.set_ticklabels(cbarticks,**labelprops)
                    ax.set_ylim(flim[0],flim[1])
                    ax.axhline(fnotch(sta.StaDepth),**notchprops)
                    ax.set_title(f'c)',loc='left', pad=0,y=1.1,x=-.1,fontsize=titlefontszie)
                    ax.set_yticklabels([])
                    ax.set_xlabel('Time after origin, hours',**labelprops)
                    ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800)/3600)
                    ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
                    ax.set_xlim(xlim)
                    # cbar2.ax.xaxis.set_ticks_position('top')
                    # cbar2.ax.xaxis.set_label_position('top')
                    ax.xaxis.set_ticks_position('top')
                    ax.xaxis.set_label_position('top')
                    ax1.set_position([ax2.get_position().x0 - ax2.get_position().width - 0.05,
                    ax2.get_position().y0,ax2.get_position().width,ax2.get_position().height])
                    # fig.tight_layout(pad=0.1)
                    fig.subplots_adjust(wspace=0.6, hspace=0.0)

                SAVE_OUTPUT = True
                if SAVE_OUTPUT: # SAVE OUTPUT
                    if sta.StaDepth<=cat.sort_values(by='StaDepth')[:15].copy().StaDepth.max():dep=f'Shallow'
                    else:dep=f'Deep'
                    file = f'{event.Name}.Mw{event.magnitudes[0].mag}_{method}_{dep}_{stanm}.{int(sta.StaDepth)}m.png'
                    if len(note)>0:file=file.replace('.png',f'_{note}.png')
                    if octave_av:file=file.replace('.png',f'_octave.png')
                    fold = dirs.Ch1/'Spectrograms'
                    note = file.split('_')[0]
                    if octave_av:file = fold/'OctaveAveraged'/magfold/dep/note/file
                    else:file = fold/'_SpecFigsRevised_420_'/'ComprehensiveSet'/subfold/magfold/note/file
                    file.parent.mkdir(parents=True,exist_ok=True)
                    title=file.name.replace('_psd.and.psd_residual.png','').split('_')[0]+'\n'+\
                    '_'.join(file.name.replace('_psd.and.psd_residual.png','').split('_')[1:])
                    fig.suptitle(title,fontweight='bold',y=1.2,fontsize=5,horizontalalignment='left',x=0)
                    save_tight(file,fig,dpi=700,margins=False)
                    # print('Output:')
                    # print(f'Folder: {file.parent}')
                    # print(f'File: {file.name}')
                plt.close('all')