from imports import *
from scipy.stats import iqr
from local_tools.math import cohstats
from scipy.ndimage import gaussian_filter
cat = catalog.copy()
octavg=lt.math.octave_average
# Figure priorities:
# 1. Transfer Functions for a deep and shallow station
# 2. Map
# 3. Single station ZZ coherence for a deep and shallow
# 4. 8x spectrograms deep/shallow x atacr/hps x pre/post


# ----------------------------------------------------------------------------------------


# -----TF-data
get_tf = lambda stanm:load_pickle(list((dirs.TransferFunctions/stanm).glob('*-*.pkl'))[0])
get_day_tfs = lambda stanm: [load_pickle(i) for i in [f for f in list((dirs.TransferFunctions/stanm).glob('*.pkl')) if str(f).find('-')<0]]
# get_day_tfs = lambda stanm,days: [load_pickle(dirs.TransferFunctions/stanm/i) for i in days]
# ----------------------------------------------------------------------------------------

minmag=6.8
maxmag=7.0
magfold = f'{int(minmag*10)}_{int(maxmag*10)}'


cat=catalog.copy()
cat = pd.concat([cat.sort_values(by='StaDepth')[:15].copy(),cat.sort_values(by='StaDepth')[-15:].copy()])
flim = [1/230,1]
# ----------------------------------------------------------------------------------------
outlierprops={'color':'r','s':10,'alpha':0.09}
inlierprops={'color':'dodgerblue','s':13,'alpha':0.2} #'#7370cb' royalblue dodgerblue
whiskerprops={'linewidth':0.5,'color':'k'}
whiskerwidth=0.002;margin=1.02
midlineprops={'linewidth':1.5,'color':'k','alpha':0.8}
loglogprops={'linewidth':1.5,'color':'k','alpha':0.8,'nonpositive':'mask'}
midscatterprops={'s':10,'color':'k'} #args.csd_pairs[pair]
notchprops = {'alpha':0.4,'linewidth':2,'color':'k','linestyle':'-.'}
# labelprops = {'fontweight':'bold','fontsize':12}
labelprops = {'fontsize':12}
methodlabelprops = {'fontweight':'bold','fontsize':15,'verticalalignment':'bottom','horizontalalignment':'right'}
# ----------------------------------------------------------------------------------------


status = lambda:print(f'S: {si+1}/{len(cat)} :: {evi+1}/{len(sta.Events)} :: {method}')
figsize=(30,7)
for si,sta in enumerate(cat.iloc):
    stanm=sta.StaName
    for evi,event in enumerate(sta.Events):
        if event.magnitudes[0].mag<minmag:continue
        if event.magnitudes[0].mag>maxmag:continue
        for method in ['NoiseCut','ATaCR']:
            status()
            if method=='NoiseCut':fig,axes=plt.subplots(1,2,figsize=figsize,width_ratios=[1,1])
            if method=='ATaCR':
                fig,axes=plt.subplots(1,3,figsize=figsize,width_ratios=[.25,1,1])
                legendcols=2
                keys=['TF_ZP-21']

                ax=axes[0]
                tfs=get_tf(stanm);f=tfs.f;ind=(f>0)&(f<=1);f=f[ind]
                # tfs.good_days
                days=get_day_tfs(stanm)
                days=np.array([abs(d.transfunc['ZP-21'][keys[0]]) for d in days])
                days=days[:,ind]

                # avgkeys = list(tfs.transfunc['ZP-21'].keys())
                avgkeys=keys
                avg = {k:np.abs(tfs.transfunc['ZP-21'][k][ind]) for k in avgkeys}
                xlim = [np.min(avg[avgkeys[0]][avg[avgkeys[0]]>0]),10*np.max(avg[avgkeys[0]][avg[avgkeys[0]]>0])]

                upper,lower,median,outliers,inliers = cohstats(days)
                # y_mid = median
                # y_mid = abs(tfs.transfunc['ZP-21'][keys[0]][ind])
                # [ax.scatter(f[i],d[i],**{'color':'r','s':.06,'alpha':0.09}) for i,d in zip(outliers,days)]
                # ax.scatter(f,y_mid,**midscatterprops)
                # ax.plot(f,y_mid,**midlineprops) #args.csd_pairs[pair]
                [ax.scatter(d[i],f[i],**inlierprops) for i,d in zip(inliers,days)]
                [ax.loglog(avg[k],f,label=k.replace('TF_',''),**loglogprops) for k in list(avg.keys())]
                ax.set_xscale('linear')
                # ax.legend(markerscale=10,ncols=legendcols,loc='lower right',columnspacing=0,prop={'weight': 'bold'})
                ax.axhline(fnotch(sta.StaDepth),**notchprops)
                ax.set_ylabel('Frequency, Hz',**labelprops)
                # ax.set_title(f'"Figure: {stanm}"',fontstyle='italic')
                ax.set_title(f'Station averaged transfer functions',**labelprops)
                ax.set_ylim(flim[0],flim[1])
                ax.set_xlabel('dB / Pa',**labelprops)
                # ----------
                # ax.set_yscale('log')
                # ax.set_xscale('linear')
                # ax.set_xlim(xlim[0],xlim[1])

            try:
                traces = get_traces(stanm,event.Name)
                specs = {tr.stats.location:spectrogram(tr.copy()) for tr in traces}
            except:
                print(f'{event.Name} : {stanm} : Data?')
                continue

            if method=='ATaCR':ax=axes[1];k='Raw'
            else:ax=axes[0];k='Raw'

            S,spec_f,spec_t=specs[k]
            trace = traces[0].copy()
            vlim=[np.min([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()]),np.max([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()])]
            specfig=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=vlim,cbar=False)
            ax.set_ylim(flim[0],flim[1])
            ax.axhline(fnotch(sta.StaDepth),**notchprops)
            ax.set_title(f'"Figure: {k} | {event.Name} Mw{event.magnitudes[0].mag}"')
            ax.set_yticklabels([])
            ax.set_xlabel('Time after origin, seconds')
            ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800))

            if method=='NoiseCut':ax.set_ylabel('Frequency, Hz',**labelprops)

            if method=='ATaCR':
                ax=axes[2];k='ATaCR'
                S,spec_f,spec_t=specs[k]
                trace = traces[0].copy()
                vlim=[np.min([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()]),np.max([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()])]
                specfig=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=vlim,cbar=True)
                ax.set_ylim(flim[0],flim[1])
                ax.axhline(fnotch(sta.StaDepth),**notchprops)
                ax.set_title(f'"Figure: {k} | {event.Name} Mw{event.magnitudes[0].mag}"')
                ax.set_yticklabels([])
                ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800))
                ax.set_xlabel('Time after origin, seconds')
            else:
                ax=axes[1];k='NoiseCut'
                S,spec_f,spec_t=specs[k]
                trace = traces[0].copy()
                vlim=[np.min([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()]),np.max([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()])]
                specfig=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=vlim,cbar=True)
                ax.set_ylim(flim[0],flim[1])
                ax.axhline(fnotch(sta.StaDepth),**notchprops)
                ax.set_title(f'"Figure: {k} | {event.Name} Mw{event.magnitudes[0].mag}"')
                ax.set_yticklabels([])
                ax.set_xlabel('Time after origin, seconds')
                ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800))
                fig.tight_layout(pad=0.01)

            fig.subplots_adjust(wspace=0.09, hspace=0.1)
            if sta.StaDepth<=cat.sort_values(by='StaDepth')[:15].copy().StaDepth.max():dep=f'Shallow'
            else:dep=f'Deep'
            file = f'{event.Name}.Mw{event.magnitudes[0].mag}_{method}_{dep}_{stanm}.{int(sta.StaDepth)}m.png'
            fold = dirs.Ch1/'Spectrograms'
            file = fold/magfold/file
            save_tight(file,fig,dpi=700)
            plt.close('all')