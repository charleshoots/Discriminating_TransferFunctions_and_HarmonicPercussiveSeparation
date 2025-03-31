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





def ax_sta_metrics(ax,report,sta,method,event_name=None,flim=[1/500,1],octave_av=False,y_faxis=True):
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = {'Z1':'#0c51a6','ZP':'#2a7e93','ZZ':'#7370cb','Z2':'#4f86c5'}
    note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    stanm = sta.StaName
    tf=''
    # stastr = ' | '.join([stanm+tf,sta.Experiment,'Depth: '+str(int(abs(sta.StaDepth)))+'m, Notch: '+str(int(1/fnotch(1000*abs(sta.StaDepth))))+'s',note])

    outlierprops={'color':'r','s':10,'alpha':0.09}
    inlierprops={'color':'dodgerblue','s':13,'alpha':0.2} #'#7370cb' royalblue dodgerblue
    whiskerprops={'linewidth':0.5,'color':'k'}
    whiskerwidth=0.002;margin=1.02
    midlineprops={'linewidth':1.0,'color':'k','alpha':0.8}
    eventcohlnie={'linewidth':1.0,'color':'r','alpha':0.8}
    midscatterprops={'s':1,'color':'k'} #args.csd_pairs[pair]
    notchprops = {'alpha':0.4,'linewidth':1,'color':'k','linestyle':'-.'}
    labelprops = {'fontsize':5} #'fontweight':'bold'
    methodlabelprops = {'fontweight':'bold','fontsize':5,'verticalalignment':'bottom','horizontalalignment':'right'}

    stastr=f'"Figure: ZZ Coherence for {stanm} | {sta.Experiment} | {str(int(abs(sta.StaDepth)))}m"'
    # columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # columns=[['NoiseCut',['Coherence','ZZ']],['ATaCR',['Coherence','ZZ']]]
    # fig,axes = plt.subplots(nrows=3,ncols=1,layout='constrained',sharex='all',squeeze=True,figsize=figsize,height_ratios=[.4,.8,.8])

    keys=['TF_ZP-21']
    fn = 1/fnotch(sta.StaDepth)
    methods=[method]
    for bi,b in enumerate(methods):
        method = methods[bi]
        metric = 'Coherence'
        pair = 'ZZ'
        f=report.f
        ax.set_title(f'i. {pair} {metric}',**labelprops)
        ind=(f>0)&(f<=1)
        f=f[ind]
        coh=report[method][stanm].coh[:,ind]
        events=report[method][stanm].events[~np.any(np.isnan(coh),axis=1)]
        coh = coh[~np.any(np.isnan(coh),axis=1)]
        if event_name:event_coh=coh[np.where(events==event_name)[0],:].reshape(-1)
        if octave_av:
            avged=np.array([octavg(c,f)[1] for c in coh])
            coh=avged
            event_coh = octavg(event_coh,f)
            f = octavg(coh[0],f)[0]

        upper,lower,median,outliers,inliers=cohstats(coh,margin=margin)

        if y_faxis:
            # [ax.scatter(x,yy,label=':'.join([pair]),s=5,facecolor=args.csd_pairs[pair],alpha=0.01) for yy in coh]
            ax.set_xlim(f[0],1)
            if metric.lower()=='phase':ax.set_ylim(0,180)
            if metric.lower()=='coherence':ax.set_ylim(0,1.02)
            y_mid = coh.mean(axis=0)
            # y_mid = gaussian_filter(median, sigma=1)
            # y_mid = median
            ax.scatter(y_mid,f,label=':'.join([pair]),**midscatterprops)
            ax.hlines(f,xmin=y_mid,xmax=upper,**{'linewidth':0.3,'color':'k'})
            ax.hlines(f,xmin=lower,xmax=y_mid,**{'linewidth':0.3,'color':'k'})
            ax.vlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.vlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.plot(y_mid,f,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]

            ax.plot(event_coh,f,**eventcohlnie)

            ax.set_ylabel('Frequency, Hz',**labelprops)
            ax.set_xlabel('Coherence',**labelprops)
            ax.axhline(1/fn,**notchprops)
            # ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            # if bi==0:ax.set_yticklabels([]);ax.set_ylabel('')
            ax.set_ylim(flim[0],flim[1])
            ax.set_xlim(0,1)
            _ = ax.set_yscale('log')
        else:
            ax.set_xlim(f[0],1)
            if metric.lower()=='phase':ax.set_ylim(0,180)
            if metric.lower()=='coherence':ax.set_ylim(0,1.02)
            y_mid = coh.mean(axis=0)
            # y_mid = gaussian_filter(median, sigma=1)
            # y_mid = median
            ax.scatter(f,y_mid,label=':'.join([pair]),**midscatterprops)
            ax.vlines(f,ymin=y_mid,ymax=upper,**{'linewidth':0.3,'color':'k'})
            ax.vlines(f,ymin=lower,ymax=y_mid,**{'linewidth':0.3,'color':'k'})
            ax.hlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.hlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
            ax.plot(f,y_mid,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]
            ax.set_xlabel('Frequency, Hz',**labelprops)
            ax.set_ylabel('Coherence',**labelprops)
            ax.axvline(1/fn,**notchprops)

            _ = [ax.set_xscale('log') for ax in axes]
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            # ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
    # ax.text(0.001+1/fn,0,str(int(np.round(fn)))+'s',verticalalignment='bottom',horizontalalignment='left',fontweight='bold',fontsize=11)
    return ax
# ----------------------------------------------------------------------------------------




# -----TF-data
report=get_reports('ZZ',catalog,dirs.Archive,dirs,AVG=True)
get_tf = lambda stanm:load_pickle(list((dirs.TransferFunctions/stanm).glob('*-*.pkl'))[0])
get_day_tfs = lambda stanm: [load_pickle(i) for i in [f for f in list((dirs.TransferFunctions/stanm).glob('*.pkl')) if str(f).find('-')<0]]
# get_day_tfs = lambda stanm,days: [load_pickle(dirs.TransferFunctions/stanm/i) for i in days]
# ----------------------------------------------------------------------------------------

minmag=6.8
maxmag=7.0

# minmag=7.4
# maxmag=8.0


magfold = f'{int(minmag*10)}_{int(maxmag*10)}'


cat=catalog.copy()
# cat = pd.concat([cat.sort_values(by='StaDepth')[:15].copy(),cat.sort_values(by='StaDepth')[-15:].copy()])
cat=cat.sort_values(by='StaDepth')[-15:].copy()

flim = [1/230,1]
# ----------------------------------------------------------------------------------------
outlierprops={'color':'r','s':10,'alpha':0.09}
inlierprops={'color':'dodgerblue','s':13,'alpha':0.2} #'#7370cb' royalblue dodgerblue
whiskerprops={'linewidth':0.5,'color':'k'}
whiskerwidth=0.002;margin=1.02
midlineprops={'linewidth':1.5,'color':'k','alpha':0.8}
loglogprops={'linewidth':1.5,'color':'k','alpha':0.8,'nonpositive':'mask'}
midscatterprops={'s':10,'color':'k'} #args.csd_pairs[pair]
notchprops = {'alpha':0.4,'linewidth':1,'color':'k','linestyle':'-.'}
# labelprops = {'fontweight':'bold','fontsize':12}
labelprops = {'fontsize':5}
methodlabelprops = {'fontweight':'bold','fontsize':5,'verticalalignment':'bottom','horizontalalignment':'right'}
# ----------------------------------------------------------------------------------------

# ind=(f>=1/500)&(f<1)
# f=f[ind]
# coh=report[method][stanm].coh[:,ind]

# ax_sta_metrics(ax,report,sta,method,event_name=None,flim=[1/500,1],octave_av=False,y_faxis=True)

status = lambda:print(f'S: {si+1}/{len(cat)} :: {evi+1}/{len(sta.Events)} :: {method}')
figsize=(8,4)
octave_av=False
y_faxis=True
for si,sta in enumerate(cat.iloc):
    stanm=sta.StaName
    for evi,event in enumerate(sta.Events):
        if event.magnitudes[0].mag<minmag:continue
        if event.magnitudes[0].mag>maxmag:continue
        for method in ['NoiseCut','ATaCR']:
            status()
            if not sta.Network=='Z6':continue
                # if method=='ATaCR':
 
            try:
                traces = get_traces(stanm,event.Name)
                specs = {tr.stats.location:spectrogram(tr.copy()) for tr in traces}
                raw=traces.select(location='Raw').copy()
                for k in ['ATaCR','NoiseCut']:
                    temp= list(specs['Raw']);temp[0]=specs[k][0] - temp[0]
                    # temp= list(specs['Raw']);temp[0]=temp[0]-specs[k][0]
                    specs[k+'_Residual'] = temp
                vlim=[np.min([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()]),np.max([[specs[k][0].min(),specs[k][0].max()] for k in specs.keys()])]
                vlim=[-100,-10]
            except:
                print(f'{event.Name} : {stanm} : Data?')
                continue

            fig = plt.figure(figsize=figsize)
            gs = gridspec.GridSpec(1, 4, width_ratios=[0.25, 1, 1, 0.05], wspace=0.09)
            # Create subplots using the gridspec
            ax1 = fig.add_subplot(gs[0, 0])  # First subplot
            ax2 = fig.add_subplot(gs[0, 1])  # Second subplot
            ax3 = fig.add_subplot(gs[0, 2])  # Third subplot
            cbar_ax = fig.add_subplot(gs[0, 3])  # Colorbar axis
            axes = [ax1,ax2,ax3,cbar_ax]

            ax=axes[0]
            ev_name = UTCDateTime.strptime(event.Name,'%Y.%j.%H.%M').strftime('%Y-%m-%d %H:%M:00')
            ax_sta_metrics(ax,report,sta,method,event_name=event.Name,flim=flim,octave_av=octave_av,y_faxis=y_faxis)
            ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
            ax.set_yticklabels(ax.get_yticklabels(),**labelprops)

            if method=='ATaCR':ax=axes[1];k='ATaCR'
            else:ax=axes[1];k='NoiseCut'
            S,spec_f,spec_t=specs[k]
            trace = traces[0].copy()
            specfig,pcm=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=vlim,cbar=False)
            cbar = fig.colorbar(pcm, cax=axes[-1],pad= 0.01)
            cbar.set_label('dB',**labelprops)
            cbar.set_ticks([vlim[0],int(vlim[0]/2),vlim[-1]])
            cbar.set_ticklabels([vlim[0],int(vlim[0]/2),vlim[-1]],**labelprops)

            ax.set_ylim(flim[0],flim[1])
            ax.axhline(fnotch(sta.StaDepth),**notchprops)
            # ax.set_title(f'ii. {k} | {ev_name} Mw{event.magnitudes[0].mag}')
            ax.set_title(f'ii. {k}, {ev_name} Mw{event.magnitudes[0].mag} observed at {f'{sta.StaName} ({int(sta.StaDepth)}m)'}',**labelprops)
            ax.set_yticklabels([])
            ax.set_xlabel('Time after origin, hours',**labelprops)
            ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800)/3600)
            ax.set_xticklabels(ax.get_xticklabels(),**labelprops)

            # if method=='NoiseCut':ax.set_ylabel('Frequency, Hz',**labelprops)
            ax=axes[2]
            k=f'{k}_Residual'
            S,spec_f,spec_t=specs[k]
            trace = traces[0].copy()

            specfig,pcm=plot_spectrogram(S,spec_f,spec_t,ax=ax,vlim=[None,None],cbar=False)
            cbar = fig.colorbar(pcm, cax=axes[-1],pad= 0.01)
            cbar.set_label('dB',**labelprops)
            # cbar.set_ticks([vlim[0],int(vlim[0]/2),vlim[-1]])
            # cbar.set_ticklabels([vlim[0],int(vlim[0]/2),vlim[-1]],**labelprops)
            ax.set_ylim(flim[0],flim[1])
            ax.axhline(fnotch(sta.StaDepth),**notchprops)
            ax.set_title(f'iii.  Residual',**labelprops)
            ax.set_yticklabels([])
            ax.set_xlabel('Time after origin, hours',**labelprops)
            ax.set_xticks(np.arange(0,int(spec_t[-1])+1800,1800)/3600)
            ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
            # ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
            fig.tight_layout(pad=0.01)

            fig.subplots_adjust(wspace=0.09, hspace=0.1)
            if sta.StaDepth<=cat.sort_values(by='StaDepth')[:15].copy().StaDepth.max():dep=f'Shallow'
            else:dep=f'Deep'
            file = f'{event.Name}.Mw{event.magnitudes[0].mag}_{method}_{dep}_{stanm}.{int(sta.StaDepth)}m.png'
            fold = dirs.Ch1/'Spectrograms'

            file = fold/(magfold)/file
            save_tight(file,fig,dpi=700)
            plt.close('all')