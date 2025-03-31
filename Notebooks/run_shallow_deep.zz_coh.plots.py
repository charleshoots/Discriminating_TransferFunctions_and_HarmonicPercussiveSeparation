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



# ---------------------------------------------------------------------------------------------------------
# -----TF-data
get_tf = lambda stanm:load_pickle(list((dirs.TransferFunctions/stanm).glob('*-*.pkl'))[0])
get_day_tfs = lambda stanm: [load_pickle(i) for i in [f for f in list((dirs.TransferFunctions/stanm).glob('*.pkl')) if str(f).find('-')<0]]
# -----ZZ Event Coherence data
report=get_reports('ZZ',catalog,dirs.Archive,dirs,AVG=True)

# ---------------------------------------------------------------------------------------------------------
cat=catalog.copy()
cat = pd.concat([cat.sort_values(by='StaDepth')[:15].copy(),cat.sort_values(by='StaDepth')[-15:].copy()])
flim = [1/230,1]
# ----------------------------------------------------------------------------------------
# outlierprops={'color':'r','s':10,'alpha':0.09}
# inlierprops={'color':'dodgerblue','s':13,'alpha':0.2} #'#7370cb' royalblue dodgerblue
# whiskerprops={'linewidth':0.5,'color':'k'}
# whiskerwidth=0.002;margin=1.02
# midlineprops={'linewidth':1.5,'color':'k','alpha':0.8}
# loglogprops={'linewidth':1.5,'color':'k','alpha':0.8,'nonpositive':'mask'}
# midscatterprops={'s':10,'color':'k'} #args.csd_pairs[pair]
# notchprops = {'alpha':0.4,'linewidth':2,'color':'k','linestyle':'-.'}
# # labelprops = {'fontweight':'bold','fontsize':12}
# labelprops = {'fontsize':12}
# methodlabelprops = {'fontweight':'bold','fontsize':15,'verticalalignment':'bottom','horizontalalignment':'right'}
# ----------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------



def sta_metrics(report,sta,flim=[1/500,1],legendcols=2,figsize=(8,7),octave_av=False):
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
    midlineprops={'linewidth':1.5,'color':'k','alpha':0.8}
    midscatterprops={'s':10,'color':'k'} #args.csd_pairs[pair]
    notchprops = {'alpha':0.4,'linewidth':1,'color':'k','linestyle':'-.'}
    labelprops = {'fontweight':'bold','fontsize':12}
    methodlabelprops = {'fontweight':'bold','fontsize':15,'verticalalignment':'bottom','horizontalalignment':'right'}

    stastr=f'"Figure: ZZ Coherence for {stanm} | {sta.Experiment} | {str(int(abs(sta.StaDepth)))}m"'
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # columns=[['NoiseCut',['Coherence','ZZ']],['ATaCR',['Coherence','ZZ']]]
    # fig,axes = plt.subplots(nrows=3,ncols=1,layout='constrained',sharex='all',squeeze=True,figsize=figsize,height_ratios=[.4,.8,.8])


    fig,axes = plt.subplots(nrows=2,ncols=1,layout='constrained',sharex='all',squeeze=True,figsize=figsize,height_ratios=[1,1])


    axes = axes.reshape(-1)

    keys=['TF_ZP-21']
    # ax=axes[0]
    # tfs=get_tf(stanm);f=tfs.f;ind=(f>=1/500)&(f<1);f=f[ind]
    # # avgkeys = list(tfs.transfunc['ZP-21'].keys())
    # avgkeys=keys
    # avg = {k:np.abs(tfs.transfunc['ZP-21'][k][ind]) for k in avgkeys}
    # days=get_day_tfs(stanm)
    # days=np.array([abs(d.transfunc['ZP-21'][keys[0]]) for d in days])
    # days=days[:,ind]

    # upper,lower,median,outliers,inliers = cohstats(days)
    # y_mid = median
    # y_mid = abs(tfs.transfunc['ZP-21'][keys[0]][ind])
    # [ax.scatter(f[i],d[i],**{'color':'r','s':.06,'alpha':0.09}) for i,d in zip(outliers,days)]
    # ax.scatter(f,y_mid,**midscatterprops)
    # ax.plot(f,y_mid,**midlineprops) #args.csd_pairs[pair]

    # [ax.scatter(f[i],d[i],**inlierprops) for i,d in zip(inliers,days)] ####-------disable
    # [ax.plot(f,avg[k],label=k.replace('TF_',''),**midlineprops) for k in list(avg.keys())]

    # ax.set_yscale('log')
    # ax.legend(markerscale=10,ncols=legendcols,loc='lower right',columnspacing=0,prop={'weight': 'bold'})
    # ax.axvline(fnotch(sta.StaDepth),**notchprops)
    # ax.set_ylabel('dB/Pa',**labelprops)
    # # ax.set_title(f'"Figure: {stanm}"',fontstyle='italic')
    # ax.set_title(f'Station averaged transfer functions',**labelprops)
    # # ----------
    # ax.set_xscale('log')
    # ax.set_xlim(flim[0],flim[1])
    # ax.set_ylim(top=10**np.log10(max(ax.get_ylim())))
    # ax.set_ylim(np.mean(avg[avgkeys[0]])/10,np.max(avg[avgkeys[0]]))
    # ax.set_yticks([10**(-9),10**(-20)])

    # fig.suptitle(stastr,fontstyle='italic')
    fn = 1/fnotch(sta.StaDepth)
    methods=['ATaCR','NoiseCut']
    # methods=['NoiseCut','ATaCR']
    # axes=axes[1:]
    for bi,b in enumerate(columns):
        method = methods[bi]
        metric = b[1][0]
        pair = b[1][1]
        ax = axes[bi]
        f=report.f
        # for tri in range(len([1])):
        ax.set_title(f'{method} {pair} {metric}',**labelprops)
        ind=(f>=1/500)&(f<1)
        f=f[ind]
        coh=report[method][stanm].coh[:,ind]
        coh = coh[~np.any(np.isnan(coh),axis=1)]
        

        if octave_av:
            avged=np.array([octavg(c,f)[1] for c in coh])
            f = octavg(coh[0],f)[0]
            coh=avged

        upper,lower,median,outliers,inliers=cohstats(coh,margin=margin)
        # [ax.scatter(x,yy,label=':'.join([pair]),s=5,facecolor=args.csd_pairs[pair],alpha=0.01) for yy in coh]
        ax.set_xlim(f[0],1)
        if metric.lower()=='phase':ax.set_ylim(0,180)
        if metric.lower()=='coherence':ax.set_ylim(0,1.02)
        # y_mid = np.mean(y,axis=0)
        # y_mid = gaussian_filter(median, sigma=1)
        y_mid = median
        ax.scatter(f,y_mid,label=':'.join([pair]),**midscatterprops)
        ax.vlines(f,ymin=y_mid,ymax=upper,**{'linewidth':0.3,'color':'k'})
        ax.vlines(f,ymin=lower,ymax=y_mid,**{'linewidth':0.3,'color':'k'})
        ax.hlines(upper,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)
        ax.hlines(lower,f-whiskerwidth/2,f+whiskerwidth/2,**whiskerprops)

        #Comment out to remove dots
        # [ax.scatter(f[i],d[i],**inlierprops) for i,d in zip(inliers,coh)]
        # [ax.scatter(f[i],d[i],**outlierprops) for i,d in zip(outliers,coh)]

        ax.plot(f,y_mid,label=':'.join([pair]),**midlineprops) #args.csd_pairs[pair]
        # ax.fill_between(x,lower,upper, alpha=0.1,color='gray')
        ax.set_xlabel('Frequency, Hz',**labelprops)
        ax.set_ylabel('Coherence',**labelprops)
        ax.axvline(1/fn,**notchprops)
        ax.set_yticklabels(ax.get_yticklabels(),**labelprops)
        ax.set_xticklabels(ax.get_xticklabels(),**labelprops)
        # ax.text(max(ax.get_xlim())-0.05,0,method,**methodlabelprops)
        if bi==0:ax.set_xticklabels([]);ax.set_xlabel('')
    _ = [ax.set_xscale('log') for ax in axes]
    # ax.text(0.001+1/fn,0,str(int(np.round(fn)))+'s',verticalalignment='bottom',horizontalalignment='left',fontweight='bold',fontsize=11)
    return fig
sta=cat.loc['7D.FS42D']
for si,sta in enumerate(cat.iloc):
    print(f'{si+1}/{len(cat)} : {sta.StaName}')
    for octave_av in [False,True]:
        fig=sta_metrics(report,sta,legendcols=1,figsize=(10,8),octave_av=octave_av)
        file = f'{sta.StaName}.{sta.StaDepth}m.{'Octave' if octave_av else 'Normal'}.png'
        fold = dirs.Ch1/'SingleStationZZCoherences'

        file = fold/file
        # file = fold/'_QC'/file

        save_tight(file,fig,600)
        plt.close('all')