def coherence_plot(coh,sta,method,**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[1,4]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0]
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pair_colors = {'ZP':'#0c51a6','ZZ':'#2a7e93','Z1':'#7370cb','Z2':'#4f86c5',
    'ZZ':'#0c51a6','11':'#2a7e93','22':'#7370cb','PP':'#4f86c5'}
    # ------------
    stanm = sta.StaName
    stastr = ' | '.join([method.upper(),stanm,sta.Experiment,'Depth: '+str(int(sta.StaDepth))+'m'])
    stastr= f'{sta.Experiment} ({sta.Network}) {sta.Station}, {int(sta.StaDepth)}m'
    nrows = len(args.csd_pairs);ncols=1
    figsize=(args.figwidth*args.figaspect[0],nrows*args.figaspect[1]) #20,16
    figsize=(10,14)
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,layout='constrained',sharex='all',squeeze=True,figsize=figsize)
    axes = axes.reshape(-1)
    context_title = f'Before/After Coherence when {method} applied to station event channel'
    fig.suptitle(f'{context_title}\n{stastr}',fontweight='bold',fontsize=10)
    fn = fnotch(sta.StaDepth)
    x=coh.f
    colors=args.csd_pair_colors
    fticks = np.array([1/500,1/300,1/100,1/50,1/10,1])
    xlim=[1/500,1]
    for pi,(pair,ax) in enumerate(zip(args.csd_pairs,axes)):
        ax.set_xlim(xlim);ax.set_xscale('log');ax.set_ylim(0,1.01)
        y_coh = coh[pair].coh
            # ------------------------------------------------------------------------------------------
        chan = pair[0].replace('P','HDH').replace('Z','HZ').replace('1','H1').replace('2','H2')
        y_avg = np.mean(y_coh,axis=0)
        y_var = np.std(y_coh,axis=0)**1
        ax.scatter(x,y_avg,s=0.5,color=colors[pair])
        ax.plot(x,y_avg,linewidth=0.8,alpha=0.8,color=colors[pair])
        [ax.scatter(x,y,s=1,color=colors[pair],alpha=0.4) for y in y_coh]
        ax.fill_between(x,y_avg-y_var,y_avg+y_var, alpha=0.2,color=colors[pair])
        ax.set_xlabel('Period, seconds ')
        ax.axvline(fn,alpha=0.4,linewidth=1,color='k',linestyle='-.')
        ax.text(fn+.001,0.01,str(np.round(1/fn,1))+'s',verticalalignment='bottom',horizontalalignment='left',
        fontweight='bold',fontsize=13)
        # ax.set_title(context_title,fontweight='bold',fontsize=10,horizontalalignment='right')
        ax.text(xlim[0]+0.00006,1,f'Event {chan}',verticalalignment='top',horizontalalignment='left',
        fontweight='bold',fontsize=13,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1))
        ax.set_xticks(fticks)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xticklabels(np.array(1/fticks,dtype=int))
        ax.set_ylabel(f'Coherence',fontweight='bold')
    return fig


# report=get_reports('ZZ',catalog,dirs.Archive,AVG=False)
Comps=['ZZ','11','22','PP']
# fig,axes = plt.subplots(nrows=4,ncols=1,figsize=(10,13));#axes=np.atleast_2d(axes).T
stanm = catalog.StaName[0]
# for stanm in catalog.StaName:
for si,stanm in enumerate(catalog.StaName):
    clear_output()
    print(f'{stanm} {si+1}/{len(catalog)}')
    sta = catalog.loc[stanm]
    coh=AttribDict()
    for comp in Comps:
        file=f'/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences/hps/complete_hps.{comp}_coh.report.pkl'
        d=unravel_report(file,AVG=False)[stanm]
        coh[comp]=d
        coh.f=d.f
    # NOISECUT HORIZONTALS COHERENCES ZZ
    fig=coherence_plot(coh,sta,'Noisecut',csd_pairs=Comps)
    save_tight(dirs.Plots/'HPS_Horizontals_SingleEvent_Pages'/f'{stanm}.png',fig,dpi=800)