from modules import *
from local_tools.math import avg_meter
import local_tools as lt
def get_gridplot():
    fig = plt.figure(figsize=(20,15),layout="constrained")
    height_ratios = [1,1,0.6,0.6,0.6,0.6]
    gs = gridspec(6, 6, figure=fig,height_ratios=height_ratios)
    ax_23_0=fig.add_subplot(gs[0,0:2]);ax_23_1=fig.add_subplot(gs[0,2:4]);ax_23_2=fig.add_subplot(gs[0,4:])
    ax_23_3=fig.add_subplot(gs[1,0:2]);ax_23_4=fig.add_subplot(gs[1,2:4]);ax_23_5=fig.add_subplot(gs[1,4:])
    ax23 = np.array([[ax_23_0,ax_23_1,ax_23_2],[ax_23_3,ax_23_4,ax_23_5]])
    ax_42_0=fig.add_subplot(gs[2,:3]);ax_42_1=fig.add_subplot(gs[2,3:])
    ax_42_2=fig.add_subplot(gs[3,:3]);ax_42_3=fig.add_subplot(gs[3,3:])
    ax_42_4=fig.add_subplot(gs[4,:3]);ax_42_5=fig.add_subplot(gs[4,3:])
    ax_42_6=fig.add_subplot(gs[5,:3]);ax_42_7=fig.add_subplot(gs[5,3:])
    ax42 = np.array([[ax_42_0,ax_42_1],[ax_42_2,ax_42_3],[ax_42_4,ax_42_5],[ax_42_6,ax_42_7]])
    return fig,ax23,ax42

def station_event_page(st_hold,sta,evmeta,method,type='stream',**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[1,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')];defargs.Noise=True
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    if args.csd_pairs[0][0][0]==args.csd_pairs[0][0][1]:args.Noise=False
    # ------------
    # ------------
    if type.lower()=='stream':
        note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:
        note = 'Noise (gray) | Raw (red)'
    staname = st_hold[0].stats.network+'.'+st_hold[0].stats.station
    stastr = ' | '.join([method.upper(),staname,sta.Experiment,
    'Depth: '+str(int(1000*abs(st_hold[0].stats.sac.stel)))+'m',
    'F-Notch: '+str(int(1/fnotch(1000*abs(st_hold[0].stats.sac.stel))))+'s',
    note])
    if not args.nev:args.nev = len(st_hold.select(location='*Raw*'))
    if type.lower()=='stream':columns = args.bands
    elif type.lower()=='metrics':columns=['Coherence','Phase']
    nev_per_plot = args.nev
    nplots = int(np.ceil(nev_per_plot/len(evmeta)))
    for ploti in range(nplots):
        fig,axes = plt.subplots(nrows=nev_per_plot,ncols=len(columns),layout='constrained',sharex='all',squeeze=True,figsize=(args.figwidth*args.figaspect[0],nev_per_plot*args.figaspect[1]))
        axes = np.atleast_2d(axes).reshape((nev_per_plot,len(columns)))
        fig.suptitle(stastr)
        fn = 1/fnotch(1000*abs(st_hold[0].stats.sac.stel))
        for bi,b in enumerate(columns):
            # print(method+'-'+type + '| Column:'+str(bi+1)+'/'+str(len(columns)))
            band_ax = axes[:,bi]
            st_band = st_hold.copy()
            correct_hold = st_band.select(location='*Correct*').copy()
            raw_hold = st_band.select(location='*Raw*').copy()
            if type.lower()=='stream':
                st_band.filter('bandpass',freqmin=1/b[1],freqmax=1/b[0],zerophase=True,corners=4)
                correct_hold = st_band.select(location='*Correct*').copy()
                raw_hold = st_band.select(location='*Raw*').copy()
                # --------------
                # This handles large scale amplitude differences between the two traces.
                # These are traces where the raw amplitudes exceed the vertical_scale threshold wrt the corrected trace amplitudes.
                ylim = np.array([np.max(np.abs(c.data))*args.vertical_scale for c in correct_hold])
                out_scaled = np.array([np.max(np.abs(c.data)) for c in raw_hold]) > ylim
                if np.any(out_scaled):
                    print('Large amplitude scale differences detected in '+staname)
                    for tr_ind,tr in enumerate(raw_hold):
                        if out_scaled[tr_ind]:tr.data = (tr.data/np.max(np.abs(tr.data)))*ylim[tr_ind]
                # --------------
            for si,s in enumerate(['Raw','Corrected']):
                st = {'Raw':raw_hold,'Corrected':correct_hold}[s]
                for tri in range(nev_per_plot):
                    tr_ind = tri*(ploti+1)
                    tr = st[tr_ind].copy()
                    tr.trim(evmeta[tr_ind].origins[0].time,evmeta[tr_ind].origins[0].time+7200,pad=True,fill_value=0).taper(0.0001)
                    colors = args.phasecolors
                    ax = band_ax[tr_ind]
                    if tr_ind==0:
                        if type.lower()=='stream':ax.set_title(''.join([str(b[0]),'s-',str(b[1]),'s']))
                        else:ax.set_title(b+':'+args.csd_pairs[0][0])
                    if type.lower()=='stream':
                        x = tr.times('relative');y = tr.data;ax.plot(x,y,color=args.linecolor[si],alpha=args.alpha[si],linewidth=args.linewidth[si])
                    elif type.lower()=='metrics':
                        x = tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[0];ind=x<=1
                        y=tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[1][ind]
                        if s=='Raw':
                            if args.Noise:[ax.scatter(xy[0],xy[1],c='darkgrey',s=0.1,label=args.csd_pairs[pi][0]+':Noise') 
                            for pi,xy in enumerate([avg_meter(st_hold.Noise,b,p[0]) for p in args.csd_pairs])]
                            [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.8,color='r') for p in args.csd_pairs]
                            [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linewidth=0.2,alpha=0.4,color='r') for p in args.csd_pairs]
                        if s=='Corrected':
                            [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.4,color=p[1]) for p in args.csd_pairs]
                            [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linestyle=':',linewidth=0.45,alpha=0.4,color=p[1]) for p in args.csd_pairs]
                    if type.lower()=='metrics':
                        ax.set_xlim(1/500,1);ax.set_xscale('log')
                        if b.lower()=='phase':ax.set_ylim(-180,180)
                        if b.lower()=='coherence':ax.set_ylim(0,1.1)
                    else:ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim[tr_ind],ylim[tr_ind]);ax.set_yticklabels('')
                    if s=='Corrected':
                        stallaz=[tr.stats.sac.stla,tr.stats.sac.stlo,tr.stats.sac.stel]
                        evllaz=[evmeta[tr_ind].origins[0].latitude,evmeta[tr_ind].origins[0].longitude,evmeta[tr_ind].origins[0].depth/1000]
                        tr.stats.sac.gcarc = lt.math.distance(sta,evmeta[tr_ind])
                        evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name,'M'+str(evmeta[tr_ind].magnitudes[0].mag),str(np.round(tr.stats.sac.gcarc,2))+'°'])
                        ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='center',fontsize=10)
                        if type.lower()=='stream':
                            if tr.stats.sac.gcarc<=100:phases=args.phases
                            else:phases=args.shadow_phases;[colors.update({p:'k'}) for p in phases]
                            ar = get_arrivals(stallaz,evllaz,model = 'iasp91',phases=phases)
                            [ax.axvline(a[1],linewidth=0.1,color=colors[a[0]]) for a in ar]
                            [ax.text(a[1],y.max(),a[0],fontsize=6,color='k',verticalalignment='bottom',horizontalalignment='center') for a in ar]
            if type.lower()=='metrics':ax.set_xlabel('frequency (hz)')
            else:ax.set_xlabel('seconds')
    return fig
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# ----------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def station_event_page_averages(st_hold,sta,evmeta,type='Metrics',raw_reference=None,**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[4,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')]
    defargs.Noise=True
    defargs.columns = [['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = {'ZP':'#0c51a6','ZZ':'#2a7e93','Z1':'#7370cb','Z2':'#4f86c5'}
    # ------------
    if type.lower()=='stream':note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:note = '\nVariance (shaded)'
    staname = st_hold[0].stats.network+'.'+st_hold[0].stats.station
    # if method.lower()=='atacr':tf='('+st_hold.select(location='*Correct*')[0].stats.location.split('.')[1]+')'
    tf=''
    stastr = ' | '.join([staname+tf,sta.Experiment,'Depth: '+str(int(1000*abs(st_hold[0].stats.sac.stel)))+'m, Notch: '+str(int(1/fnotch(1000*abs(st_hold[0].stats.sac.stel))))+'s',note])
    if not args.nev:args.nev = len(st_hold.select(location='*Raw*'))
    nrows = 1;ncols=2
    columns=[['ATaCR',['Coherence','ZZ']],['NoiseCut',['Coherence','ZZ']]]
    # if type.lower()=='metrics':columns=['Coherence']
    # if type.lower()=='metrics':columns=['Coherence','Phase','Admittance']
    fig,axes = plt.subplots(nrows=ncols,ncols=nrows,layout='constrained',sharex='all',squeeze=True,figsize=(8,7))
    axes = axes.reshape(-1)
    fig.suptitle(stastr)
    fn = 1/fnotch(1000*abs(st_hold[0].stats.sac.stel))
    methods=['ATaCR','NoiseCut']
    for bi,b in enumerate(columns):
        method = methods[bi]
        metric = b[1][0]
        pair = b[1][1]
        # for pi,pair in enumerate(args.csd_pairs):
        channels = pair
        if channels[0]==channels[1]:args.Noise=False
        else: args.Noise=True
        # print(method+'-'+type + '| Column:'+str(bi+1)+'/'+str(len(columns)))
        ax = axes[bi]
        st_band = st_hold.copy()
        correct_hold = st_band.select(location=f'*{method}*').copy()
        raw_hold = st_band.select(location='*Raw*').copy()
        if raw_reference:raw_hold = raw_reference
        for si,s in enumerate(['Corrected']):
            st = {'Raw':raw_hold,'Corrected':correct_hold}[s].copy()
            for tri in range(len(correct_hold)):
                tr_ind = tri
                tr = st[tr_ind].copy()
                if tr_ind==0:
                    if type.lower()=='stream':ax.set_title(''.join([str(metric),'s-',str(pair),'s']))
                    else:ax.set_title(f'{method} {pair} {metric}')
                if type.lower()=='metrics':
                    x = tr.Metrics.__getattribute__(metric)(pair)[0];ind=x<=1
                    y=tr.Metrics.__getattribute__(metric)(pair)[1][ind]
                    if s=='Raw':
                        if args.Noise:[
                        ax.scatter(xy[0],np.abs(xy[1]),c='darkgrey',s=0.1,label=pair+':Noise')
                        for ni,xy in enumerate([avg_meter(st_hold.Noise,metric,pair) for p in [1]])]
                        [ax.scatter(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),s=0.4,color='r',alpha=0.1) for p in [1]]
                        [ax.plot(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),linewidth=0.05,alpha=0.05,color='r') for p in [1]]
                    if s=='Corrected':
                        [ax.scatter(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),s=0.2,color=args.csd_pairs[pair],alpha=0.1) for p in [1]]
                        [ax.plot(x[ind],np.abs(tr.Metrics.__getattribute__(metric)(pair)[1][ind]),label=':'.join([s,pair]),linestyle=':',linewidth=0.05,alpha=0.05,color=args.csd_pairs[pair]) for p in [1]]
                if type.lower()=='metrics':
                    ax.set_xlim(1/500,1);ax.set_xscale('log')
                    # if metric.lower()=='phase':ax.set_ylim(-180,180)
                    if metric.lower()=='phase':ax.set_ylim(0,180)
                    if metric.lower()=='coherence':ax.set_ylim(0,1.01)
                else:ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim[tr_ind],ylim[tr_ind]);ax.set_yticklabels('')
                # ------------------------------------------------------------------------------------------
            y_avg = np.mean(([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)
            y_var = np.std(([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)**2

            # Should I set Phase domain to [0,180] instead of [-180,180]? Would be easier to read...
            # y_avg = np.mean(np.abs([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)
            # y_var = np.std(np.abs([rt.Metrics.__getattribute__(metric)(pair)[1][ind] for rt in st]),axis=0)**2
            if s=='Corrected':
                stallaz=[tr.stats.sac.stla,tr.stats.sac.stlo,tr.stats.sac.stel]
                evllaz=[evmeta[tr_ind].origins[0].latitude,evmeta[tr_ind].origins[0].longitude,evmeta[tr_ind].origins[0].depth/1000]
                tr.stats.sac.gcarc = lt.math.distance(sta,evmeta[tr_ind])
                # if type.lower()=='metrics':
                    # evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name])
                # el:
                if not type.lower()=='metrics':
                    evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name,'M'+str(evmeta[tr_ind].magnitudes[0].mag),str(np.round(tr.stats.sac.gcarc,2))+'°'])
                    ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='center',fontsize=10)
                # ------
                ax.scatter(x[ind],y_avg,label=':'.join([s,pair]),s=0.5,color=args.csd_pairs[pair])
                ax.plot(x[ind],y_avg,label=':'.join([s,pair]),linewidth=0.8,alpha=0.8,color=args.csd_pairs[pair])
                ax.fill_between(x[ind],y_avg-y_var,y_avg+y_var, alpha=0.2,color=args.csd_pairs[pair])
            elif s=='Raw':
                ax.scatter(x[ind],y_avg,label=':'.join([s+'-Average',pair]),s=0.5,color='r')
                ax.plot(x[ind],y_avg,label=':'.join([s+'-Average',pair]),linewidth=0.8,alpha=0.5,color='r')
        if type.lower()=='metrics':ax.set_xlabel('frequency (hz)')
        else:ax.set_xlabel('seconds')
        ax.axvline(1/fn,alpha=0.4,linewidth=1,color='k',linestyle='-.')
        ax.text(1/fn,1,str(int(fn))+'s',verticalalignment='top',horizontalalignment='right')
    return fig


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ----------------------------------------------------------------------------------------------------------------
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# =XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX==XX=
# ----------------------------------------------------------------------------------------------------------------
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
def station_event_single(st_hold,sta,evmeta,type='stream',**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[1,1]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')];defargs.Noise=True
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    if args.csd_pairs[0][0][0]==args.csd_pairs[0][0][1]:args.Noise=False
    # ------------# ------------
    method=''
    if type.lower()=='stream':
        note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:
        note = 'Noise (gray) | Raw (red)'
    staname = st_hold.Raw.stats.network+'.'+st_hold.Raw.stats.station
    stastr = ' | '.join([method.upper(),staname,sta.Experiment,
    'Depth: '+str(int(1000*abs(st_hold.Raw.stats.sac.stel)))+'m',
    'F-Notch: '+str(int(1/fnotch(1000*abs(st_hold.Raw.stats.sac.stel))))+'s',note])
    if not args.nev:args.nev = len([st_hold.Raw])
    if type.lower()=='stream':columns = args.bands
    elif type.lower()=='metrics':columns=['Coherence','Phase']
    nev_per_plot = args.nev
    nplots = int(np.ceil(nev_per_plot/len(evmeta)))

    channel_color=AttribDict()
    channel_color.ZP='#0c51a6';channel_color.ZZ='#2a7e93'
    channel_color.Z1='#7370cb';channel_color.Z2='#4f86c5'

    # XX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XX
    fig,ax23,ax42 = get_gridplot()
    nax = ax42.size + ax23.size
    plt.tight_layout()
    axis = [0,1,2,3,4,5,0,1,2,3,4,5,6,7]
    fig.suptitle(stastr,y=1.025)
    fn = 1/fnotch(1000*abs(st_hold.Raw.stats.sac.stel))
    meters = ['Coherence','Phase','Coherence','Phase','Coherence','Phase','Coherence','Phase']
    channels = ['ZP','ZP','ZZ','ZZ','Z1','Z1','Z2','Z2']
    bands = [[1,10],[10,30],[30,100],[1,10],[10,30],[30,100]]
    methods = ['ATaCR','ATaCR','ATaCR','Noisecut','Noisecut','Noisecut']
    # XX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XXXX-------XX

    # x = tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[0];ind=x<=1
    # y=tr.Metrics.__getattribute__(b)(args.csd_pairs[0][0])[1][ind]
    # if s=='Raw':
    #     if args.Noise:[ax.scatter(xy[0],xy[1],c='darkgrey',s=0.1,label=args.csd_pairs[pi][0]+':Noise') 
    #     for pi,xy in enumerate([avg_meter(st_hold.Noise,b,p[0]) for p in args.csd_pairs])]
    #     [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.8,color='r') for p in args.csd_pairs]
    #     [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linewidth=0.2,alpha=0.4,color='r') for p in args.csd_pairs]
    # if s=='Corrected':
    #     [ax.scatter(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),s=0.4,color=p[1]) for p in args.csd_pairs]
    #     [ax.plot(x[ind],tr.Metrics.__getattribute__(b)(p[0])[1][ind],label=':'.join([s,p[0]]),linestyle=':',linewidth=0.45,alpha=0.4,color=p[1]) for p in args.csd_pairs]
    # if type.lower()=='metrics':
    # ax.set_xlim(1/500,1);ax.set_xscale('log')
    # if b.lower()=='phase':ax.set_ylim(-180,180)
    # if b.lower()=='coherence':ax.set_ylim(0,1.1)

    # if args.Noise:
    #     [ax.scatter(xy[0],np.abs(xy[1]),c='darkgrey',s=0.1,label=channels+':Noise') 
    #     for ni,xy in enumerate([avg_meter(st_hold.Noise,b,channels) for p in [pair]])]
    for axi in range(nax):
        print(str(evmeta.indi)+'-'+staname+'- Single event page- '+str(axi+1)+'/'+str(nax))
        if axi<ax23.size:
            ax=ax23.reshape(-1)[axis[axi]];c_band=bands[axi];c_method=methods[axi]
            # ------------
            # Trace section
            # ------------
            raw = st_hold.Raw.copy()
            atacr = st_hold.ATaCR.copy()
            hps = st_hold.Noisecut.copy()

            raw.taper(0.0001)
            raw.filter('bandpass',freqmin=1/c_band[1],freqmax=1/c_band[0],zerophase=True,corners=4)
            raw.trim(evmeta[0].origins[0].time,evmeta[0].origins[0].time+7200,pad=True,fill_value=0)

            hps.taper(0.0001)
            hps.filter('bandpass',freqmin=1/c_band[1],freqmax=1/c_band[0],zerophase=True,corners=4)
            hps.trim(evmeta[0].origins[0].time,evmeta[0].origins[0].time+7200,pad=True,fill_value=0)

            atacr.taper(0.0001)
            atacr.filter('bandpass',freqmin=1/c_band[1],freqmax=1/c_band[0],zerophase=True,corners=4)
            atacr.trim(evmeta[0].origins[0].time,evmeta[0].origins[0].time+7200,pad=True,fill_value=0)

            if c_method=='ATaCR':
                method_tr=atacr
            else:
                method_tr=hps
            ylim = np.max(np.array(np.max([np.max(np.abs(c.data))*args.vertical_scale for c in [atacr]])))
            out_scaled = np.max(np.abs(raw.data)) > ylim
            if np.any(out_scaled):
                print('Large amplitude scale differences detected in '+staname)
                raw.data = (raw.data/np.max(np.abs(raw.data)))*ylim
            x = raw.times()
            ax.plot(raw.times(),raw.data,color='red',alpha=1.0,linewidth=0.1)
            ax.plot(method_tr.times(),method_tr.data,color='black',alpha=1.0,linewidth=0.2)
            ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim,ylim);ax.set_yticklabels('')
            ax.set_title(''.join([str(c_band[0]),'s-',str(c_band[1]),'s']),bbox=dict(facecolor='white', alpha=1),pad=0,verticalalignment='baseline',loc='left',y=0)
            if axi>2:
                ax.set_xlabel('seconds')
            else:
                ax.set_xticklabels('')

            # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            stallaz=[raw.stats.sac.stla,raw.stats.sac.stlo,raw.stats.sac.stel]
            evllaz=[evmeta[0].origins[0].latitude,evmeta[0].origins[0].longitude,evmeta[0].origins[0].depth/1000]
            raw.stats.sac.gcarc = locations2degrees(sta,evmeta)
            evstr = '|'.join(['['+str(evmeta.indi)+'] '+c_method,evmeta[0].Name,'M'+str(evmeta[0].magnitudes[0].mag),str(np.round(raw.stats.sac.gcarc,2))+'°'])
            ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='bottom',fontsize=10)
            if type.lower()=='stream':
                colors = args.phasecolors
                if raw.stats.sac.gcarc<=100:phases=args.phases
                else:phases=args.shadow_phases;[colors.update({p:'k'}) for p in phases]
                ar = get_arrivals(stallaz,evllaz,model = 'iasp91',phases=phases)
                [ax.axvline(a[1],linewidth=0.1,color=colors[a[0]]) for a in ar]
                [ax.text(a[1],ylim,a[0],fontsize=6,color='k',verticalalignment='bottom',horizontalalignment='center') for a in ar]
            # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        else:
            ax=ax42.reshape(-1)[axis[axi]];c_chan=channels[axis[axi]];c_meter=meters[axis[axi]]
            # ------------
            # Metrics section
            # ------------
            method_color=dict()
            method_color['ATaCR']='magenta'
            method_color['Noisecut']='blue'
            for ci,c_method in enumerate(['ATaCR','Noisecut']):
                tr = st_hold[c_method]
                x = tr.Metrics.__getattribute__(c_meter)(c_chan)[0]
                ind=x<=1
                y=tr.Metrics.__getattribute__(c_meter)(c_chan)[1][ind]
                if ci==0:
                    if not c_chan=='ZZ':
                        if args.Noise:
                            [ax.scatter(xy[0],xy[1],c='darkgrey',edgecolor='darkgrey',facecolor='darkgrey',s=0.1,label='Noise') 
                            for pi,xy in enumerate([avg_meter(st_hold.Noise,c_meter,p) for p in [c_chan]])]    
                            # [ax.plot(xy[0],xy[1],c='darkgrey',linewidth=0.3,label='Noise') 
                            # for pi,xy in enumerate([avg_meter(st_hold.Noise,c_meter,p) for p in [c_chan]])]
                    [ax.scatter(x[ind],st_hold.Raw.Metrics.__getattribute__(c_meter)(p)[1][ind],s=0.8,color='r') for p in [c_chan]]
                    [ax.plot(x[ind],st_hold.Raw.Metrics.__getattribute__(c_meter)(p)[1][ind],label='Raw',linewidth=0.8,alpha=0.4,color='r') for p in [c_chan]]
                [ax.scatter(x[ind],tr.Metrics.__getattribute__(c_meter)(p)[1][ind],s=2,facecolor=channel_color[c_chan],edgecolors=method_color[c_method]) for p in [c_chan]]
                [ax.plot(x[ind],tr.Metrics.__getattribute__(c_meter)(p)[1][ind],label=c_method,linewidth=0.8,alpha=0.4,color=method_color[c_method]) for p in [c_chan]]
                ax.set_xlim(1/500,1);ax.set_xscale('log')
                if c_meter.lower()=='phase':ax.set_ylim(-180,180)
                if c_meter.lower()=='coherence':ax.set_ylim(0,1.1)
                ax.set_title(c_meter+':'+c_chan,bbox=dict(facecolor='white', alpha=1),pad=0,verticalalignment='baseline',loc='left',y=0)
                if axi>=12:
                    ax.set_xlabel('frequency (hz)')
                if axi==6:
                    lg = ax.legend(ncols=2,loc='upper left')
                    [l.set_linewidth(4) for l in lg.legendHandles]
                    [l.set_sizes([5,5,5]) for l in [lg.legendHandles[0]]]
                if ci==0:
                    ax.axvline(1/fn,alpha=0.4,linewidth=1,color='k',linestyle='-.')

    return fig

def dataset_averaged_coherence_plot(cat,dirs,
    title='Station Averaged Noise',
    meter='Coherence',
    chans='ZP',
    figsize=[10,4],fontsize=12):
    font = {'weight':'bold','size':fontsize};matplotlib.rc('font', **font)
    f = avg_meter(get_noise(dirs,cat.StaName[0]),meter,chans)[0]
    coh = np.array([avg_meter(get_noise(dirs,stanm),meter,chans)[1] 
    for stanm in cat.StaName[np.argsort(cat.StaDepth)]])
    x = np.round(np.sort(cat.StaDepth))
    # ---------------------------------------------------------------
    # More spatially accurate approach but matrix is sparse.
    # xx = np.arange(int(x.min()),int(x.max())+1,int(np.diff(x).min()))
    # zz = np.empty((len(xx),len(f)));zz[:] = 0
    # for ii,xi in enumerate(x):zz[np.where(xi==xx)[0][0],:] = coh_plottable[ii,:]
    # ---------------------------------------------------------------
    fig,ax=plt.subplots(figsize=figsize);x = np.round(np.sort(cat.StaDepth))
    # coh_plottable = np.array([smooth(d,k=3) for d in coh]);coh_plottable = gaussian_filter(coh,.5)
    coh_plottable = coh
    ax.contourf(f,x,coh_plottable,linewidth=2,levels=10,extend='max',vmin=coh.min(),vmax=coh.max())
    fn = [fnotch(fq) for fq in x]
    ax.plot(fn,x,linestyle='dashed',color='w',linewidth=0.4*fontsize)
    ax.set_xlim(1/500,1);ax.set_xscale('log')
    fticks = np.array([1/500,1/300,1/100,1/50,1/10,1])
    ax.set_xticks(fticks)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticklabels(np.array(1/fticks,dtype=int))
    ax.set_ylabel('Station depth (m)',fontweight='bold')
    ax.set_xlabel('Period (s)',fontweight='bold')
    ax.set_title(f'{coh.shape[0]} {title}')
    plt.tight_layout()
    return fig
def catalog_map_plot(sta_inv,ev_cat,
    mode='Depth',
    water_fill_color='lightblue',
    continent_fill_color='lightgrey',
    projection = 'ortho',
    figsize = (15,8),
    cmap=None):
    # projection = ['ortho','global']
    inv = sta_inv.copy()
    evm_plottable = ev_cat.copy()
    if cmap==None:cmap={'mag':'bwr','depth':'tab20c'}[mode.lower()]
    if mode.lower()=='mag':
        for e in evm_plottable:e.origins[0].depth = e.magnitudes[0].mag*1000
    clear_output(wait=False)
    nev = len(evm_plottable)
    nsta = len(np.unique(list(itertools.chain.from_iterable([e.Stations for e in evm_plottable]))))
    fig = evm_plottable.plot(resolution='f',
    water_fill_color=water_fill_color,continent_fill_color=continent_fill_color,
    projection=projection,
    color='depth',label=None)
    settings = AttribDict();settings.map_ortho = AttribDict();settings.map_global=AttribDict()
    settings.map_ortho.shrink=0.4;settings.map_global.shrink=0.7
    title=' | '.join([str(nev)+' Events',str(nsta)+ ' Stations'])
    clear_output(wait=False)
    events_plot = fig.get_axes()[0].get_children()[-14]
    events_plot.set_cmap(cmap)
    events_plot.set_edgecolor('k');events_plot.set_linewidth(0.3)
    sizes = events_plot.get_sizes();events_plot.set_sizes(sizes/3)
    fig.set_size_inches(figsize);fig.set_tight_layout('tight')
    norm = fig.get_axes()[1]._colorbar.norm
    if mode.lower()=='mag':norm._vmin=6.0;norm._vmax=8.0
    fig.get_axes()[1].remove()
    if len(np.unique([e.magnitudes[0].magnitude_type for e in evm_plottable]))>1:mtype='M'
    else:mtype=np.unique([e.magnitudes[0].magnitude_type for e in evm_plottable])[0]
    if mode.lower()=='mag':colorbarlabel = 'Magnitude ('+mtype+')'
    else:colorbarlabel = 'Depth (km)'
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    ax=fig.get_axes()[0], orientation='horizontal',
    label=colorbarlabel,
    shrink=settings['map_'+projection].shrink,aspect=75,pad=0.01)
    plt.draw()
    fig.suptitle(title,y=0.93)
    # ------------
    # -- Station plot
    fig = inv.plot(size=50,label=False,fig=fig,color_per_network=True,marker='v')
    stations_plot = fig.get_axes()[0].get_children()[-13]
    stations_plot.set_edgecolor('k')
    stations_plot.set_linewidth(0.5)
    fig.show()
    lg = fig.get_axes()[0].get_legend()
    labels=lg.texts;handles=lg.legend_handles;lg.remove()
    labels = [l._text for l in labels]
    lg = fig.get_axes()[0].legend(handles, labels, ncols=11,loc='lower left',
    labelspacing=0.0,columnspacing=0.1,
    handletextpad=-.5,borderpad=0.2,edgecolor='k',fontsize=11)
    lg.set_zorder(1000)
    fig.show()
    return fig


def plot_spec_coh_adm_ph(Metrics):
    pairs = ['ZP']
    meters = ['psd','Coherence','Admittance','Phase']
    fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(8,10),layout='constrained',squeeze=False,sharey='row',sharex='all')
    axes = axes.reshape(-1)
    stam = Metrics['ATaCR'].traces[0].stats.network + '.' + Metrics['ATaCR'].traces[0].stats.station
    label = Metrics['ATaCR'].traces.select(channel='*Z')[0].stats.location
    Pre = Metrics['Raw']
    Post = Metrics['ATaCR']
    Noise = Metrics['Noise']
    fn = fnotch(abs(Post.traces[0].stats.sac.stel*1000))
    tstamp = Pre.traces[0].stats.starttime.strftime('%Y.%j.%H.%M')
    for pi,(ax,m) in enumerate(zip(axes,meters)):
        if m=='psd':
            p = 'Z'
        else:
            p = pairs[0]
        evf,prey = Pre.__getattribute__(m)(p)
        evf,posty = Post.__getattribute__(m)(p)
        if m=='psd':
            noisef,noisey = Noise.f,Noise.StaNoise.power.__dict__['cZZ']
        else:
            noisef,noisey = Noise.__getattribute__(m)(p)
        noisey = noisey[noisef>0]
        noisef = noisef[noisef>0]
        if m=='psd':
            noisey = 10*np.log10(noisey)
            prey = 10*np.log10(prey)
            posty = 10*np.log10(posty)
        ax.scatter(noisef,noisey,s=0.5,c='gray',label='Noise')
        if pi==0:
            lbl = stam + ' | ' + label + ' ' + tstamp + ' | '
        else:
            lbl = ''
        ax.set_xlabel('Frequency')
        ax.set_ylabel(m.replace('psd','Power Density'))
        ax.set_title(lbl + p + '-' + m.replace('psd','PSD'),fontweight='bold')
        ax.scatter(evf,prey,c='k',label='PRE',marker='o',s=1)
        ax.scatter(evf,posty,c='m',label='POST',marker='o',s=0.5)
        ax.axvline(fn,linewidth=0.2,color='k')
        if pi==0:
            ax.text(fn*1.05,0.99*min(ax.get_ylim()),'Fn:' + str(round(1/fn*100)/100) + 's',alpha=0.4)
            ax.set_xscale('log')
            ax.set_xlim(evf[1],evf[-1])
            ax.legend(markerscale=10,ncols=len(meters))
    plt.tight_layout()
    return fig
def dataset_averaged_coherence_plot(f,z,coh,
    title='Station Averaged Noise',
    figsize=[15,6],fontsize=12,vlim=[0,1],levels=None,fig=None,ax=None):
    font = {'weight':'bold','size':fontsize};matplotlib.rc('font', **font)
    z = np.round(z)
    i = np.argsort(z)
    z,coh = z[i],coh[i,:]
    if levels is None:levels=np.linspace(np.min(coh),np.max(coh),20)
    # if vlim is None:vlim=[np.min(coh),np.max(coh)]
    if ax is None:fig,ax=plt.subplots(figsize=figsize)
    # coh_plottable = np.array([smooth(d,k=3) for d in coh]);coh_plottable = gaussian_filter(coh,.5)
    coh_plottable = coh
    cnt = ax.contourf(f,z,coh_plottable,vmin=vlim[0],vmax=vlim[1],extend="both",levels=levels)
    fn = [fnotch(fqz) for fqz in z]
    ax.plot(fn,z,linestyle='dashed',color='w',linewidth=0.4*fontsize)
    ax.set_xlim(1/500,1);ax.set_xscale('log')
    fticks = np.array([1/500,1/300,1/100,1/50,1/10,1])
    ax.set_xticks(fticks)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticklabels(np.array(1/fticks,dtype=int))
    ax.set_ylabel('Station depth (m)',fontweight='bold')
    ax.set_xlabel('Period (s)',fontweight='bold')
    # ax.set_title(f'{coh.shape[0]} {title}')
    ax.set_facecolor('k')
    if fig is not None:fig.suptitle(title)
    # plt.colorbar(cnt)
    plt.tight_layout()
    if fig is not None:return fig
    if ax is not None:return ax

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| HPS Spectrogram Plots (Replicates Zali) |||||||||||||||||||||||||
# _________________________________________________________________________________________________________________

def hps_spectrograms(S_full, S_background, S_hps, frequencies, times,figsize=(11,12)):
        t = times
        f = frequencies
        s = S_full
        cmap = 'magma'
        xlabel = True
        yscale = 'log'
        ax = None
        vmin,vmax = None,None
        fig, axes = plt.subplots(nrows=3, ncols=1,figsize=figsize,height_ratios=[1,1,1],width_ratios=[1],layout='constrained',squeeze=False,sharey='col',sharex='row')
        titles = ['Raw','Noise','Noise Removed']
        for r,s in enumerate([S_full, S_background, S_hps]):
                gax = axes[r,0]
                pc = gax.pcolormesh(t,f, 10*np.log10(s), cmap = cmap, shading= 'auto')
                if (vmin is not None) and (vmax is not None):
                        pc.set_clim(vmin, vmax)
                else:
                        vmin,vmax = pc.get_clim()
                if ylabel:
                        gax.set_ylabel('Frequency (Hz)',fontweight='bold')
                if xlabel:
                        gax.set_xlabel('Time (s)',fontweight='bold')
                gax.set_yscale(yscale)
                gax.set_ylim(f[1],f[-1])
                gax.set_xlim(t[0],t[-1])
                gax.set_title(titles[r],fontweight='bold')
                if fig is not None:
                        fig.colorbar(pc, ax=gax, pad=0.01, label='dB')

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| EVENT RECORD FUNCTION |||||||||||||||||||||||||||||||||||||||||||
# _________________________________________________________________________________________________________________
def event_record_plot(evstream,evstream_back=None,linewidth=0.2,trim = (None,None),scales = [1,1],band = None,facecolor=('b','r'),norm = 'trace',figsize = (20,13),phases = ('P','S','SKS','PKiKP','SKiKS','SKSSKS',),evdepth=None,title='',sortindex=None,ax=None,normscale=1.0,residual_fraction=1.0):
        # linewidth=0.2
        # prepare traces
        if evstream_back:
            sets = [evstream_back,evstream] # [UNCORRECT_SET , CORRECT_SET] 
        else:
            sets = [evstream]
        sets = [preparetraces(stream,trim=trim,band=band,sortindex=sortindex) for stream in sets]
        if len(sets)==2:
            residuals = [ev0.data - ev.data for ev0,ev in zip(sets[0],sets[1])]
            residuals = [normscale * residual_fraction * np.array(res) / np.max(np.abs(np.array(res))) for res in residuals]
        fig = None
        if ax is None:
            fig, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
            ax = axes[0,0]
        # norming for plots
        postsetmax = [np.max(ev.data) for ev in  sets[-1]]
        for si in range(len(sets)):
            for i in range(len(sets[0])): #norm the uncorrected set to the max in the corrected set
                # -----------
                if isinstance(norm,list):
                        sets[si][i].data = sets[si][i].data/norm[i]
                elif norm.lower()=='postset':
                        sets[si][i].data = sets[si][i].data/postsetmax[i]
                elif norm.lower()=='trace':
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                elif norm.lower()=='col':
                        # print('Trace norm scaling by r')
                        dist = [s.stats.sac.dist for s in sets[si]]
                        norms = [d/np.max(dist) for d in dist]
                        norms = [d**(1) for d in norms]
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                        sets[si][i].data = sets[si][i].data / norms[i]
                        colmax = np.max([abs(d.data.max()) for d in sets[si]])
                        sets[si][i].data = sets[si][i].data / colmax
                # -----------
        correctset = sets[-1]
        [ax.plot(tr.times(),tr.data + ysep,linewidth=linewidth,color='k') for ysep,tr in enumerate(correctset)]
        if len(sets)>1:
            for si,s in enumerate(sets):
                [ax.fill_between(tr.times(),tr.data*scales[si] + ysep,tr.data*0 + ysep, where=np.abs(tr.data)>=0, facecolor=facecolor[si]) for ysep,tr in enumerate(s)]
        [ax.plot(tr.times(),tr.data*0 + ysep,linewidth=0.4,color='k') for ysep,tr in enumerate(correctset)]
        if evdepth is not None:
            arrivals = [ObsQA.io.get_arrivals(sta_llaz=(sta.stats.sac.stla,sta.stats.sac.stlo,sta.stats.sac.stel),ev_llaz=(sta.stats.sac.evla,sta.stats.sac.evlo,evdepth),phases=phases) for sta in correctset]
            ardict = dict()
            corephase_dy = 0.2
            direcphase_dy = 0.03
            [[ardict.update({ph[0]:[]}) for ph in a] for a in arrivals]
            [[ardict[ph[0]].extend([ph[1]]) for ph in a] for a in arrivals]
            [ardict.update({k:np.max(ardict[k])}) for k in list(ardict.keys())]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='b',linewidth=0.2) for ph in a] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='P'] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='S'] for ysep,a in enumerate(arrivals)]
            [ax.text(ardict[k], len(correctset) - corephase_dy, k, color='b',horizontalalignment='center') for k in list(ardict.keys()) if (k!='P') and (k!='S')]
            [ax.text(ardict[k], len(correctset) + direcphase_dy, k, color='r',horizontalalignment='center') for k in list(ardict.keys()) if ((k=='P') or (k=='S'))]
        yl = (-1,len(correctset))
        ax.set_ylim(yl)
        ax.set_xlim(correctset[0].times()[0],correctset[0].times()[-1])
        ax.set_yticks([i for i in range(len(correctset))])
        labels = [str(int(ev.stats.sac.dist)) +'km' + ' [' + ev.stats.network + '] ' + ev.stats.station  + '\n depth:' + str(int(abs(ev.stats.sac.stel*1000))) + 'm' for ev in correctset]
        ax.set_yticklabels(labels)
        ax.set_xlabel('Time(s)')
        if fig is not None:
            fig.suptitle(title,fontweight='bold',fontsize=15)
        return ax

# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # ++++++++++++++++++++++++++ CONSTRUCTOR AREA ++++++++++++++++++++++++++++++
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# methods = ['PostATACR']
# atacrdatafolder = archive / 'ATaCR_Data' / 'ATaCR_Python'
# for correction_method in methods:
#     coh_comp = correction_method.replace('PostHPS','HPS').replace('PostATACR','ATaCR')
#     if correction_method=='PostHPS':
#       return_hps = True
#     else:
#       return_hps = False
#     OutFolder = Path(plotfolder)
#     SubFolders = Path('EventRecords') / correction_method / 'coherence'
#     OutFolder = OutFolder / SubFolders
#     OutFolder.mkdir(parents=True,exist_ok=True)
#     for station in catalog.iloc:
#       stations = [station.Station]
#       networks = [station.Network]
#       events = station.Events
#       for i,(net,sta) in enumerate(zip(networks,stations)):
#         Metrics = []
#         for evi,event in enumerate(events):
#           depth = round(station.Metadata[evi].origins[0].depth/1000)
#           mag = station.Metadata[evi].magnitudes[0].mag
#           File = '.'.join([net,sta]) + '.m' + str(mag) + '.z' + str(depth) + 'km' + '.' + event + '.' + correction_method.replace('Post','') + '_SPECCOHPHADM.png'
#           title = File.replace('_',' | ').replace('z','z: ').replace('m','mag: m')
#           print('[' + str(evi) + '/' + str(len(events)) + '] ' + File)
#           post_record = Stream()
#           pre_record = Stream()
#           M,Comp = get_metrics_comp(net,sta,atacrdatafolder,event,return_hps=return_hps,events_folder='EVENTS')
#           # M['Noise'] = get_Noise(atacrdatafolder,net,sta,'sta')['Noise']
#           Metrics.append(M.copy())
#           fig = plot_spec_coh_adm_ph(M)
#           save_tight(str(plotfolder / 'MeetingFigs' / 'SPECCOHPHADM' / File),fig,dpi=600)
