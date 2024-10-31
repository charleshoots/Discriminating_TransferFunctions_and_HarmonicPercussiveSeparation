from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from obspy import read
from obspy.geodetics import locations2degrees
# ========================================================================================================================================================
from IPython.display import clear_output
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
from local_tools import *

def get_gridplot():
    fig = plt.figure(figsize=(20,15),layout="constrained")
    height_ratios = [1,1,0.6,0.6,0.6,0.6]
    gs = GridSpec(6, 6, figure=fig,height_ratios=height_ratios)
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
            print(method+'-'+type + '| Column:'+str(bi+1)+'/'+str(len(columns)))
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
                        tr.stats.sac.gcarc = locations2degrees(stallaz[0],stallaz[1],evllaz[0],evllaz[1])
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
def station_event_page_averages(st_hold,sta,evmeta,method,type='stream',raw_reference=None,**args):
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2;defargs.figwidth=20;defargs.figaspect=[1,4]
    defargs.linewidth=[.1,.2];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')];defargs.Noise=True
    defargs.columns = ['Coherence']
    [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    args.csd_pairs = [('ZP','#0c51a6'),('ZZ','#2a7e93'),('Z1','#7370cb'),('Z2','#4f86c5')]
    # ------------
    if type.lower()=='stream':note = 'Corrected ('+args.linecolor[1]+') | Raw ('+args.linecolor[0]+')'
    else:note = 'Noise (gray) | Raw (red) | Variance (shaded)'
    staname = st_hold[0].stats.network+'.'+st_hold[0].stats.station
    if method.lower()=='atacr':tf='('+st_hold.select(location='*Correct*')[0].stats.location.split('.')[1]+')'
    else:tf=''
    stastr = ' | '.join([method.upper(),staname+tf,sta.Experiment,'Depth: '+str(int(1000*abs(st_hold[0].stats.sac.stel)))+'m',
    'F-Notch: '+str(int(1/fnotch(1000*abs(st_hold[0].stats.sac.stel))))+'s',note])
    if not args.nev:args.nev = len(st_hold.select(location='*Raw*'))
    nrows = 2;ncols=2
    columns = args.columns
    # if type.lower()=='metrics':columns=['Coherence']
    # if type.lower()=='metrics':columns=['Coherence','Phase','Admittance']
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,layout='constrained',sharex='all',squeeze=True,figsize=(args.figwidth*args.figaspect[0],nrows*args.figaspect[1]))
    axes = axes.reshape(-1)
    fig.suptitle(stastr)
    fn = 1/fnotch(1000*abs(st_hold[0].stats.sac.stel))
    for bi,b in enumerate(columns):
        for pi,pair in enumerate(args.csd_pairs):
            channels = pair[0]
            if channels[0]==channels[1]:args.Noise=False
            else: args.Noise=True
            print(method+'-'+type + '| Column:'+str(bi+1)+'/'+str(len(columns)))
            ax = axes[pi]
            st_band = st_hold.copy()
            correct_hold = st_band.select(location='*Correct*').copy()
            raw_hold = st_band.select(location='*Raw*').copy()
            if raw_reference:raw_hold = raw_reference
            for si,s in enumerate(['Raw','Corrected']):
                st = {'Raw':raw_hold,'Corrected':correct_hold}[s].copy()
                for tri in range(len(correct_hold)):
                    tr_ind = tri;tr = st[tr_ind].copy()
                    if tr_ind==0:
                        if type.lower()=='stream':ax.set_title(''.join([str(b[0]),'s-',str(b[1]),'s']))
                        else:ax.set_title(channels +' '+b)
                    if type.lower()=='metrics':
                        x = tr.Metrics.__getattribute__(b)(channels)[0];ind=x<=1
                        y=tr.Metrics.__getattribute__(b)(channels)[1][ind]
                        if s=='Raw':
                            if args.Noise:[
                            ax.scatter(xy[0],np.abs(xy[1]),c='darkgrey',s=0.1,label=channels+':Noise')
                            for ni,xy in enumerate([avg_meter(st_hold.Noise,b,channels) for p in [pair]])]
                            [ax.scatter(x[ind],np.abs(tr.Metrics.__getattribute__(b)(p[0])[1][ind]),label=':'.join([s,p[0]]),s=0.4,color='r',alpha=0.1) for p in [pair]]
                            [ax.plot(x[ind],np.abs(tr.Metrics.__getattribute__(b)(p[0])[1][ind]),label=':'.join([s,p[0]]),linewidth=0.05,alpha=0.05,color='r') for p in [pair]]
                        if s=='Corrected':
                            [ax.scatter(x[ind],np.abs(tr.Metrics.__getattribute__(b)(p[0])[1][ind]),label=':'.join([s,p[0]]),s=0.2,color=p[1],alpha=0.1) for p in [pair]]
                            [ax.plot(x[ind],np.abs(tr.Metrics.__getattribute__(b)(p[0])[1][ind]),label=':'.join([s,p[0]]),linestyle=':',linewidth=0.05,alpha=0.05,color=p[1]) for p in [pair]]
                    if type.lower()=='metrics':
                        ax.set_xlim(1/500,1);ax.set_xscale('log')
                        # if b.lower()=='phase':ax.set_ylim(-180,180)
                        if b.lower()=='phase':ax.set_ylim(0,180)
                        if b.lower()=='coherence':ax.set_ylim(0,1.01)
                    else:ax.set_xlim(x[0],x[-1]);ax.set_ylim(-1*ylim[tr_ind],ylim[tr_ind]);ax.set_yticklabels('')
                    # ------------------------------------------------------------------------------------------
                y_avg = np.mean(([rt.Metrics.__getattribute__(b)(channels)[1][ind] for rt in st]),axis=0)
                y_var = np.std(([rt.Metrics.__getattribute__(b)(channels)[1][ind] for rt in st]),axis=0)**2

                # Should I set Phase domain to [0,180] instead of [-180,180]? Would be easier to read...
                # y_avg = np.mean(np.abs([rt.Metrics.__getattribute__(b)(channels)[1][ind] for rt in st]),axis=0)
                # y_var = np.std(np.abs([rt.Metrics.__getattribute__(b)(channels)[1][ind] for rt in st]),axis=0)**2
                if s=='Corrected':
                    stallaz=[tr.stats.sac.stla,tr.stats.sac.stlo,tr.stats.sac.stel]
                    evllaz=[evmeta[tr_ind].origins[0].latitude,evmeta[tr_ind].origins[0].longitude,evmeta[tr_ind].origins[0].depth/1000]
                    tr.stats.sac.gcarc = locations2degrees(stallaz[0],stallaz[1],evllaz[0],evllaz[1])
                    evstr = '|'.join(['['+str(tr_ind+1)+'] ',evmeta[tr_ind].Name,'M'+str(evmeta[tr_ind].magnitudes[0].mag),str(np.round(tr.stats.sac.gcarc,2))+'°'])
                    ax.text(np.max(ax.get_xlim()),np.min(ax.get_ylim()),evstr,bbox=dict(facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='center',fontsize=10)
                    # ------
                    ax.scatter(x[ind],y_avg,label=':'.join([s,channels]),s=0.5,color=pair[1])
                    ax.plot(x[ind],y_avg,label=':'.join([s,channels]),linewidth=0.8,alpha=0.8,color=pair[1])
                    ax.fill_between(x[ind],y_avg-y_var,y_avg+y_var, alpha=0.2,color=pair[1])
                elif s=='Raw':
                    ax.scatter(x[ind],y_avg,label=':'.join([s+'-Average',channels]),s=0.5,color='r')
                    ax.plot(x[ind],y_avg,label=':'.join([s+'-Average',channels]),linewidth=0.8,alpha=0.5,color='r')
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
            raw.stats.sac.gcarc = locations2degrees(stallaz[0],stallaz[1],evllaz[0],evllaz[1])
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