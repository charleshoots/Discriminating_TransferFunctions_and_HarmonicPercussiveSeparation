from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from obspy import read
from obspy.geodetics import locations2degrees
from IPython.display import clear_output
from matplotlib.gridspec import GridSpec
def run_sections(st_hold,mdi,mode,mthdi,method,evi,event,cat,dirs):
    stations=np.array([stanm for stanm in event.Stations if mirror(dirs.Events,dirs.Events_HPS,stanm,event.Name)])
    defargs = AttribDict();defargs.bands=[(1,10),(10,30),(30,100)];defargs.vertical_scale=1.2
    defargs.figwidth=20;defargs.figaspect=[1,1]
    defargs.linewidth=[.4,.5];defargs.linecolor=['red','black'];defargs.alpha=[1.0,1.0];defargs.nev=None
    defargs.phases=('P','S');defargs.shadow_phases=('PKIKP','SKS','SKIKSSKIKS');defargs.phasecolors = {'P':'r','S':'b'}
    defargs.csd_pairs=[('ZP','blue'),('ZZ','gray')];defargs.Noise=True
    defargs.prepostclr={'Raw':'r','Corrected':'k'}
    defargs.prepostsz={'Raw':0.8,'Corrected':0.4}
    defargs.prepostlw={'Raw':0.2,'Corrected':0.45}
    defargs.prepostls={'Raw':'-','Corrected':':'}
    defargs.prepostalpha={'Raw':0.7,'Corrected':0.9}
    defargs.height_per_sta=1.0
    # [defargs.update({k:args[k]}) for k in list(args.keys())]
    args = defargs
    # ___________________________________________________________________________________________________________________
    vertical_scale=1.05
    figsize=(8,10) #width x height
    bands=[(1,10),(30,100)]
    colors={'Raw':'red','Corrected':'black'}
    st_raw = st_hold.select(location='*Raw*')
    if isinstance(mode,str):
        if mode.lower()=='traces':columns=bands
        args.figwidth=18
    else:
        columns=mode;mode='Metrics'
        args.figwidth=18
    nsta = len(st_raw)
    clear_output(wait=False);os.system('cls' if os.name == 'nt' else 'clear')
    print(' | '.join([event.Name,str(evi+1)+'/'+str(len(evs)),method,str(mthdi+1)+'/'+str(len(methods)),mode,str(mdi+1)+'/'+str(len(modes))]))

    # ___________________________________________________________________________________________________________________
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # ___________________________________________________________________________________________________________________

    ncols = len(columns)
    figsize=(args.figwidth*args.figaspect[0],args.height_per_sta*nsta*args.figaspect[1])
    fig, axes = plt.subplots(nrows=nsta, ncols=ncols,figsize=figsize,layout='constrained',squeeze=False)
    evstr = '|'.join([method.replace('HPS','Noisecut') +' ',event.Name,str(event.magnitudes[0].mag) + event.magnitudes[0].magnitude_type,str(int(event.origins[0].depth/1000))+'km'])
    # .__getattribute__(b)(p[0])[1][ind]
    note=''
    notches=np.array([fnotch(d) for d in [cat.loc[sta].StaDepth for sta in [f'{s.stats.network}.{s.stats.station}' for s in st_hold]]])
    notches=np.round(notches,4)
    if mode.lower()=='traces':
        for stai in range(nsta):
            arrivals = [event_stream_arrivals(tr,event) for tr in Stream(st_hold.select(location='*Raw*')[stai])][0]
            sta=cat[cat.StaName==stations[stai]].iloc[0]
            statr=st_hold.select(location='*Raw*')[stai]
            stastr = '|'.join([sta.StaName+' ('+sta.Experiment+')',
            'Depth: '+str(int(1000*abs(statr.stats.sac.stel)))+'m',
            'F-Notch: '+str(int(1/fnotch(1000*abs(statr.stats.sac.stel))))+'s',
            note])

            for bi in range(len(columns)):
                notch_filt=['Left','Right'][bi]
                # _______________________________________________________________
                # Trace Plots
                # ===============================================================
                cur_band = bands[bi]
                # -------- Filter and rel. amplitudes

                sta_tr_filt=Stream(st_hold.select(location='*Raw*')[stai].copy())+Stream(st_hold.select(location='*Correct*')[stai].copy())
                sta_tr_filt.taper(.05)
                # if notch_filt=='Left':sta_tr_filt.filter('lowpass',freq=notches[stai],zerophase=True,corners=4)
                # if notch_filt=='Right':sta_tr_filt.filter('highpass',freq=notches[stai],zerophase=True,corners=4)

                if notch_filt=='Left':sta_tr_filt.filter('bandpass',freqmin=1/100,freqmax=notches[stai],zerophase=True,corners=4)
                if notch_filt=='Right':sta_tr_filt.filter('bandpass',freqmin=notches[stai],freqmax=1,zerophase=True,corners=4)
                # sta_tr_filt.taper(.001)
                st_raw = sta_tr_filt.select(location='*Raw*').copy()
                st_corrected = sta_tr_filt.select(location='*Corrected*').copy()
                # st_raw=detect_outscale(st_raw,st_corrected,vertical_scale=vertical_scale,suppress=True)
                del sta_tr_filt;sta_tr_filt=st_raw+st_corrected
                # -------- Plotting
                ax = axes[stai,bi]
                tr = sta_tr_filt
                x=tr.select(location='*Raw*')[0].times()
                ylim = vertical_scale*abs(tr.select(location='*Raw*')[0].data).max()
                [ax.plot(x,tr.select(location='*'+m+'*')[0].data,
                color=args.linecolor[si],
                alpha=args.alpha[si],
                linewidth=args.linewidth[si]) for si,m in enumerate(['Raw','Corrected'])]
                ax.set_xlim(x[0],x[-1]);ax.set_ylim(-ylim,ylim)
                ax.set_yticks([])
                if stai<(nsta-1):
                    ax.set_xticks([])
                    ax.set_xticklabels('')
                else:ax.set_xlabel('Time (s)')
                # if stai==0:ax.set_title(''.join([str(cur_band[0]),'-',str(cur_band[1]),'s']),fontweight='bold',pad=15)
                if stai==0:ax.set_title(f'{notch_filt.replace('Left','Longer').replace('Right','Shorter')} period than noise notch',fontweight='bold',pad=15,fontsize=16)
                # ===============================================================
                if bi==0:
                    colors = args.phasecolors
                    stallaz=[tr[0].stats.sac.stla,tr[0].stats.sac.stlo,tr[0].stats.sac.stel]
                    evllaz=[event.origins[0].latitude,event.origins[0].longitude,event.origins[0].depth/1000]
                    tr[0].stats.sac.gcarc = locations2degrees(stallaz[0],stallaz[1],evllaz[0],evllaz[1])
                    if tr[0].stats.sac.gcarc<=100:phases=args.phases
                    else:phases=args.shadow_phases;[colors.update({p:'k'}) for p in phases]
                    arrivals=[event_stream_arrivals(tr,event) for tr in Stream(st_hold.select(location='*Raw*')[stai])][0]
                    arrivals=[[n,t] for n,t in zip(list(arrivals.keys()),list(arrivals.values()))]
                [ax.axvline(a[1],
                linewidth=0.1,color=colors[a[0]]) for a in arrivals]
                [ax.text(a[1],ylim,a[0],
                fontsize=6,color='k',verticalalignment='bottom',horizontalalignment='center') for a in arrivals]
                if bi==0:
                    ax.text(np.max(ax.get_xlim())*0.995,np.max(ax.get_ylim()),stastr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='top',fontsize=9)
    else:
        for stai in range(nsta):
            tr = Stream(st_hold.select(location='*Raw*')[stai].copy())+Stream(st_hold.select(location='*Correct*')[stai].copy())
            sta=cat[cat.StaName==stations[stai]].iloc[0]
            statr=st_hold.select(location='*Raw*')[stai]
            stastr = '|'.join([sta.StaName+' ('+sta.Experiment+')',
            'Depth: '+str(int(1000*abs(statr.stats.sac.stel)))+'m',
            'F-Notch: '+str(int(1/fnotch(1000*abs(statr.stats.sac.stel))))+'s',
            note])
            for bi in range(len(columns)):
                ax = axes[stai,bi]
                b,p = columns[bi]
                if stai==0:ax.set_title(b+':'+p,pad=15,fontweight='bold')
                if not (p[0]==p[1]):
                    if args.Noise:
                        noiseind = avg_meter(tr[0].Noise,b,p)[0]<=1.0
                        [ax.scatter(xy[0][noiseind],xy[1][noiseind],c='darkgrey',s=0.1,label=p+':Noise') 
                        for xy in [avg_meter(tr[0].Noise,b,p)]]
                # ------
                f=tr.select(location='*Raw*')[0].Metrics.__getattribute__(b)(p)[0]
                ind=f<=1.0;f=f[ind]
                # ------
                [ax.scatter(f,tr.select(location='*'+m+'*')[0].Metrics.__getattribute__(b)(p)[1][ind],
                label=m,
                s=args.prepostsz[m],
                color=args.prepostclr[m]) for m in ['Raw','Corrected']]
                # ------
                [ax.plot(f,tr.select(location='*'+m+'*')[0].Metrics.__getattribute__(b)(p)[1][ind],
                label=m,
                color=args.prepostclr[m],
                alpha=args.prepostalpha[m],
                linestyle=args.prepostls[m],
                linewidth=args.prepostlw[m]) for m in ['Raw','Corrected']]
                ax.set_xlim(1/500,f[-1])
                # ylim={'Coherence':[0,1],'Phase':[-180,180],'Admittance':None}[b]
                # if ylim:ax.set_ylim(ylim[0],ylim[1])
                ax.set_yticklabels('')
                ax.set_xscale('log')
                ax.axvline(fnotch(1000*abs(statr.stats.sac.stel)),linestyle=':',alpha=0.4,color='k')
                if stai<(nsta-1):
                    ax.set_xticks([])
                    ax.set_xticklabels('')
                    plt.setp( ax.get_xticklabels(), visible=False)
                if bi==0:
                    ax.text(np.max(ax.get_xlim())*0.995,np.max(ax.get_ylim()),stastr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='top',fontsize=9)
            del tr
    plt.tight_layout(h_pad=0.00001)
    plt.subplots_adjust(hspace=0)
    fig.suptitle(evstr,fontweight='bold',y=1.01,fontsize=13)
    evstr = '|'.join([event.Name,str(event.magnitudes[0].mag) + event.magnitudes[0].magnitude_type,
    str(int(event.origins[0].depth/1000))+'km'])
    file = (evstr.replace(' ','').replace('|','.')+'_'+mode.lower()+'.png')
    if mode=='Metrics':
        file = file.replace('.png','_cohph.'+''.join([c[0][:2]+c[1] for c in columns])+'.png')
    file = file.replace('.png','__'+method.replace('HPS','Noisecut').upper()+'.png')
    save_tight(dirs.Plots/'_Plots'/'RecordSections'/file,fig,dpi=800)
    plt.close('all')

def distance(sta,ev,unit='deg'):
    origins=ev.origins[0]
    stalla,evlla=[sta.Latitude,sta.Longitude],[origins.latitude,origins.longitude]
    dist=locations2degrees(stalla[0],stalla[1],evlla[0],evlla[1])
    if unit.lower()=='km':dist=degrees2kilometers(dist)
    return dist
def mirror(Afold,Bfold,stanm,evname):
    Afold=dirs.Events/'corrected';Bfold=dirs.Events_HPS/'corrected'
    return np.all([len(list((f/stanm).glob(f'*{evname}*')))>0 for f in [Afold,Bfold]])


cat = catalog.copy()
evcat=Catalog(unravel([e.copy() for e in cat.Events]))
evcat=Catalog([evcat[i].copy() for i in np.unique([e.Name for e in evcat],return_index=True)[1]])
evcat=Catalog([evcat[i].copy() for i in [15]])
# evcat=Catalog(evcat[47].copy()).copy()
# evcat[0].Stations
evs=evcat
# evs=Catalog([e for e in evcat if (len(e.Stations)<=24) & (len(e.Stations)>=20)])

# stations=np.array([stanm for stanm in ev.Stations if mirror(dirs.Events,dirs.Events_HPS,stanm,ev.Name)])

min_sta=20
Events = [c.Events for c in cat.iloc]
EvHold = Events[0]
for e in Events:EvHold+=e
Events = [e for e in EvHold if len(e.Stations)>=min_sta]
Events = [Events[i] for i in np.unique([e.Name for e in Events],return_index=True)[1]]
# evs=Events
# method = 'HPS'
methods = ['HPS','ATaCR']
modes = ['Traces']
for evi,event in enumerate(evs):
    for mthdi,method in enumerate(methods):
        # ___________________________________________________________________________________________________________________
        # ___________________________________________________________________________________________________________________
        # ------------------------------------------------------------
        if method.lower()=='hps':
            evdir=[dirs.Events_HPS,dirs.Events]
        else:
            evdir=dirs.Events
        st_hold = Stream()
        stations=np.array([stanm for stanm in event.Stations if mirror(dirs.Events,dirs.Events_HPS,stanm,event.Name)])
        for si,s in enumerate(stations):
            try:
                if method.lower()=='hps':
                    st_hold+=get_station_events_hps(s,evdir,evmeta=Catalog([event]),type='metrics',tf='HZ.SAC')[0]
                else:
                    st_hold+=get_station_events(s,evdir,evmeta=Catalog([event]),type='metrics')[0]
            except:
                stations.pop(si)
                continue
        # ------------------------------------------------------------
        stasort = np.argsort([abs(s.stats.sac.stel*1000) for s in st_hold.select(location='*Raw*')])
        # stasort = np.argsort([distance(catalog.loc[f'{st.stats.network}.{st.stats.station}'],event) for st in st_hold])
        stasort = np.argsort([distance(catalog.loc[s],event) for s in stations])
        stations=list(np.array(stations)[stasort])
        st_hold = Stream([st_hold.select(location='*Raw*')[i] for i in stasort])+Stream([st_hold.select(location='*Corrected*')[i] for i in stasort])
        # ___________________________________________________________________________________________________________________
        # ___________________________________________________________________________________________________________________
        for mdi,mode in enumerate(modes):
            # xxxxxxx
            run_sections(st_hold,mdi,mode,mthdi,method,evi,event,cat,dirs)
            # xxxxxxx
        del st_hold
