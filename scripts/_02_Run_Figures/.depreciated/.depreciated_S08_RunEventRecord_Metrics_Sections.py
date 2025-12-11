# Author: Charles Hoots
# This code was developed as part of my PhD research in the
# Department of Earth Sciences, University of Hawai‘i at Mānoa.
# Unless otherwise noted, the code is my own original work.
# External libraries and standard research software packages are used as cited.

### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

from local_tools.io import *
from local_tools.quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from obspy import read
from obspy.geodetics import locations2degrees
from IPython.display import clear_output
from matplotlib.gridspec import GridSpec
# run sections
def run_sections(st_hold,mdi,mode,mthdi,method,evi,event,cat,plotfold):
    # stations=np.array([stanm for stanm in event.Stations if mirror(dirs.Events,dirs.Events_HPS,stanm,event.Name)])
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
    # figsize value
    figsize=(8,10) #width x height
    # bands value
    bands=[(1,10),(30,100)]
    # colors value
    colors={'Raw':'red','Corrected':'black'}
    # st raw value
    st_raw = st_hold.select(location='*Raw*')
    if isinstance(mode,str):
        if mode.lower()=='traces':columns=bands
        else:columns=[['Coherence','ZZ'],['Phase','ZZ']]
        args.figwidth=18
    else:
        # columns value
        columns=mode;mode='Metrics'
        args.figwidth=18
    # nsta value
    nsta = len(st_raw)
    # successfully loaded value
    successfully_loaded = [f'{s.stats.network}.{s.stats.station}' for s in st_raw]
    # stations value
    stations=successfully_loaded
    # ev distances value
    ev_distances = [distance(cat.loc[sta],event) for sta in stations]
    clear_output(wait=False);os.system('cls' if os.name == 'nt' else 'clear')
    print(' | '.join([event.Name,str(evi+1)+'/'+str(len(evs)),method,str(mthdi+1)+'/'+str(len(methods)),mode,str(mdi+1)+'/'+str(len(modes))]))

    # ___________________________________________________________________________________________________________________
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # ___________________________________________________________________________________________________________________

    # ncols value
    ncols = len(columns)
    # figsize value
    figsize=(args.figwidth*args.figaspect[0],args.height_per_sta*nsta*args.figaspect[1])
    fig, axes = plt.subplots(nrows=nsta, ncols=ncols,figsize=figsize,layout='constrained',squeeze=False)
    # file evstr value
    file_evstr='.'.join([event.Name,str(event.magnitudes[0].mag) + str(event.magnitudes[0].magnitude_type).replace('None','M'),str(int(event.origins[0].depth/1000))+'km depth'])
    # .__getattribute__(b)(p[0])[1][ind]
    note=None
    # notches value
    notches=np.array([fnotch(d) for d in [cat.loc[sta].StaDepth for sta in [f'{s.stats.network}.{s.stats.station}' for s in st_hold]]])
    # notches value
    notches=np.round(notches,4)
    # file value
    file = file_evstr+'_'+mode.lower()+'.png'
    if mode=='Metrics':file = file.replace('.png','_cohph.'+''.join([c[0][:2]+c[1] for c in columns])+'.png')
    # file value
    file = file.replace('.png','__'+method.replace('HPS','Noisecut')+'.png')
    # file value
    file=plotfold/file
    # evstr value
    evstr = ' | '.join(([method.replace('HPS','Noisecut') +' ',event.Name,str(event.magnitudes[0].mag) + str(event.magnitudes[0].magnitude_type).replace('None','M'),str(int(event.origins[0].depth/1000))+'km']))
    if mode.lower()=='traces':evstr = evstr + '\n Sorted by distance'
    else:evstr = evstr + '\n Sorted by depth'
    # if file.exists()&(not ovr):return
    if mode.lower()=='traces':
        for stai in range(nsta):
            net,sta=stations[stai].split('.')
            # sta value
            sta=cat[cat.StaName==stations[stai]].iloc[0]
            # arrivals value
            arrivals = [event_stream_arrivals(tr,event) for tr in st_hold.select(network=sta.Network,station=sta.Station,location='*Raw*')][0]
            # statr value
            statr=st_hold.select(network=sta.Network,station=sta.Station,location='*Raw*')[0]
            # stastr value
            stastr = f'{sta.StaName} ({sta.Experiment}) |  Depth: {str(int(1000*abs(statr.stats.sac.stel)))}m, Notch: {str(int(1/fnotch(1000*abs(statr.stats.sac.stel))))}s |  {int(np.round(ev_distances[stai]))}°'

            for bi in range(len(columns)):
                # notch filt value
                notch_filt=['Left','Right'][bi]
                # -------- Filter and rel. amplitudes
                sta_tr_filt=st_hold.select(network=sta.Network,station=sta.Station,location='*Raw*').copy()+st_hold.select(network=sta.Network,station=sta.Station,location='*Correct*').copy()
                sta_tr_filt.taper(.01)
                # sta_tr_filt.normalize(global_max=True)
                # if notch_filt=='Left':sta_tr_filt.filter('lowpass',freq=notches[stai],zerophase=True,corners=4)
                # if notch_filt=='Right':sta_tr_filt.filter('highpass',freq=notches[stai],zerophase=True,corners=4)
                if notch_filt=='Left':sta_tr_filt.filter('bandpass',freqmin=1/200,freqmax=fnotch(sta.StaDepth),zerophase=True,corners=9)
                if notch_filt=='Right':sta_tr_filt.filter('bandpass',freqmin=fnotch(sta.StaDepth),freqmax=1,zerophase=True,corners=9)
                sta_tr_filt.taper(.01)
                # st_raw=detect_outscale(st_raw,st_corrected,vertical_scale=vertical_scale,suppress=True)
                # -------- Plotting
                ax = axes[stai,bi]
                # tr value
                tr = sta_tr_filt
                x=tr.select(network=sta.Network,station=sta.Station,location='*Raw*')[0].times()
                # ylim value
                ylim = vertical_scale*abs(tr.select(network=sta.Network,station=sta.Station,location='*Raw*')[0].data).max()
                [ax.plot(x,tr.select(location='*'+m+'*')[0].data,
                # color value
                color=args.linecolor[si],
                # alpha value
                alpha=args.alpha[si],
                # linewidth value
                linewidth=args.linewidth[si]) for si,m in enumerate(['Raw','Corrected'])]
                ax.set_xlim(x[0],x[-1]);ax.set_ylim(-ylim,ylim)
                ax.set_yticks([])
                if stai<(nsta-1):
                    ax.set_xticks([])
                    ax.set_xticklabels('')
                else:ax.set_xlabel('Time (s)')
                # if stai==0:ax.set_title(''.join([str(cur_band[0]),'-',str(cur_band[1]),'s']),fontweight='bold',pad=15)
                if stai==0:ax.set_title(f'{notch_filt} of infragravity limit',fontweight='bold',pad=15,fontsize=16)
                # ===============================================================
                if bi==0:
                    # colors value
                    colors = args.phasecolors
                    # stallaz value
                    stallaz=[tr[0].stats.sac.stla,tr[0].stats.sac.stlo,tr[0].stats.sac.stel]
                    # evllaz value
                    evllaz=[event.origins[0].latitude,event.origins[0].longitude,event.origins[0].depth/1000]
                    tr[0].stats.sac.gcarc = locations2degrees(stallaz[0],stallaz[1],evllaz[0],evllaz[1])
                    if tr[0].stats.sac.gcarc<=100:phases=args.phases
                    else:phases=args.shadow_phases;[colors.update({p:'k'}) for p in phases]
                    arrivals=[event_stream_arrivals(tr,event) for tr in st_hold.select(network=sta.Network,station=sta.Station,location='*Raw*')][0]
                    arrivals=[[n,t] for n,t in zip(list(arrivals.keys()),list(arrivals.values()))]
                [ax.axvline(a[1],
                linewidth=0.1,color=colors[a[0]]) for a in arrivals]
                [ax.text(a[1],ylim,a[0],
                fontsize=6,color='k',verticalalignment='bottom',horizontalalignment='center') for a in arrivals]
                if bi==0:
                    ax.text(np.max(ax.get_xlim())*0.995,np.max(ax.get_ylim()),stastr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='top',fontsize=9)
    else:
        for stai in range(nsta):
            net,sta=stations[stai].split('.')
            tr = st_hold.select(network=net,station=sta,location='*Raw*').copy()+st_hold.select(network=net,station=sta,location='*Correct*').copy()
            sta=cat[cat.StaName==stations[stai]].iloc[0]
            statr=st_hold.select(network=net,station=sta.Station,location='*Raw*')[0]
            stastr = f'{sta.StaName} ({sta.Experiment}) |  Depth: {str(int(1000*abs(statr.stats.sac.stel)))}m, Notch: {str(int(1/fnotch(1000*abs(statr.stats.sac.stel))))}s |  {int(np.round(ev_distances[stai]))}°'
            if note:stastr+=f'{note} |'
            for bi in range(len(columns)):
                ax = axes[stai,bi]
                b,p = columns[bi]
                if stai==0:ax.set_title(b+': '+p,pad=15,fontweight='bold')
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
                ylim={'Coherence':[0,1],'Phase':[-180,180],'Admittance':None}
                ax.set_ylim(ylim[b][0],ylim[b][1])
                ax.set_yticklabels('')
                ax.set_xscale('log')
                ax.axvline(fnotch(1000*abs(statr.stats.sac.stel)),linestyle=':',alpha=0.4,color='k')
                if stai<(nsta-1):
                    ax.set_xticks([])
                    ax.set_xticklabels('')
                    plt.setp( ax.get_xticklabels(), visible=False)
                if bi==0:
                    ax.text(np.max(ax.get_xlim())*0.995,min(ylim[b]),stastr,bbox=dict(boxstyle="square,pad=0.1",facecolor='white', alpha=1),horizontalalignment='right',verticalalignment='bottom',fontsize=9)
            del tr
    plt.tight_layout(h_pad=0.00001)
    plt.subplots_adjust(hspace=0)
    fig.suptitle(evstr,fontweight='bold',y=1.01,fontsize=13)
    evstr = '|'.join([event.Name,str(event.magnitudes[0].mag) + str(event.magnitudes[0].magnitude_type).replace('None','M'),
    str(int(event.origins[0].depth/1000))+'km'])
    save_tight(file,fig,dpi=800)
    plt.close('all')

def distance(sta,ev,unit='deg'):
    origins=ev.origins[0]
    stalla,evlla=[sta.Latitude,sta.Longitude],[origins.latitude,origins.longitude]
    dist=locations2degrees(stalla[0],stalla[1],evlla[0],evlla[1])
    if unit.lower()=='km':dist=degrees2kilometers(dist)
    return dist
def mirror(Afold,Bfold,stanm,evname):
    if Afold is None:Afold=dirs.Events/'corrected'
    if Bfold is None:Bfold=dirs.Events_HPS/'corrected'
    return np.all([len(list((f/stanm).glob(f'*{evname}*')))>0 for f in [Afold,Bfold]])


cat = catalog.copy()
evcat=lt.cat.unravel_cat(cat)
evs=evcat
min_sta=10
# methods = ['HPS','ATaCR']
methods = ['ATaCR']
modes = ['Metrics']
columns = [['Coherence','ZZ'],['Phase','ZZ']]
ovr=False
plotfold = dirs.P01.S08;plotfold.mkdir(exist_ok=True,parents=True)
for evi,event in enumerate(evs):
    clear_output()
    print(f'{evi+1}/{len(evs)} | {event.Name} | Stations:{len(event.Stations)}')
    if len(event.Stations)<min_sta:
        continue #This should never happen since the minimum density is now set to atleast 10 stations per event
    for mthdi,method in enumerate(methods):
        mode=modes[0]
        evstr='.'.join([event.Name,str(event.magnitudes[0].mag) + str(event.magnitudes[0].magnitude_type).replace('None','M'),str(int(event.origins[0].depth/1000))+'km'])
        file = evstr+'_'+mode.lower()+'.png'
        if mode=='Metrics':file = file.replace('.png','_cohph.'+''.join([c[0][:2]+c[1] for c in columns])+'.png')
        file = file.replace('.png','__'+method.replace('HPS','Noisecut').upper()+'.png')
        file=plotfold/file
        fold=plotfold/mode
        if file.exists()&(not ovr):continue
        # ___________________________________________________________________________________________________________________
        # ___________________________________________________________________________________________________________________
        # ------------------------------------------------------------
        if method.lower()=='hps':
            evdir=[dirs.Events_HPS,dirs.Events]
        else:
            evdir=dirs.Events
        st_hold = Stream()
        stations=np.array([stanm for stanm in event.Stations if mirror(dirs.Events/'corrected',dirs.Events_HPS/'raw',stanm,event.Name)])
        if mode=='Metrics':
            stasort=np.argsort([catalog.loc[s].StaDepth for s in stations])
        else:stasort = np.argsort([distance(catalog.loc[s],event) for s in stations])
        stations=np.array(stations)[stasort]
        for si,s in enumerate(stations):
            try:
                if method.lower()=='hps':
                    st_hold+=get_station_events_hps(s,evdir,evmeta=Catalog([event]),type='metrics',tf='HZ.SAC')[0]
                else:
                    e=Catalog([event])
                    st_hold+=get_station_events(s,dirs.Events,evmeta=Catalog([event]),type='metrics')[0]
            except:
                # stations.pop(si)
                continue
        # ------------------------------------------------------------
        # stasort = np.argsort([abs(s.stats.sac.stel*1000) for s in st_hold.select(location='*Raw*')])
        # stasort = np.argsort([distance(catalog.loc[f'{st.stats.network}.{st.stats.station}'],event) for st in st_hold])
        raw=st_hold.select(location='*Raw*')
        corrected=st_hold.select(location='*Corrected*')
        if len(corrected)<min_sta:continue
        if len(raw)<min_sta:continue
        st_hold=raw+corrected
        for mdi,mode in enumerate(modes):
            # xxxxxxx
            # (st_hold,mdi,mode,mthdi,method,evi,event,cat,plotfold)
            run_sections(st_hold,mdi,mode,mthdi,method,evi,event,cat,fold)
            # xxxxxxx
        del st_hold