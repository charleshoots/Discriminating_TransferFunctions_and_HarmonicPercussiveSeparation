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

from local_tools.quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
import operator
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
# dirs value
dirs=io.dir_libraries()
# ========================================================================================================================================================
# ========================================================================================================================================================
# ========================================================================================================================================================
# ========================================================================================================================================================
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
# variable
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
# seismometer marker value
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
# write pickle
def write_pickle(file,var):
    import pickle
    with open(str(file), 'wb') as handle:
        pickle.dump(var, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('Saved to :' + str(file))
# load pickle
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        # b value
        b = pickle.load(handle)
    return b
# function mirror events
def mirror_events(reports):
    # nkeys value
    nkeys = [n for n in list(reports[0].__dict__.keys()) if not n=='f']
    # mirror value
    mirror = dict()
    for ni,n in enumerate(nkeys):
        # skeys value
        skeys = list(reports[0][n].__dict__.keys())
        for si,s in enumerate(skeys):
            # stanm value
            stanm = '.'.join([n,s]).replace('n','')
            # ev0 value
            ev0 = [k.replace('.','') for k in reports[0][n][s].events]
            # ev1 value
            ev1 = [k.replace('.','') for k in reports[1][n][s].events]
            mirror[stanm] = np.intersect1d(ev0,ev1)
    return mirror
# ========================================================================================================================================================
# ========================================================================================================================================================
# ========================================================================================================================================================
# ========================================================================================================================================================
reportfolder = dirs.Data/'Analysis'/'NetworkCoherences'
# bands value
bands = ['1-10','10-30','30-100']
# cat value
cat = catalog.copy()
# nets value
nets = list(cat.sr.Network.unique())
# method value
method = 'HPS'
# file value
file = str(reportfolder / method.lower() / 'Complete' / ('complete_' + method.lower() + '.ZZ_coh.report.pkl'))
# hps report value
hps_report = load_pickle(file)
# method value
method = 'ATaCR'
# file value
file = str(reportfolder / method.lower() / 'Complete' /('complete_' + 'ATaCR'.lower() + '.ZZ_coh.report.pkl'))
# atacr report value
atacr_report = load_pickle(file)
reports,methods = [atacr_report,hps_report],['ATaCR','HPS']
# reports,methods = [atacr_report,atacr_report],['ATaCR','ATaCR']
mirror = mirror_events(reports)


def report_parser(cat,report,sort='StaDepth',network=None,band=[],average=False,slice=None):
    # icat value
    icat = cat.copy()
    # f value
    f = report['f'];f_ind = np.isfinite(f)
    # ireport value
    ireport = report.copy()
    if network:ireport = ireport['n'+network];icat = icat[icat.Network==network]
    if not isinstance(band,str):
        if len(band)>0:f_ind = (f>=band[0]) & (f<=band[1])
    else:
        # f ind value
        f_ind = f>-1
    # f ind notch value
    f_ind_notch = lambda f,z,lr: f<=fnotch(z) if lr=='left' else f>fnotch(z)
    if slice:
        # slice bools value
        slice_bools = np.array([[d[key]<slice[key] or d[key]==slice[key]for d in icat.Deployment.iloc]for key in list(slice.keys())])
        # slice filter value
        slice_filter = np.sum(slice_bools,axis=0)==slice_bools.shape[0]
        # icat value
        icat = icat[slice_filter]
    if network==None:
        # events value
        events=[ireport['n'+sta.Network][sta.Station].events for sta in icat.iloc]
        # inds value
        inds = [np.intersect1d([e.Name for e in s.Events],events[si],return_indices=True)[2] for si,s in enumerate(icat.iloc)]
        coh = [ireport['n'+sta.Network][sta.Station].coh[inds[si],:] for si,sta in enumerate(icat.iloc)]
    else:
        events=[ireport[sta.Station].events for sta in icat.iloc]
        inds = [np.intersect1d([e.Name for e in s.Events],events[si],return_indices=True)[2] for si,s in enumerate(icat.iloc)]
        coh = [ireport[sta.Station].coh[inds[si],:] for si,sta in enumerate(icat.iloc)]
    zzcoh_xf_ysta = np.array([np.mean(np.array(c)[:,f_ind if not isinstance(band,str) else f_ind_notch(f,z,band)],axis=1) for c,z in zip(coh,icat.StaDepth)],dtype=object)
    zzcoh_xsta_yevent = [np.mean(np.array(c)[:,f_ind if not isinstance(band,str) else f_ind_notch(f,z,band)],axis=1) for c,z in zip(coh,icat.StaDepth)]
    if sort=='Magnitude':
        x = []
        for si,(s,e) in enumerate(zip(icat.iloc,events)):
            ets = np.intersect1d(s.Events,e,return_indices=True)[2];evs = [s.EventMeta[eid].magnitudes[0].mag for eid in ets]
            x.append(evs)
        y = zzcoh_xsta_yevent
    if sort=='f':x = f;y = zzcoh_xf_ysta
    if sort=='stalta':x=[];y=[]
    if sort=='StaDepth':x = list(icat.StaDepth);y = zzcoh_xf_ysta
    return x,y
slice = {
'Water_Depth_m':5000,'Pressure_Gauge': 'DPG','Environment': 'North Pacific',
'Distance_from_Land_km': 100000,'Distance_to_Plate_Boundary_km': 50000,
'Sediment_Thickness_m': 5000,'Deployment_Length_days': 10000,'Seismometer':'Trillium Compact'}
# -----

for method_report,method in zip(reports,methods):
    cat = catalog.copy()
    f = method_report['f']
    Instrument_Design = [d.Instrument_Design for d in cat.Deployment]
    Seismometer = np.array([d.Seismometer for d in cat.Deployment])
    networkaverage = dict();[networkaverage.update({n:[0,0,0]}) for n in cat.Network]
    instrumentaverage = dict();[instrumentaverage.update({n:[0,0,0]}) for n in np.unique(Instrument_Design)]
    seismometeraverage = dict();[seismometeraverage.update({n:[0,0,0]}) for n in np.unique(Seismometer)]
    for n in cat.Network:
        for bi,b in enumerate(bands):
            band_hz = [1/int(a) for a in b.split('-')];ind = (f<band_hz[0]) & (f>=band_hz[1])
            x,y = report_parser(cat,method_report,sort='StaDepth',band=np.flip(band_hz),network=n)
            y=[np.array([j for j in a]) for a in y]
            y=[a[~np.isnan(a)] for a in y]
            networkaverage[n][bi] = np.mean([np.mean(a) for a in y])
    for d in np.unique(Instrument_Design):
        di = np.array(Instrument_Design)==d
        for bi,b in enumerate(bands):
            band_hz = [1/int(a) for a in b.split('-')];ind = (f<band_hz[0]) & (f>=band_hz[1])
            x,y = report_parser(cat,method_report,sort='StaDepth',band=np.flip(band_hz))
            y=[np.array([j for j in a]) for a in y]
            y=[a[~np.isnan(a)] for a in y]
            instrumentaverage[d][bi] = np.mean([np.mean(a) for i,a in zip(di,y) if i])
            # np.sum([np.sum(k) for k in y[di]]) / np.sum([len(k) for k in y[di]])
    for d in np.unique(Seismometer):
        di = np.array(Seismometer)==d
        for bi,b in enumerate(bands):
            band_hz = [1/int(a) for a in b.split('-')];ind = (f<band_hz[0]) & (f>=band_hz[1])
            x,y = report_parser(cat,method_report,sort='StaDepth',band=np.flip(band_hz))
            y=[np.array([j for j in a]) for a in y]
            y=[a[~np.isnan(a)] for a in y]
            seismometeraverage[d][bi] = np.mean([np.mean(a) for i,a in zip(di,y) if i])
            # np.sum([np.sum(k) for k in y[di]]) / np.sum([len(k) for k in y[di]])
    # ======XXXXX===========XXXXX===========XXXXX===========XXXXX===========XXXXX=====

    # ======XXXXX===========XXXXX===========XXXXX===========XXXXX===========XXXXX=====
    for ni,n in enumerate(nets):
        if (method=='HPS') & (n=='7D'):
            k=1
        icat = cat[cat.Network==n].copy()
        icat = icat.sort_values(by='StaDepth',ascending=True)
        if len(icat)>11:fig,axes = plt.subplots(nrows=3,ncols=1,squeeze=True,figsize=(23,7),sharey='all')
        else:fig,axes = plt.subplots(nrows=3,ncols=1,squeeze=True,figsize=(14,8.66),sharey='all')
        depth =[int(s.StaDepth) for s in icat.iloc]
        fn = [round(100/fnotch(s.StaDepth))/100 for s in icat.iloc]
        stanm = [s.StaName for s in icat.iloc]
        seismometer = [s.Deployment.Seismometer for s in icat.iloc]
        note = f'(UPDATED 2.14.25) ZZ Coherences| '
        fig.suptitle(note + icat.iloc[0].Experiment + ' (' + n + ')' + '| '+method.replace('HPS','Noisecut')+'' + '\n Checkers: Station bands not corrected in ATaCR'
        + '\n' + 'Seismometer average:  Guralp CMG3T 120=●, Trillium 240=X , Trillium Compact=▲'
        '\nInstrument average: Square  ,  Network average: Dashed line'
        ,y=0.98)
        print(method,n,'['+str(ni+1)+'/'+str(len(nets))+']')
        # labels = [s +',n=' + '\n[' + str(round(10*(d/1000))/10) + 'km] \nNotch:' + str(int(notch*10)/10) + 's' for seis,s,d,notch in zip(stanm,depth,fn)]
        f = method_report['f']
        inst_hold = []


        for axi,(ax,b )in enumerate(zip(axes,bands)):
            band_hz = [1/int(a) for a in b.split('-')];ind = (f<band_hz[0]) & (f>=band_hz[1])
            outside_band = [si for si,sta in enumerate(icat.iloc) if len(np.where(f[ind] < fnotch(sta.StaDepth))[0])==0]
            band = np.flip(band_hz)
            x,y = report_parser(cat,method_report,sort='StaDepth',band=np.flip(band_hz),network=n)
            # x,y = report_parser(cat,method_report,sort='Magnitude',band=np.flip(band_hz),network=None)
            labels,outside_band,yy = [],[],[]
            for si,sta in enumerate(icat.iloc):
                    Instrument_Design = sta.Deployment.Instrument_Design
                    sta_report = method_report['n'+sta.Network][sta.Station]
                    notch = round(100/fnotch(sta.StaDepth))/100
                    coh = sta_report.coh
                    events = sta_report.events;events = [e.replace('.','') for e in events]
                    evind = list(np.sort(np.intersect1d(events,mirror[sta.StaName],return_indices=True)[1]))
                    nev = len(evind)
                    fband_coh = coh[evind][:,ind]
                    fband_coh = np.mean(abs(fband_coh),axis=1)
                    if len(np.where(f[ind] < fnotch(sta.StaDepth))[0])==0:
                        if method=='ATaCR':outside_band.append(si)
                            # coh = coh*0+1
                    yy.append(fband_coh)
                    labels.append(sta.StaName+'\n n='+str(nev)+'\n' + str(round(10*(int(sta.StaDepth)/1000))/10) + 'km,' + str(int(int(notch*10)/10)) + 's')
            if len(outside_band)<len(icat):
                pass
                k=1
            ax.set_title(b + 's')
            positions = np.arange(len(yy))+1

            yy=[i[~np.isnan(i)] for i in yy]

            if axi==2:
                bplot = ax.boxplot(yy,patch_artist=True,labels=labels,positions=positions)
            else:
                bplot = ax.boxplot(yy,patch_artist=True,labels=None,positions=positions)
                ax.set_xticklabels([])
            if (len(outside_band)>0) & (method=='ATaCR'):
                halfwdth=0.25
                errorboxes = [ax.add_patch(Rectangle((positions[si]-halfwdth,0),width=2*halfwdth,height=1.1,hatch='//',fill=False,alpha=0.2,edgecolor='grey')) for si in outside_band]
                [ax.text(positions[si],0.5,'no corrections',rotation='vertical',horizontalalignment='center',alpha=0.2,verticalalignment='center') for si in outside_band]
            inst_hold = []
            for ii,patch in enumerate(bplot['boxes']):
                marker = seismometer_marker[icat.iloc[ii].Deployment.Seismometer]
                fliercolor = instrument_colors[icat.iloc[ii].Deployment.Instrument_Design]
                patch.set_facecolor(fliercolor)
                tracker = icat.iloc[ii].Deployment.Instrument_Design
                if tracker not in inst_hold:
                    patch.set_label(icat.iloc[ii].Deployment.Instrument_Design)
                    inst_hold.append(tracker)
                if (method=='ATaCR') & (ii in outside_band):
                        alpha = 0.5
                        bplot['fliers'][ii].set_alpha(alpha);bplot['medians'][ii].set_alpha(alpha)
                        bplot['caps'][2*ii].set_alpha(alpha);bplot['caps'][2*ii+1].set_alpha(alpha)
                        bplot['whiskers'][2*ii].set_alpha(alpha);bplot['whiskers'][2*ii+1].set_alpha(alpha)
                        bplot['boxes'][ii].set_alpha(alpha)
                        del alpha
                patch.set_linewidth(0.1)
                patch.set_edgecolor(fliercolor);patch.set_linewidth(0.1)
                bplot['fliers'][ii].set_markersize(20)
                bplot['fliers'][ii].set_marker('_')
                bplot['fliers'][ii].set_markerfacecolor(fliercolor)
                bplot['fliers'][ii].set_markeredgecolor(fliercolor)
                bplot['medians'][ii].set_color('k')
                bplot['caps'][2*ii].set_color(fliercolor);bplot['caps'][2*ii+1].set_color(fliercolor)
                bplot['whiskers'][2*ii].set_color(fliercolor);bplot['whiskers'][2*ii+1].set_color(fliercolor)
            ax.set_ylim(0,1.05) # ----# seismometeraverage# ----
            # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            ax.axhline(networkaverage[n][axi],color='k',linestyle=':')
            [ax.scatter(positions[ri]+0.25,instrumentaverage[r][axi],color=instrument_colors[r],marker='s',edgecolors='k',s=35) for ri,r in enumerate([d.Instrument_Design for d in icat.Deployment])]
            [ax.scatter(positions[ri]-0.25,seismometeraverage[r][axi],marker=seismometer_marker[r],color='k',s=35) for ri,r in enumerate([d.Seismometer for d in icat.Deployment])]
            axes[2].legend(ncols=len(inst_hold),loc='upper right', bbox_to_anchor=(1, 1.35))
            if len(icat)>37:
                ax.set_xticklabels(ax.get_xticklabels(),fontsize=6)
                fig.set_figwidth(23)
            # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        # handles = axes[2].get_legend_handles_labels()[0]
        handles,labels = axes[2].get_legend_handles_labels()
        axes[2].get_legend().remove()
        # labels_tmp = [c.Deployment.Instrument_Design for c in icat.iloc]
        # hl = sorted(zip(handles, labels_tmp), key=operator.itemgetter(1))
        # handles2, labels2 = zip(*hl)
        axes[0].legend(handles, labels,ncols=len(inst_hold),loc='upper right', bbox_to_anchor=(1, 1.35))
        plt.tight_layout()
        figfold = dirs.P01.S06/'bands'
        figfold.mkdir(parents=True,exist_ok=True)
        figfile = figfold / ('banded_'+icat.iloc[0].Experiment.replace(' ','_') + '.' + n + '.' + method.lower() + '.png')
        save_tight(figfile,dpi = 600)
        plt.close('all')