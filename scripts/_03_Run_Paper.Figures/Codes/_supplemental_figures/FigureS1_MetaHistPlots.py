### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *





# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_supplemental_figures'/'FigureS1_meta_plot';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'

# cat value
cat = catalog.copy()
# icat value
icat=cat.sr.copy()

# full_event_catalog = lt.cat.unravel_cat(cat)
MetaHold = dict()
MetaHold['Event_depths_km'] = icat.EvDepth.to_numpy()
MetaHold['Event_Magnitude_M'] = icat.Magnitude.to_numpy()
MetaHold['Event_Distances'] = icat.Distance.to_numpy()
MetaHold['Events_per_Station'] = np.array([sum(icat.StaName==s) for s in icat.StaName.unique()])
MetaHold['Stations_per_Event'] = np.array([len(icat.loc[e].iloc[0].Stations) for e in icat.Name.unique()])
MetaHold['Deployment_Depths'] = cat.r.StaDepth.to_numpy()
MetaHold['Pressure_Gauges'] = cat.r.Pressure_Gauge.to_numpy()
MetaHold['Experiment'] = '('+cat.r.Network.to_numpy()+')'+cat.r.Experiment.to_numpy()
MetaHold['Noise_Notch_period'] = 1/fnotch(cat.r.StaDepth)
MetaHold['Regions'] = icat.Environment
# deployment metas value
deployment_metas = [
 'Distance_from_Land_km','Distance_to_Plate_Boundary_km',
 'Sediment_Thickness_m','Surface_Current_ms','Crustal_Age_Myr',
 'Deployment_Length_days','Instrument_Design','Seismometer']
# _ = [MetaHold.update({m:icat[m]}) for m in deployment_metas]
_ = [MetaHold.update({m:np.array([icat[icat.StaName==s].iloc[0][m] for s in icat.StaName.unique()])}) for m in deployment_metas]
# labels value
labels=list(MetaHold.keys())
# n ax value
n_ax=len(list(MetaHold.keys()))

# rows value
rows=int(np.ceil(n_ax/3))
# cols value
cols=3

fig,axes=plt.subplots(ncols=cols,nrows=rows,figsize=(20,8),layout='constrained',squeeze=True)
# axes value
axes=axes.reshape(-1)
# bins value
bins={'Event_depths_km':[0,100,200,300,400,500,600,700],
 'Stations_per_Event':[0,5,10,15,20,25,30,35,40],
 'Events_per_Station':[20,30,40,50,60,70,80],
 'Event_Magnitude_M':[6,6.5,7.0,7.5,8.0],
 'Event_Distances':np.arange(0,180+45,45),
 'Deployment_Depths':np.arange(0,7000,1000),
 'Pressure_Gauges':2,
 'Noise_Notch_period':np.arange(5,70,5),
 'Experiment':11,
 'Regions':10,
 'Distance_from_Land_km':np.arange(0,2750,250)[::2],
 'Distance_to_Plate_Boundary_km':np.arange(0,4250,250)[::2],
 'Sediment_Thickness_m':np.arange(0,9000,1000),
 'Surface_Current_ms':(np.arange(0,1,.1)[::2]*10).astype(int)/10,
 'Crustal_Age_Myr':np.arange(0,160,30),
 'Deployment_Length_days':np.arange(180,450,60),
 'Instrument_Design':8,
 'Seismometer':3}


# centers value
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
# categorical labels value
categorical_labels=['Seismometer','Instrument_Design', 'Regions','Experiment','Pressure_Gauges']
txi=0
for axi,(ax,label) in enumerate(zip(axes,labels)):
    if 'Event' in label:ylabel='Events'
    if 'Station' in label:ylabel='Stations'
    if 'Stations_per_Event' in label:ylabel='Events'
    else:ylabel='Stations'
    xlabel = label.replace('_period',' period, s').replace('_km',', km').replace('_ms',' m/s').replace('_m',', m').replace('_Myr',', Myr').replace('_days',', days').replace('_',' ')
    xlabel= xlabel[0].upper() +xlabel[1:].lower()
    xlabel = xlabel.replace('Deployment depths','Deployment depths, m')
    if xlabel=='Event magnitude m':xlabel=('Event magnitude, M')
    if xlabel=='Event distances':xlabel=f'{xlabel}, °'
    isnumber= not isinstance(MetaHold[label][0],str)
    if not isnumber:bn=len(np.unique(MetaHold[label]))
    else:bn=bins[label]
    n,b,p=ax.hist(MetaHold[label],bins=bn,edgecolor='k',rwidth=0.9,facecolor='darkgrey')
    # ax.set_xticks(b[:-1] + np.diff(b)[0]/2)
    if isnumber:
        if label=='Surface_Current_ms':ax.set_xticklabels([f for f in [tx._text for tx in ax.get_xticklabels()]])
        elif not label=='Event_Magnitude_M':ax.set_xticklabels([f for f in [tx._text for tx in ax.get_xticklabels()]])
    else:ax.set_xticklabels(ax.get_xticklabels(),rotation=[0,32,15,0,0,0][txi]);txi+=1
    ax.set_xlabel(xlabel.lower().replace('mw','Mw'),fontweight='bold')
    if np.any(np.array(['Event_Magnitude_M','Event_Distances','Event_depths_km'])==label):ax.set_ylabel('Events'.lower())
    else:ax.set_ylabel(ylabel.lower())

    if label in categorical_labels:ax.set_xticks(centers(b))
    else:ax.set_xticks(bn[1:]);ax.set_xticklabels(bn[1:])
[ax.axis('off') for axi,ax in enumerate(axes) if (axi+1)>len(labels)]
file=f'meta_plot.{save_format}'
save_tight(plotfolder/file,fig,dpi=700)
