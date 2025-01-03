from imports import *
full_event_catalog=Catalog(unravel([s.Events for s in catalog.iloc]))
full_event_catalog=Catalog([full_event_catalog[i] for i in np.unique([e.Name for e in full_event_catalog],return_index=True)[1]])
ev=full_event_catalog[0]
MetaHold = dict()
MetaHold['Event_depths_km'] = [ev.origins[0].depth/1000 for ev in full_event_catalog]
MetaHold['Event_Magnitude_M'] = [ev.magnitudes[0].mag for ev in full_event_catalog]
MetaHold['Event_Distances'] = unravel([[distance(catalog.loc[s],ev) for s in ev.Stations if s in catalog.StaName] for ev in full_event_catalog])
MetaHold['Stations_per_Event'] = [len(ev.Stations) for ev in full_event_catalog]
MetaHold['Deployment_Depths'] = catalog.StaDepth
MetaHold['Pressure_Gauges'] = catalog.Pressure_Gauge
MetaHold['Experiment'] = '('+catalog.Network+')'+catalog.Experiment
MetaHold['Noise_Notch_period'] = 1/fnotch(catalog.StaDepth)
MetaHold['Regions'] = catalog.Environment
deployment_metas = [
 'Distance_from_Land_km','Distance_to_Plate_Boundary_km',
 'Sediment_Thickness_m','Surface_Current_ms','Crustal_Age_Myr',
 'Deployment_Length_days','Instrument_Design','Seismometer']
_ = [MetaHold.update({m:unravel([[s[m] for s in catalog.Deployment]])}) for m in deployment_metas]
labels=list(MetaHold.keys())
n_ax=len(list(MetaHold.keys()))
fig,axes=plt.subplots(ncols=3,nrows=int(np.ceil(n_ax/3)),figsize=(20,8),layout='constrained',squeeze=True)
axes=axes.reshape(-1)
bins={'Event_depths_km':[0,100,200,300,400,500,600,700],
 'Stations_per_Event':None,
 'Event_Magnitude_M':[6,6.5,7.0,7.5,7.9],
 'Event_Distances':np.arange(0,180,30),
 'Deployment_Depths':np.arange(0,6000,1000),
 'Pressure_Gauges':2,
 'Noise_Notch_period':11,
 'Experiment':11,
 'Regions':10,
 'Distance_from_Land_km':None,
 'Distance_to_Plate_Boundary_km':None,
 'Sediment_Thickness_m':np.arange(0,5000,500),
 'Surface_Current_ms':None,
 'Crustal_Age_Myr':np.arange(0,160,30),
 'Deployment_Length_days':None,
 'Instrument_Design':8,
 'Seismometer':3}

txi=0
for axi,(ax,label) in enumerate(zip(axes,labels)):
    if 'Event' in label:ylabel='Events'
    if 'Station' in label:ylabel='Stations'
    else:ylabel='Stations'
    xlabel = label.replace('_period',', period').replace('_km',', km').replace('_m',', m').replace('_ms',', m/s').replace('_Myr',', Myr').replace('_days',', days').replace('_',' ')
    xlabel= xlabel[0].upper() +xlabel[1:].lower()
    xlabel = xlabel.replace('Deployment depths','Deployment depths, m')
    if xlabel=='Event magnitude m':xlabel=('Event magnitude, M')
    if xlabel=='Event distances':xlabel=f'{xlabel}, °'
    isnumber= not isinstance(MetaHold[label][0],str)
    if not isnumber:bn=len(np.unique(MetaHold[label]))
    else:bn=bins[label]
    n,b,p=ax.hist(MetaHold[label],bins=bn,edgecolor='k',rwidth=0.9,facecolor='darkgrey')
    ax.set_xticks(b[:-1] + np.diff(b)[0]/2)
    if isnumber:
        if label=='Surface_Current_ms':ax.set_xticklabels([str(np.round(float(f),2)) for f in [tx._text for tx in ax.get_xticklabels()]])
        elif not label=='Event_Magnitude_M':ax.set_xticklabels([str(int(float(f))) for f in [tx._text for tx in ax.get_xticklabels()]])
    else:ax.set_xticklabels(ax.get_xticklabels(),rotation=[0,32,15,0,0,0][txi]);txi+=1
    ax.set_xlabel(xlabel,fontweight='bold')
    if np.any(np.array(['Event_Magnitude_M','Event_Distances','Event_depths_km'])==label):ax.set_ylabel('Events')
    else:ax.set_ylabel(ylabel)
[ax.axis('off') for axi,ax in enumerate(axes) if (axi+1)>len(labels)]
save_tight(dirs.Plots/'_Plots'/'meta_plot.png',fig,dpi=700)