from imports import *
from modules import *
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
colors=np.array([instrument_colors[sta.Deployment.Instrument_Design] for sta in catalog.iloc])
markers=np.array([seismometer_marker[sta.Deployment.Seismometer] for sta in catalog.iloc])
Colors_InstrumentDesigns=np.unique([sta.Deployment.Instrument_Design for sta in catalog.iloc])
Markers_Seismometers=np.unique([sta.Deployment.Seismometer for sta in catalog.iloc])

Archive=dirs.Archive

ymethod='ATaCR'
xmethod='Noisecut'

colors=np.array([instrument_colors[sta.Deployment.Instrument_Design] for sta in catalog.iloc])
markers=np.array([seismometer_marker[sta.Deployment.Seismometer] for sta in catalog.iloc])
Colors_InstrumentDesigns=np.unique([sta.Deployment.Instrument_Design for sta in catalog.iloc])
Markers_Seismometers=np.unique([sta.Deployment.Seismometer for sta in catalog.iloc])
leg=0
legmarkersize=400
markersize=300 #40

fbands = [[1,10],[30,100]]
# Comps=['ZP','ZP','ZP','ZP']
# Comps=['Z1','Z1','Z1','Z1']
# Comps=['ZP','Z1','Z2']
Comps=['ZZ','ZZ','ZZ','ZZ']

c,r=len(fbands),len(Comps)
fig,axes = plt.subplots(nrows=len(Comps),ncols=len(fbands),figsize=(1.3*8*c,1.3*6*r))
axes=np.atleast_2d(axes)
AVG=False
for ci,Comp in enumerate(Comps):
    cat = catalog.copy()
    Report = get_reports(Comp,cat,Archive,AVG=AVG)
    f=Report.f
    X=Report[xmethod]
    Y=Report[ymethod]
    cat=cat[[np.isin(stanm,list(X.keys()))==True for stanm in cat.StaName]]
    cat=cat[[np.isin(stanm,list(Y.keys()))==True for stanm in cat.StaName]]
    # if ymethod=='ATaCR':
        # if Comp=='ZP':
        #     for stanm,stadepth in zip(cat.StaName,cat.StaDepth):
        #         Y[stanm].coh[:,f>fnotch(stadepth)]=1.0
    if AVG:
        X=np.array([X[stanm].coh for stanm in cat.StaName])
        Y=np.array([Y[stanm].coh for stanm in cat.StaName])
    for bi,band in enumerate(fbands):
        ax = axes[ci,bi]
        if bi==0:
                axtitle='Below infragravity limit frequency'
                Xband=[xx.mean(axis=1) for xx in [X[stanm].coh[:,f<=fnotch(stadepth)] for stanm,stadepth in zip(cat.StaName,cat.StaDepth)]]
                Yband=[yy.mean(axis=1) for yy in [Y[stanm].coh[:,f<=fnotch(stadepth)] for stanm,stadepth in zip(cat.StaName,cat.StaDepth)]]
        else:
                axtitle='Above infragravity limit frequency'
                Xband=[xx.mean(axis=1) for xx in [X[stanm].coh[:,(f>(20*fnotch(stadepth))) & (f<1)] for stanm,stadepth in zip(cat.StaName,cat.StaDepth)]]
                Yband=[yy.mean(axis=1) for yy in [Y[stanm].coh[:,(f>(20*fnotch(stadepth))) & (f<1)] for stanm,stadepth in zip(cat.StaName,cat.StaDepth)]]
        if AVG:
            # Takes station average from the event averages
            Xband=np.array([np.nanmean(i) for i in Xband])
            Yband=np.array([np.nanmean(i) for i in Yband])

        ev_catalogs=dict()
        for i,(evs,evm,stanm) in enumerate(zip([X[stanm].events for stanm in cat.StaName],cat.Events,cat.StaName)):
            evs=evs[:Xband[i].shape[0]]
            c,ia,ib=np.intersect1d(evs,[e.Name for e in evm],return_indices=True)
            Xband[i]=Xband[i][ia]
            Yband[i]=Yband[i][ia]
            X[stanm].events=X[stanm].events[ia]
            Y[stanm].events=Y[stanm].events[ia]
            ev_catalogs[stanm]=Catalog([evm[k] for k in ib])

        magwin=[[6.0,6.5],[6.5,7.0],[7.0,7.5],[7.5,8.0]][ci]
        distances=[[0,40],[40,60],[60,100],[100,150]][ci]
        depths=[[0,2000],[2000,3000],[3000,4000],[4000,6000]][ci]
        # magwin=[[0.0,9.0],[0.0,9.0],[0.0,9.0],[0.0,9.0]][ci]

        # for i,stadepth in enumerate(cat.StaDepth):
        #     if (stadepth>=min(depths)) & (stadepth<max(depths)):
        #         ax.scatter(Xband[i],Yband[i],color=colors[i],s=markersize,marker=markers[i],edgecolor='k',alpha=1.0)

        # sc = [ax.scatter(
        # xx[(np.array([e.Distance for e in ev_catalogs[stanm]])>=min(distances)) & (np.array([e.Distance for e in ev_catalogs[stanm]])<max(distances))],
        # yy[(np.array([e.Distance for e in ev_catalogs[stanm]])>=min(distances)) & (np.array([e.Distance for e in ev_catalogs[stanm]])<max(distances))],
        # color=clr,s=markersize,marker=mkr,edgecolor='k',alpha=1.0)
        # for xx,yy,clr,mkr,stanm in zip(Xband,Yband,colors,markers,cat.StaName)]

        sc = [ax.scatter(
        xx[(np.array([e.magnitudes[0].mag for e in ev_catalogs[stanm]])>=min(magwin)) & (np.array([e.magnitudes[0].mag for e in ev_catalogs[stanm]])<max(magwin))],
        yy[(np.array([e.magnitudes[0].mag for e in ev_catalogs[stanm]])>=min(magwin)) & (np.array([e.magnitudes[0].mag for e in ev_catalogs[stanm]])<max(magwin))],
        color=clr,s=markersize,marker=mkr,edgecolor='k',alpha=1.0)
        for xx,yy,clr,mkr,stanm in zip(Xband,Yband,colors,markers,cat.StaName)]
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)

        [ax.scatter(-100,-100,s=legmarkersize,marker='s',color=instrument_colors[c],label=c) for c in Colors_InstrumentDesigns]
        [ax.scatter(-100,-100,s=legmarkersize,marker=seismometer_marker[c],color='k',label=c) for c in Markers_Seismometers]

        if (bi==1) and (ci==0) and (leg==0):
            lg = ax.legend(
            ncols=len(Colors_InstrumentDesigns)+len(Markers_Seismometers),loc='upper center',
            bbox_to_anchor=(0.5, 1.20),fontsize=20);leg=1

        font_props = {'weight' : 'bold','size': 18}
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):_=label.set_font(font_props)

        if bi==0:ax.set_ylabel(f'{ymethod}\nEvent {Comp} Coherence',fontweight='bold',fontsize=18)
        if (ci+1)==len(Comps):ax.set_xlabel(f'Event {Comp} Coherence\n{xmethod}',fontweight='bold',fontsize=18)
        # ax.legend(loc='lower right')
        ax.plot([0,1],[0,1],linestyle=(0, (3, 1, 1, 1)),c='k',linewidth=1.5,alpha=0.9)
        ax.set_title(axtitle,fontsize=20,fontweight='bold')

        # ========================
        # f'm{depths[0]}-{depths[1]}'
        # ax.text(0.99,0.01,f'{depths[0]}-{depths[1]}m',verticalalignment='bottom',fontsize=28,horizontalalignment='right',
        # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        # # ax.text(.01,0.01,Comp,verticalalignment='bottom',fontweight='bold',fontsize=40,horizontalalignment='left',
        # # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        # ax.text(0,1.0,Comp,verticalalignment='top',fontweight='bold',fontsize=40,horizontalalignment='left',
        # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        # ========================
        # f'm{distances[0]}-{distances[1]}'
        # ax.text(0.99,0.01,f'{distances[0]}-{distances[1]}°',verticalalignment='bottom',fontsize=28,horizontalalignment='right',
        # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        # # ax.text(.01,0.01,Comp,verticalalignment='bottom',fontweight='bold',fontsize=40,horizontalalignment='left',
        # # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        # ax.text(0,1.0,Comp,verticalalignment='top',fontweight='bold',fontsize=40,horizontalalignment='left',
        # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        f'm{magwin[0]}-{magwin[1]}'
        ax.text(0.99,0.01,f'M{magwin[0]}-{magwin[1]}',verticalalignment='bottom',fontsize=28,horizontalalignment='right',
        bbox=dict(facecolor='w', edgecolor='k', pad=3))
        # ax.text(.01,0.01,Comp,verticalalignment='bottom',fontweight='bold',fontsize=40,horizontalalignment='left',
        # bbox=dict(facecolor='w', edgecolor='k', pad=3))
        ax.text(0,1.0,Comp,verticalalignment='top',fontweight='bold',fontsize=40,horizontalalignment='left',
        bbox=dict(facecolor='w', edgecolor='k', pad=3))

        # ================================================
        # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        # ================================================
c='.'.join(Comps)
if AVG:c=c+'.AVG'
fold = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/RunningNotebooks/znb_images/UpdatedData_1.31.25')
# save_tight(fold/f'Scatter.{xmethod}.{ymethod}.{c}.bymag.notch.bands.png',fig,dpi=600)
save_tight(fold/f'Scatter.{xmethod}.{ymethod}.{c}.bymag.and.notch.png',fig,dpi=600)
# save_tight(fold/f'Scatter.{xmethod}.{ymethod}.{c}.bydist.notch.bands.png',fig,dpi=600)
# save_tight(fold/f'Scatter.{xmethod}.{ymethod}.{c}.bydepth.notch.bands.png',fig,dpi=600)