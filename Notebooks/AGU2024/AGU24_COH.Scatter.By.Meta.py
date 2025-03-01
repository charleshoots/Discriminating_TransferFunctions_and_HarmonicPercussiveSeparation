import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
from imports import *
from modules import *
instrument_colors = {'B2':[227,26,28], 'KE':[178,223,138], 'AB':[166,206,227], 'BA':[202,178,214], 'AR':[255,127,0], 'TRM':[31,120,180], 'BG':[51,160,44], 'BD':[106,61,154]}
_ = [instrument_colors.update({k:list(np.array(instrument_colors[k])/255)}) for k in list(instrument_colors.keys())]
seismometer_marker = {'Guralp CMG3T 120':'o','Trillium 240':'x','Trillium Compact':'^'}
colors=np.array([instrument_colors[sta.Deployment.Instrument_Design] for sta in catalog.iloc])
markers=np.array([seismometer_marker[sta.Deployment.Seismometer] for sta in catalog.iloc])
Colors_InstrumentDesigns=np.unique([sta.Deployment.Instrument_Design for sta in catalog.iloc])
Markers_Seismometers=np.unique([sta.Deployment.Seismometer for sta in catalog.iloc])

Comp='ZZ'
Archive=dirs.Archive
AVG=False
xmethod = 'NoiseCut'
ymethod = 'ATaCR'

magwin=[6.0,8.0]
fbands = [[30,100],[30,100],[30,100]]
# fbands = [[30,100]]

# fig,axes = plt.subplots(ncols=len(fbands),nrows=2,figsize=(15+5,8));axes=np.atleast_2d(axes)
fig,axes = plt.subplots(nrows=len(fbands),ncols=1,figsize=(11,13))#;axes=np.atleast_2d(axes)


# axes = axes.T
# axes=axes.reshape(-1)
# ax=axes[0]
bi=0

get_reports(Comp,catalog,Archive,dirs,AVG=True,methods=['ATaCR','NoiseCut'])

colors=np.array([instrument_colors[sta.Deployment.Instrument_Design] for sta in catalog.iloc])
markers=np.array([seismometer_marker[sta.Deployment.Seismometer] for sta in catalog.iloc])
# # -----
ylabel='Station depth';meta_title = 'by deployment depth';meta_i=0;ylim=[-5,6200];yunit='Station depth, m'
# ylabel='Event magnitude';meta_title = 'by event magnitude';meta_i=1;ylim=[5.8,8.1];yunit='Magnitude, M'
# ylabel='Event distance';meta_title = 'by event distance';meta_i=2;ylim=[28,123];yunit='Event distance, degrees'
Report = get_reports(Comp,catalog,Archive,dirs,AVG=True,methods=['ATaCR','NoiseCut'])

f=Report.f
# if ymethod=='ATaCR':
    # if Comp=='ZP':
    #     for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth):
    #         Y[stanm].coh[:,f>fnotch(stadepth)]=1.0
R=AttribDict()
for m in [i for i in list(Report.keys()) if not i=='f']:
    R[m]=AttribDict()
    for n in [i for i in list(Report[m].keys()) if not i=='f']:
        for s in list(Report[m][n].keys()):
            stanm=f'{n[1:]}.{s}'
            R[m][stanm]=AttribDict()
            # if AVG:R[m][stanm].coh=np.mean(Report[m][n][s].coh,axis=0)
            # else:
            R[m][stanm].coh=Report[m][n][s].coh
            R[m][stanm].events=Report[m][n][s].events

Report = R
Report.f = f
f=Report.f
X=Report[xmethod]
Y=Report[ymethod]
catalog=catalog[[np.isin(stanm,list(X.keys()))==True for stanm in catalog.StaName]]
catalog=catalog[[np.isin(stanm,list(Y.keys()))==True for stanm in catalog.StaName]]
if ymethod=='ATaCR':
    for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth):
        if np.isin(stanm,list(Y.keys())):Y[stanm].coh[:,f>fnotch(stadepth)]=1.0
if AVG:
    X=np.array([X[stanm].coh for stanm in catalog.StaName])
    Y=np.array([Y[stanm].coh for stanm in catalog.StaName])
# -----
for bi,band in enumerate(fbands):
    ax=axes[bi]
    if bi<2:
            axtitle='Below infragravity limit frequency'
            Xband=[xx.mean(axis=1) for xx in [X[stanm].coh[:,f<=fnotch(stadepth)] for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth)]]
            Yband=[yy.mean(axis=1) for yy in [Y[stanm].coh[:,f<=fnotch(stadepth)] for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth)]]
    else:
            axtitle='Above infragravity limit frequency'
            Xband=[xx.mean(axis=1) for xx in [X[stanm].coh[:,(f>fnotch(stadepth)) & (f<=1)] for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth)]]
            Yband=[yy.mean(axis=1) for yy in [Y[stanm].coh[:,(f>fnotch(stadepth)) & (f<=1)] for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth)]]
    # ------------------
    # f_ind=(f<=(1/band[0])) & (f>=(1/band[1]))
    # axtitle=f'{min(band)}-{max(band)}s'
    # Xband=[xx.mean(axis=1) for xx in [X[stanm].coh[:,f_ind] for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth)]]
    # Yband=[yy.mean(axis=1) for yy in [Y[stanm].coh[:,f_ind] for stanm,stadepth in zip(catalog.StaName,catalog.StaDepth)]]
    # ------------------
    if AVG:Xband=np.array([np.nanmean(i) for i in Xband]);Yband=np.array([np.nanmean(i) for i in Yband])

    ev_catalogs=dict()
    for i,(evs,evm,stanm) in enumerate(zip([X[stanm].events for stanm in catalog.StaName],catalog.Events,catalog.StaName)):
        evs=evs[:Xband[i].shape[0]]
        c,ia,ib=np.intersect1d(evs,[e.Name for e in evm],return_indices=True)
        Xband[i]=Xband[i][ia];Yband[i]=Yband[i][ia]
        X[stanm].events=X[stanm].events[ia];Y[stanm].events=Y[stanm].events[ia]
        ev_catalogs[stanm]=Catalog([evm[k] for k in ib])

    for si,(xx,yy,clr,mkr,stanm) in enumerate(zip(Xband,Yband,colors,markers,catalog.StaName)):
        filter = np.ones(yy.shape[0],dtype=bool)
        # filter=(np.array([e.magnitudes[0].mag for e in ev_catalogs[stanm]])>=min(magwin)) & (np.array([e.magnitudes[0].mag for e in ev_catalogs[stanm]])<=max(magwin))
        sta=catalog.loc[stanm];events=ev_catalogs[stanm]
        mags = np.array([e.magnitudes[0].mag for e in events])
        dist = np.array([distance(sta,ev) for ev in events])
        water_depth=catalog.loc[stanm].StaDepth
        plotx=xx[filter]
        ploty=yy[filter]
        meta_y = [water_depth+plotx*0,mags,dist][meta_i]

        if bi>0:
            method='NoiseCut'
            ax.scatter(plotx,meta_y,color=clr,s=30,marker=mkr,edgecolor='k',alpha=1.0) #Noisecut
        else:
            method='ATaCR'
            ax.scatter(ploty,meta_y,color=clr,s=30,marker=mkr,edgecolor='k',alpha=1.0) #ATaCR
        # for ax,method in zip(ax,[xmethod,ymethod]):
        ax.set_xlim(0,1.04)
        ax.set_ylim(ylim)
        # ax.set_title(method,fontweight='bold',fontsize=25)
        # if bi==0:ax.set_xticks([])
        corner=min(ax.get_xlim()),max(ax.get_ylim())
        ax.text(corner[0]+0.00,corner[1]*1.005,f'{method} by {ylabel.lower()} : {axtitle}',verticalalignment='bottom',fontsize=18,fontweight='bold',horizontalalignment='left',
        bbox=dict(facecolor='w', edgecolor='k', pad=3))
        if bi==2:ax.set_xlabel(f'Pre-Post ({Comp}) Coherence',fontweight='bold',fontsize=19)
        # ax.set_yticklabels(ax.get_yticklabels(),fontweight='bold',fontsize=19)
        ax.set_xticklabels(ax.get_xticklabels(),fontweight='bold',fontsize=23)
        ax.set_yticklabels(ax.get_yticklabels(),fontweight='bold',fontsize=24)
        if bi==1:ax.set_ylabel(yunit,fontweight='bold',fontsize=23)
fig.suptitle(f'{ylabel}',fontweight='bold',fontsize=30,y=1.01)

plt.tight_layout()
metaname = ylabel.replace(' ','_').replace('(','').replace(')','').replace(',','_').replace('°','')
save_tight(dirs.P01.S05/f'02.28.25.Meta.{metaname}.Scatter.{xmethod}.{ymethod}.{Comp}.NotchFilters.png',fig,dpi=700)
# [e.magnitudes[0].mag for e in catalog.loc[stanm].Events]