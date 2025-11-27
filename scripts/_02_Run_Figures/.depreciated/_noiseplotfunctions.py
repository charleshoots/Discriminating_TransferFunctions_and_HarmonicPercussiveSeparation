### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.modules import *

import matplotlib.cm as cm
from modules import *
# function ph adm
def ph_adm(noise,sta):
    # colors value
    colors = [cm.Categorical.__dict__[e].resampled(4).colors for e in ['devonS']][0]
    # stanm value
    stanm = sta.StaName
    fig = plt.figure(figsize=(25, 20))
    # gs value
    gs = gridspec.GridSpec(4, 3)
    # topaxes value
    topaxes=[]
    topaxes.append(fig.add_subplot(gs[0, 0])) 
    topaxes.append(fig.add_subplot(gs[0, 1]))  
    topaxes.append(fig.add_subplot(gs[0, 2]))  
    topaxes.append(fig.add_subplot(gs[1, 0])) 
    topaxes.append(fig.add_subplot(gs[1, 1]))  
    topaxes.append(fig.add_subplot(gs[1, 2]))  
    # bottomaxes value
    bottomaxes=[]
    bottomaxes.append(fig.add_subplot(gs[2, 0])) 
    bottomaxes.append(fig.add_subplot(gs[2, 1]))  
    bottomaxes.append(fig.add_subplot(gs[2, 2]))  
    bottomaxes.append(fig.add_subplot(gs[3, 0])) 
    bottomaxes.append(fig.add_subplot(gs[3, 1]))  
    bottomaxes.append(fig.add_subplot(gs[3, 2]))
    # f value
    f=noise.f
    # m value
    m = 'ph';mname = 'Phase'
    # keys value
    keys = ['ZP','Z1','Z2','P1','P2','21']
    # daycolor value
    daycolor = 'darkgray';stacolor = 'b'
    # daysize value
    daysize=.5;stasize=2
    # gooddays value
    gooddays = load_sta_atacrnoise(stanm).gooddays
    # dayplot value
    dayplot=[[ax.scatter(f[f>0],e[f>0],s=daysize,color=daycolor if good else 'r',label='Day average') for e,good in zip(noise.day[m][a],gooddays)] for ax,a in zip(topaxes,keys)][4][0]
    # staplot value
    staplot=[ax.scatter(f[f>0],noise.sta[m][a][f>0],s=stasize,color=colors[0] if m is not 'ph' else 'b',label='Station average') for ax,a in zip(topaxes,keys)][4]
    # variable
    _=[ax.set_title(f'{a} {mname}',fontweight='bold') for ax,a in zip(topaxes,keys)]
    # variable
    _=[ax.set_xscale('log')for ax in topaxes]
    # variable
    _=[ax.set_xlim(1/500,1)for ax in topaxes]
    # variable
    _=[ax.axvline(fnotch(sta.StaDepth),linewidth=1,color='k',linestyle=':') for ax in topaxes]
    # variable
    _=topaxes[4].legend(handles=[staplot,dayplot],markerscale=5,ncols=2,fontsize=12,loc='upper center', prop={'weight': 'bold'})
    m = 'adm';mname = 'Admittance'
    keys = ['ZP','Z1','Z2','P1','P2','21']
    daycolor = 'darkgray';stacolor = 'b'
    dayplot=[[ax.scatter(f[f>0],e[f>0],s=1,color=daycolor if good else 'r',label='Day average') for e,good in zip(noise.day[m][a],gooddays)] for ax,a in zip(bottomaxes,keys)][4][0]
    staplot=[ax.scatter(f[f>0],noise.sta[m][a][f>0],s=3,color=colors[0] if m is not 'ph' else 'b',label='Station average') for ax,a in zip(bottomaxes,keys)][4]
    _=[ax.set_title(f'{a} {mname}',fontweight='bold') for ax,a in zip(bottomaxes,keys)]
    _=[ax.set_xscale('log')for ax in bottomaxes]
    _=[ax.set_xlim(1/500,1)for ax in bottomaxes]
    _=[ax.axvline(fnotch(sta.StaDepth),linewidth=1,color='k',linestyle=':') for ax in bottomaxes]
    _=bottomaxes[4].legend(handles=[staplot,dayplot],markerscale=5,ncols=2,fontsize=12,loc='upper center', prop={'weight': 'bold'})
    # disptoaccel = lambda f,s: 10*np.log10(s) + 40*np.log10(2*np.pi*f,where=f>0.) - np.mean(abs(10*np.log10(s)+40*np.log10(2*np.pi*f,where=f>0.)))
    keys = ['Z','1','2','P']
    fig.suptitle(f'Noise | {stanm} ({sta.Experiment}) | {np.round(sta.StaDepth)}m',fontsize=20,fontweight='bold',y=0.91)
    return fig
def coh_psd(noise,sta):
    colors = [cm.Categorical.__dict__[e].resampled(4).colors for e in ['devonS']][0]
    stanm = sta.StaName
    fig = plt.figure(figsize=(25, 20))
    meter = 'coh'
    # gs = gridspec.GridSpec(4, 6)
    # Create a main GridSpec with 4 rows and 6 columns
    gs = gridspec.GridSpec(4, 6, figure=fig)
    # Create sub-gridspecs for the top and bottom rows
    gs_top = gs[:2, :]  # First two rows
    gs_bottom = gs[2:, :]  # Last two rows
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 6,subplot_spec=gs_top)
    w=0.2
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 6, width_ratios=[2,2,w,w,2,2],subplot_spec=gs_bottom)
    if meter=='coh':
        topaxes=[]
        topaxes.append(fig.add_subplot(gs1[0, 0:2])) 
        topaxes.append(fig.add_subplot(gs1[0, 2:4]))  
        topaxes.append(fig.add_subplot(gs1[0, 4:]))  
        topaxes.append(fig.add_subplot(gs1[1, 0:2])) 
        topaxes.append(fig.add_subplot(gs1[1, 2:4]))  
        topaxes.append(fig.add_subplot(gs1[1, 4:]))  
        bottomaxes=[]
        bottomaxes.append(fig.add_subplot(gs2[0, 0:2]))
        bottomaxes.append(fig.add_subplot(gs2[0, 4:]))
        bottomaxes.append(fig.add_subplot(gs2[1, 0:2]))
        bottomaxes.append(fig.add_subplot(gs2[1, 4:]))
        goodax = fig.add_subplot(gs2[0:,2:4])
        # badax = fig.add_subplot(gs[3, 3])
    f=noise.f
    m = 'coh';mname = 'Coherence'
    keys = ['ZP','Z1','Z2','P1','P2','21']
    daycolor = 'darkgray';grayalpha=0.4
    stacolor = 'b';daysize=.5;stasize=2
    gooddays = load_sta_atacrnoise(stanm).gooddays
    days = [f.replace('.spectra.pkl','') for f in load_sta_atacrnoise(stanm).day_files]
    dayplot=[[ax.scatter(f[f>0],e[f>0],s=daysize,color=daycolor if good else 'r',label='Day average',alpha=grayalpha if good else 1) for e,good in zip(noise.day.coh[a],gooddays)] for ax,a in zip(topaxes,keys)][4][0]
    staplot=[ax.scatter(f[f>0],noise.sta.coh[a][f>0],s=stasize,color='k' if m is not 'psd' else 'k',label='Station average') for ax,a in zip(topaxes,keys)][4]
    # _=[[ax.scatter(f[f>0],e[f>0],s=daysize,color=daycolor if good else 'r',label='Day average',alpha=grayalpha if good else 1) for e,good in zip(noise.day.coh[a],gooddays) if not good] for ax,a in zip(topaxes,keys)]
    _=[ax.set_title(f'{a} {mname}',fontweight='bold') for ax,a in zip(topaxes,keys)]
    # _=[ax.set_xscale('log')for ax in topaxes]
    _=[ax.set_xlim(1/500,1)for ax in topaxes]
    _=[ax.axvline(fnotch(sta.StaDepth),linewidth=1,color='k',linestyle=':') for ax in topaxes]
    _=topaxes[4].legend(handles=[staplot,dayplot],markerscale=5,ncols=2,fontsize=12,loc='upper center', prop={'weight': 'bold'})
    disptoaccel = lambda f,s: 10*np.log10(s) + 40*np.log10(2*np.pi*f,where=f>0.) - np.mean(abs(10*np.log10(s)+40*np.log10(2*np.pi*f,where=f>0.)))
    keys = ['Z','1','2','P']
    daycolor = 'darkgray'
    stacolor = 'b'
    m = 'psd';mname = 'Power Spectral Density'

    dayplot=[[ax.scatter(f[f>0],disptoaccel(f[f>0],e[f>0]) if a is not 'P' else 10*np.log10(e[f>0]),s=daysize,color=daycolor if good else 'r',label='Day average',alpha=grayalpha if good else 1) for e,good in zip(noise.day.psd[a],gooddays)] for ax,a in zip(bottomaxes,keys)][0][0]
    staplot=[ax.scatter(f[f>0],disptoaccel(f[f>0],noise.sta.psd[a][f>0]) if a is not 'P' else 10*np.log10(noise.sta.psd[a][f>0]),s=stasize,color='k' if m is not 'psd' else 'k',label='Station average') for ax,a in zip(bottomaxes,keys)][0]
    # _=[[ax.scatter(f[f>0],disptoaccel(f[f>0],e[f>0]) if a is not 'P' else 10*np.log10(e[f>0]),s=daysize,color=daycolor if good else 'r',label='Day average',alpha=grayalpha if good else 1) for e,good in zip(noise.day.psd[a],gooddays) if not good] for ax,a in zip(bottomaxes,keys)]
    _=[ax.set_title(f'{a} {mname}, dB',fontweight='bold') for ax,a in zip(bottomaxes,keys)]
    # _=[ax.set_xscale('log')for ax in bottomaxes]
    _=[ax.set_xlim(1/500,1)for ax in bottomaxes]
    _=[ax.axvline(fnotch(sta.StaDepth),linewidth=1,color='k',linestyle=':') for ax in bottomaxes]
    lg=bottomaxes[0].legend(handles=[staplot,dayplot],markerscale=5,ncols=2,fontsize=12,loc='upper center', prop={'weight': 'bold'})
    fig.suptitle(f'Noise | {stanm} ({sta.Experiment}) | {np.round(sta.StaDepth)}m',fontsize=20,fontweight='bold',y=0.91)
    [ax.set_xticks([1/500,1/100,1/10,1]) for ax in topaxes+bottomaxes]
    _=[ax.set_xscale('log')for ax in bottomaxes]
    _=[ax.set_xscale('log')for ax in topaxes]
    [a.set_ylim(-75,-275) for a in bottomaxes[:2]]
    [a.invert_yaxis() for a in bottomaxes[:2]]
    good = [days[i] for i,good in enumerate(gooddays) if good]
    bad = [days[i] for i,good in enumerate(gooddays) if not good]
    goodax.text(0.5,1.1,f'Good days ({len(good)})',color='b',horizontalalignment='center',verticalalignment='top',fontweight='bold')
    goodax.text(0.5,1,'\n'.join(good) if len(good) else 'None',color='b',verticalalignment='top',horizontalalignment='center')
    goodax.text(0.5,-.1,f'Bad days ({len(bad)})',color='r',horizontalalignment='center',verticalalignment='bottom',fontweight='bold')
    goodax.text(0.5,0,'\n'.join(bad) if len(bad) else 'None',color='r',verticalalignment='bottom',horizontalalignment='center')
    goodax.plot(0.5,0);goodax.plot(0.5,1)
    goodax.set_xticks([]);goodax.set_yticks([])
    goodax.set_xlim(0.3,0.7)
    return fig


get_statf = lambda stanm: load_pickle(list((dirs.TransferFunctions/stanm).glob('*-*.pkl'))[0])
get_daytfs = lambda stanm: [load_pickle(e) for e in [f for f in list((dirs.TransferFunctions/stanm).glob('*.pkl')) if not str(f).find('-')>=0]]
load_sta_atacrnoise = lambda stanm: load_pickle(list((dirs.SpectraAvg/stanm).glob('*.pkl'))[0])
def plot_TFs(stanm):
    statf = get_statf(stanm)
    daytf=get_daytfs(stanm)
    # daytf = [d for i,d in zip(load_sta_atacrnoise(stanm).gooddays,daytf) if i]
    statf.__dict__.keys()
    f = statf.f
    fn = statf.frequency_notch
    z = statf.station_depth
    statransfunc = statf.transfunc['ZP-21']
    daytransfunc = [tf.transfunc['ZP-21'] for tf in daytf]
    tfnames = list(statransfunc.keys())
    fig,axes = plt.subplots(len(tfnames),1,figsize=(16,(20/6)*len(tfnames)))
    grayalpha=0.4
    _=[[ax.loglog(f[(f>0)&(f<=1)],np.abs(dtf[k][(f>0)&(f<=1)]),color='lightgray' if good else 'r',alpha=grayalpha if good else 1) for dtf,good in zip(daytransfunc,load_sta_atacrnoise(stanm).gooddays)] for ax,k in zip(axes,tfnames)]
    _=[ax.loglog(f[(f>0)&(f<=1)],np.abs(statransfunc[k][(f>0)&(f<=1)]),color='blue') for ax,k in zip(axes,tfnames)]
    _=[ax.set_xlim(f[(f>0)&(f<=1)][0],f[(f>0)&(f<=1)][-1]) for ax in axes]
    _=[ax.set_title(f'{stanm} | {int(np.round(z))}m | Transfer Function: {k}') for k,ax in zip(tfnames,axes)]
    _=[ax.axvline(fn,linestyle=':',linewidth=0.4,color='k') for ax in axes]
    return fig
