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

import time;start=time.time()
# runtime value
runtime=lambda:int(time.time()-start)

# ======================================================================================================
# ======================================================================================================
# ---Options:
sigma = False
if sigma:stat=np.std
else:stat=np.mean
# oct av value
oct_av = True # ~3.4 minutes
# oct_av = False # ~2.19 minutes
notched = True
# modes = ['Z','Mag']
# modes = ['Mag']
# modes = ['Z','Noise','Mag']

# modes value
modes = ['Z']
# modes value
modes = ['Mag']
# modes = ['Noise']

# modes = ['Noise.vs.Z'] #do not combine this mode with others. code breaks.
noiselims=[-200,-20] #default: noiselims=[-180,-30]
# tallfigure=False;figsize=(6.5,6)
tallfigure=True;figsize=(4.5,6.5)

# pairs value
pairs = False
# ======================================================================================================
# ======================================================================================================
# 
if len(modes)==1:
    nrows,ncols = [3,1] if tallfigure else [1,3]
    if 'Noise.vs.Z' in modes:nrows,ncols=[2,1];figsize=[4,4]
    fig,axes=plt.subplots(nrows,ncols,figsize=figsize,sharex='all',sharey='row',constrained_layout=True)
else:
    nrows,ncols = [3,len(modes)] if tallfigure else [len(modes),3]
    fig,axes=plt.subplots(nrows,ncols,figsize=figsize,sharex='col' if tallfigure else 'row',sharey='col',constrained_layout=True)

# if len(modes)==1:axes=axes.T

# axes value
axes=np.atleast_2d(axes);col=-1
# axes value
axes=axes.reshape(3,-1)

# fold value
fold = dirs.P01.S12
if sigma:fold=fold/'sigma'

# Basic variables
cat = catalog.copy()
# SR value
SR = cat.sr.copy() #Source-receiver pairs
# s value
s=cat.r.iloc[0];f=s.Data.Coherence().f
# foct value
foct=octavg(s.Data.Coherence().ATaCR.zp_21.coh,f)[0]
# nets value
nets = cat.r.Network.unique()
# snm value
snm=cat.r.StaName.unique()

# Plot vars
sz = 80;alpha=.02
# psd ttl value
psd_ttl=lambda:rf'10log10$(m^{2}/s^{4}/Hz)$ dB'
# yttl value
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# sigma yttl value
sigma_yttl = lambda c:fr"$\sigma(\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}})$"
if sigma:yttl=sigma_yttl
# lambda fnotch value
lambda_fnotch = fnotch if notched else lambda z:1
# opts value
opts = '.'.join(['octav' if oct_av else '','notched' if notched else '','tall' if tallfigure else 'wide']).replace('..','.')
# Noise
s=SR.iloc[0];c='cZZ';f=s.Data.Noise.Averaged().f;faxis=f>0;f=f[faxis]
# noise f value
noise_f=f
# c value
c='cZZ';Znoise={s.StaName:np.atleast_2d( PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__[c][faxis]) ) for s in cat.r.iloc}
# c value
c='H';Hnoise={s.StaName:np.atleast_2d( np.mean([PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__[c][faxis]) for c in ['c11','c22']],axis=0) ) for s in cat.r.iloc}
# f value
f=SR.iloc[0].Data.Coherence().f
# noise xf value
noise_xf=noise_f;xf=f

if oct_av: #Adds about 1-min more to the run time.
    print('Running octave averaging...')
    # new noise xf value
    new_noise_xf=octavg(Znoise[list(Znoise.keys())[0]],noise_xf)[0]
    # new xf value
    new_xf=octavg(SR.iloc[0].Data.Coherence().ATaCR.zp_21.coh,xf)[0]
    # loop over iloc
    for s in SR.iloc:
        usnr=unpack_metrics(SR)
        xf=1/usnr.coh.TF_Z.bands.copy()
        coh=usnr.coh.TF_Z.D
        s.Coherence.TF=octavg(usnr.coh.TF_Z.D,xf)[1]
        s.Coherence.HPS_Z=octavg(usnr.coh.HPS_Z.D,xf)[1]
        s.Coherence.HPS_H=octavg(usnr.coh.HPS_1.D,xf)[1]
    for k in Znoise.keys():
        Znoise[k]=octavg(Znoise[k],noise_xf)[1]
        Hnoise[k]=octavg(Hnoise[k],noise_xf)[1]
    print(f'...complete ({np.round(start-time.time(),1)}s)')
if oct_av:
    # noise xf value
    noise_xf=new_noise_xf
    # xf value
    xf=new_xf


# yl value
yl=[-0.02,1.02] if not sigma else [-.02,.5]
if np.isin('Noise',modes):
    col+=1;caxes=axes[:,col]
    # ROW-1 ||||||||||||||||||||||||||||||||||||||||
    # ---------
    # # ------------------------------------------------------------------------------------------------
    # # [Noise vs TF-Z]
    xl=noiselims

    print('|| Noise vs Coherence ||' + ' (Octave averaged)' if oct_av else '')
    ax=caxes[0]
    print(f'AXES:{0}')
    # variable
    _=[ax.scatter(
    Znoise[s.StaName][:,noise_xf<=lambda_fnotch(s.StaDepth)].mean(),
    stat(s.Coherence.TF[:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0),
    # c value
    c=ColorStandard.instrument[s.Instrument_Design],marker=ColorStandard.seismometer_marker[s.Seismometer],alpha=alpha) for s in SR.iloc]
    # variable
    _=[ax.scatter(
    Znoise[n][:,noise_xf<=lambda_fnotch(cat.r.loc[n].StaDepth)].mean(),
    stat(np.vstack([i.Coherence.TF[:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean(),
    # c value
    c=ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]],
    # marker value
    marker=ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for n in snm]
    if not tallfigure:ax.set_xlabel('Z',fontweight='bold',fontsize=11)
    ax.set_ylabel(yttl('TF.Z') if not sigma else sigma_yttl('TF.Z'),fontweight='bold',fontsize=11)
    ax.set_xlim(xl);ax.set_ylim(yl)
    # # # ------------------------------------------------------------------------------------------------
    # # # [Noise vs HPS-Z]
    ax=caxes[1]
    print(f'AXES:{1}')
    # variable
    _=[ax.scatter(
    Znoise[s.StaName][:,noise_xf<=lambda_fnotch(s.StaDepth)].mean(),
    stat(s.Coherence.HPS_Z[:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0),
    # c value
    c=ColorStandard.instrument[s.Instrument_Design],marker=ColorStandard.seismometer_marker[s.Seismometer],alpha=alpha) for s in SR.iloc]
    # variable
    _=[ax.scatter(
    Znoise[n][:,noise_xf<=lambda_fnotch(cat.r.loc[n].StaDepth)].mean(),
    stat(np.vstack([i.Coherence.HPS_Z[:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean(),
    # c value
    c=ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]],
    # marker value
    marker=ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for n in snm]
    ax.set_ylabel(yttl('HPS.Z') if not sigma else sigma_yttl('HPS.Z'),fontweight='bold',fontsize=11)
    # xlbl value
    xlbl = 'Noise Power spectral density,\n'+psd_ttl() if tallfigure else 'Z\nNoise Power spectral density'+psd_ttl()
    if not tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    ax.set_xlim(xl);ax.set_ylim(yl)
    # # ------------------------------------------------------------------------------------------------
    # # [Noise vs HPS-H]
    ax=caxes[2]
    print(f'AXES:{2}')
    # variable
    _=[ax.scatter(
    Hnoise[s.StaName][:,noise_xf<=lambda_fnotch(s.StaDepth)].mean(),
    stat(s.Coherence.HPS_H[:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0),
    # c value
    c=ColorStandard.instrument[s.Instrument_Design],marker=ColorStandard.seismometer_marker[s.Seismometer],alpha=alpha) for s in SR.iloc]
    # variable
    _=[ax.scatter(
    Hnoise[n][:,noise_xf<=lambda_fnotch(cat.r.loc[n].StaDepth)].mean(),
    stat(np.vstack([i.Coherence.HPS_H[:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean(),
    # c value
    c=ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]],
    # marker value
    marker=ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for n in snm] 
    if tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    else:ax.set_xlabel('H',fontweight='bold',fontsize=11)
    ax.set_ylabel(yttl('HPS.H'),fontweight='bold',fontsize=11)
    ax.set_xlim(xl);ax.set_ylim(yl)

    print(f'....Done ({(runtime())}s)')
    if len(modes)==1:
        # file value
        file=fold/f'Noise.vs.Coherence.{opts}.png'
        save_tight(file,fig,dpi=600);print('-Saved-')


if 'Noise.vs.Z' in modes:
    col+=1;caxes=axes[:,col]
    # ROW-1 ||||||||||||||||||||||||||||||||||||||||
    # ---------
    # # ------------------------------------------------------------------------------------------------
    # # [Noise vs TF-Z]
    ylbl='Water depth, m'
    # xlbl value
    xlbl = 'Noise Power spectral density,\n'+psd_ttl() if tallfigure else 'Z\nNoise Power spectral density'+psd_ttl()
    # swap axes value
    swap_axes = True
    # yl value
    yl = [-200,6200]
    # xl value
    xl=noiselims
    if swap_axes:xl,yl = yl,xl;xlbl,ylbl = ylbl,xlbl
    print('|| Noise vs Z ||' + ' (Octave averaged)' if oct_av else '')
    ax=caxes[0]
    print(f'AXES:{0}')
    x=[Znoise[n][:,noise_xf<=lambda_fnotch(cat.r.loc[n].StaDepth)].mean() for n in snm];y=[cat.r.loc[n].StaDepth for n in snm]
    if swap_axes:x,y=y,x
    _=[ax.scatter(xx,yy,
    c=ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]],
    marker=ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for xx,yy,n in zip(x,y,snm)]
    ax.set_ylabel(ylbl,fontweight='bold',fontsize=11,y=0 if swap_axes else 0.5)
    ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    # # ------------------------------------------------------------------------------------------------
    # # [Noise vs HPS-H]
    ax=caxes[1]
    print(f'AXES:{1}')
    x=[Hnoise[n][:,noise_xf<=lambda_fnotch(cat.r.loc[n].StaDepth)].mean() for n in snm];y=[cat.r.loc[n].StaDepth for n in snm]
    if swap_axes:x,y=y,x
    _=[ax.scatter(xx,yy,
    c=ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]],
    marker=ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for xx,yy,n in zip(x,y,snm)]
    ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    if not swap_axes:ax.set_ylabel(ylbl,fontweight='bold',fontsize=11)
    ax.set_xlim(xl);ax.set_ylim(yl)
    print(f'....Done ({(runtime())}s)')
    if len(modes)==1:
        file=fold/f'Noise.vs.Z.{opts + ('.axesswapped' if swap_axes else '')}.png'
        save_tight(file,fig,dpi=600);print('-Saved-')


req='StaDepth';xlbl='Water depth, m';xl=[0,SR[req].max()+300];ftitle=f'Z.vs.Coherence.{opts}.png'
if np.isin('Z',modes):
    col+=1;caxes=axes[:,col]
    yl=[-0.02,1.02] if not sigma else [-.02,.5]
    # ---------
    # # ------------------------------------------------------------------------------------------------
    # print(f'|| {req} vs Coherence ||' + ' (Octave averaged)' if oct_av else '')
    # for ci,(method,tl) in enumerate(zip(['TF','HPS_Z','HPS_H'],['TF.Z','HPS.Z','HPS.H'])):
    #     print(f'AXES:{ci}')
    #     ax=caxes[ci]
    #     c=[ColorStandard.instrument[s.Instrument_Design] for s in SR.iloc]
    #     m=[ColorStandard.seismometer_marker[s.Seismometer] for s in SR.iloc]
    #     x=SR[req]
    #     y=[stat(s.Coherence[method][:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0) for s in SR.iloc]
    #     _=[ax.scatter(xx,yy,c=cc,alpha=alpha,marker=mm) for xx,yy,cc,mm in zip(x,y,c,m)]
    #     # ----
    #     c=[ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]] for n in snm]
    #     m=[ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]] for n in snm]
    #     x=[SR[SR.StaName==n][req].mean() for n in snm]
    #     y=[stat(np.vstack([i.Coherence[method][:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean() for n in snm]
    #     _=[ax.scatter(xx,yy,c=cc,marker=mm,alpha=1,s=sz,edgecolor='k') for xx,yy,cc,mm in zip(x,y,c,m)]
    #     # ---
    #     temp_ttl = yttl(tl)
    #     ax.set_ylabel(temp_ttl,fontweight='bold',fontsize=11)
    #     ax.set_xlim(xl);ax.set_ylim(yl)
    #     if not tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    # print(f'....Done ({(runtime())}s)')
    # if len(modes)==1:file=fold/ftitle;save_tight(file,fig,dpi=600);print('-Saved-')

    for ci,(method,tl) in enumerate(zip(['TF','HPS_Z','HPS_H'],['TF.Z','HPS.Z','HPS.H'])):
        print(f'AXES:{ci}')
        ax=caxes[ci]
        c=[ColorStandard.instrument[s.Instrument_Design] for s in SR.iloc]
        m=[ColorStandard.seismometer_marker[s.Seismometer] for s in SR.iloc]
        x=SR[req]
        y=np.array([stat(s.Coherence[method][:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0) for s in SR.iloc])
        if ci==0:b=np.array([stat(s.Coherence['HPS_Z'][:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0)  for s in SR.iloc]);y=(y-b)/b

        _=[ax.scatter(xx,yy,c=cc,alpha=alpha,marker=mm) for xx,yy,cc,mm in zip(x,y,c,m)]
        # ----
        c=[ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]] for n in snm]
        m=[ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]] for n in snm]
        x=[SR[SR.StaName==n][req].mean() for n in snm]
        y=np.array([stat(np.vstack([i.Coherence[method][:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean() for n in snm])
        if ci==0:b=np.array([stat(np.vstack([i.Coherence['HPS_Z'][:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean() for n in snm]);y=(y-b)/b
        _=[ax.scatter(xx,yy,c=cc,marker=mm,alpha=1,s=sz,edgecolor='k') for xx,yy,cc,mm in zip(x,y,c,m)]
        # ---
        temp_ttl = yttl(tl)
        if ci==0:temp_ttl=yttl('HPS.Z')+'/'+temp_ttl
        ax.set_ylabel(temp_ttl,fontweight='bold',fontsize=11)
        if ci==0:ax.set_xlim(xl);ax.set_ylim([-1,1]);ax.set_yticks(np.arange(-1,1.1,.1))
        else:
            ax.set_xlim(xl);ax.set_ylim(yl)
        if not tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    print(f'....Done ({(runtime())}s)')
    if len(modes)==1:file=fold/('RelRatios.'+ftitle);save_tight(file,fig,dpi=600);print('-Saved-')



# ------------------------------------------------------------------------------------------------
req='Magnitude';xlbl='Magnitude, Mw';xl=[5.9,8.0];ftitle=f'Magnitudes.vs.Coherence.{opts}.png'
if np.isin('Mag',modes):
    col+=1;caxes=axes[:,col]
    yl=[-0.02,1.02] if not sigma else [-.02,.5]
    # ---------
    # # ------------------------------------------------------------------------------------------------
    print(f'|| {req} vs Coherence ||' + ' (Octave averaged)' if oct_av else '')
    # for ci,(method,tl) in enumerate(zip(['TF','HPS_Z','HPS_H'],['TF.Z','HPS.Z','HPS.H'])):
    #     print(f'AXES:{ci}')
    #     ax=caxes[ci]
    #     c=[ColorStandard.instrument[s.Instrument_Design] for s in SR.iloc]
    #     m=[ColorStandard.seismometer_marker[s.Seismometer] for s in SR.iloc]
    #     x=SR[req]
    #     y=np.array([stat(s.Coherence[method][:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0) for s in SR.iloc])
    #     _=[ax.scatter(xx,yy,c=cc,alpha=alpha,marker=mm) for xx,yy,cc,mm in zip(x,y,c,m)]
    #     # ----
    #     c=[ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]] for n in snm]
    #     m=[ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]] for n in snm]
    #     x=[SR[SR.StaName==n][req].mean() for n in snm]
    #     y=[stat(np.vstack([i.Coherence[method][:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean() for n in snm]
    #     _=[ax.scatter(xx,yy,c=cc,marker=mm,alpha=1,s=sz,edgecolor='k') for xx,yy,cc,mm in zip(x,y,c,m)]
    #     # --
    #     temp_ttl = yttl(tl)
    #     if ci==0:temp_ttl=yttl('HPS.Z')+'/'+temp_ttl
    #     ax.set_ylabel(temp_ttl,fontweight='bold',fontsize=11)
    #     ax.set_xlim(xl);ax.set_ylim(yl)
    #     if not tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    # print(f'....Done ({(runtime())}s)')
    # if len(modes)==1:file=fold/ftitle;save_tight(file,fig,dpi=600);print('-Saved-')

    for ci,(method,tl) in enumerate(zip(['TF','HPS_Z','HPS_H'],['TF.Z','HPS.Z','HPS.H'])):
        print(f'AXES:{ci}')
        ax=caxes[ci]
        c=[ColorStandard.instrument[s.Instrument_Design] for s in SR.iloc]
        m=[ColorStandard.seismometer_marker[s.Seismometer] for s in SR.iloc]
        x=SR[req]
        y=np.array([stat(s.Coherence[method][:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0) for s in SR.iloc])
        if ci==0:b=np.array([stat(s.Coherence['HPS_Z'][:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0)  for s in SR.iloc]);y=(y-b)/b
        _=[ax.scatter(xx,yy,c=cc,alpha=alpha,marker=mm) for xx,yy,cc,mm in zip(x,y,c,m)]
        # ------------------------------------------------------------------------------------------------------------------------------------
        c=[ColorStandard.instrument[SR[SR.StaName==n].Instrument_Design[0]] for n in snm]
        m=[ColorStandard.seismometer_marker[SR[SR.StaName==n].Seismometer[0]] for n in snm]
        x=[SR[SR.StaName==n][req].mean() for n in snm]
        y=np.array([stat(np.vstack([i.Coherence[method][:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean() for n in snm])
        if ci==0:b=np.array([stat(np.vstack([i.Coherence['HPS_Z'][:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in SR[SR.StaName==n].iloc]),axis=0).mean() for n in snm]);y=(y-b)/b
        _=[ax.scatter(xx,yy,c=cc,marker=mm,alpha=1,s=sz,edgecolor='k') for xx,yy,cc,mm in zip(x,y,c,m)]
        # ------------------------------------------------------------------------------------------------------------------------------------
        temp_ttl = yttl(tl)
        if ci==0:temp_ttl=yttl('HPS.Z')+'/'+temp_ttl
        ax.set_ylabel(temp_ttl,fontweight='bold',fontsize=11)
        if ci==0:ax.set_xlim(xl);ax.set_ylim([-1,1])
        else:ax.set_xlim(xl);ax.set_ylim(yl)
        if not tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
    print(f'....Done ({(runtime())}s)')
    if len(modes)==1:file=fold/('RelRatios.'+ftitle);save_tight(file,fig,dpi=600);print('-Saved-')




# --------------- A BINNED version of scatter plots by event magnitude
# if 'Mag' in modes:
#     col+=1;caxes=axes[:,col]
#     dmag=.3
#     # dmag=0.5
#     mags=np.array([[i,i+dmag] for i in np.arange(6,8,dmag)])
#     yl=[-0.02,1.02] if not sigma else [-.02,.5]
#     xl,yl=yl,xl

#     for m in mags:
#         print(f'mags: {m}')
#         # iSR=SR[SR[req]==m].copy()
#         iSR=SR[(SR[req]>=min(m))&(SR[req]<=max(m))].copy()
#         snm = iSR.StaName.unique()
#         # ROW-1 ||||||||||||||||||||||||||||||||||||||||
#         # ---------
#         # # ------------------------------------------------------------------------------------------------
#         # # [Z vs TF-Z]
#         # # ------------------------------------------------------------------------------------------------
#         print('|| Z vs Coherence ||' + ' (Octave averaged)' if oct_av else '')
#         ax=caxes[0];print(f'AXES:{0}')
#         if pairs:
#             _=[ax.scatter(
#             y=s[req],
#             x=stat(s.Coherence.TF[:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0),
#             c=ColorStandard.instrument[s.Instrument_Design],marker=ColorStandard.seismometer_marker[s.Seismometer],alpha=alpha) for s in iSR.iloc]
#         _=[ax.scatter(
#         y=iSR[iSR.StaName==n][req].mean(),
#         x=stat(np.vstack([i.Coherence.TF[:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in iSR[iSR.StaName==n].iloc]),axis=0).mean(),
#         c=ColorStandard.instrument[iSR[iSR.StaName==n].Instrument_Design[0]],
#         marker=ColorStandard.seismometer_marker[iSR[iSR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for n in snm]
#         # ax.set_xlabel('Z',fontweight='bold',fontsize=11)
#         ax.set_xlabel(yttl('TF.Z') if not sigma else sigma_yttl('TF.Z'),fontweight='bold',fontsize=11)
#         ax.set_xlim(xl);ax.set_ylim(yl)
#         # # # ------------------------------------------------------------------------------------------------
#         # # # [Z vs HPS-Z]
#         # # # ------------------------------------------------------------------------------------------------

#         ax=caxes[1];print(f'AXES:{1}')
#         if pairs:
#             _=[ax.scatter(
#             y=s[req],
#             x=stat(s.Coherence.HPS_Z[:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0),
#             c=ColorStandard.instrument[s.Instrument_Design],marker=ColorStandard.seismometer_marker[s.Seismometer],alpha=alpha) for s in iSR.iloc]
#         _=[ax.scatter(
#         y=iSR[iSR.StaName==n][req].mean(),
#         x=stat(np.vstack([i.Coherence.HPS_Z[:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in iSR[iSR.StaName==n].iloc]),axis=0).mean(),
#         c=ColorStandard.instrument[iSR[iSR.StaName==n].Instrument_Design[0]],
#         marker=ColorStandard.seismometer_marker[iSR[iSR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for n in snm]
#         ax.set_xlabel(yttl('HPS.Z') if not sigma else sigma_yttl('HPS.Z'),fontweight='bold',fontsize=11)
#         if not tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
#         else:ax.set_ylabel(xlbl,fontweight='bold',fontsize=11)
#         ax.set_xlim(xl);ax.set_ylim(yl)
#         # # ------------------------------------------------------------------------------------------------
#         # # [Z vs HPS-H] | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | | |
#         # # ------------------------------------------------------------------------------------------------
#         ax=caxes[2];print(f'AXES:{2}')
#         if pairs:
#             _=[ax.scatter(
#             y=s[req],
#             x=stat(s.Coherence.HPS_H[:,xf<=lambda_fnotch(s.StaDepth)],axis=1).mean(axis=0),
#             c=ColorStandard.instrument[s.Instrument_Design],marker=ColorStandard.seismometer_marker[s.Seismometer],alpha=alpha) for si,s in enumerate(iSR.iloc)]
#         _=[ax.scatter(
#         y=iSR[iSR.StaName==n][req].mean(),
#         x=stat(np.vstack([i.Coherence.HPS_H[:,xf<=lambda_fnotch(cat.r.loc[n].StaDepth)] for i in iSR[iSR.StaName==n].iloc]),axis=0).mean(),
#         c=ColorStandard.instrument[iSR[iSR.StaName==n].Instrument_Design[0]],
#         marker=ColorStandard.seismometer_marker[iSR[iSR.StaName==n].Seismometer[0]],alpha=1,s=sz,edgecolor='k') for n in snm]
#         ax.set_xlabel(yttl('HPS.H') if not sigma else sigma_yttl('HPS.H'),fontweight='bold',fontsize=11)
#         # if tallfigure:ax.set_xlabel(xlbl,fontweight='bold',fontsize=11)
#         ax.set_xlim(xl);ax.set_ylim(yl)
#     print(f'....Done ({(runtime())}s)')
#     if len(modes)==1:
#         file=fold/ftitle;save_tight(file,fig,dpi=600);print('-Saved-')



if len(modes)>1:
    file=fold/f'Combined.Noise.{'.'.join(modes)}.Coherence.{opts}.png'
    save_tight(file,fig,dpi=600);print('-Saved-')
print(f"Elapsed time: {(runtime())/60:.2f} minutes")