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


# cat value
cat = catalog.r.copy()
# cat value
# cat = cat.loc[['7D.G17B','Z6.11']]
# cat = cat.loc[['7D.G30B']]
Stream_LogPSD = lt.math.Stream_LogPSD
# plotfold value
plotfold = dirs.P01.S09
# daysinfold value
daysinfold = lambda fo: ['.'.join(f.name.split('.')[:2]) for f in list(fo.glob('*HZ.SAC'))]
for si,sta in enumerate(cat.iloc):
    print(f'{si+1}/{len(cat)} : {sta.StaName}')
    # -----Define
    stanm=sta.StaName
    # q value
    q=daysinfold(dirs.Noise/'rmresp'/stanm/'_quarantine')
    # if len(q)==0:continue
    qfold = dirs.Noise/'rmresp'/stanm/'_quarantine'
    # efold value
    efold = dirs.Noise/'rmresp'/stanm/'extra.days'
    # gfold value
    gfold = dirs.Noise/'rmresp'/stanm
    # quarantined value
    quarantined = daysinfold(qfold)
    # extradays value
    extradays = daysinfold(efold)
    # gooddays value
    gooddays = daysinfold(gfold)
    # qdata value
    qdata = lambda quarantined=quarantined,qfold=qfold: Stream([load_sac(i) for i in lt.cat.unravel([fo for fo in [list(qfold.glob(f'*{f}*')) for f in quarantined] if len(fo)==4])])
    # edata value
    edata = lambda extradays=extradays,efold=efold: Stream([load_sac(i) for i in lt.cat.unravel([fo for fo in [list(efold.glob(f'*{f}*')) for f in extradays] if len(fo)==4])])
    # gdata value
    gdata = lambda gooddays=gooddays,gfold=gfold: Stream([load_sac(i) for i in lt.cat.unravel([fo for fo in [list(gfold.glob(f'*{f}*')) for f in gooddays] if len(fo)==4])])
    # -------
    # -------
    fig,axes = plt.subplots(4,3,figsize=(8,7))
    # col titles value
    col_titles = ['Good days','Quarantined days','Extra days']
    # sets value
    sets = [gdata,qdata,edata]
    for col,data in enumerate(sets):
        # -----Load
        st=data()
        if len(st)==0:
            for ci,ax in enumerate(axes[:,col]):
                # ndays value
                ndays=len(st)
                # comp value
                comp = comps[ci]
                if comp=='P':ax.set_ylim([-60,60]);ticks=[-50,-25,0,25,50]
                elif np.isin(comp,['1','2']):ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
                else:ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
                # ticks value
                ticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],5,dtype=int)
                ax.set_yticks(ticks)
                ax.set_yticklabels(ax.get_yticklabels())
                if ci==0:ax.set_title(f'{col_titles[col]} ({ndays})')
                if ci<3:ax.set_xticks([])
                if ci==3:ax.set_xlabel('Frequency (Hz)')
                if (col==0)&(ci==1):ax.set_ylabel('dB Power Density 10log10 $m^{2}/s^{4}$ / Hz')
                if (col==0)&(ci==3):ax.set_ylabel('dB Power Density \n 10log10 $Pa^{2}$ / Hz')
            continue
        ndays=int(len(st)/4)
        # smoothed value
        smoothed = True
        comps,days,f,t,psd1=Stream_LogPSD(st.select(channel='*1*'),smoothed=smoothed) # (comp, day, time, freq)
        # psd1 value
        psd1=psd1.squeeze().reshape((-1,len(f))).T
        comps,days,f,t,psd2=Stream_LogPSD(st.select(channel='*2*'),smoothed=smoothed) # (comp, day, time, freq)
        # psd2 value
        psd2=psd2.squeeze().reshape((-1,len(f))).T
        comps,days,f,t,psdZ=Stream_LogPSD(st.select(channel='*Z*'),smoothed=smoothed) # (comp, day, time, freq)
        # psdZ value
        psdZ=psdZ.squeeze().reshape((-1,len(f))).T
        comps,days,f,t,psdP=Stream_LogPSD(st.select(channel='*DH*'),smoothed=smoothed) # (comp, day, time, freq)
        # psdP value
        psdP=psdP.squeeze().reshape((-1,len(f))).T
        logpsds = [psd1,psd2,psdZ,psdP];comps=['1','2','Z','P']
        del st
        # -----
        for ci,(ax,comppsd) in enumerate(zip(axes[:,col],logpsds)):
            # -----Plot
            comp = comps[ci]
            comppsd=comppsd[:,~np.isinf(comppsd.sum(axis=0))]
            _=ax.semilogx(f,comppsd, 'k', lw=0.5)
            ax.set_xlim(1/150,5)
            if comp=='P':ax.set_ylim([-60,60]);ticks=[-50,-25,0,25,50]
            elif np.isin(comp,['1','2']):ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
            else:ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
            # ticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],5,dtype=int)
            ax.set_yticks(ticks)
            ax.set_yticklabels(ax.get_yticklabels())
            if ci==0:ax.set_title(f'{col_titles[col]} ({ndays})')
            if ci<3:ax.set_xticks([])
            if ci==3:ax.set_xlabel('Frequency (Hz)')
            if (col==0)&(ci==1):ax.set_ylabel('dB Power Density 10log10 $m^{2}/s^{4}$ / Hz')
            if (col==0)&(ci==3):ax.set_ylabel('dB Power Density \n 10log10 $Pa^{2}$ / Hz')
        if col==0:[ax.text(min(ax.get_xlim())+.001,min(ax.get_ylim()),c,verticalalignment='bottom',horizontalalignment='left') for c,ax in zip(comps,axes[:,col])]
        fig.suptitle(f'{sta.Experiment} {sta.StaName},{int(sta.StaDepth)}m ({np.round(fnotch(sta.StaDepth),2)}Hz)',y=0.94)
        file = f'{stanm}.Noise.QC.png'
        # file = file.replace('.png',f'.{note}.png')
    save_tight(plotfold/file)
    plt.close('all')



# ===========
# ===========


# for si,sta in enumerate(cat.iloc):
#     print(f'{si+1}/{len(cat)} : {sta.StaName}')
#     # -----Define
#     stanm=sta.StaName
#     q=daysinfold(dirs.Noise/'rmresp'/stanm/'_quarantine')
#     # if len(q)==0:continue
#     qfold = dirs.Noise/'rmresp'/stanm/'_quarantine'
#     efold = dirs.Noise/'rmresp'/stanm/'extra.days'
#     gfold = dirs.Noise/'rmresp'/stanm
#     quarantined = daysinfold(qfold)
#     extradays = daysinfold(efold)
#     gooddays = daysinfold(gfold)
#     qdata = lambda quarantined=quarantined,qfold=qfold: Stream([load_sac(i) for i in lt.cat.unravel([fo for fo in [list(qfold.glob(f'*{f}*')) for f in quarantined] if len(fo)==4])])
#     edata = lambda extradays=extradays,efold=efold: Stream([load_sac(i) for i in lt.cat.unravel([fo for fo in [list(efold.glob(f'*{f}*')) for f in extradays] if len(fo)==4])])
#     gdata = lambda gooddays=gooddays,gfold=gfold: Stream([load_sac(i) for i in lt.cat.unravel([fo for fo in [list(gfold.glob(f'*{f}*')) for f in gooddays] if len(fo)==4])])
#     # -------
#     # -------
#     fig,axes = plt.subplots(4,3,figsize=(8,7))
#     col_titles = ['Good days','Quarantined days','Extra days']
#     sets = [gdata,qdata,edata]
#     for col,data in enumerate(sets):
#         # -----Load
#         st=data()
#         if len(st)==0:
#             for ci,ax in enumerate(axes[:,col]):
#                 ndays=len(st)
#                 comp = comps[ci]
#                 if comp=='P':ax.set_ylim([-60,60]);ticks=[-50,-25,0,25,50]
#                 elif np.isin(comp,['1','2']):ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
#                 else:ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
#                 # ticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],5,dtype=int)
#                 ax.set_yticks(ticks)
#                 ax.set_yticklabels(ax.get_yticklabels(),fontsize=6)
#                 if ci==0:ax.set_title(f'{col_titles[col]} ({ndays})', fontdict={'fontsize': 9})
#                 if ci<3:ax.set_xticks([])
#                 if ci==3:ax.set_xlabel('Frequency (Hz)', fontdict={'fontsize': 9})
#                 if (col==0)&(ci==1):ax.set_ylabel('dB Power Density 10log10 $m^{2}/s^{4}$ / Hz')
#                 if (col==0)&(ci==3):ax.set_ylabel('dB Power Density \n 10log10 $Pa^{2}$ / Hz')
#             continue
#         times = [tr.stats.starttime for tr in st.select(channel='*Z')]
#         st0 = st.copy()
#         for tt in times:
#             st = st0.copy()
#             st = Stream([st[i] for i in np.where(np.array([tr.stats.starttime for tr in st])==tt)[0]]).copy()
#             note = st[0].stats.starttime.strftime('%Y.%j')
#             ndays=int(len(st)/4)
#             smoothed = True
#             comps,days,f,t,psd1=Stream_LogPSD(st.select(channel='*1*'),smoothed=smoothed) # (comp, day, time, freq)
#             psd1=psd1.squeeze().reshape((-1,len(f))).T
#             comps,days,f,t,psd2=Stream_LogPSD(st.select(channel='*2*'),smoothed=smoothed) # (comp, day, time, freq)
#             psd2=psd2.squeeze().reshape((-1,len(f))).T
#             comps,days,f,t,psdZ=Stream_LogPSD(st.select(channel='*Z*'),smoothed=smoothed) # (comp, day, time, freq)
#             psdZ=psdZ.squeeze().reshape((-1,len(f))).T
#             comps,days,f,t,psdP=Stream_LogPSD(st.select(channel='*DH*'),smoothed=smoothed) # (comp, day, time, freq)
#             psdP=psdP.squeeze().reshape((-1,len(f))).T
#             logpsds = [psd1,psd2,psdZ,psdP];comps=['1','2','Z','P']
#             del st
#             # -----
#             for ci,(ax,comppsd) in enumerate(zip(axes[:,col],logpsds)):
#                 # -----Plot
#                 comp = comps[ci]
#                 comppsd=comppsd[:,~np.isinf(comppsd.sum(axis=0))]
#                 _=ax.semilogx(f,comppsd, 'k', lw=0.5)
#                 ax.set_xlim(1/150,5)
#                 if comp=='P':ax.set_ylim([-60,60]);ticks=[-50,-25,0,25,50]
#                 elif np.isin(comp,['1','2']):ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
#                 else:ax.set_ylim([-250,10]);ticks=[-250,-150,-50,0]
#                 # ticks=np.linspace(ax.get_ylim()[0],ax.get_ylim()[1],5,dtype=int)
#                 ax.set_yticks(ticks)
#                 ax.set_yticklabels(ax.get_yticklabels(),fontsize=6)
#                 if ci==0:ax.set_title(f'{col_titles[col]} ({ndays})', fontdict={'fontsize': 9})
#                 if ci<3:ax.set_xticks([])
#                 if ci==3:ax.set_xlabel('Frequency (Hz)', fontdict={'fontsize': 9})
#                 if (col==0)&(ci==1):ax.set_ylabel('dB Power Density 10log10 $m^{2}/s^{4}$ / Hz')
#                 if (col==0)&(ci==3):ax.set_ylabel('dB Power Density \n 10log10 $Pa^{2}$ / Hz')
#             if col==0:[ax.text(min(ax.get_xlim())+.001,min(ax.get_ylim()),c,verticalalignment='bottom',horizontalalignment='left',fontsize=11) for c,ax in zip(comps,axes[:,col])]
#             fig.suptitle(f'{sta.Experiment} {sta.StaName},{int(sta.StaDepth)}m ({np.round(fnotch(sta.StaDepth),2)}Hz)',y=0.94)
#             file = f'{stanm}.Noise.QC.png'
#             file = file.replace('.png',f'.{note}.png')
#             save_tight(plotfold/file)
#     plt.close('all')