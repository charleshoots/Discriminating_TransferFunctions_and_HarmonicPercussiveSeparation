from imports import *



# # Z-Detector:
# Swindell, W. H., and N. S. Snell. "Station processor automatic signal detection system, 
# phase I: Final report, station processor software development." 
# Texas Instruments Report 
# AFTAC Contract Number FO8606-76-C-0025, 
# Texas Instruments Incorporated, Dallas, Texas (1977).

# Recursive STA/LTA:
# [Withers1998]
# Withers, M., Aster, R., Young, C., Beiriger, J., Harris, M., Moore, S., and Trujillo, J. (1998),
# A comparison of select trigger algorithms for automated global seismic phase and event detection,
# Bulletin of the Seismological Society of America, 88 (1), 95-106.
# http://www.bssaonline.org/content/88/1/95.abstract

# Delayed STA/LTA:
# [Trnkoczy2012]
# Trnkoczy, A. (2012),Understanding and parameter setting of STA/LTA trigger algorithm,
# in New Manual of Seismological Observatory Practice 2 (NMSOP-2), IS 8.1, 20 pp.
# http://nmsop.gfz-potsdam.de

# ------------------------------------------------------------------------------------------------------------------------



# def EnsembleSNR(A,B,C,nsta=60,nlta=300,ratio=1,quiet=1,fband=[1,500]):
#     A.filter('bandpass',freqmin=1/fband[1],freqmax=1/fband[0],zerophase=True,corners=4)
#     B.filter('bandpass',freqmin=1/fband[1],freqmax=1/fband[0],zerophase=True,corners=4)
#     C.filter('bandpass',freqmin=1/fband[1],freqmax=1/fband[0],zerophase=True,corners=4)
#     n=min([A[0].data.shape[0],B[0].data.shape[0],C[0].data.shape[0]])
#     A[0].data=A[0].data[:n];B[0].data=B[0].data[:n];C[0].data=C[0].data[:n]
#     fs=A[0].stats.sampling_rate
#     nsta,nlta=int(nsta*fs),int(nlta*fs)
#     algorithms = ['classicstalta','delayedstalta','recstalta','zdetect',]
#     Ascores=[];Bscores=[];Cscores=[]
#     for algo in algorithms:
#         if algo=='carlstatrig':Ascores.append(A.copy().trigger('carlstatrig',nsta=nsta,nlta=nlta,ratio=ratio,quiet=quiet)[0].data)
#         if algo=='classicstalta':Ascores.append(A.copy().trigger('classicstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='delayedstalta':Ascores.append(A.copy().trigger('delayedstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='recstalta':Ascores.append(A.copy().trigger('recstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='zdetect':Ascores.append(A.copy().trigger('zdetect',nsta=nsta)[0].data)
#         if algo=='carlstatrig':Bscores.append(B.copy().trigger('carlstatrig',nsta=nsta,nlta=nlta,ratio=ratio,quiet=quiet)[0].data)
#         if algo=='classicstalta':Bscores.append(B.copy().trigger('classicstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='delayedstalta':Bscores.append(B.copy().trigger('delayedstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='recstalta':Bscores.append(B.copy().trigger('recstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='zdetect':Bscores.append(B.copy().trigger('zdetect',nsta=nsta)[0].data)
#         if algo=='carlstatrig':Cscores.append(C.copy().trigger('carlstatrig',nsta=nsta,nlta=nlta,ratio=ratio,quiet=quiet)[0].data)
#         if algo=='classicstalta':Cscores.append(C.copy().trigger('classicstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='delayedstalta':Cscores.append(C.copy().trigger('delayedstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='recstalta':Cscores.append(C.copy().trigger('recstalta',nsta=nsta,nlta=nlta)[0].data)
#         if algo=='zdetect':Cscores.append(C.copy().trigger('zdetect',nsta=nsta)[0].data)
#     Ascores=np.array(Ascores)
#     Bscores=np.array(Bscores)
#     Cscores=np.array(Cscores)
#     peaks=[(Ascores[i,:].argmax(),Bscores[i,:].argmax(),Cscores[i,:].argmax()) for i in range(Ascores.shape[0])]

#     # Each traces average of the two peak positions found in each test.
#     Apeaks=np.array([Ascores[pi,p] for pi,p in enumerate(peaks)]).mean(axis=1)
#     Bpeaks=np.array([Bscores[pi,p] for pi,p in enumerate(peaks)]).mean(axis=1)
#     Cpeaks=np.array([Cscores[pi,p] for pi,p in enumerate(peaks)]).mean(axis=1)
#     # Just the magnitudes to keep it simple
#     # Apeaks = 10*np.log10(Apeaks);Bpeaks = 10*np.log10(Bpeaks)
#     Apeaks=(Apeaks**2).mean()**0.5
#     Bpeaks=(Bpeaks**2).mean()**0.5
#     Cpeaks=(Cpeaks**2).mean()**0.5
#     score = [Apeaks,Bpeaks,Cpeaks]
#     #positive values indicate better snr in trace A.
#     #negative values indicate better snr in trace B.
#     # score=(Apeaks/Bpeaks)-1
#     # score=10*np.log10(Apeaks/Bpeaks)
#     # score=Apeaks/Bpeaks
#     # score=[score.max(),score.min()] #The overall best A and B scores, respectively.
#     # score = score.mean()
#     #---Final scores represents how many TIMES greater SNR across the ensemble models 
#     # ----A is then B (positives) or B is than A (negatives).
#     # -----The best A and the best B scores are returned.
#     return score
# def mirror(afold,bfold,events,comp='HZ'):
#     mirrored=[ev for ev in events if 
#     (len(list(afold.glob(f'*{ev.Name}*{comp}.SAC')))>0) 
#     and (len(list(bfold.glob(f'*{ev.Name}*{comp}.SAC')))>0)]
#     return Catalog(mirrored)
# cat = catalog.copy()
# # cat=cat.iloc[13:].copy()
# score_hold=AttribDict()
# for si,stanm in enumerate(cat.StaName):
#     # ---------------------------------------------------------------------------------------
#     # ---------------------------------------------------------------------------------------
#     events=cat.loc[stanm].Events
#     afold = dirs.Events/'corrected'/stanm
#     bfold = dirs.Events_HPS/'corrected'/stanm
#     cfold = dirs.Events/'rmresp'/stanm
#     events=mirror(afold,bfold,events)

#     scores=[]
#     for evi,ev in enumerate(events):
#         clear_output()
#         print(f'{stanm} | {si+1}/{len(cat)} : '+'|'*(si+1)+'<')
#         print(f'{evi+1}/{len(events)}'+'-'*(evi+1)+'<')
#         A=load_sac(afold/f'{stanm}.{ev.Name}.day.ZP-21.HZ.SAC')[0]
#         B=load_sac(bfold/f'{stanm}.{ev.Name}.HZ.SAC')[0]
#         C=load_sac(cfold/f'{ev.Name}.HZ.SAC')[0]
#         scores.append([EnsembleSNR(A,B,C,nsta=60,nlta=300,ratio=1,quiet=1)])
#         del A,B
#     score_hold[stanm]=np.array(scores)

# write_pickle(dirs.Analysis/'scores_1_500s_simple.pkl',score_hold)


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------

file='scores_1_500s_simple.pkl'
# file='scores_1_500s_linear_ratio.pkl'
# file='scores_1_500s_log_ratio_algo_avg.pkl'
score_hold=load_pickle(dirs.Analysis/file)
# ---------------------------------------------------------------
# ---------------------------------------------------------------
cat = catalog.copy()
cat=cat.sort_values(by='StaDepth')
nvariance=1
snr_stds=nvariance*np.array([score_hold[stanm].std(axis=0) for stanm in cat.StaName]).T**0.5
snr_means=np.array([score_hold[stanm].mean(axis=0) for stanm in cat.StaName]).T
fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(25,5));axes=np.atleast_1d(axes).T
ax=axes[0]
x = np.array([i for i in range(len(cat))])
k=2
colors =np.array(['r','b','k']);labels=['ATaCR','Noisecut','Raw']
for si,stanm in enumerate(cat.StaName):
    atacr_snr=score_hold[stanm][:,0]
    noisecut_snr=score_hold[stanm][:,1]
    raw_snr=score_hold[stanm][:,2]
    ax.scatter(raw_snr*0+si,raw_snr,s=2,c='k')
    ax.scatter(atacr_snr*0+si,atacr_snr,s=1,c='r')
    ax.scatter(noisecut_snr*0+si,noisecut_snr,s=0.5,c='b')
alphas=[0.3,0.3,0.1];zorder=[100,50,1]
_=[ax.fill_between(x,smooth(m,k=k),smooth(m-s,k=k), alpha=a,color=c,zorder=z) for z,a,m,s,c,l in zip(zorder,alphas,snr_means,snr_stds,colors,labels)]
_=[ax.fill_between(x,smooth(m,k=k),smooth(m+s,k=k), alpha=a,color=c,zorder=z) for z,a,m,s,c,l in zip(zorder,alphas,snr_means,snr_stds,colors,labels)]
_=[ax.plot(x,smooth(m,k=k),linestyle='dashdot',linewidth=lw,color=c,alpha=0.5,zorder=z) for z,lw,m,s,c,l in zip(zorder,[1,1,2],snr_means,snr_stds,colors,labels)]
# This is just for the legend marker
_=[ax.plot(x,m-1000,linestyle='dashdot',linewidth=lw,label=l,color=c,alpha=1,zorder=z) for z,lw,m,s,c,l in zip(zorder,[4,4,4],snr_means,snr_stds,colors,labels)]
ax.set_xlim(x[0],x[-1])
ax.set_ylim(bottom=1,top=20)
ax.set_ylabel('SNR',fontweight='bold',fontsize=13)
ax.legend(loc='lower left',ncols=3)
labels=[int(a) for a in np.round([cat.iloc[int(i)].StaDepth for i in np.linspace(0,len(cat)-1,ax.get_xticks().shape[0],dtype=int)])]
ax.set_xticklabels(labels,fontsize=11)
ax.set_xlabel('Station Depth, meters',fontweight='bold',fontsize=11)
# ---------------------------------------------------------------
save_tight(dirs.Plots/'simple_snr.png',fig=fig,dpi=800)
# ---------------------------------------------------------------