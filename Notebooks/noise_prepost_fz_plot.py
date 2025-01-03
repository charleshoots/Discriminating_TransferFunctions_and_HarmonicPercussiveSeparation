k=1
from imports import *
from modules import *
def unravel(lst):return list(itertools.chain.from_iterable(lst))
def get_noise(dirs,stanm):return load_pickle(list((dirs.SpectraAvg/stanm).glob('*sta.pkl'))[0])
def dataset_averaged_coherence_plot(f,z,coh,
    title='Station Averaged Noise',
    figsize=[15,6],fontsize=12,vlim=None,levels=None):
    font = {'weight':'bold','size':fontsize};matplotlib.rc('font', **font)
    z = np.round(z)
    i = np.argsort(z)
    z,coh = z[i],coh[i,:]
    if levels is None:levels=np.linspace(np.min(coh),np.max(coh),20)
    if vlim is None:vlim=[np.min(coh),np.max(coh)]
    fig,ax=plt.subplots(figsize=figsize)
    # coh_plottable = np.array([smooth(d,k=3) for d in coh]);coh_plottable = gaussian_filter(coh,.5)
    coh_plottable = coh
    cnt = ax.contourf(f,z,coh_plottable,vmin=vlim[0],vmax=vlim[1],extend="both",levels=levels)
    fn = [fnotch(fqz) for fqz in z]
    ax.plot(fn,z,linestyle='dashed',color='w',linewidth=0.4*fontsize)
    ax.set_xlim(1/500,1);ax.set_xscale('log')
    fticks = np.array([1/500,1/300,1/100,1/50,1/10,1])
    ax.set_xticks(fticks)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xticklabels(np.array(1/fticks,dtype=int))
    ax.set_ylabel('Station depth (m)',fontweight='bold')
    ax.set_xlabel('Period (s)',fontweight='bold')
    # ax.set_title(f'{coh.shape[0]} {title}')
    ax.set_facecolor('k')
    fig.suptitle(title)
    plt.colorbar(cnt)
    plt.tight_layout()
    return fig

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# ====================================================================================================================================
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

cat = catalog.copy();all_noise=[]
# for si,stanm in enumerate(cat.StaName):
#      print(f'{si+1}/{len(cat)} {stanm}')
#      all_noise.append(get_noise_spectra(dirs,stanm))

def slice_meter(all_noise,M='Coherence',R='ZP'):
    R = ''.join(sorted(R))
    RawCoh=[];stage='Raw'
    for a in all_noise:
        R=''.join(sorted(R))
        AB=a[stage][f'c{R[0]}{R[1]}']
        AA=a[stage][f'c{R[0]}{R[0]}']
        BB=a[stage][f'c{R[1]}{R[1]}']
        RawCoh.append(Meters(M,AB,AA,BB))
    RawCoh = np.array(RawCoh)
    CorrectedCoh=[];stage='Corrected'
    for a in all_noise:
        R=''.join(sorted(R))
        AB=a[stage][f'c{R[0]}{R[1]}']
        AA=a[stage][f'c{R[0]}{R[0]}']
        BB=a[stage][f'c{R[1]}{R[1]}']
        CorrectedCoh.append(Meters(M,AB,AA,BB))
    CorrectedCoh = np.array(CorrectedCoh)
    return RawCoh,CorrectedCoh

def _calc_coherence(ab,aa,bb):coh = ((abs(ab)**2)/abs(aa*bb));return coh
def _calc_phase(ab):ph = np.angle(ab,deg=True);return ph
def _calc_admittance(ab,bb):ad = np.abs(ab)/bb;return ad
def Meters(M,AB,AA,BB):
    if M=='Coherence':return _calc_coherence(AB,AA,BB)
    if M=='Phase':return _calc_phase(AB)
    if M=='Admittance':return _calc_admittance(AB,BB)
def NoiseMeter(all_noise,R='ZP',M='Coherence'):
    R = ''.join(sorted(R))
    Output=AttribDict()
    stage='Raw'
    Output.f=all_noise[0].f
    Output.f=Output.f[Output.f>=0]
    for a in all_noise:Output[a.stanm]=AttribDict()
    for a in all_noise:
        AB=a[stage][f'c{R[0]}{R[1]}']
        AA=a[stage][f'c{R[0]}{R[0]}']
        BB=a[stage][f'c{R[1]}{R[1]}']
        Output[a.stanm].f=a.f[a.f>=0]
        Output[a.stanm]['Uncorrected']=Meters(M,AB,AA,BB)[a.f>=0]
    stage='Corrected'
    for a in all_noise:
        AB=a[stage][f'c{R[0]}{R[1]}']
        AA=a[stage][f'c{R[0]}{R[0]}']
        BB=a[stage][f'c{R[1]}{R[1]}']
        Output[a.stanm]['Corrected']=Meters(M,AB,AA,BB)[a.f>=0]
    return Output
# ============
# XXXXXXXXXXXX
# ============


file = '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/Analysis/NetworkCoherences/AllNoise.pkl'
all_noise = load_pickle(file)

f = all_noise[0].f;depths=[a.depth for a in all_noise]
M = 'Coherence';vlim = [0,1];levels=np.linspace(0,1,40)

Comps = ['ZP','Z1','Z2']
# Comps = ['ZP']
for Comp in Comps:
    R = Comp
    RawCoh,CorrectedCoh = slice_meter(all_noise,R=R)
    Z = RawCoh;Z[Z<levels[1]]=0
    print(f'{Comp} |--Begin First--|')
    stage='Uncorrected';title = f'{stage} Station Averaged Noise {Comp} Coherence'
    fig=dataset_averaged_coherence_plot(f,depths,Z,vlim=vlim,title=title,levels=levels)
    save_tight(dirs.Plots/'_Plots'/f'Contour.{Comp}.{stage}.NoiseAverge.png',fig)
    print(f'{Comp} |--First Finished, Begin Second--|')
    Z=CorrectedCoh;Z[Z<levels[1]]=0
    stage='Corrected';title = f'{stage} Station Averaged Noise {Comp} Coherence'
    fig=dataset_averaged_coherence_plot(f,depths,Z,vlim=vlim,title=title,levels=levels)
    save_tight(dirs.Plots/'_Plots'/f'Contour.{Comp}.{stage}.NoiseAverge.png',fig)
    print(f'{Comp} |--Second Finished--|')
