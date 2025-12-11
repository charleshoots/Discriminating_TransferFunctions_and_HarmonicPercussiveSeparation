### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

# --Imports
import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

from matplotlib.colors import LinearSegmentedColormap
# Noise Spectra
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
# f value
f=cat.r.iloc[0].Data.Noise.Averaged().f;faxis=(f>0)&(f<=1);f=f[faxis];noise_f=f
cat.r['Noise']=[AttribDict({'f':f,
'Z':PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__['cZZ'][faxis]),
'P':PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__['cPP'][faxis]),
'H':np.mean([PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__[c][faxis]) for c in ['c11','c22']],axis=0)
}) for s in cat.r.iloc]
# figs value
figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 #Darkness on a scale from 0 (black) to 1 (white)
# centers value
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
# dbin value
dbin = lambda x:np.array([x[:-1],x[1:]]).T
# yttl value
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# yttl eta value
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"



# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure10_Coherence.with.HPS.Traces.Magnitude';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'

# ---Data
magwins=dbin(np.arange(6,8.25,.25))
# icat value
icat=catalog.sr.copy()
usnr=unpack_metrics(icat)

# data_audit = np.array([sr.Traces()==None for sr in icat.iloc])


sr=icat.iloc[0];st=sr.Traces()
sn=Signal(st.select(location='*NoiseCut')[0],st.select(location='*ATaCR')[0])
sn_f=sn.coherence()[0];faxis=(f>=(1/200))&(f<=1);f=f[faxis]
xf=sn_f
state=lambda i,n,mod=10:print('%...'.join(np.round(100*((np.arange(0,n,np.round(n/mod))[np.arange(0,n,np.round(n/mod))<=i])/n)).astype(str))) if np.isin(i,np.arange(0,n,np.round(n/mod))) else None
print('Calculating coherences between TF-Z and HPS-Z. \nThis is not a metric normally stored. \n-Please wait.')
for sri,sr in enumerate(icat.iloc):
    state(sri+1,len(icat)) if sri>0 else print('0.0%')
    st=sr.Traces()
    sn=None if st is None else Signal(st.select(location='*NoiseCut')[0],st.select(location='*ATaCR')[0])
    faxis=(xf>=(1/200))&(xf<=1)
    sr.Coherence.update({'TFZ_HPSZ':None if sn is None else sn.coherence()[1][faxis]})
TFZ_HPSZ=np.array([sr.Coherence['TFZ_HPSZ'] for sr in icat.iloc])
mags=icat.Magnitude.to_numpy()


# --Plotting
fig, axes = figs(3,1,f=(6,6))
mthds = ['TF_Z','HPS_Z','TFZ_HPSZ']
obs_colors=['#001325','#073958','#425978','#6a5e74','#ed7968']
obs_cmap=LinearSegmentedColormap.from_list('obs_cmap', obs_colors)
cmap=obs_cmap
cmap=cmap.resampled(magwins.shape[0])
# --- magnitudes → colors ---
magvals = np.array([mg.mean() for mg in magwins])
vmin,vmax=6.0,7.75
den = vmax - vmin if vmax > vmin else 1.0
frac = (magvals - vmin) / den # scale to [0,1]
colors = cmap(frac) # one color per magwin
# cmap = ListedColormap(colors)
# ScalarMappable just for the colorbar
norm=mpl.colors.BoundaryNorm(np.unique(magwins.reshape(-1)),len(magwins),clip=True)
sm = mpl.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])
mids = magwins.mean(axis=1)  # 1 value per bin
colors = cmap(norm(mids))    # 8 RGBA colors
for axi, (ax, mthd) in enumerate(zip(axes, mthds)):
    if axi == 2:xf=sn_f;y=np.array([TFZ_HPSZ[(mags>=mg[0]) & (mags<=mg[1]), :].mean(axis=0) for mg in magwins]);ylbl=yttl('TFZ  HPSZ')
    else:xf=(1/usnr.coh.bands);y=np.array([usnr.coh.__dict__[mthd].D[(mags>=mg[0]) & (mags<=mg[1]), :].mean(axis=0) for mg in magwins]);ylbl = yttl(mthd.replace('_',' '))
    faxis=(xf>=(1/200))&(xf<=1)
    for yy, c in zip(y, colors):
        ax.plot(xf[faxis], yy, c=c,lw=3)
    ax.set_ylabel(ylbl);ax.set_yticks([.4,.6,.8,1.0])
ax.set_xlabel('frequency (Hz)')
ax.set_xscale('log')
ax.set_xlim([1/100,1])

# cax: same height as middle axis, just to the right of it
pos=axes[1].get_position()  # middle subplot
# cax = fig.add_axes([pos.x1 + 0.02, pos.y0, 0.01, .3])
cax = fig.add_axes([1.00,.4, 0.01, .28])
cbar = fig.colorbar(sm, cax=cax,pad=1000,extend='max')
cbar.set_label('earthquake magnitude, Mw')   # or 'Mw', etc.
cbar.set_ticks(np.unique(magwins.reshape(-1))[:-1])

file=f"S03.Figure10_TF.Coherence.with.HPS.Traces.Magnitude.{save_format}"
_=save_tight(plotfolder/file,fig,dpi=900)