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

from imports import * #Standard imports for doing anything in this project. Approx. ~19 seconds.
# Spectra scatters for Deployment Parameters vs SNR and Coherence
import matplotlib as mpl
# figs value
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)


# dirs value
dirs = dir_libraries()
# OUT CSV value
OUT_CSV=dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx' #Tilt orientations and coherence from Janiszewski et al. (2023)
df=pd.read_excel(OUT_CSV)
# theta deg value
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in catalog.sr.StaName])
catalog.sr['TiltDirection']=theta_deg
# oriencohere value
oriencohere=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].oriencohere for stnm in catalog.sr.StaName])
catalog.sr['TiltCoherence']=oriencohere

# icat value
icat=catalog.sr.copy()
# usnr value
usnr=unpack_metrics(icat)


# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure8_SNR.SpectraBands.by.DeploymentParam';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'


# --- helpers ---
def _listed_cmap(base_cmap, n):
    # base value
    base = cm.cmaps[base_cmap] if isinstance(base_cmap, str) else base_cmap
    return mpl.colors.ListedColormap(base(np.linspace(0, 1, n)))
# function add discrete colorbar
def add_discrete_colorbar(xcat,ax, fig, labels, cmap='glasgow', edges=None,
    # orientation value
    orientation='vertical', fraction=0.035, pad=0.04, tickfmt=None):
    """
    If edges is None -> categorical strings in `labels`.
    If edges is array-like -> numeric bins with edges defining bins; `labels` are strings to show.
    Returns a ScalarMappable you can use to map values to colors consistently.
    """
    if edges is None:
        # K value
        K = len(labels)
        # bounds value
        bounds = np.arange(-0.5, K + 0.5, 1.0)
        # norm value
        norm = mpl.colors.BoundaryNorm(bounds, K)
        # ticks value
        ticks = np.arange(K)
        # ticklabels value
        ticklabels = labels
        # cmap value
        cmap = _listed_cmap(cmap, K)
    else:
        # edges value
        edges = np.asarray(edges)
        # K value
        K = len(edges) - 1
        # centers value
        centers = 0.5 * (edges[:-1] + edges[1:])
        # bounds value
        bounds = edges
        # norm value
        norm = mpl.colors.BoundaryNorm(bounds, K)
        # ticks value
        ticks = centers
        if tickfmt is None:
            # ticklabels = [f"{hi:.2g}" for lo, hi in zip(edges[:-1], edges[1:])]
            ticklabels = labels
        else:
            ticklabels = [tickfmt(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])]
        cmap = _listed_cmap(cmap, K)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation=orientation, fraction=fraction, pad=pad)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticklabels,rotation=90 if xcat=='Seismometer' else 0,ha='center')
    cbar.ax.tick_params(pad=12)  # try 4–8
    if xcat=='TiltCoherence':cbar.set_label('ZH coherence (tilt)'.lower())
    elif xcat=='StaDepth':cbar.set_label('Water depth, m'.lower())
    elif xcat=='EvDepth':cbar.set_label('Earthquake depth, km'.lower())
    else:
        ttl=xcat.replace('TiltCoherence','ZH coherence (Tilt)').replace('_',' ').lower()
        cbar.set_label(f'{ttl[0].upper()}{ttl[1:]}'.lower().replace('mw','Mw'))
    return sm, norm
def cbarlam(axes,fig,cmap=cm.cmaps['glasgow'],vmin=0,vmax=1,orientation='vertical',fraction=0.025,pad=0.04):
    # (kept for compatibility if you still want a continuous bar somewhere)
    norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cbar=fig.colorbar(sm, ax=axes, orientation=orientation, fraction=fraction, pad=pad)
    return cbar


# --- your plotting ---


Horiz = False # whether to make horizontal channel plots (True) or vertical plots (False)


cmap = 'glasgow'
cmap = 'grayC'
# for mthd in ['TF_Z','HPS_Z']:
sets = ['Instrument_Design','Pressure_Gauge','TiltCoherence','Seismometer']
# sets = ['StaDepth','Distance','Magnitude','EvDepth']
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"
# for mtr in ['snr','coh']:
for mtr in ['snr']:
    fig,axes = figs(2,2,f=(9,4),x=True,y=True)
    if mtr=='coh':f = 1 / usnr.coh.TF_Z.bands
    else:f = 1 / usnr.snr.TF_Z.bands
    nanidx = np.array([f > fnotch(sr.StaDepth) for sr in icat.iloc])  # shape: (Nstations, Nfreq)
    for axi,(ax, xcat) in enumerate(zip(axes.reshape(-1), sets)):
        xc_raw = icat[xcat].to_numpy()

        # build color mapping + colorbar for this subplot
        if xcat == 'TiltCoherence':
            dx = 0.1
            edges = np.arange(0, 1 + dx, dx)
            bins = list(zip(edges[:-1], edges[1:]))                  # [(lo,hi), ...]
            labels = [f"{hi:.1f}" for lo,hi in bins]
            sm, norm = add_discrete_colorbar(xcat,ax, fig, labels, cmap=cmap, edges=edges,
            orientation='vertical', fraction=0.04, pad=0.03)
            def color_for_bin(lo, hi):
                mid = 0.5 * (lo + hi)
                return sm.to_rgba(mid)
            u = bins
        elif xcat in ['StaDepth','Magnitude','Distance','EvDepth']:

            dx = {'StaDepth':500,'Magnitude':.5,'Distance':20,'EvDepth':50}[xcat]
            dmax = {'StaDepth':6000,'Magnitude':8.0,'Distance':120,'EvDepth':700}[xcat]
            edges = np.arange(0, dmax + dx, dx)
            bins = list(zip(edges[:-1], edges[1:]))                  # [(lo,hi), ...]
            labels = [f"{hi}" for lo,hi in bins]
            sm, norm = add_discrete_colorbar(xcat,ax, fig, labels, cmap=cmap, edges=edges,
            orientation='vertical', fraction=0.04, pad=0.03)
            def color_for_bin(lo, hi):
                mid = 0.5 * (lo + hi)
                return sm.to_rgba(mid)
            u = bins
        else:
            # categorical strings
            u_vals = np.unique(xc_raw)
            u_labels = [str(v) for v in u_vals]
            if xcat=='Seismometer':u_labels=[f'{u.split(' ')[0]}\n{' '.join(u.split(' ')[1:])}' for u in u_labels]
            sm, norm = add_discrete_colorbar(xcat,ax, fig, u_labels, cmap=cmap,
            orientation='vertical', fraction=0.04, pad=0.03)

            idx_map = {val: i for i, val in enumerate(u_vals)}
            def color_for_val(val):
                return sm.to_rgba(idx_map[val])

            u = u_vals

        # plot TF_Z and HPS_Z grouped by xci with colors tied to the colorbar
        # for mthd, mkr in [('TF_Z','s'), ('HPS_Z','o')]:
        if Horiz:mthds_mkr = [('HPS_1','s'), ('HPS_2','o')]
        else:mthds_mkr = [('HPS_Z','s'), ('TF_Z','o')]
        for mthd, mkr in mthds_mkr:
            if mtr=='coh':Y = usnr.__dict__[mtr].__dict__[mthd].D.copy()  # (Nstations, Nfreq)
            else:Y=usnr.__dict__[mtr].__dict__[mthd].R().D.copy();Y=np.nanmax(Y,axis=2)  # (Nstations, Nfreq)
            Y[nanidx] = np.nan

            for xci in u:
                if xcat in ['TiltCoherence','StaDepth','Magnitude','Distance','EvDepth']:
                    lo, hi = xci
                    idx = (xc_raw >= lo) & (xc_raw <= hi)
                    c = color_for_bin(lo, hi)
                else:
                    idx = (xc_raw == xci)
                    c = color_for_val(xci)

                if not np.any(idx): 
                    continue
                ymean = np.nanmean(Y[idx, :], axis=0)
                ax.scatter(f, ymean, s=10, marker=mkr, c=[c],ec='k' if mthd=='TF_Z' else 'r',lw=0.5 if mthd=='TF_Z' else 0.3)
                ax.plot(f, ymean, color='k' if mthd=='TF_Z' else 'r',lw=1,zorder=-1e5,alpha=0.3)
                ax.plot(f, ymean, color=mcolors.to_hex(c),lw=1,zorder=-1e5,alpha=0.3)

        ax.set_xscale('log')
        ax.set_xlim(1/100, 1/7)
        ax.set_ylim(0, 1.05) if mtr=='coh' else None
        ttl=xcat.replace('TiltCoherence','Tilt').replace('_',' ').lower()
        # ax.set_title(f'{ttl[0].upper()}{ttl[1:]}')
        ax.grid(True, zorder=0, alpha=0.25)
        ax.set_ylabel(r'$\eta$' if mtr=='snr' else r'$\gamma$')

    if axi==3:
        ax=axes[0,0];ax.scatter(np.nan,0,label=yttl(mthds_mkr[0][0]) if mtr=='coh' else yttl_eta(mthds_mkr[0][0]),marker='s',c='None',ec='r',lw=.5,s=12);ax.scatter(np.nan,0,label=yttl(mthds_mkr[1][0]) if mtr=='coh' else yttl_eta(mthds_mkr[1][0]),c='None',marker='s',ec='k',lw=2,s=12)
        lg=ax.legend(fancybox=False,frameon=False,ncols=2,labelspacing=0.0,columnspacing=0.5)

    ax=axes[1,0];ax.set_xlabel(r'$\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;\;$frequency, Hz',x=.5,ha='left')
    fig.subplots_adjust(wspace=0)
    file = f'S03.Figure08_{mthd}.{mtr}.by.DeploymentParameters.{save_format}'
    if Horiz:file = f'Horiz.{mtr.upper()}.{'.'.join(sets)}.{save_format}'
    save_tight(plotfolder/file,fig,dpi=700)