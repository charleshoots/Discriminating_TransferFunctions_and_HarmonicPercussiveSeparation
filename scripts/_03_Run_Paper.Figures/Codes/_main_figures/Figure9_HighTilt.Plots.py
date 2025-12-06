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

from scipy.stats import iqr
from local_tools.quick_class import *
from local_tools.math import spectra
from obspy.geodetics import kilometers2degrees
# yttl value
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# yttl eta value
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"
# cat value
cat = catalog.copy()
# icat value
icat = cat.sr.copy()
# usnr value
usnr=unpack_metrics(icat)
# dirs value
dirs = dir_libraries()
# OUT CSV value
OUT_CSV=dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx'
df=pd.read_excel(OUT_CSV)
# theta deg value
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in cat.sr.StaName])
cat.sr['TiltDirection']=theta_deg
# oriencohere value
oriencohere=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].oriencohere for stnm in cat.sr.StaName])
cat.sr['TiltCoherence']=oriencohere




# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure9_HighTilt.Plots';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'

# octavg value
octavg=lt.math.octave_average
# --- Figure: coherence spectra by tilt quantile bins, with stems for IQR variability ---
tilt_bin_edges = [0.08267542, 0.22943596, 0.50754827, 0.88716358, 0.9808266 ]


# ============================== user params ==============================
# Band for the "band-limited averages" objects at the end (not used for plotting)
tiltband = (1, 100)  # seconds

# How many tilt bins? (4 => quartiles)
n_tilt_bins = 4

# Frequency range to display (1–100 s)
period_min, period_max = 1.0, 100.0  # seconds
# reverse x value
reverse_x = True   # show long periods on the left (100 → 1 s)
# n stems value
n_stems   = 18     # number of stem positions across the band
# cap scale value
cap_scale = 1.06   # multiplicative cap width on log-x
# =======================================================================

# ---------- required from your environment ----------
# cat.sr.TiltCoherence : vector (Npairs,)
# usnr.coh.<METHOD>.D  : array (Npairs × Nfreq)
# usnr.coh.<METHOD>.bands : period vector (Nfreq,) in seconds
# fnotch(depth) exists; cat.sr has .StaDepth
# ----------------------------------------------------

# --- inputs ---
tilt_coherences = cat.sr.TiltCoherence.to_numpy()
# mtr value
mtr='coh'
# mtr value
mtr = 'snr'
if mtr=='coh':
    # M value
    M = {
        "TF_Z":  usnr.__dict__[mtr].TF_Z,
        "HPS_Z": usnr.__dict__[mtr].HPS_Z,
        "HPS_1": usnr.__dict__[mtr].HPS_1,
        "HPS_2": usnr.__dict__[mtr].HPS_2,
    }
else:
    # M value
    M = {
    "TF_Z":  usnr.__dict__[mtr].TF_Z.R(),
    "HPS_Z": usnr.__dict__[mtr].HPS_Z.R(),
    "HPS_1": usnr.__dict__[mtr].HPS_1.R(),
    "HPS_2": usnr.__dict__[mtr].HPS_2.R(),
    }

# --- period/frequency axis and display mask ---
P_full = np.asarray(M["TF_Z"].bands, float)  # seconds, shape (Nfreq,)
if P_full.ndim != 1:
    raise ValueError("bands must be a 1-D period vector in seconds")
# msk disp value
msk_disp = (P_full >= period_min) & (P_full <= period_max)
if not np.any(msk_disp):
    raise ValueError("Display mask empty; check period_min/period_max vs bands")

# P value
P = P_full[msk_disp]                    # (Nf_disp,)
# logP value
logP = np.log(P)

# --- infragravity notch mask per station on FULL axis, then slice to display ---
f_full = 1.0 / P_full
# nanidx full value
nanidx_full = np.array([f_full > fnotch(sr.StaDepth) for sr in cat.sr.iloc], dtype=bool)  # (Npairs × Nfreq)
# nanidx value
nanidx = nanidx_full[:, msk_disp]                                                          # (Npairs × Nf_disp)



# --- gather per-method spectra and apply masks ---
D = {k: v.D.copy() for k, v in M.items()}    # (Npairs × Nfreq_full)

if mtr=='snr':
    for k in D.keys():D[k]=np.nanmax(D[k],axis=2)


for k in D.keys():
    D[k] = D[k][:, msk_disp]                 # (Npairs × Nf_disp)
    D[k][nanidx] = np.nan

# ================== tilt quantile bins (to match "qsets" style margins) ==================
# Compute n_tilt_bins+1 quantiles from 0..1 (e.g., quartiles with 0, .25, .5, .75, 1)
q = np.linspace(0, 1, n_tilt_bins + 1)
# raw edges value
raw_edges = np.nanquantile(tilt_coherences, q)

# Guard against duplicate edges (e.g., many identical tilt values); expand slightly if needed
edges = raw_edges.copy()
for i in range(1, len(edges)):
    if edges[i] <= edges[i-1]:
        edges[i] = np.nextafter(edges[i-1], 1.0)  # bump minimally toward 1.0
# Build masks and labels per bin
idx_bins = []
# tilt labels value
tilt_labels = []
for i in range(len(edges) - 1):
    lo, hi = edges[i], edges[i+1]
    # include left, exclude right except the final bin which includes right
    if i < len(edges) - 2:
        # msk pairs value
        msk_pairs = (tilt_coherences >= lo) & (tilt_coherences < hi)
    else:
        msk_pairs = (tilt_coherences >= lo) & (tilt_coherences <= hi)
    idx_bins.append(msk_pairs)
    tilt_labels.append(f"{lo:.2f}–{hi:.2f}")

# Colors for up to 6 bins; truncate/extend as needed
palette = ["#9e9e9e", "#7cb342", "#42a5f5", "#ef5350", "#8d6e63", "#ab47bc"]
tilt_colors = palette[:len(idx_bins)]

# --- helper: robust percentiles across station dimension ---
def qstats(A, q=(25, 50, 75)):
    if A.size == 0:
        n = A.shape[1] if A.ndim == 2 else 0
        return (np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan))
    return np.nanpercentile(A, q, axis=0)  # (len(q) × Nf_disp)

# --- choose stem positions and map to nearest indices in log-period space ---
def nearest_log_indices(P_vec, targets, base=None):
    # nearest in log-space, robust to uneven grids and axis reversal
    logx = np.log(P_vec)
    logt = np.log(targets)
    idx = np.array([np.argmin(np.abs(logx - lt)) for lt in logt], dtype=int)
    return np.unique(np.clip(idx, 0, len(P_vec) - 1))

stem_P_targets = np.geomspace(max(period_min, P.min()), min(period_max, P.max()), n_stems)
stem_idx = nearest_log_indices(P, stem_P_targets)

# ============================== plotting ==============================
fig, axes = plt.subplots(1, 4, figsize=(6, 3), sharey=True)
methods = ["TF_Z", "HPS_Z", "HPS_1", "HPS_2"]

for axi,( ax, m )in enumerate(zip(axes, methods)):
    ax.set_axisbelow(True)
    for lab, col, msk_pairs in zip(tilt_labels, tilt_colors, idx_bins):
        if not np.any(msk_pairs):
            continue
        A = D[m][msk_pairs, :]  # (Npairs_in_bin × Nf_disp)
        q25, q50, q75 = qstats(A)

        # median line
        ax.plot(P, q50, lw=.9, label=lab, color=col, alpha=0.95, zorder=6)

        # stems at sparse positions (nearest on log-axis)
        x = P[stem_idx]; y1 = q25[stem_idx]; y2 = q75[stem_idx]
        valid = np.isfinite(y1) & np.isfinite(y2)

        if np.any(valid):
            nz  = valid & (y2 > y1)
            ziq = valid & ~(y2 > y1)  # zero IQR

            # stems with nonzero height
            if np.any(nz):
                xi  = x[nz]; yi1 = y1[nz]; yi2 = y2[nz]
                ax.vlines(xi, yi1, yi2, lw=.6, color=col, alpha=0.95, zorder=7, clip_on=False)
                ax.hlines(yi1, xi/cap_scale, xi*cap_scale, lw=.6, color=col, alpha=0.95, zorder=7, clip_on=False)
                ax.hlines(yi2, xi/cap_scale, xi*cap_scale, lw=.6, color=col, alpha=0.95, zorder=7, clip_on=False)

            # zero-IQR markers so they’re visible
            if np.any(ziq):
                ax.plot(x[ziq], y1[ziq], marker="_", ms=9, color=col, lw=0, alpha=0.95, zorder=7, clip_on=False)

    ax.set_xscale('log')
    ax.set_xlim(period_min, period_max)
    if reverse_x:
        lo, hi = ax.get_xlim()
        ax.set_xlim(hi, lo)
    if mtr=='coh':ax.set_ylim(0, 1)
    ax.set_title(yttl(m.replace('_',' ')) if mtr=='coh' else yttl_eta(m.replace('_',' ')), pad=12)
    ax.grid(True, which='both', ls=':', lw=0.4, alpha=0.6)
    if ax is axes[0]:
        ax.set_ylabel(r'$\gamma$' if mtr=='coh' else r'$\eta$')
    if axi==1:spaces='\;'*19;ax.set_xlabel(f'${spaces}$period, s',ha='left')

for ax in axes:_=ax.set_xlim(100,8)

# # unified legend
_=fig.tight_layout(rect=1.3*np.array([0, .3, 1,1]))
handles, labels = axes[0].get_legend_handles_labels()
lg=fig.legend(handles, labels, ncol=min(len(idx_bins), 4), loc='upper center',
bbox_to_anchor=(0.69, 1.34), frameon=False)
for h in lg.legend_handles:h.set_linewidth(2)


file=f'{mtr.upper()}.HighTiltPlot.{save_format}'
_=save_tight(plotfolder/file,fig,dpi=800)