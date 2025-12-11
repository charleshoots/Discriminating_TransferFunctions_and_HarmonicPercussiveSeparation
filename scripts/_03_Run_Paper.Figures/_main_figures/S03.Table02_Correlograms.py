# Author: Charles Hoots
# This code was developed as part of my PhD research in the
# Department of Earth Sciences, University of Hawai‘i at Mānoa.
# Unless otherwise noted, the code is my own original work.
# External libraries and standard research software packages are used as cited.

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

# sys.path.append(str(Path(__file__).parent.parent))
from imports import * #Standard imports for doing anything in this project. Approx. ~19 seconds.
# cat value
cat = catalog.copy()
# octavg value
octavg=lt.math.octave_average
# Noise Spectra
f=cat.r.iloc[0].Data.Noise.Averaged().f
# faxis value
faxis=(f>0)&(f<=1)
# f value
f=f[faxis]
# noise f value
noise_f=f
cat.r['Noise']=[AttribDict({'f':f,
'Z':PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__['cZZ'][faxis]),
'P':PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__['cPP'][faxis]),
'H':np.mean([PowDisp_to_AcceldB(f,s.Data.Noise.Averaged().power.__dict__[c][faxis]) for c in ['c11','c22']],axis=0)
}) for s in cat.r.iloc]
# rmse value
rmse=lambda y:( (  ( abs(y)-abs(y).mean() )**2  ).mean())**.5 
# rms value
rms=lambda y:np.mean(y**2)**0.5
# s value
s=cat.r.iloc[0];faxis=(s.Noise.f>(1/100) )& (s.Noise.f<=1)
# f value
f=s.Noise.f[faxis]
cat.r['NoiseAverage']=[{f'{b[0]}_{b[1]}':-rms(s.Noise.Z[faxis][(f<=(1/b[0]))&(f>=(1/b[1]))]) for b in [[1,10],[10,30],[30,100]]} for s in cat.r.iloc]
cat.sr['NoiseAverage']=[cat.r.loc[sr.StaName].NoiseAverage[0] for sr in cat.sr.iloc]
# function custom cmap
def custom_cmap(ind=0,nbins=5):
    if ind==0:cmap = cm.cmaps['glasgow'].reversed().resampled(nbins)
    if ind==1:cmap = cm.cmaps['batlow'].reversed().resampled(nbins)
    return cmap
# figs value
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
from obspy.signal.trigger import classic_sta_lta,carl_sta_trig,recursive_sta_lta
# stalta methods value
stalta_methods={'classic':classic_sta_lta,'carl':carl_sta_trig,'recurssive':recursive_sta_lta}
# darken value
darken=lambda cmap,frac=0.8:ListedColormap([cmap(i) for i in np.arange(0,frac,0.01)]).resampled(100)
# luminance value
luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098 
# suspect stations value
suspect_stations=np.array(['ZA.B02','YL.C09W','7D.G25B','7D.FS08D','7D.G17B','YL.A14W'])
# baz value
baz=lambda s:obspy.geodetics.base.gps2dist_azimuth(s.Latitude,s.Longitude,s.LaLo[0],s.LaLo[1])[1]
# bootstrap value
bootstrap = lambda y,nruns=10000,nchoose=100,aggregate=np.mean: np.mean([aggregate(np.random.choice(y[~np.isnan(y)],nchoose)) for _ in range(nruns)])
# centers value
centers=lambda x:np.array((x[:-1] + x[1:]) / 2)
# dbin value
dbin = lambda x:np.array([x[:-1],x[1:]]).T
# phases value
phases=['P','S','Rg'];preferred_pbands={'P':'1_10','S':'10_30','Pdiff':'1_10','Sdiff':'10_30','Rg':'30_100'}
# methods value
methods=['NoiseCut','ATaCR'];mnames={'NoiseCut':'HPS','ATaCR':'TF'}
# mnames r value
mnames_r={mnames[k]:k for k in mnames.keys()};mname_comp={f'HPS_Z':'NoiseCut','TF':'ATaCR','Original':'Original'}
# mnames comp r value
mnames_comp_r={mname_comp[k]:k for k in mname_comp.keys()}
# cohnames2snrnames value
cohnames2snrnames=c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_1':'HPS.1','HPS_2':'HPS.2'}
# c2s value
c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_1':'HPS.1','HPS_2':'HPS.2'}
# c2s r value
c2s_r={c2s[k]:k for k in c2s.keys()}
# make bands
def make_bands(N, width=None,line=np.linspace, lo=1.0, hi=100.0):
    if width is None:width=(hi-lo)/N
    if width<=0 or width>(hi-lo): raise ValueError("width must be in (0, hi-lo]")
    # s value
    s=line(lo, hi-width, N)
    return np.c_[s, s+width]


# dirs value
dirs=io.dir_libraries()


# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Table2_Correlograms';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'png'

# OUT CSV value
OUT_CSV=dirs.Catalogs/'Janiszewski_etal_2023_StationAverages.xlsx'
df=pd.read_excel(OUT_CSV)
# theta deg value
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in cat.sr.StaName])
cat.sr['TiltDirection']=theta_deg
# oriencohere value
oriencohere=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].oriencohere for stnm in cat.sr.StaName])
cat.sr['TiltCoherence']=oriencohere


# icat value
icat=cat.sr.copy()
# usnr value
usnr=unpack_metrics(icat)

# function argsort luminence
def argsort_luminence(cc,cmap):
    # cc value
    cc=np.array(cc).ravel()
    # luminance value
    luminance=lambda rgb:np.sum([scl*(x / 255.0) for x,scl in zip(rgb[:-1],[0.2126,0.7152,0.0722])])/0.00392156862745098
    # color value
    color=[cmap(c) for c in enumerate(cc)]
    # zorder value
    zorder=np.argsort(np.array([luminance(c) for c in color]))
    return zorder
from matplotlib.colors import LinearSegmentedColormap

df=pd.read_excel(OUT_CSV)
# theta deg value
theta_deg=np.array([df[(df.network==stnm.split('.')[0])&((df.station==stnm.split('.')[1]))].iloc[0].orientation for stnm in icat.StaName])
cat.sr['Tilt']=theta_deg

from scipy.stats import spearmanr, kendalltau, rankdata
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rankdata, kendalltau
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ---------- categorical helpers ----------
def _is_multilabel_cell(v): return isinstance(v, (list, tuple, set))

# function to columns categorical
def _to_columns_categorical(vec, max_levels=None):
    # v value
    v = np.asarray(vec, dtype=object); N = len(v)
    labels_per_row, all_labels = [], set()
    # loop over v
    for x in v:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            labels_per_row.append(None); continue
        # s value
        s = set(map(str, x)) if _is_multilabel_cell(x) else {str(x)}
        labels_per_row.append(s); all_labels |= s
    # levels value
    levels = sorted(all_labels)
    if (max_levels is not None) and (len(levels) > max_levels): levels = levels[:max_levels]
    # K value
    K = len(levels)
    if K == 0: return np.full((N,1), np.nan, float)
    # idx value
    idx = {lab:j for j,lab in enumerate(levels)}
    # M value
    M = np.full((N,K), np.nan, float)
    for i, s in enumerate(labels_per_row):
        if s is None: continue
        # row value
        row = np.zeros(K, float)
        # row = np.ones(K, float)

        # loop over s
        for lab in s:
            j = idx.get(lab); 
            if j is not None: row[j] = 1.0
            # if j is not None: row[j] = 0.0
        M[i,:] = row
    return M

# ---------- SAFE column access (never coerce whole table) ----------
def _as_columns(X):
    if isinstance(X, np.ndarray):
        if X.ndim == 1: return [X]
        return [X[:, j] for j in range(X.shape[1])]
    return list(X)

# function try float column
def _try_float_column(col):
    try:
        # arr value
        arr = np.asarray(col, dtype=float)
        return arr, True
    except Exception:
        return np.asarray(col, dtype=object), False

def _is_categorical_vector(vec):
    arr, numeric_ok = _try_float_column(vec)
    if numeric_ok: return False
    arr = np.asarray(vec, dtype=object)
    for k in range(min(32, len(arr))):
        v = arr[k]
        if isinstance(v, str) or _is_multilabel_cell(v): return True
    return True

# ---------- pairwise scalar associations for numeric-numeric ----------
def _pearson_pair(a,b):
    m = np.isfinite(a)&np.isfinite(b)
    return np.corrcoef(a[m],b[m])[0,1] if m.sum()>=2 else np.nan

def _spearman_pair(a,b):
    m = np.isfinite(a)&np.isfinite(b)
    if m.sum()<2: return np.nan
    ar = rankdata(a[m], nan_policy='omit'); br = rankdata(b[m], nan_policy='omit')
    return np.corrcoef(ar,br)[0,1]

def _kendall_pair(a,b):
    m = np.isfinite(a)&np.isfinite(b)
    if m.sum()<2: return np.nan
    tau,_ = kendalltau(a[m], b[m], nan_policy='omit')
    return float(tau) if tau is not None else np.nan

# ---------- blockwise reducers (numeric-categorical & categorical-categorical) ----------
def _multiple_corr_signed(y, X):
    # y: (N,), X: (N,k) one-hot (or numeric block)
    if X.ndim==1: X = X[:,None]
    m = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if m.sum() < max(3, X.shape[1]+1): return np.nan
    y = y[m]; X = X[m]
    X = np.column_stack([np.ones(len(X)), X])
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None); yhat = X @ beta
    except np.linalg.LinAlgError:
        XtX = X.T@X + 1e-8*np.eye(X.shape[1]); beta = np.linalg.solve(XtX, X.T@y); yhat = X@beta
    r = np.corrcoef(y,yhat)[0,1]
    return float(r)

def _rv_coeff(X,Z):
    if X.ndim==1: X=X[:,None]
    if Z.ndim==1: Z=Z[:,None]
    m = np.all(np.isfinite(X),axis=1) & np.all(np.isfinite(Z),axis=1)
    if m.sum()<3: return np.nan
    Xc = X[m]-X[m].mean(0,keepdims=True); Zc = Z[m]-Z[m].mean(0,keepdims=True)
    A = Xc@Xc.T; B = Zc@Zc.T
    num = np.trace(A@B); den = np.sqrt(np.trace(A@A)*np.trace(B@B))
    return float(num/den) if den>0 else np.nan

# ---------- expand only categoricals; keep per-variable blocks ----------
def _expand_categoricals(X, max_levels=500):
    cols = _as_columns(X)
    mats, groups, is_cat, blocks = [], [], [], []
    col_start = 0
    for col in cols:
        if _is_categorical_vector(col):
            M = _to_columns_categorical(col, max_levels=max_levels)  # (N, K)
            mats.append(M); blocks.append(M)
            groups.append(list(range(col_start, col_start + M.shape[1])))
            col_start += M.shape[1]; is_cat.append(True)
        else:
            v, _ = _try_float_column(col)
            M = v.reshape(-1, 1)
            mats.append(M); blocks.append(M)
            groups.append([col_start]); col_start += 1; is_cat.append(False)
    Xexp = np.concatenate(mats, axis=1)
    return Xexp, groups, is_cat, blocks

# ---------- correlation using blockwise association ----------
def _corr_matrix(X, method="pearson"):
    # expand only categoricals; retain per-variable blocks
    _, groups, is_cat, blocks = _expand_categoricals(X, max_levels=500)
    p = len(blocks)
    C = np.eye(p, dtype=float)

    for i in range(p):
        for j in range(i+1, p):
            if not is_cat[i] and not is_cat[j]:
                a = blocks[i][:,0]; b = blocks[j][:,0]
                if method == "pearson":   r = _pearson_pair(a,b)
                elif method == "spearman": r = _spearman_pair(a,b)
                elif method == "kendall":  r = _kendall_pair(a,b)
                else: raise ValueError("method must be 'pearson','spearman','kendall'")
            elif is_cat[i] and is_cat[j]:
                r = _rv_coeff(blocks[i], blocks[j])
            else:
                # numeric ↔ categorical
                if is_cat[i]:
                    y = blocks[j][:,0]; Xblk = blocks[i]
                else:
                    y = blocks[i][:,0]; Xblk = blocks[j]
                r = _multiple_corr_signed(y, Xblk)
            C[i,j] = C[j,i] = r if np.isfinite(r) else np.nan
    return C

# ---------- plotting helpers (safe numeric view) ----------
def _encode_categorical_column(col):
    arr = np.asarray(col, dtype=object)
    out = np.full(len(arr), np.nan, float)
    def key(v):
        if v is None or (isinstance(v, float) and np.isnan(v)): return None
        if isinstance(v, (list, tuple, set)): return tuple(sorted(map(str, v)))
        return str(v)
    keys, seen = [], set()
    for v in arr:
        k = key(v)
        if k is not None and k not in seen:
            seen.add(k); keys.append(k)
    lut = {k: i for i, k in enumerate(keys)}
    for i, v in enumerate(arr):
        k = key(v)
        if k is not None: out[i] = float(lut[k])
    return out

def _numeric_view_for_plot(X, encode_categoricals=False):
    cols = X if isinstance(X, (list, tuple)) else ([X] if getattr(X, "ndim", 1) == 1 else [X[:, j] for j in range(X.shape[1])])
    out = []
    for col in cols:
        try:
            v = np.asarray(col, dtype=float)
            out.append(v)
        except Exception:
            if encode_categoricals:
                out.append(_encode_categorical_column(col))
            else:
                out.append(np.full(len(col), np.nan, float))
    return np.column_stack(out)

# ---------- plotting (unchanged layout/sizing) ----------
def tri_correlogram(X, names=None, method="spearman", bins=24, 
    s=15, a=.35, cmap=None, block_fnt=40, fnt=20, figscaling=2.2):
    n = len(X)
    C = _corr_matrix(X, method=method)

    if cmap is None:
        cmap = cm.cmaps['vik']

    fig, axes = plt.subplots(n, n, figsize=(figscaling*n, figscaling*n))
    im_for_cb = None
    Xarr = _numeric_view_for_plot(X, encode_categoricals=True)
    for i in range(n):
        for j in range(n):
            ax = axes[i,j]
            xi = Xarr[:, j]; yi = Xarr[:, i]
            m  = np.isfinite(xi) & np.isfinite(yi)

            if i == j:
                ax.hist(xi[m], bins=bins)
            elif i > j:
                ax.scatter(xi[m], yi[m], s=s, alpha=a)
            else:
                r = C[i,j]
                im = ax.imshow([[r]], vmin=-1, vmax=1, interpolation="nearest", cmap=cmap)
                im_for_cb = im
                ax.text(0, 0, f"{r:.2f}" if np.isfinite(r) else "NaN",
                        ha="center", va="center", fontsize=block_fnt, fontweight='bold')

            ax.set_xticks([]); ax.set_yticks([])
            if names is not None:
                if i == n-1: ax.set_xlabel(names[j], ha='center', va='center', labelpad=45, )
                if j == 0:   ax.set_ylabel(names[i], ha='center', va='center', labelpad=30, )
            if i >= j and np.any(m):
                xmn,xmx = np.nanmin(xi[m]), np.nanmax(xi[m])
                ymn,ymx = np.nanmin(yi[m]), np.nanmax(yi[m])
                if i>j: ax.set_xlim(xmn,xmx); ax.set_ylim(ymn,ymx)
            if i != j:
                for spine in ax.spines.values(): spine.set_visible(False)

    if im_for_cb is not None:
        ax = axes.reshape(-1)[-1]
        cax = inset_axes(ax, width='40%', height=f'{n*110}%',
                         loc="lower left", bbox_to_anchor=(1.08, 0.0, 1, 1.0),
                         bbox_transform=ax.transAxes, borderpad=0)
        cb = fig.colorbar(im_for_cb, cax=cax, cmap=cmap)
        cb.set_label({"pearson":"Pearson r","spearman":"Spearman ρ","kendall":"Kendall τ"}[method], )
        cb.ax.tick_params(labelsize=fnt)
    fig.tight_layout(rect=[0,0,1.0,1])
    return C, fig
yttl = lambda c: fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
yttl_eta = lambda c:fr"$\underset{{{c}}}{{\eta\;\;\;\;\;\;\;}}$"


# ====================================================================================================================
# ====================================================================================================================

# ====================================================================================================================
# Setup: "Robust" - The most detailed layout of correlation measurements.
# ====================================================================================================================
props=AttribDict()
props.X=[lambda:usnr.coh.TF_Z.Average(b,fn=fn), #band averaged coherence
lambda:usnr.coh.HPS_Z.Average(b,fn=fn), #----- 
lambda:usnr.coh.HPS_1.Average(b,fn=fn), #-----
lambda:usnr.coh.HPS_2.Average(b,fn=fn), #-----
lambda:usnr.snr.TF_Z.R().P.Average(b,fn=fn), #----- #band averaged SNR for P/Pdiff
lambda:usnr.snr.HPS_Z.R().P.Average(b,fn=fn), #-----
lambda:usnr.snr.HPS_1.R().P.Average(b,fn=fn), #-----
lambda:usnr.snr.HPS_2.R().P.Average(b,fn=fn), #-----
lambda:usnr.snr.TF_Z.R().S.Average(b,fn=fn), #----- #band averaged SNR for S/Sdiff
lambda:usnr.snr.HPS_Z.R().S.Average(b,fn=fn), #-----
lambda:usnr.snr.HPS_1.R().S.Average(b,fn=fn), #-----
lambda:usnr.snr.HPS_2.R().S.Average(b,fn=fn), #-----
lambda:usnr.snr.TF_Z.R().Rg.Average(b,fn=fn), #----- #band averaged SNR for Rayleigh
lambda:usnr.snr.HPS_Z.R().Rg.Average(b,fn=fn), #-----
lambda:usnr.snr.HPS_1.R().Rg.Average(b,fn=fn), #-----
lambda:usnr.snr.HPS_2.R().Rg.Average(b,fn=fn), #-----
lambda:usnr.ST.TF_Z.R().P.Average(b,fn=fn), #----- #band averaged signal window reduction for P/Pdiff
lambda:usnr.ST.HPS_Z.R().P.Average(b,fn=fn), #-----
lambda:usnr.ST.HPS_1.R().P.Average(b,fn=fn), #-----
lambda:usnr.ST.HPS_2.R().P.Average(b,fn=fn), #-----
lambda:usnr.ST.TF_Z.R().S.Average(b,fn=fn), #----- #band averaged signal window reduction for S/Sdiff
lambda:usnr.ST.HPS_Z.R().S.Average(b,fn=fn), #-----
lambda:usnr.ST.HPS_1.R().S.Average(b,fn=fn), #-----
lambda:usnr.ST.HPS_2.R().S.Average(b,fn=fn), #-----
lambda:usnr.ST.TF_Z.R().Rg.Average(b,fn=fn), #----- #band averaged signal window reduction for Rayleigh
lambda:usnr.ST.HPS_Z.R().Rg.Average(b,fn=fn), #-----
lambda:usnr.ST.HPS_1.R().Rg.Average(b,fn=fn), #-----
lambda:usnr.ST.HPS_2.R().Rg.Average(b,fn=fn), #-----
lambda:usnr.LT.TF_Z.R().Average(b,fn=fn), #----- #band averaged noise window reduction
lambda:usnr.LT.HPS_Z.R().Average(b,fn=fn), #-----
lambda:usnr.LT.HPS_1.R().Average(b,fn=fn), #-----
lambda:usnr.LT.HPS_2.R().Average(b,fn=fn), #-----
lambda:np.array(list(icat.StaDepth)), #----- #water depth
lambda:np.array([i[bn] for i in np.array(list(icat.NoiseAverage))]) if bn in ['1_10','10_30','30_100'] else np.array([np.mean(list(i.values())) for i in np.array(list(icat.NoiseAverage))]), # band averaged station averaged noise levels
lambda:np.array(list(icat.TiltCoherence)), #----- Station average ZH coherence
lambda:np.array(list(icat.Instrument_Design)), #-----
lambda:np.array(list(icat.Seismometer)), #-----
lambda:np.array(list(icat.Pressure_Gauge)), #-----
lambda:np.array(list(icat.Sediment_Thickness_m)), #-----
lambda:np.array(list(icat.Magnitude)), #-----
lambda:np.array(list(icat.Distance)), #-----
lambda:np.array(list(icat.EvDepth)), #-----
]
props.names=[
yttl('TFZ'), yttl('HPSZ'), yttl('HPS1'), yttl('HPS2'), #band averaged coherence
yttl_eta('TFZ') + '\n(P/Pdiff)',yttl_eta('HPSZ') + '\n(P/Pdiff)',yttl_eta('HPS1') + '\n(P/Pdiff)',yttl_eta('HPS2') + '\n(P/Pdiff)', #band averaged SNR for P/Pdiff
yttl_eta('TFZ') + '\n(S/Sdiff)',yttl_eta('HPSZ') + '\n(S/Sdiff)',yttl_eta('HPS1') + '\n(S/Sdiff)',yttl_eta('HPS2') + '\n(S/Sdiff)', #band averaged SNR for S/Sdiff
yttl_eta('TFZ') + '\n(Rayleigh)',yttl_eta('HPSZ') + '\n(Rayleigh)',yttl_eta('HPS1') + '\n(Rayleigh)',yttl_eta('HPS2') + '\n(Rayleigh)', #band averaged SNR for Rayleigh
'TFZ\nsignal\nratio (P/Pdiff)','HPSZ\nsignal\nratio (P/Pdiff)','HPS1\nsignal\nratio (P/Pdiff)','HPS2\nsignal\nratio (P/Pdiff)',#band averaged signal ratio for P/Pdiff
'TFZ\nsignal\nratio (S/Sdiff)','HPSZ\nsignal\nratio (S/Sdiff)','HPS1\nsignal\nratio (S/Sdiff)','HPS2\nsignal\nratio (S/Sdiff)',#band averaged signal ratio for P/Pdiff
'TFZ\nsignal\nratio (Rayleigh)','HPSZ\nsignal\nratio (Rayleigh)','HPS1\nsignal\nratio (Rayleigh)','HPS2\nsignal\nratio (Rayleigh)',#band averaged signal ratio for P/Pdiff
'TF Z\nnoise\nratio','HPS Z\nnoise\nratio','HPS 1\nnoise\nratio','HPS 2\nnoise\nratio',#band averaged noise window reduction
'Water\ndepth',
'Station\naverage\nnoise',
'Tilt (ZH)',
'Instrument\ndesign',
'Seismometer',
'Pressure\ngauge',
'Sediment\nThickness',
'Magnitude',
'Distance',
'Quake\ndepth'
]


# ====================================================================================================================
# Setup: "Peak metrics" - A more condensed layout of correlation measurements.
# ====================================================================================================================
# Option for how we aggregate the SNR average across all three phases in a trace (P/Pdiff, S/Sdiff, and Rayleigh).
# Options are either by the peak (max) or mean.
# snravg = np.nanmean #Nan version ensures protection from values that may get masked either by notch sensitivity or filter band.
# snravg = np.nanmax
# props=AttribDict()
# props.snravg = np.nanmax
# props.X=[lambda:usnr.coh.TF_Z.Average(b,fn=fn), #band averaged coherence
# lambda:usnr.coh.HPS_Z.Average(b,fn=fn), #----- 
# lambda:usnr.coh.HPS_1.Average(b,fn=fn), #-----
# lambda:usnr.coh.HPS_2.Average(b,fn=fn), #-----
# lambda:props.snravg(usnr.snr.TF_Z.R().Average(b,fn=fn),axis=1), #----- Maximum (peak) band averaged SNR from all three phases in each source-receiver.
# lambda:props.snravg(usnr.snr.HPS_Z.R().Average(b,fn=fn),axis=1), #-----
# lambda:props.snravg(usnr.snr.HPS_1.R().Average(b,fn=fn),axis=1), #-----
# lambda:props.snravg(usnr.snr.HPS_2.R().Average(b,fn=fn),axis=1), #-----
# lambda:props.snravg(usnr.ST.TF_Z.R().Average(b,fn=fn),axis=1), #----- #Maximum (peak, so the least) band averaged signal window reduction.
# lambda:props.snravg(usnr.ST.HPS_Z.R().Average(b,fn=fn),axis=1), #-----
# lambda:props.snravg(usnr.ST.HPS_1.R().Average(b,fn=fn),axis=1), #-----
# lambda:props.snravg(usnr.ST.HPS_2.R().Average(b,fn=fn),axis=1), #-----
# lambda:usnr.LT.TF_Z.R().Average(b,fn=fn), #----- #band averaged noise window reduction
# lambda:usnr.LT.HPS_Z.R().Average(b,fn=fn), #-----
# lambda:usnr.LT.HPS_1.R().Average(b,fn=fn), #-----
# lambda:usnr.LT.HPS_2.R().Average(b,fn=fn), #-----
# lambda:np.array(list(icat.StaDepth)), #----- #water depth
# lambda:np.array([i[bn] for i in np.array(list(icat.NoiseAverage))]) if bn in ['1_10','10_30','30_100'] else np.array([np.mean(list(i.values())) for i in np.array(list(icat.NoiseAverage))]), # band averaged station averaged noise levels
# lambda:np.array(list(icat.TiltCoherence)), #----- Station average ZH coherence
# lambda:np.array(list(icat.Instrument_Design)), #-----
# lambda:np.array(list(icat.Seismometer)), #-----
# lambda:np.array(list(icat.Pressure_Gauge)), #-----
# lambda:np.array(list(icat.Sediment_Thickness_m)), #-----
# lambda:np.array(list(icat.Magnitude)), #-----
# lambda:np.array(list(icat.Distance)), #-----
# lambda:np.array(list(icat.EvDepth)), #-----
# ]
# props.names=[
# yttl('TF Z'), yttl('HPS Z'), yttl('HPS 1'), yttl('HPS 2'),
# r'$R_{\mathrm{TF}\ Z}$',r'$R_{\mathrm{HPS}\ Z}$' ,r'$R_{\mathrm{HPS}\ 1}$' ,r'$R_{\mathrm{HPS}\ 2}$',
# 'TF Z\nsignal \npeak','HPS Z\nsignal \npeak','HPS 1\nsignal \npeak','HPS 2\nsignal \npeak',
# 'TF Z\nnoise \nwindow','HPS Z\nnoise \nwindow','HPS 1\nnoise \nwindow','HPS 2\nnoise \nwindow',
# 'Water\ndepth','Station\naverage\nnoise',f'Tilt\n({yttl('ZH')})','Instrument\ndesign','Seismometer','Pressure\ngauge','Sediment\nthickness','Magnitude','Distance','Earthquake\ndepth']


#------Example set : WaterDepth, TiltCoherence, Instrument_Design -- Comment out to disable
# N=10 #How many values (source-receivers) to use in the example
# instlist=['AB','TRM','KE']
# inst_inds=unravel([[np.where(np.array(list(icat.Instrument_Design))==i)[0][:2],np.where(np.array(list(icat.Instrument_Design))==i)[0][-2:]] for i in instlist])[:N]
# props=AttribDict()
# props.X=[lambda:np.array(list(icat.StaDepth))[inst_inds],
# lambda:np.array(list(icat.TiltCoherence))[inst_inds],
# lambda:np.array(list(icat.Instrument_Design))[inst_inds]]
# props.names=['Water\ndepth',f'Tilt\n({yttl('ZH')})','Instrument\ndesign']







names=props.names
X=props.X
assert len(X)==len(names), 'missing properties'
icat=cat.sr.copy()
bands=[[1,100]] #period bands to use
corrmethods=['pearson'] # correlation function used. can be: 'pearson', 'spearman', or 'kendall'
# fns = [None,'IG','MS']
fns = [None] # List of Notch sensitivity option to run. 
# None completely ignores the notch when aggregating data. 
# IG forces aggregation to only periods within a given band that are sensitive to the IG wave.
# MS runs the code only at periods outside the IG sensitivity.
for fn in fns:
    for corrmethod in corrmethods:
        for b in bands:
            print(f'{fn} : {corrmethod.upper()} : {b[0]}-{b[1]}s')
            bn='_'.join(tuple(map(str,b)))
            yttl = lambda c: fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"

            Xin = [xi() for xi in X]
            C,fig = tri_correlogram(Xin, names, method=corrmethod, figscaling=2.2)

            if fn=='IG':fig.suptitle(f'{b[0]} to {b[1]}s below the notch',y=1.01,fontweight='bold')
            elif fn=='MS':fig.suptitle(f'{b[0]} to {b[1]}s above the notch',y=1.01,fontweight='bold')
            else:fig.suptitle(f'{b[0]} to {b[1]}s regardless of notch',y=1.1,fontweight='bold')

            xcnames = names.copy()
            for i,j in enumerate(xcnames):
                for k in ['$','mathrm','underset','{','}',';','\n','\;']:
                    xcnames[i]=xcnames[i].replace(k,'')
            A=np.asarray(C);rows=list(xcnames);cols=list(xcnames)

            df=pd.DataFrame(A, index=rows, columns=cols)
            df.index.name="row_cat"; df.columns.name="col_cat"

            fold=plotfolder;fold.mkdir(parents=True,exist_ok=True)

            fold = fold/f'_{corrmethod}'
            if fn is not None:fold=fold/f'_{fn}'
            fold.mkdir(parents=True,exist_ok=True)
            (fold/'_xls').mkdir(parents=True,exist_ok=True)

            prefix = (fn or "").replace("IG","Notched.").replace("MS","MS.")
            file = f"S03.Table02_{prefix}{corrmethod}.correlogram.{len(X)}x{len(X)}.{bn}.{save_format}"
            df.to_excel(fold/'_xls'/file.replace(f'.{save_format}','.xlsx'))
            _=save_tight(fold/file,dpi=50)
            plt.close('all')
