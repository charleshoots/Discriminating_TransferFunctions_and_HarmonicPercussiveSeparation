import numpy as np
from typing import Optional, Sequence, Tuple, Union, Literal
import copy
AggFn = Literal["mean","median","sum","min","max","std"]
BandRange = Union[Tuple[float,float], Sequence[float], None]
try:
    import bottleneck as bn
    _fast_nanmedian=lambda a,axis: bn.nanmedian(a,axis=axis)
except Exception:
    _fast_nanmedian=lambda a,axis: np.nanmedian(a,axis=axis)

def _as_pair(b:BandRange)->Optional[Tuple[float,float]]:
    if b is None:return None
    lo,hi=float(b[0]),float(b[1]);return (lo,hi) if lo<=hi else (hi,lo)

def _safe_attr(name:str)->Optional[str]:
    s=str(name).replace(" ","_").replace(".","_")
    if not s or not (s[0].isalpha() or s[0]=="_"):return None
    return s if all(ch.isalnum() or ch=="_" for ch in s) else None

def fnotch(d,n=1):
        '''The frequency knotch root function described in Crawford et al., 1998.
        depth (d) is in meters. Returned (f) is in Hz.'''
        assert (0.5<=n) & (n<=2), f'n={n} violates 0.5<=n<=2'
        g = 9.80665
        f = (g/(2*np.pi*n*d))**0.5
        return f

def _preferred_nominal_centers(fraction: int):
    if fraction == 3:
        return 10
    return 3 * int(fraction)

def _build_centers_ansi(fmin, fmax, fraction=8, ref=1000.0, centers="nominal"):
    if centers not in ("exact", "nominal"):
        raise ValueError("centers must be 'exact' or 'nominal'")
    if centers == "exact":
        k_min = int(np.ceil(fraction * np.log2(fmin / ref)))
        k_max = int(np.floor(fraction * np.log2(fmax / ref)))
        k = np.arange(k_min, k_max + 1, dtype=int)
        fc = ref * (2.0 ** (k.astype(float) / float(fraction)))
        return fc
    d = _preferred_nominal_centers(fraction)
    n_min = int(np.ceil(d * np.log10(fmin)))
    n_max = int(np.floor(d * np.log10(fmax)))
    n = np.arange(n_min, n_max + 1, dtype=int)
    fc = 10.0 ** (n.astype(float) / float(d))
    return fc

def _band_edges_from_centers_ansi(fc, fraction):
    q = 2.0 ** (0.5 / float(fraction))
    f1 = fc / q
    f2 = fc * q
    return f1, f2

def octave_average_ansi(Y, f, fraction=8, fmin=None, fmax=None, agg="mean", nan_ok=True,
                        centers="nominal", ref=1000.0, return_edges=False):
    Y = np.asarray(Y)
    f = np.asarray(f, float)
    if f.ndim != 1 or f.size == 0:
        raise ValueError("f must be a non-empty 1-D array.")
    if np.any(f <= 0):
        raise ValueError("f must be strictly positive.")
    if Y.shape[-1] != f.size:
        raise ValueError("The last dimension of Y must match len(f).")
    order = np.argsort(f)
    f = f[order]
    Y = Y[..., order]
    if fmin is None: fmin = float(f[0])
    if fmax is None: fmax = float(f[-1])
    if not (fmin > 0 and fmax > fmin):
        raise ValueError("Invalid fmin/fmax.")
    fc = _build_centers_ansi(fmin, fmax, fraction=fraction, ref=ref, centers=centers)
    if fc.size == 0:
        raise ValueError("No band centers within [fmin, fmax].")
    f1, f2 = _band_edges_from_centers_ansi(fc, fraction)
    keep = (f2 >= f[0]) & (f1 <= f[-1])
    fc, f1, f2 = fc[keep], f1[keep], f2[keep]
    if fc.size == 0:
        raise ValueError("All bands fell outside data frequency span.")
    if   agg == "mean":   red = (np.nanmean if nan_ok else np.mean)
    elif agg == "median": red = (np.nanmedian if nan_ok else np.median)
    elif agg == "sum":    red = (np.nansum  if nan_ok else np.sum)
    elif agg == "min":    red = (np.nanmin  if nan_ok else np.min)
    elif agg == "max":    red = (np.nanmax  if nan_ok else np.max)
    elif agg == "std":    red = (np.nanstd  if nan_ok else np.std)
    else:
        raise ValueError("agg must be one of mean, median, sum, min, max, std")
    lead = Y.shape[:-1]
    Yb = np.empty(lead + (fc.size,), dtype=float)
    for i in range(fc.size):
        lo, hi = f1[i], f2[i]
        sel = (f >= lo) & (f < hi) if i < fc.size - 1 else (f >= lo) & (f <= hi)
        if not np.any(sel):
            Yb[..., i] = np.nan
        else:
            Yb[..., i] = red(Y[..., sel], axis=-1)
    if return_edges:
        return fc, Yb, f1, f2
    return fc, Yb

def octave_average(d, f, fmin=1/100, fmax=1, fraction=8, domain="geo", N=291,
                   centers="nominal", ref=1000.0, agg="mean", nan_ok=True):
    d = np.asarray(d)
    f = np.asarray(f, float)
    if domain == "db":      x = 10.0**(d/10.0)
    elif domain == "log":   x = np.exp(d)
    else:                   x = d
    fc, xb = octave_average_ansi(x, f, fraction=fraction, fmin=fmin, fmax=fmax,
                                 agg=agg, nan_ok=nan_ok, centers=centers, ref=ref)
    if domain == "db":      yb = 10.0*np.log10(xb)
    elif domain == "log":   yb = np.log(xb)
    else:                   yb = xb
    return fc, yb
class AggregateMeasurements:
    def __init__(self,D:Optional[np.ndarray]=None,bands:Optional[Sequence[float]]=None,
        phases:Optional[Sequence[Union[str,int]]]=None,**kwargs):
        self.D=None if D is None else np.asarray(D)
        self.bands=None if bands is None else np.asarray(bands)
        self.phases=phases if phases is not None else None
        self.comp = kwargs.get('comp','Original_Z')
        self.depth = None if ('depth' not in kwargs or kwargs['depth'] is None) else np.asarray(kwargs['depth'])
        self._parent=None; self._name=None; self._idx_cache={}
        for k,v in kwargs.items():
            if isinstance(v,AggregateMeasurements): v._parent=self; v._name=k
            setattr(self,k,v)
        if self.D is not None:
            if self.D.ndim not in (2,3): raise ValueError(f"`D` must be 2D or 3D; got {self.D.shape}.")
            if self.bands is None: raise ValueError("`bands` is required.")
            nb=self.D.shape[1]
            if self.bands.ndim!=1: raise ValueError("`bands` must be 1D.")
            if len(self.bands) not in (nb,nb+1):
                raise ValueError(f"`bands` len must be {nb} or {nb+1}; got {len(self.bands)}.")
            if self.depth is not None and len(self.depth)!=self.D.shape[0]:
                raise ValueError("`depth` length must match number of source-receivers (D.shape[0]).")
        self._phase_children={}
        if self.D is not None and self.D.ndim==3 and self.phases is not None:
            if len(self.phases)!=self.D.shape[2]: raise ValueError("`phases` must match D.shape[2].")
            for i,name in enumerate(self.phases):
                child=AggregateMeasurements(D=self.D[:,:,i],bands=self.bands,comp=self.comp,depth=self.depth)
                child._parent=self; child._name=str(name)
                self._phase_children[str(name)]=child
                safe=_safe_attr(str(name))
                if safe and not hasattr(self,safe): setattr(self,safe,child)


    def _band_idx(self,p:Optional[Tuple[float,float]])->np.ndarray:
        if self.D is None or self.bands is None: raise ValueError("No data/bands.")
        nb=self.D.shape[1]
        if p is None: return np.arange(nb)
        key=(p[0],p[1])
        hit=self._idx_cache.get(key)
        if hit is not None: return hit
        lo,hi=p
        if len(self.bands)==nb:
            idx=np.flatnonzero((self.bands>=lo)&(self.bands<=hi))
        else:
            L,R=self.bands[:-1],self.bands[1:]
            idx=np.flatnonzero((R>=lo)&(L<=hi))
        self._idx_cache[key]=idx
        return idx

    def Average(self, band: BandRange = (1, 100), agg: AggFn = "mean", nan_ok: bool = True,
                log10: bool = False, fn: Optional[str] = None, octave: bool = True
                ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        if self.D is None: raise ValueError("No data in `D`.")
        p = _as_pair(band); idx = self._band_idx(p)
        if idx.size == 0: raise ValueError(f"No bands in range {band}.")
        a = self.D.take(idx, axis=1)

        # period centers for all bands (used by gating and/or octave averaging)
        nb = self.D.shape[1]
        if len(self.bands) == nb:
            centers = np.asarray(self.bands, dtype=float)
        else:
            L, Rb = np.asarray(self.bands[:-1], float), np.asarray(self.bands[1:], float)
            centers = 0.5 * (L + Rb)
        csel = centers[idx]  # (B_sel,)

        # Optional MS/IG gating by per-SR fnotch (period = 1/fnotch(depth))
        if fn is not None:
            if self.depth is None: raise ValueError("`fn` requires `depth` array on the instance.")
            Tn = 1.0 / fnotch(self.depth)  # (SR,)
            if   fn == 'MS': mask = (csel[None, :] <  Tn[:, None])   # periods shorter than notch
            elif fn == 'IG': mask = (csel[None, :] >= Tn[:, None])   # periods longer  than notch
            else: raise ValueError("fn must be None, 'MS', or 'IG'")
            a = np.where(mask[:, :, None] if a.ndim == 3 else mask, a, np.nan)

        T_new = None

        # Optional octave-band averaging BEFORE the final aggregation
        if octave:
            fvec = 1.0 / csel  # period → frequency
            if a.ndim == 3:
                SR, B, Pn = a.shape
                a2 = a.transpose(0, 2, 1).reshape(SR * Pn, B)   # (SR*P, B)
                fc, ao = octave_average(a2, fvec, fraction=8, centers="nominal")   # ANSI path
                Bn = ao.shape[1]
                a = ao.reshape(SR, Pn, Bn).transpose(0, 2, 1)   # (SR, B_new, P)
            else:
                SR, B = a.shape
                fc, a = octave_average(a, fvec, fraction=8, centers="nominal")      # ANSI path
            T_new = 1.0 / fc                                    # back to period

        # Final aggregation across bands (axis=1)
        if   agg == "mean":   f = (np.nanmean if nan_ok else np.mean)
        elif agg == "median": f = (lambda x, axis: _fast_nanmedian(x, axis)) if nan_ok else np.median
        elif agg == "sum":    f = (np.nansum  if nan_ok else np.sum)
        elif agg == "min":    f = (np.nanmin  if nan_ok else np.min)
        elif agg == "max":    f = (np.nanmax  if nan_ok else np.max)
        elif agg == "std":    f = (np.nanstd  if nan_ok else np.std)
        else: raise ValueError("agg must be one of mean, median, sum, min, max, std")

        out = f(a, axis=1)
        if log10: out = np.where(out > 0, np.log10(out), np.nan)

        # return (out, T_new) if octave else out
        return out
    def Octave(self, band: BandRange = None, fraction: int = 8, centers: str = "nominal",
               fmin: float = None, fmax: float = None, domain: str = "geo",
               agg: AggFn = "mean", nan_ok: bool = True):
        """
        Octave-band the container's data using ANSI/IEC fractional-octave bands.

        Parameters
        ----------
        band : (Tmin, Tmax) in seconds (period). If None, use full band extent.
        fraction : int
            Fractional octave denominator (e.g., 8 for 1/8-oct).
        centers : {"nominal","exact"}
            ANSI/IEC nominal preferred numbers (default) or exact base-2 centers.
        fmin, fmax : float, optional
            Band-limiting in Hz (frequency domain). If None, inferred from data.
        domain : {"geo","db","log"}
            Pre/post domain handling:
              - "geo": linear averaging (default)
              - "db" : average in linear power, return in dB
              - "log": average in linear amplitude, return in log
        agg : AggFn
            Reducer within each band (mean/median/sum/min/max/std).
        nan_ok : bool
            If True, use NaN-aware reducers.

        Returns
        -------
        Yb : np.ndarray
            Octave-banded data with shape:
              - 2D input: (SR, B_new)
              - 3D input: (SR, B_new, P)  (phases preserved)
        T_new : np.ndarray
            New period-band centers (seconds), ascending.
        """
        if self.D is None: raise ValueError("No data in `D`.")
        if self.bands is None: raise ValueError("`bands` required.")

        # Select band indices
        idx = self._band_idx(_as_pair(band)) if band is not None else self._band_idx((min(self.bands), max(self.bands)))
        if idx.size == 0: raise ValueError("No bands selected.")

        # Slice data along band axis
        A = self.D.take(idx, axis=1)

        # Period centers for selection
        nb = self.D.shape[1]
        if len(self.bands) == nb:
            centers_T = np.asarray(self.bands, dtype=float)
        else:
            L, Rb = np.asarray(self.bands[:-1], float), np.asarray(self.bands[1:], float)
            centers_T = 0.5 * (L + Rb)
        T_sel = centers_T[idx]  # (B_sel,)
        fvec = 1.0 / T_sel      # Hz

        # Prepare input for octave_average (expects last axis = frequency)
        if A.ndim == 3:
            SR, B, Pn = A.shape
            A2 = A.transpose(0, 2, 1).reshape(SR * Pn, B)  # (SR*P, B)
            fc, Yb2 = octave_average(A2, fvec, fraction=fraction, centers=centers,
                                      fmin=fmin, fmax=fmax, domain=domain, agg=agg, nan_ok=nan_ok)
            Bn = Yb2.shape[1]
            Yb = Yb2.reshape(SR, Pn, Bn).transpose(0, 2, 1)  # (SR, B_new, P)
        else:
            fc, Yb = octave_average(A, fvec, fraction=fraction, centers=centers,
                                     fmin=fmin, fmax=fmax, domain=domain, agg=agg, nan_ok=nan_ok)

        T_new = 1.0 / fc  # back to period (s)
        return Yb, T_new


    def R(self, comp=None, eps:float=0.0, log10:bool=True)->"AggregateMeasurements":
        if comp is None: comp=self.comp
        if isinstance(comp,str):
            if self._parent is None: raise ValueError("No parent to resolve comparison.")
            if not hasattr(self._parent,comp): raise AttributeError(f"Comparison '{comp}' not found on parent.")
            other=getattr(self._parent,comp)
        else:
            other=comp
        if not isinstance(other,AggregateMeasurements) or other.D is None: raise ValueError("Invalid comparison target.")
        if self.D is None: raise ValueError("No data to ratio.")
        if self.D.shape!=other.D.shape: raise ValueError(f"Shape mismatch: {self.D.shape} vs {other.D.shape}")

        # optional: depth sanity check (keeps behavior minimal)
        if getattr(other,'depth',None) is not None and self.depth is not None:
            if len(other.depth)!=len(self.depth): raise ValueError("`depth` length mismatch between operands.")

        with np.errstate(divide="ignore", invalid="ignore"):
            den = other.D if eps==0 else other.D+eps
            out = np.divide(self.D, den)
            if log10: out = np.where(out>0, np.log10(out), np.nan)
            out[~np.isfinite(out)] = np.nan

        return AggregateMeasurements(D=out, bands=self.bands, phases=self.phases, comp=comp, depth=self.depth)


    def copy(self, deep=True):
        return copy.deepcopy(self) if deep else copy.copy(self)
    def __copy__(self):
        out = self.__class__.__new__(self.__class__)
        out.__dict__ = self.__dict__.copy()   # views for np arrays; shared refs for others
        return out
    def __deepcopy__(self, memo):
        out = self.__class__.__new__(self.__class__)
        memo[id(self)] = out
        d = out.__dict__
        for k, v in self.__dict__.items():
            if isinstance(v, np.ndarray):
                d[k] = v.copy()
            elif hasattr(v, "copy") and not isinstance(v, (str, bytes)):
                try: d[k] = v.copy()
                except Exception: d[k] = copy.deepcopy(v, memo)
            else:
                d[k] = copy.deepcopy(v, memo)
        return out

    def __getitem__(self,key:Union[int,str])->"AggregateMeasurements":
        if self.D is None or self.D.ndim!=3: raise TypeError("Phase indexing only on 3D containers.")
        if isinstance(key,int):
            if self.phases is None: return AggregateMeasurements(D=self.D[:,:,key],bands=self.bands,comp=self.comp)
            key=str(self.phases[key])
        try: return self._phase_children[str(key)]
        except KeyError: raise KeyError(f"Phase {key!r} not found. Available: {list(self._phase_children)}")

    @property
    def shape(self): return None if self.D is None else self.D.shape
    def __repr__(self):
        dims="None" if self.D is None else f"{self.D.shape}"
        ph=None if self.phases is None else list(self._phase_children.keys())
        return f"AggregateMeasurements(shape={dims}, bands_len={None if self.bands is None else len(self.bands)}, phases={ph})"
