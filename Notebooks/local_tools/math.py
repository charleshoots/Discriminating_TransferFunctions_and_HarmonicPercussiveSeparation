from modules import *

def fnotch(d):
        '''The frequency knotch root function described in Crawford et al., 1998.
        depth (d) is in meters. Returned (f) is in Hz.'''
        g = 9.80665
        f = (g/(2*np.pi*d))**0.5
        return f
def smooth(d,k=10):return np.convolve(d, np.ones(k) / k, mode='same')
def distance(sta,ev,unit='deg'):
    origins=ev.origins[0]
    stalla,evlla=[sta.Latitude,sta.Longitude],[origins.latitude,origins.longitude]
    dist=locations2degrees(stalla[0],stalla[1],evlla[0],evlla[1])
    if unit.lower()=='km':dist=degrees2kilometers(dist)
    return dist
def detect_outscale(raw,correct,vertical_scale=1.2,ylim=None,suppress=False):
    # Occasionally the ampltide changes (ie noise reduction) after correction is so 
    # significant that it makes it challenging to plot the raw and corrected on the 
    # same plot.This function detects when the scale differences exceed vertical_scale 
    # (120% by default) of the maximum value in the corrected trace. When it occurs, 
    # the distance in amplitudes the raw is from the corrected is shrunk to within this margin.
    if not ylim:ylim = np.array([np.max(np.abs(c.data))*vertical_scale for c in correct])
    out_scaled = np.array([np.max(np.abs(c.data),
    where=~(np.isinf(c.data)+np.isnan(c.data)),
    initial=0) for c in raw]) > ylim
    if np.any(out_scaled):
        if not suppress:print('Large amplitude scale differences detected in '+raw[0].id)
        for tr_ind,tr in enumerate(raw):
            if out_scaled[tr_ind]:tr.data = (tr.data/np.max(np.abs(tr.data)))*ylim[tr_ind]
    return raw
def _calc_phase(ab,**args):
        ph = np.angle(ab,deg=True)
        return ph
def _calc_admittance(ab,bb,**args):
        ad = np.abs(ab)/bb
        return ad
def _calc_coherence(ab,aa,bb,**args):
        coh = ((np.abs(ab)**2)/(aa*bb))
        return coh
def avg_meter(avg,m,r):
    Meters={'Coherence':_calc_coherence,'Admittance':_calc_admittance,'Phase':_calc_phase}
    if not r=='ZP':r=''.join(sorted(r))
    AA=avg.power.__dict__['c'+r[0]+r[0]]
    BB=avg.power.__dict__['c'+r[1]+r[1]]
    AB=avg.cross.__dict__['c'+r[0]+r[1]]
    x=avg.f
    y=Meters[m](ab=AB,aa=AA,bb=BB)
    return x[x>=0],y[x>=0]

def octave_average(d,f,fmin=.002,fmax=1.025,fraction=3):
    # 1/3 octave (fraction=3) → Standard in engineering seismology and ground motion studies.
    # 1/6 octave (fraction=6) → Provides better resolution while still smoothing noise.
    # 1/12 octave (fraction=12) → High resolution, often used for detailed spectral analysis.
    # 1/8 octave (fraction=8) → Sometimes used in geophysical applications for a balance between smoothing and resolution.
    k=np.arange(np.log2(fmin), np.log2(fmax), 1/fraction)
    freq_centers = 2**k
    out=[]
    for i, fc in enumerate(freq_centers):
        f_low = fc / 2**(1/(2*fraction))
        f_high = fc * 2**(1/(2*fraction))
        # Find indices within this band
        indices = np.where((f >= f_low) & (f <= f_high))[0]
        if len(indices) > 0:out.append(np.mean(d[indices]))
    return out