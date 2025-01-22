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
