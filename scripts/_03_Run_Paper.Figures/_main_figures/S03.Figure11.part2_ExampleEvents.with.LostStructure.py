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

from imports import * #Standard imports for doing anything in this project. Approx. ~19 seconds.
from local_tools.quick_class import *
# octavg value
octavg=lt.math.octave_average
# figs value
figs = lambda r=4,c=1,f=(5,6),x='all',y='all',w=None,layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout,width_ratios=np.ones(c) if w is None else w)

# plotfolder value
plotfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure11part2_ExampleRecord.of.LostStructure';plotfolder.mkdir(parents=True,exist_ok=True)
save_format = 'pdf'

# cat value
cat=catalog.copy()
# icat value
icat=cat.sr.copy()
# usnr value
usnr=unpack_metrics(icat)

# apparent killed structure value
apparent_killed_structure=[
['ZA.B01', '2012.302.03.04'],
['ZA.B24', '2012.302.03.04'],
['XF.B19', '2012.302.03.04'],
['ZA.B06', '2012.304.02.49'],
['7D.J26C', '2013.298.17.10'],
['7D.FN12C', '2014.108.14.27'],
]

# idx value
idx=[np.where((cat.sr.StaName==i[0]) & (cat.sr.Name==i[1]))[0][0] for i in apparent_killed_structure]
# icat value
icat=cat.sr.copy()
icat=icat.iloc[idx]
icat.sort_values(by=['StaDepth','Distance'],inplace=True)

bands = [[1,10],[10,30],[30,100],[1,100]]
# bands=[[1,10]]
tmax=3600
mthds=['NoiseCut','ATaCR']
# [2.0,4.2]
# rgvel = 2000
# rgvel = 3100


trace_linewidth = 0.1 #works for png, not pdf
trace_linewidth = 0.2
text_fontsize = 6
for b in bands:
    fig,axes=figs(1,2,f=(6,.4*len(icat)),x=True,y=True)
    for si,sr in enumerate(icat.iloc):
        st=sr.Traces()
        ar = sr.Phases()
        arphases=[ph if ph in list(ar.keys()) else f'{ph}diff' for ph in ['P','S']]
        times=[ar[ph][0] for ph in arphases]
        # dist=distance(sr,sr.Event,'km')*1000
        # r1=dist/rgvel
        for ax,mthd in zip(axes,mthds):
            step=si*2.5
            st.taper(0.001);st.filter('bandpass',freqmin=1/b[1],freqmax=1/b[0]);st.taper(0.001)
            t=st[0].times()
            tidx=t<=tmax
            y=st.select(location='Original')[0].data;norm=np.max(np.abs(y));y=y/norm
            ax.plot(t[tidx],y[tidx]+step,lw=trace_linewidth,c='r')
            y=st.select(location=mthd)[0].data;y=y/norm
            ax.plot(t[tidx],y[tidx]+step,lw=trace_linewidth,c='k')
            ax.text(t[tidx][-1]-30,step+1,rf'${sr.StaName}$ (${sr.StaDepth}m$), $Mw{sr.Magnitude}$',ha='right',fontsize=text_fontsize,)
            ax.text(0.02,0.95,f'{b[0]}-{b[1]}s',ha='left',transform=ax.transAxes,fontsize=text_fontsize,)
            [ax.vlines(t,step-1,step+1,color='k',ls=':',lw=.2) for ph,t in zip(arphases,times)]
            [ax.text(t-10,step+1-.1,ph,c='k',ha='right',va='top',
            fontsize=text_fontsize,
            ) for ph,t in zip(arphases,times)]
            # ax.vlines(r1,step-1,step+1,color='k',ls=':',lw=.2);ax.text(r1-10,step+1-.1,'R1',c='k',ha='right',va='top')
            ax.tick_params(size=2,labelsize=3)
            ax.set_xlabel('time after origin, s',)
    for ax in axes:ax.set_xlim(0,tmax)
    file=f'S03.Figure11.part2_{len(icat)}.pairs.{b[0]}-{b[1]}s.{save_format}'
    save_tight(plotfolder/file,fig,dpi=800)
