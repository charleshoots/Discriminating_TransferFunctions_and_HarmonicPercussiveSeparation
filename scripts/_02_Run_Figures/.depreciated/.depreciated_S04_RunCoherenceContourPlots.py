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

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

# figs value
figs = lambda r=3,c=1,f=(5,6),x='all',y='all',layout='constrained':plt.subplots(r,c,figsize=f,sharex=x,sharey=y,layout=layout)
# yttl value
yttl = lambda c:fr"$\underset{{{c}}}{{\gamma\;\;\;\;\;\;\;}}$"
# cohnames2snrnames value
cohnames2snrnames=c2s={'TF':'TF.Z','HPS_Z':'HPS.Z','HPS_H':'HPS.H'}


# icat value
icat=catalog.sr.copy();note=''



# plotfold value
plotfold = dirs.P01.S04
icat.sort_values(by='StaDepth',inplace=True)


# icat=icat[icat.Magnitude>=7.0];note='7to8'
# icat=icat[icat.Magnitude<=7.0];note='6to7'


fig,axes=figs(f=[6,6])

# usnr value
usnr=unpack_metrics(icat)
# f value
f=1/usnr.coh.__dict__['TF_Z'].bands
z=np.array(icat.StaDepth.tolist())
for m,ax in zip(['TF_Z','HPS_Z','HPS_H'],axes):
    if m=='HPS_H':coh=np.array([usnr.coh.__dict__['HPS_1'].D,usnr.coh.__dict__['HPS_2'].D]).mean(axis=0)
    else:coh=usnr.coh.__dict__[m].D #data

    dataset_averaged_coherence_plot(f,z,coh,ax=ax) #plot

    axr = ax.twinx();axr.set_ylim(ax.get_ylim());axr.set_yticks([])
    axr.yaxis.set_label_position("right")
    axr.set_ylabel(f'{yttl(m.replace('_',' '))}', labelpad=16,fontsize=9,rotation=-90)


[ax.set_xlabel(None) for ax in axes[:-1]]

cbar_axes = [a for a in fig.axes if a.get_label() == '<colorbar>']
cbar_axes[0].remove()
cbar_axes[2].remove()
cax=cbar_axes[1]
pos = cax.get_position() # [x0, y0, width, height] in figure coords
new = [pos.x0+0.0,0.02,pos.width*0.21,pos.height*3 + .15]
cax.set_position(new)
cax.tick_params(labelsize=8)
cax.set_ylabel("Coherence", fontsize=9)

plotfold=dirs.Plots/'_Papers'/'Source';plotfold.mkdir(parents=True,exist_ok=True)
file=plotfold/f'{f'{note}_' if len(note)>0 else ''}Coherence.StationEventAverage.png'
save_tight(file,fig,dpi=800)
plt.close('all')
print('DONE')