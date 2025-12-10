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


# update rmresp folder
def update_rmresp_folder(datafold,reg='*SAC',ovr=False):
    # rawfold value
    rawfold=Path(datafold)/'raw'
    (Path(datafold)/'rmresp').mkdir(exist_ok=True)
    # files value
    files=list(rawfold.rglob(reg))
    # files value
    files=[f for f in files if (not f.parent.name=='_depreciated')]
    # files value
    files=[f for f in files if (not f.parent.name=='_quarantine')]
    files=[f for f in files if (f.parent.parent.name=='raw')]
    for fi,f in enumerate(files):
        state=f'{np.round(100*((fi+1)/len(files)),3)}% | {str(fi+1).zfill(len(str(len(files))))}/{len(files)}'
        dest=Path(str(f).replace('/raw/','/rmresp/'))
        dest.parent.mkdir(exist_ok=True,parents=True)
        # clear_output()
        if dest.exists():
            if not ovr:print(f'{state} || Already exists');continue
        try:
            tr=load_sac(f,rmresp=True)
            if len(tr.data)>0:
                tr.write(str(dest))
                del tr
                print(f'{state} || Updating.')
            else:print(f'{state} || Failed IO')
        except:
            print(f'{state} || Failed IO')


folders = [dirs.Events, #ATaCR event data directory
dirs.Noise, #ATaCR noise data directory
dirs.Events_HPS, #NoiseCut event (24 hr traces) data directory
]

for datafold in folders:
    update_rmresp_folder(datafold,reg='*SAC',ovr=False)
    k=1