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
from source.imports import *
from source.modules import *
from local_tools.quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
from obspy.core.util.attribdict import AttribDict as attr
# function unravel
def unravel(lst):return list(itertools.chain.from_iterable(lst))
# function trim spectrograms
def trim_spectrograms(out,endlen=7200):
    # S value
    S = out.Spectrograms
    # t value
    t = S[0].t
    # twind value
    twind = t>=(t[-1]-endlen)
    # td twind value
    td_twind = (S[0].Original.times()>=(S[0].Original.times()[-1]-endlen))
    for i in range(len(S)):
        S[i].t = S[i].t[twind]
        S[i].Original.data = S[i].Original.data[td_twind]
        S[i].Removed = S[i].Removed[td_twind]
        S[i].Background = S[i].Background[:,twind]
        S[i].Phase = S[i].Phase[:,twind]
        S[i].HPS = S[i].HPS[:,twind]
        S[i].Full = S[i].Full[:,twind]
    return S



# =====OPTIONS=====

# verbose=True
# verbose=False
verbose=[0,5,6,8] #Verbose can be set to True/False for all/none or to a list of specific stages in the NoiseCut script to plot.

# verbose argument is a list of stages in the NoiseCut method to pause the execution for a plot output.clear_output.clear_output
# stages are as follows:
# 0 : Lists all arguments used in calculating the NoiseCut spectrogram.
# 1 : Plots spectrogram of input data before any processing.
# 2 : Same spectrogram of input data but isolated to the band MED algorithm is applied.
# 3 : Same spectrogram of input data but isolated to the band SIM algorithm is applied.
# 4 : W-hat from Zali '23
# 5 : Repeating spectrogram, Eq. 5 from Zali '23.
# 6 : Harmonic (MED output), H from Zali '23.
# 7 : Percussive component from the MED output.
# 8 : N = R + H in Zali '23.
# 9 : S_full - S_background



# Directories for input/output. Must be Path objects.
noisecut_input_dir=dirs.Events_HPS/'rmresp' #Input directory for NoiseCut (after removing instrument response)
# noisecut output dir value
noisecut_output_dir = dirs.Events_HPS/'corrected' #Output directory for NoiseCut results

# ovr value
ovr=True #Skips if output already exists.
# save SAC enabled value
save_SAC_enabled = True #Save outputs to SAC

# channel value
channel=['*Z','*1','*2','*H'] #Channels to process
# save spectrograms value
save_spectrograms=False

# icat value
icat=catalog.sr.copy()
icat.sort_values(by=['Magnitude','StaName'],ascending=False,inplace=True)

# If you are to subset the data catalog to be processed, do it here.
# icat=icat[icat.StaName.isin(['7D.M08A'])]

# ===============



# status value
status = lambda : print(f'SR-Pair | {str(sri+1)}/{str(len(icat.StaName.unique()))} | Event {Event.Name} | Station {stanm}')
for sri,sr in enumerate(icat.iloc):
    # - Source-receiver info
    stanm = sr.StaName
    # Event value
    Event = sr.Event
    status()

    # - Paths and existence checks
    sacoutfold = noisecut_output_dir / stanm
    sacoutfold.mkdir(exist_ok=True,parents=True)
    (sacoutfold/'Plots').mkdir(exist_ok=True,parents=True)
    sacfile ='.'.join([stanm,Event.Name,'SAC'])
    all_exist=np.all([len(list(sacoutfold.glob(sacfile.replace('.SAC',c+'.SAC'))))>0 for c in channel])
    if (not ovr) & all_exist:print(f'EXISTS: {sacoutfold.name}/{sacfile}');continue
    fin=noisecut_input_dir/stanm/f'{Event.Name}.HZ.SAC'
    if not fin.exists():print(f'FILE: {fin} | Does not exist. Continuing');continue

    # - Run noisecut
    out=get_noisecut_event(noisecut_input_dir,stanm,Event.Name,channel=channel,verbose=verbose)[0]

    # - Save spectrograms (optional)
    #   If saving spectrograms: 24hr spectrograms are huge files, so only save last 2hrs.
    if save_spectrograms:
        S=trim_spectrograms(out,7200) #Trim to last 2hrs
        specfold = sacoutfold/'Spectrograms'
        specfold.mkdir(exist_ok=True,parents=True)
        S.to_pickle(str(specfold / sacfile.replace('.SAC','.pkl')))

    # -- NoiseCut output in a single ObsPy Stream object. 
    #       If processed over multiple channels (e.g. Z and H), this Stream will contain all of them.
    st = Stream([c for c in out.Corrected])

    # if len(out)==0:badcut+=1;continue
    assert len(out), 'No data returned after NoiseCut processing.'
    #A basic AllClose check. Three times I notice the end time downloaded was off by as much as ~90s.
    tlen=(st[0].stats.endtime - st[0].stats.starttime)
    assert tlen>(7200-300),'Unexpected length after NoiseCut processing.'
    print(stanm,' ------ | Saving output:',sacfile)

    # Save NoiseCut output trace to SAC
    if save_SAC_enabled:_=[tr.write(str(sacoutfold/sacfile.replace('.SAC','.'+tr.stats.channel+'.SAC')),format='SAC') for tr in st]

    # Verbose plotting of raw vs corrected traces
    if verbose:
        raw = Stream([r for r in out.Raw])
        fig,ax=plt.subplots(nrows=len(raw),ncols=1,figsize=(10,2*len(raw)),sharex='all',sharey=False,squeeze=True)
        [axi.plot(r.times(),r.data,c='r',linewidth=0.3) for r,axi in zip(raw,ax)]
        [axi.plot(r.times(),r.data,c='k',linewidth=0.3) for r,axi in zip(st,ax)]
        [axi.set_title(r.stats.channel) for r,axi in zip(st,ax)]
        fig.suptitle(Event.Name,fontweight='bold')
        save_tight(str(sacoutfold /'Plots'/ sacfile.replace('.SAC','.png')),dpi=600)
        plt.close('all')

    del st,out #Clean up memory and move on to next pair
print('--Complete--')