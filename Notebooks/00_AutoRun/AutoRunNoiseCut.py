import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory
from obspy.core.util.attribdict import AttribDict as attr

HPSDataFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data')
HPSDataFolder=dirs.Events_HPS
# get_noisecut_event(parentfolder,staname,event,channel=['*1','*2','*Z','*H'],len_hrs=2,pre_trim=False,post_trim=True,win_length=163.84,width=None)
ovr = False
cat = catalog.copy()
def unravel(lst):return list(itertools.chain.from_iterable(lst))
# cat = cat.iloc[np.min(np.where(catalog.Network=='XE')):]
channel=['*Z','*1','*2','*H']
# channel=['*Z','*1','*2','*H']
baddata=0;badcut=0;unexpectedlen=0
for stai,Station in enumerate(cat.iloc):
    _=os.system('cls' if os.name == 'nt' else 'clear')
    for evi,Event in enumerate(Station.Events):
        stanm = Station.StaName
        sacoutfold = HPSDataFolder  / 'corrected' / stanm
        sacoutfold.mkdir(exist_ok=True,parents=True)
        (sacoutfold / 'Plots').mkdir(exist_ok=True,parents=True)
        sacfile = '.'.join([stanm,Event.Name,'SAC'])
        all_exist = np.all([len(list(sacoutfold.glob(sacfile.replace('.SAC',c+'.SAC'))))>0 for c in channel])
        state = f'S {str(stai+1)}/{str(len(cat))} | E {str(evi+1)}/{str(len(Station.Events))} | {stanm} | Dead trace : {baddata} | Unexpected length : {unexpectedlen} | Bad cut : {badcut}'
        # state = f'S {str(stai+1)}/{str(len(cat))} | E {str(evi+1)}/{str(len(Station.Events))} | {stanm} | Bad data : {baddata}'
        if (not ovr) & all_exist:
            # print(state + '| Exists')
            continue
        print(state)
        fin = HPSDataFolder/'rmresp'/stanm/f'{Event.Name}.HZ.SAC'
        if not fin.exists():print(f'FILE: {fin} | Does not exist. Continuing');continue
        try:out = (get_noisecut_event(HPSDataFolder/'rmresp',stanm,Event.Name,channel=channel))[0]
        except:baddata+=1;continue
        if len(out)==0:badcut+=1;continue
        st = out.Corrected[0].copy()
        if (st[0].stats.endtime - st[0].stats.starttime)<(7200-500):unexpectedlen+=1;continue
        print(stanm,' ------ | Saving output:',sacfile)
        _=[tr.write(str(sacoutfold/sacfile.replace('.SAC','.'+tr.stats.channel+'.SAC')),format='SAC') for tr in st]
        raw = out.Raw[0].copy()
        # os.rename(str(f),str(f).replace(f'{c}.SAC',f'.{c}.SAC'))
        for r,s in zip(raw,st):
            fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(13,7),sharex='all',sharey='all',squeeze=True)
            ax[0].plot(r.times(),r.data,c='k',linewidth=0.3)
            ax[0].set_title(r.stats.channel + 'Raw | ' + stanm + ' | ' + Event.Name)
            ax[1].plot(s.times(),s.data,c='k',linewidth=0.3)
            ax[1].set_title(s.stats.channel + '| Corrected | ' + stanm + ' | ' + Event.Name)
            ax[0].set_xlim(s.times()[0],s.times()[-1])
            save_tight(str(sacoutfold / 'Plots' / sacfile.replace('.SAC','.'+s.stats.channel+'.png')),dpi=600)
            plt.close('all')
        #|||| Spectrograms are huge files at 24hrs, dont save
        # specoutfold = HPSDataFolder / 'Output' / stanm / 'Spectrograms'
        # specoutfold.mkdir(exist_ok=True,parents=True)
        # specfile = sacfile.replace('.SAC','.pkl')
        # print(stanm,'| Saving spectrograms:',specfile)
        # out.to_pickle(str(specoutfold / specfile))
        del st,out
print('--Complete--')