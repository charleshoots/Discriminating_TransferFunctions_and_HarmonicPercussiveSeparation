from imports import *
from quick_class import *
import warnings
import fnmatch
from obspy.core.inventory.inventory import read_inventory


HPSDataFolder = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp/_DataArchive/HPS_Data')

# get_noisecut_event(parentfolder,staname,event,channel=['*1','*2','*Z','*H'],len_hrs=2,pre_trim=False,post_trim=True,win_length=163.84,width=None)
ovr = False
cat = catalog.copy()
for stai,Station in enumerate(cat.iloc):
    for evi,Event in enumerate(Station.Events):
        stanm = Station.StaName
        print('||'*30)
        print('||'*5)
        print('S ' + str(stai+1) + '/' + str(len(cat)) + ' | E ' + str(evi+1) + '/' + str(len(Station.Events)) + ' | ' + stanm)
        print('||'*5)
        print('||'*30)
        sacoutfold = HPSDataFolder / 'Output' / stanm
        sacoutfold.mkdir(exist_ok=True,parents=True)
        (sacoutfold / 'Plots').mkdir(exist_ok=True,parents=True)
        sacfile = '.'.join([stanm,Event,'HZ.SAC'])
        if not ovr:
            if (sacoutfold/sacfile).exists():
                print(stanm + '| OVR Disabled. Output already exists. Skipping')
                continue
        out = get_noisecut_event(HPSDataFolder,stanm,Event)
        if len(out)==0:
            print(stanm+' | '+Event+' | File cut incorrectly. Skipping')
            continue
        st = out.Corrected[0].copy()
        if (st[0].stats.endtime - st[0].stats.starttime)<7200:
            print(stanm+' | '+Event+' | Unexpected trace length. Skipping')
            continue
        # if st[0].
        raw = out.Raw[0].copy()
        print(stanm,'| Saving output:',sacfile)
        st.write(str(sacoutfold/sacfile), format='SAC')
        fig,ax = plt.subplots(nrows=2,ncols=1,figsize=(13,7),sharex='all',sharey='all',squeeze=True)
        ax[0].plot(raw[0].times(),raw[0].data,c='k',linewidth=0.3)
        ax[0].set_title('Raw | ' + stanm + ' | ' + Event)
        ax[1].plot(st[0].times(),st[0].data,c='k',linewidth=0.3)
        ax[1].set_title('Corrected | ' + stanm + ' | ' + Event)
        ax[0].set_xlim(st[0].times()[0],st[0].times()[-1])
        save_tight(str(sacoutfold / 'Plots' / sacfile.replace('.SAC','.png')),dpi=600)
        plt.close('all')
        #|||| Spectrograms are huge files at 24hrs, dont save
        # specoutfold = HPSDataFolder / 'Output' / stanm / 'Spectrograms'
        # specoutfold.mkdir(exist_ok=True,parents=True)
        # specfile = sacfile.replace('.SAC','.pkl')
        # print(stanm,'| Saving spectrograms:',specfile)
        # out.to_pickle(str(specoutfold / specfile))
        k=1
        del st,out
print('--Complete--')