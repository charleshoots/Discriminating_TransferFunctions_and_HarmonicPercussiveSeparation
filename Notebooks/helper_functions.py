import numpy as np
from imports import *
from modules import *
def save_tight(filename,fig=None,dpi=200,format=None):
        # Saves figure to PDF with no margins. Do not modify
        # plt.gca().set_axis_off()
        # plt.subplots_adjust(top = 1, bottom = 0.0, right = 1, left = 0,hspace = 0.07, wspace = 0.03)
        plt.margins(0.1,0.1)
        # plt.gca().xaxis.set_major_locator(plt.NullLocator())
        # plt.gca().yaxis.set_major_locator(plt.NullLocator())
        if fig is None:
                plt.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
                print('Complete')
        else:
                fig.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
                print('Complete')
def fnotch(d):
        '''The frequency knotch root function described in Crawford et al., 1998.
        depth (d) is in meters. Returned (f) is in Hz.'''
        g = 9.80665
        f = (g/(2*np.pi*d))**0.5
        return f
# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| HPS Spectrogram Plots (Replicates Zali) |||||||||||||||||||||||||
# _________________________________________________________________________________________________________________
def hps_spectrograms(S_full, S_background, S_hps, frequencies, times,figsize=(11,12)):
        t = times
        f = frequencies
        s = S_full
        cmap = 'magma'
        xlabel = True
        yscale = 'log'
        ax = None
        vmin,vmax = None,None
        fig, axes = plt.subplots(nrows=3, ncols=1,figsize=figsize,height_ratios=[1,1,1],width_ratios=[1],layout='constrained',squeeze=False,sharey='col',sharex='row')
        titles = ['Raw','Noise','Noise Removed']
        for r,s in enumerate([S_full, S_background, S_hps]):
                gax = axes[r,0]
                pc = gax.pcolormesh(t,f, 10*np.log10(s), cmap = cmap, shading= 'auto')
                if (vmin is not None) and (vmax is not None):
                        pc.set_clim(vmin, vmax)
                else:
                        vmin,vmax = pc.get_clim()
                if ylabel:
                        gax.set_ylabel('Frequency (Hz)',fontweight='bold')
                if xlabel:
                        gax.set_xlabel('Time (s)',fontweight='bold')
                gax.set_yscale(yscale)
                gax.set_ylim(f[1],f[-1])
                gax.set_xlim(t[0],t[-1])
                gax.set_title(titles[r],fontweight='bold')
                if fig is not None:
                        fig.colorbar(pc, ax=gax, pad=0.01, label='dB')

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| EVENT RECORD FUNCTION |||||||||||||||||||||||||||||||||||||||||||
# _________________________________________________________________________________________________________________
def event_record_plot(evstream,evstream_back=None,linewidth=0.2,trim = (None,None),scales = [1,1],band = None,facecolor=('b','r'),norm = 'trace',figsize = (20,13),phases = ('P','S','SKS','PKiKP','SKiKS','SKSSKS',),evdepth=None,title='',sortindex=None,ax=None,normscale=1.0,residual_fraction=1.0):
        # linewidth=0.2
        # prepare traces
        if evstream_back:
            sets = [evstream_back,evstream] # [UNCORRECT_SET , CORRECT_SET] 
        else:
            sets = [evstream]
        sets = [preparetraces(stream,trim=trim,band=band,sortindex=sortindex) for stream in sets]
        if len(sets)==2:
            residuals = [ev0.data - ev.data for ev0,ev in zip(sets[0],sets[1])]
            residuals = [normscale * residual_fraction * np.array(res) / np.max(np.abs(np.array(res))) for res in residuals]
        fig = None
        if ax is None:
            fig, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
            ax = axes[0,0]
        # norming for plots
        postsetmax = [np.max(ev.data) for ev in  sets[-1]]
        for si in range(len(sets)):
            for i in range(len(sets[0])): #norm the uncorrected set to the max in the corrected set
                # -----------
                if isinstance(norm,list):
                        sets[si][i].data = sets[si][i].data/norm[i]
                elif norm.lower()=='postset':
                        sets[si][i].data = sets[si][i].data/postsetmax[i]
                elif norm.lower()=='trace':
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                elif norm.lower()=='col':
                        # print('Trace norm scaling by r')
                        dist = [s.stats.sac.dist for s in sets[si]]
                        norms = [d/np.max(dist) for d in dist]
                        norms = [d**(1) for d in norms]
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                        sets[si][i].data = sets[si][i].data / norms[i]
                        colmax = np.max([abs(d.data.max()) for d in sets[si]])
                        sets[si][i].data = sets[si][i].data / colmax
                # -----------
        correctset = sets[-1]
        [ax.plot(tr.times(),tr.data + ysep,linewidth=linewidth,color='k') for ysep,tr in enumerate(correctset)]
        if len(sets)>1:
            for si,s in enumerate(sets):
                [ax.fill_between(tr.times(),tr.data*scales[si] + ysep,tr.data*0 + ysep, where=np.abs(tr.data)>=0, facecolor=facecolor[si]) for ysep,tr in enumerate(s)]
        [ax.plot(tr.times(),tr.data*0 + ysep,linewidth=0.4,color='k') for ysep,tr in enumerate(correctset)]
        if evdepth is not None:
            arrivals = [ObsQA.io.get_arrivals(sta_llaz=(sta.stats.sac.stla,sta.stats.sac.stlo,sta.stats.sac.stel),ev_llaz=(sta.stats.sac.evla,sta.stats.sac.evlo,evdepth),phases=phases) for sta in correctset]
            ardict = dict()
            corephase_dy = 0.2
            direcphase_dy = 0.03
            [[ardict.update({ph[0]:[]}) for ph in a] for a in arrivals]
            [[ardict[ph[0]].extend([ph[1]]) for ph in a] for a in arrivals]
            [ardict.update({k:np.max(ardict[k])}) for k in list(ardict.keys())]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='b',linewidth=0.2) for ph in a] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='P'] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='S'] for ysep,a in enumerate(arrivals)]
            [ax.text(ardict[k], len(correctset) - corephase_dy, k, color='b',horizontalalignment='center') for k in list(ardict.keys()) if (k!='P') and (k!='S')]
            [ax.text(ardict[k], len(correctset) + direcphase_dy, k, color='r',horizontalalignment='center') for k in list(ardict.keys()) if ((k=='P') or (k=='S'))]
        yl = (-1,len(correctset))
        ax.set_ylim(yl)
        ax.set_xlim(correctset[0].times()[0],correctset[0].times()[-1])
        ax.set_yticks([i for i in range(len(correctset))])
        labels = [str(int(ev.stats.sac.dist)) +'km' + ' [' + ev.stats.network + '] ' + ev.stats.station  + '\n depth:' + str(int(abs(ev.stats.sac.stel*1000))) + 'm' for ev in correctset]
        ax.set_yticklabels(labels)
        ax.set_xlabel('Time(s)')
        if fig is not None:
            fig.suptitle(title,fontweight='bold',fontsize=15)
        return ax