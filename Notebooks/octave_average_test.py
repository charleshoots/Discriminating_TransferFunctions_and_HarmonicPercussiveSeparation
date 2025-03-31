# load_pickle(list((noisedir.parent/'AVG_STA'/stakey).glob('*sta.pkl'))[0])
get_sta_noise = lambda stanm:load_pickle(list((dirs.SpectraAvg/stanm).glob('*avg_sta.pkl'))[0])
# [f,w in zip(noise.day_files[noise.gooddays],noise.day_goodwins[noise.gooddays])]
get_day_noise = lambda stanm,days:[load_pickle(dirs.Spectra/stanm/f) for f in days]


sta_noise=get_sta_noise(stanm)
day_noise=get_day_noise(stanm,sta_noise.day_files[sta_noise.gooddays])
wins=sta_noise.day_goodwins[sta_noise.gooddays]


f=sta_noise.f
ind = (f>0) & (f<=1)
f = f[ind]
# sta_noise.power.cZZ

fig,axes=plt.subplots(1,2,figsize=(13,3))

# unit_scale=-((2*np.pi*f)**2)
unit_scale=40*np.log10(2*np.pi*f,where=f>0.)

# d=np.abs(sta_noise.power.cZZ[ind])
d=sta_noise.power.cZZ[ind]
# d=d-np.mean(d)
d=10*np.log10(d)
d = d+np.mean(d)
d=d+unit_scale
# d = d+np.mean(d)
of,od = octavg(d,f)


ax = axes[0]
ylim=[-150,-70]
# d=sta_noise.power.cZZ[ind]
ax.plot(1/of,od,linewidth=1)
# ax.invert_xaxis()
# ax.set_xscale('log')
ax.set_xlabel('Period',**labelprops)
ax.set_xlim(.01,1/of[0])
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_title('1/8th Octave Average of Station Average Z PSD',**labelprops)
# ax.set_ylim(-130,-60)
ax.set_ylim(ylim);ax.set_yticks(np.arange(-150,-40,10));ax.grid('log')
ax = axes[1]
ax.plot(1/f,d,linewidth=1)
ax.set_xscale('log')
# ax.set_yscale('log')
# ax.invert_xaxis()
ax.set_xlabel('Period',**labelprops)
ax.set_xlim(.01,1/of[0])
ax.set_title('Normal Average of Station Average Z PSD',**labelprops)
# ax.set_ylim(-130,-60)
ax.set_ylim(ylim);ax.set_yticks(np.arange(-150,-40,10));ax.grid('log')
