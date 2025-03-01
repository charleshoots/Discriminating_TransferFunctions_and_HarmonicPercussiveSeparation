from imports import *
from modules import *
# from imports import *
# from modules import *
def catalog_map_plot(sta_inv,ev_cat,
    mode='Depth',
    water_fill_color='verylightblue',
    continent_fill_color='black',
    projection = 'ortho',
    figsize = (50,8),
    cmap=None,content='both',file=None,markers=None,colors=None):
    n=3
    # figsize = (18*n,2.5*n)
    figsize = (9.6,5) #(19.2,10) The size of my screen
    if projection=='local':figsize = (24,7);station_marker_size=20;event_size_downscale=5;latlon_lw=0.2
    if projection=='global':figsize = (24,7);station_marker_size=20;event_size_downscale=5;latlon_lw=0.2
    if projection=='ortho':figsize = (20,20);station_marker_size=300;event_size_downscale=0.5;latlon_lw=0.2
    inv = sta_inv.copy()
    evm_plottable = ev_cat.copy()
    if cmap==None:cmap={'mag':'bwr','depth':'tab20c'}[mode.lower()]
    if mode.lower()=='mag':
        for e in evm_plottable:e.origins[0].depth = e.magnitudes[0].mag*1000
    clear_output(wait=False)
    nev = len(evm_plottable)
    nsta = len(np.unique(list(itertools.chain.from_iterable([e.Stations for e in evm_plottable]))))
    titles = []

    # water_fill_color='darkturquoise'
    # water_fill_color='lightcyan';continent_fill_color='darkslategrey'
    water_fill_color='lightcyan';continent_fill_color='dimgrey'
    cmap='tab20c'
    cmap=mcolors.ListedColormap([mpl.colormaps['tab20b'].resampled(100)(i) for i in np.linspace(0,1,100)])

    if (content=='both') or (content=='events'):titles.append(str(nev)+' Events')
    if (content=='both') or (content=='stations'):titles.append(str(nsta)+ ' Stations')
    title=' | '.join(titles)
    # plt.ioff();plt.show();plt.pause(0.001)
    if (content=='both') or (content=='events'):
        fig = evm_plottable.plot(resolution='f',
        water_fill_color=water_fill_color,continent_fill_color=continent_fill_color,
        projection=projection,
        color='depth',label=None)
        # plt.ion();plt.show(block=False);plt.pause(0.001)
        clear_output(wait=False)
        # events_plot = fig.get_axes()[0].collections[0]
        # events_plot.set_cmap(cmap)
        # events_plot.set_edgecolor('k');events_plot.set_linewidth(0.05)
        # [p.set_edgecolor('k') for p in fig.get_axes()[0].collections]
        # [p.set_linewidth(0.05) for p in fig.get_axes()[0].collections]
        # sizes = events_plot.get_sizes();events_plot.set_sizes(sizes/event_size_downscale)
        # events_plot.set_alpha(0.001)
        # events_plot.set_linewidths(0.001)
        # events_plot.set_linewidth(0.001)

        # stations_paths=[d for d in fig.get_axes()[0].collections if isinstance(d,plt.matplotlib.collections.PathCollection)]
        # other_paths=[d for d in fig.get_axes()[0].collections if not isinstance(d,plt.matplotlib.collections.PathCollection)]
        # [p.set_edgecolor('k') for p in stations_paths]
        # [p.set_linewidth(0) for p in stations_paths]
        # [p.set_edgecolor('k') for p in stations_paths]
        # [p.set_linewidths(0) for p in stations_paths]
        # [p.set_alpha(1) for p in stations_paths]
        # # -----------------------------------------
        # [p.set_edgecolor('k') for p in other_paths]
        # [p.set_linewidth(latlon_lw) for p in other_paths]
        # [p.set_edgecolor('k') for p in other_paths]
        # [p.set_linewidths(latlon_lw) for p in other_paths]
        # [p.set_alpha(0.5) for p in other_paths]

        # fig.set_size_inches(figsize)
    # else:
    #     fig = None; #plt.figure(figsize=figsize)
    # -- Station plot
    # fig,ax = plt.subplots(figsize=figsize,subplot_kw=dict(projection="geo"))
    # if (content=='both') or (content=='stations'):
    #     if markers is None:markers='v'
    #     fig = inv.plot(size=station_marker_size,label=False,fig=fig,color_per_network=colors,marker='v',
    #     water_fill_color=water_fill_color,continent_fill_color=continent_fill_color,projection=projection)
    #     plt.show(block=False);plt.pause(1)
    #     stations_plot = fig.get_axes()[0].get_children()#[i]
    #     stations_plot.set_edgecolor('k')
    #     stations_paths=[d for d in fig.get_axes()[0].collections if isinstance(d,plt.matplotlib.collections.PathCollection)]
    #     other_paths=[d for d in fig.get_axes()[0].collections if not isinstance(d,plt.matplotlib.collections.PathCollection)]
    #     [p.set_edgecolor('k') for p in stations_paths]
    #     [p.set_linewidth(1) for p in stations_paths]
    #     [p.set_edgecolor('k') for p in stations_paths]
    #     [p.set_linewidths(1) for p in stations_paths]
    #     # [p.set_sizes(0) for p in stations_paths]
    #     # -----------------------------------------
    #     [p.set_edgecolor('k') for p in other_paths]
    #     [p.set_linewidth(latlon_lw) for p in other_paths]
    #     [p.set_edgecolor('k') for p in other_paths]
    #     [p.set_linewidths(latlon_lw) for p in other_paths]
    # if projection=='ortho':bbox_to_anchor=[0.5,-0.1]
    # if projection=='local':bbox_to_anchor=[0.5,-0.1]
    # if projection=='global':bbox_to_anchor=[0.5,-0.1]
    # if (content=='both') or (content=='stations'):
    #     lg = fig.get_axes()[0].get_legend()
    #     labels=lg.texts;handles=lg.legend_handles;lg.remove()
    #     labels = [l._text for l in labels]
    #     lg = fig.get_axes()[0].legend(handles, labels, ncols=8,loc='lower center',
    #     labelspacing=0.0,columnspacing=0.1,
    #     handletextpad=-.5,borderpad=0.2,edgecolor='k',fontsize=13,markerscale=5,bbox_to_anchor=bbox_to_anchor)
    #     lg.set_zorder(10000000)
    # fig.canvas.draw_idle()
    # fig.set_size_inches(figsize)
    # # fig.set_tight_layout('tight')
    # # ============================================================================================
    # # ============================================================================================
    # settings = AttribDict()
    # settings.map_ortho = AttribDict();settings.map_ortho.shrink=0.1
    # settings.map_global=AttribDict();settings.map_global.shrink=0.1
    # settings.map_local=AttribDict();settings.map_local.shrink=0.1
    # if (content=='both') or (content=='events'):
    #     fig.get_axes()[-1].remove()
    #     depths=[d.origins[0].depth for d in evm_plottable]
    #     # norm = fig.get_axes()[1]._colorbar.norm
    #     norm=mpl.colors.Normalize(6.0,8.0)
    #     if mode.lower()=='mag':norm._vmin=6.0;norm._vmax=8.0
    #     # if mode.lower()=='depth':norm._vmin=np.min(depths)/1000;norm._vmax=np.max(depths)/1000
    #     if mode.lower()=='depth':norm._vmin=0;norm._vmax=700
    #     if len(np.unique([e.magnitudes[0].magnitude_type for e in evm_plottable]))>1:mtype='M'
    #     else:mtype=np.unique([e.magnitudes[0].magnitude_type for e in evm_plottable])[0]
    #     if mode.lower()=='mag':colorbarlabel = 'Magnitude ('+mtype+')'
    #     else:colorbarlabel = 'Depth (km)'
    #     if projection=='local':cbaxes = fig.add_axes([0.29, 0.32, 0.41, 0.02])  
    #     if projection=='global':cbaxes = fig.add_axes([0.29, 0.32, 0.41, 0.02])  
    #     if projection=='ortho':cbaxes = fig.add_axes([0.445, 0.032, 0.11, 0.01])  
    #     cbr=fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    #     ax=fig.get_axes()[0], orientation='horizontal',
    #     label=colorbarlabel,
    #     shrink=settings['map_'+projection].shrink,aspect=75,pad=0.01,cax = cbaxes)
    #     cbr.set_alpha(0.7) #75
    # ============================================================================================
    # ============================================================================================
    # fig.suptitle(title,y=0.93)
    # fig.set_tight_layout('tight')
    # fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
    # fig.canvas.draw_idle()
    # fig.set_size_inches(figsize)
    # fig.canvas.draw_idle()
    # ------------
    # if file is not None:
    #     plt.margins(0.1,0.1)
    #     save_tight(file,fig,dpi=700)
    return fig
def unravel(lst):return list(itertools.chain.from_iterable(lst))

# ----------------------------------------------------------------------------------------
bbox=[[-65.444,74.2055],[705.44,106.65]]
evcat=unravel([catalog.loc[stanm].Events for stanm in catalog.StaName])
c,i=np.unique([e.Name for e in evcat],return_index=True)
evcat=Catalog([evcat[ii] for ii in i])
inv=Inventory()
for stanm in catalog.StaName:
    i_inv=catalog.loc[stanm].Inventory.copy()
    for i in i_inv:i.code=f'{catalog.loc[stanm].Deployment.Seismometer}.{catalog.loc[stanm].Deployment.Instrument_Design}'
    # i.code = f'{catalog.loc[stanm].Experiment}.{catalog.loc[stanm].Network}'
    inv+=i_inv
color_per_network = ColorStandard.instrument.copy().__dict__;color_per_network.update(ColorStandard.seismometer_marker.copy())
# seismometers,instrumentdesign = [d.Seismometer for d in catalog.Deployment],[d.Instrument_Design for d in catalog.Deployment]
# [[colors.update({f'{s}.{n}':ColorStandard.instrument[n]}) for n in instrumentdesign] for s in seismometers]

# xxxxxxxxxxxxxxxxxxxxxxxxx
# xxxxxxxxxxxxxxxxxxxxxxxxx
# projection='ortho'
projection='local'
# projection='global'
# ----------------------------------------
content = 'events'
# content = 'both'
# content = 'stations'
# ----------------------------------------
file=dirs.Plots/'_Plots'/'_Maps'/f'Map.{content}.{projection}_test.png'

        # - interactive backends:
        #   GTK3Agg, GTK3Cairo, GTK4Agg, GTK4Cairo, MacOSX, nbAgg, QtAgg,
        #   QtCairo, TkAgg, TkCairo, WebAgg, WX, WXAgg, WXCairo, Qt5Agg, Qt5Cairo

        # - non-interactive backends:
        #   agg, cairo, pdf, pgf, ps, svg, template

# matplotlib.use('TkAgg')

matplotlib.use('MacOSX')

fig=catalog_map_plot(inv,evcat,content=content,file=file,projection=projection,colors=color_per_network)

# ----------------------------------------
plt.close('all')
# fig.set_size_inches(30,10)
# display(fig)

save_tight(file,fig,dpi=800)

fig.savefig(str(file.parent/'png'/file.name).replace('.png','_b.png'),format='png',dpi=500)
fig.savefig(str(file.parent/'svg'/file.name).replace('.png','.svg'),format='svg')
fig.savefig(str(file.parent/'eps'/file.name).replace('.png','.eps'),format='eps')

k=1
# ----------------------------------------
# ----------------------------------------