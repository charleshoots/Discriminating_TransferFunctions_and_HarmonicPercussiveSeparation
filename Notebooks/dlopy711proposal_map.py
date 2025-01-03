from imports import *
import pygmt
cat=catalog.copy()
cat=catalog[catalog.Network.isin(['X9','7A','7D'])].copy()


doran_laske_calcs=Path('/Users/charlesh/Desktop/711_Proposal/figures/2016165_esupp_Table_S1.csv')
janiszewski23_meta=Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/ATaCR_Data/ATaCR_Python/Catalogs/Janiszewski_etal_2023_StationList.xlsx')
dl_data = pd.read_csv(doran_laske_calcs)
haj_dat= pd.read_excel(janiszewski23_meta)
haj_dat=haj_dat[haj_dat.Station.isin(dl_data.Station)]
dl_data=dl_data[dl_data.Station.isin(haj_dat.Station)]
dl_data=dl_data.sort_values(by='Station')
haj_dat=haj_dat.sort_values(by='Station')
haj_dat=haj_dat.rename(columns={
'Latitude (deg)':'Latitude','Longitude (deg)':'Longitude',
'Pressure Gauge':'Pressure_Gauge','Water Depth (m)':'StaDepth','Instrument Design':'Instrument_Design'})[['Station', 'Network',
'Latitude','Longitude','Experiment',
'Environment','Pressure_Gauge','StaDepth','Start','End','Instrument_Design','Seismometer']]
dl_data=dl_data[['Station','Orientation','Uncertainty','Total Phases','Unique Events','Include High Fqs','CC used']]
dl_data.reset_index(drop=True)
haj_dat.reset_index(drop=True)
dl_data=pd.merge(dl_data,haj_dat)
cat = dl_data.copy()

dl_orients_atacr=load_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/Research/_DataArchive/DLOPY_Data/dl_orients_formap.pkl')
# -------
TRN=PyGMT_PLT_Scatter_Translator
# evcat=unravel([e for e in cat.Events])
# evcat=Catalog([evcat[i] for i in np.unique([e.Name for e in evcat],return_index=True)[1]])
# evcat=pd.DataFrame(evcat)
# for i in range(len(evcat)):
#     evcat.at[i,'Longitude']=evcat.iloc[i].origins[0].longitude;evcat.at[i,'Latitude']=evcat.iloc[i].origins[0].latitude
#     evcat.at[i,'Type']=evcat.magnitudes[i][0].magnitude_type
# evcat[['Name','Longitude','Latitude','Depth', 'Mag','Type', 'Stations']]
# evcat.Depth=evcat.Depth/1000

### evcat=evcat.iloc[[np.where(evcat.Mag==evcat.Mag.min())[0][0],np.where(evcat.Mag==evcat.Mag.max())[0][0]]]

# quake_cmap='wysiwyg'
quake_cmap='romaO'
elevation_cmap="davos"
# elevation_cmap="imola"
# elevation_cmap="lapaz"
#k=1#######---------------------------------------------------------------------------------------------------

# ----------
# ----------
# ----------
# projection="Cyl_stere/150/-20/25i"
projection='G210/25/90/20i'
# dip=90;scale='20i';projection=f'B{int(lons.mean())}/{int(lats.mean())}/{int(lats.min())}/{int(lats.max())}/{scale}'
# projection='B-55/-15/-25/0/12c'
# ----------
# ----------
# ----------
ind=0
resolution='15s'
# resolution='03m'
ortho_region=[-180, 180, -64, 70]
stereo_region=[130, 300, -30, 70]
stereo_region=ortho_region=[-131, -122, 39, 50]
Q=0
for file in ['Stations']:
    for proji,projection in enumerate(['Cyl_stere/150/-20/25i']):
        projname = ['Stereo','Ortho'][proji]
        for quake_cmap in ['phase']:
            for elevation_cmap in ['lapaz']:
                ind+=1
                status = lambda:f'{ind}/2 : {file} : {projection} : {quake_cmap} : {elevation_cmap}'
                print(status())
                # -----BASEMAP
                # region="g"
                # region=[-180, 180, -64, 70]
                pygmt.config(FONT_ANNOT="40p")
                if (projname=='Stereo') and (file=='Stations'):region=stereo_region;gridregion=region
                else:region=ortho_region;gridregion="g"

                grid = pygmt.datasets.load_earth_relief(data_source='synbath',resolution=resolution,region=gridregion)
                grid.data[grid.data>0]=0;grid.data=grid.data/1000;#grid.data=-1*grid.data
                fig = pygmt.Figure()
                
                fig.basemap(region=region, projection=projection, frame=True)
                fig.grdimage(grid=grid, cmap=elevation_cmap)
                fig.coast(land="dimgrey")
                fig.basemap(frame="g01")
                fig.colorbar(frame=["a2f1", "xaf+lElevation (km)"],position="JMR+o1c/0c+w7c/0.5c+n+mc")

                cmap=cm.cmaps['hawaii_r']
                norm = mpl.colors.Normalize(vmin=cat.Orientation.min(), vmax=cat.Orientation.max())

                # -----STATIONS
                if (file=='Stations') or (file=='Both'):
                    if projection=='Cyl_stere/150/-20/25i':xmkrscl=2;markersize=1
                    # else:xmkrscl=1.3;markersize=2
                    data=cat;alpha=1.0;lons,lats=data.Longitude,data.Latitude
                    TRN=PyGMT_PLT_Scatter_Translator
                    fills=[ColorStandard.instrument[s.Instrument_Design] for s in cat.iloc]
                    markers = [ColorStandard.seismometer_marker[s.Seismometer] for s in cat.iloc]
                    styles=[f'{TRN[i]}{markersize}c' for i in markers]
                    styles=[s.replace(f'x{markersize}',f'x{xmkrscl*markersize}') for s in styles]
                    Q=0
                    for si,(sta,lon,lat,style,fill) in enumerate(zip(data.iloc,lons,lats,styles,fills)):

                        if not sta.Station[-1]=='D':
                            if np.any(sta.Station[-1]==np.array(['A','B','C','D'])):
                                nltr=['A','B','C','D'][np.where(sta.Station[-1]==np.array(['A','B','C','D']))[0][0]+1:]
                                if sum([np.sum(data.Station==(sta.Station[:-1]+l)) for l in nltr])>0:continue

                        length=4
                        clrhex=mcolors.to_hex(cmap(norm(sta.Orientation)))

                        # if np.any(np.array(list(dl_orients_atacr.keys()))==(sta.Network+'.'+sta.Station)):
                        #     angle = dl_orients_atacr[sta.Network+'.'+sta.Station].atacr
                        #     fig.plot(x=lon, y=lat, style='v0c', direction=([angle], [length-1]), pen="6p,magenta", fill=clrhex)
                        #     fig.plot(x=lon, y=lat, style='v0c', direction=([angle], [length-1]), pen="1p,black", fill=clrhex)
                        angle=sta.Orientation
                        fig.plot(x=lon, y=lat, style='v0c', direction=([angle], [length]), pen="4p,black", fill=clrhex)
                        fig.plot(x=lon, y=lat, style='v0c', direction=([angle], [length]), pen="1p,black", fill=clrhex)

                        if np.any(np.array(list(dl_orients_atacr.keys()))==(sta.Network+'.'+sta.Station)):
                            angle = dl_orients_atacr[sta.Network+'.'+sta.Station].atacr
                            fig.plot(x=lon, y=lat, style='v0c', direction=([angle], [length-1]), pen="4p,magenta", fill=clrhex)
                            fig.plot(x=lon, y=lat, style='v0c', direction=([angle], [length-1]), pen="1p,magenta", fill=clrhex)

                        if 'x' in style:
                            fig.plot(x=lon,y=lat,
                            style=style.replace(str(xmkrscl*markersize),str((xmkrscl-0.1)*1)),
                            fill='black',pen='10p,black' if 'x' in style else '10p,black')

                        fig.plot(x=lon,y=lat,style=f'x{str((xmkrscl)*1)}' if 'x' in style else style,fill=fill,
                        pen='7p,black' if 'x' in style else '2p,black')
                        Q+=1

                        if sta.Station=='J41C':
                            fig.text(x=lon,y=lat-.13,text=sta.Station,font='35p,Courier-Bold,white',fill='black')


                # ---------------------------------------------------------------
                # ---------------------------------------------------------------
                # ---------------------------------------------------------------
                # -----SAVE
                fig.show()
                folder=Path('/Users/charlesh/Desktop/711_Proposal/figures')
                fig.savefig(folder/f'{file}_{projname}.{quake_cmap}.{elevation_cmap}.png',dpi=700)
                clear_output()
k=1