import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent))
from imports import *
import pygmt
cat=catalog.copy()
# cat=catalog[catalog.Network.isin(['X9','7A','7D'])].copy()

# -------
TRN=PyGMT_PLT_Scatter_Translator
evcat=unravel([e for e in cat.Events])
evcat=Catalog([evcat[i] for i in np.unique([e.Name for e in evcat],return_index=True)[1]])
evcat=pd.DataFrame(evcat)

# evcat=lt.cat.unravel_cat(cat)

for i in range(len(evcat)):
    evcat.at[i,'Name'] = evcat.iloc[i].Name
    evcat.at[i,'Longitude']=evcat.iloc[i].origins[0].longitude;evcat.at[i,'Latitude']=evcat.iloc[i].origins[0].latitude
    evcat.at[i,'Mag'] = evcat.magnitudes[i][0].mag
    evcat.at[i,'Depth'] = evcat.iloc[i].origins[0].depth/1000
    evcat.at[i,'Type'] = evcat.magnitudes[i][0].magnitude_type
    evcat.at[i,'Stations'] = evcat.iloc[i].Stations
evcat=evcat[['Name','Longitude','Latitude','Depth', 'Mag','Type', 'Stations']]

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
resolution='03m'
ortho_region=[-180, 180, -64, 70];stereo_region=[130, 300, -30, 70]
# stereo_region=ortho_region=[-131, -122, 39, 50];resolution='15s'
for file in ['Both']:
    for proji,projection in enumerate(['Cyl_stere/150/-20/25i','G210/25/90/20i']):
        projname = ['Stereo','Ortho'][proji]
        # for quake_cmap in ['romaO','wysiwyg','phase']:
        #     for elevation_cmap in ['davos','imola','lapaz']:
        for quake_cmap,elevation_cmap in zip(['phase','romaO','wysiwyg'],['lapaz','lapaz','lapaz']):
            ind+=1
            # status = lambda:f'{ind} {file} : {projection} : {quake_cmap} : {elevation_cmap}'
            # print(status())
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
            fig.basemap(frame="g15")
            fig.colorbar(frame=["a2f1", "xaf+lElevation (km)"],position="JMR+o1c/0c+w7c/0.5c+n+mc")
            # ---------------------------------------------------------------
            # ---------------------------------------------------------------
            # ---------------------------------------------------------------
            # -----EVENTS
            if (file=='Events') or (file=='Both'):
                if file=='Both':alpha=0.30
                else:alpha=1.0
                data=evcat
                sizes = 0.01*2**data.Mag;fills=data.Depth;depths=data.Depth;markers='c'
                lons,lats=data.Longitude,data.Latitude
                pygmt.makecpt(cmap=quake_cmap, series=[data.Depth.min(), data.Depth.max()])
                fig.plot(x=lons,y=lats,size=sizes,fill=data.Depth,
                cmap=True,style='cc',pen="1p,black",transparency=100-(alpha*100))
                fig.colorbar(frame="xaf+lDepth (km)",position="JMR+o1c/10c+w7c/0.5c+n+mc")
            # ---------------------------------------------------------------
            # ---------------------------------------------------------------
            # ---------------------------------------------------------------
            # -----STATIONS
            if (file=='Stations') or (file=='Both'):
                if projection=='Cyl_stere/150/-20/25i':xmkrscl=.5;markersize=.5
                else:xmkrscl=2;markersize=2
                data=cat;alpha=1.0;lons,lats=data.Longitude,data.Latitude
                TRN=PyGMT_PLT_Scatter_Translator
                fills=[ColorStandard.instrument[s.Deployment.Instrument_Design] for s in cat.iloc]
                markers = [ColorStandard.seismometer_marker[s.Deployment.Seismometer] for s in cat.iloc]
                styles=[f'{TRN[i]}{markersize}c' for i in markers]
                styles=[s.replace(f'x{markersize}',f'x{xmkrscl*markersize}') for s in styles]
                for lon,lat,style,fill in zip(lons,lats,styles,fills):
                    if 'x' in style:
                        fig.plot(x=lon,y=lat,
                        style=style.replace(str(xmkrscl*markersize),str((xmkrscl-0.2)*markersize)),
                        fill='black',pen='25p,black' if 'x' in style else '0.5p,black')
                    fig.plot(x=lon,y=lat,
                    style=f'x{str((xmkrscl-0.5)*markersize)}' if 'x' in style else style,fill=fill,
                    pen='25p,black' if 'x' in style else '0.5p,black')
            # ---------------------------------------------------------------
            # ---------------------------------------------------------------
            # ---------------------------------------------------------------
            # -----SAVE
            fig.show()
            folder=dirs.Ch1/'Maps'
            fig.savefig(folder/f'{file}_{projname}.{quake_cmap}.{elevation_cmap}.png',dpi=700)
            clear_output()
