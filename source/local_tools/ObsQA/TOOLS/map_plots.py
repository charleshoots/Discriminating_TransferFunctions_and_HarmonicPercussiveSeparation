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

from modules import *
from imports import *
def map(color_by='Instrument',figsize=None,mapkey=0,legend_title='Legend'):
        # Data--------------------------------------------------------
        if color_by=='Station':
            key='StaName';bools=np.array([d[key].split('.')[0] for d in catalog.Deployment]);filts=np.unique(bools)
            hexclrs=[ColorStandard.network[k] for k in filts]
        if color_by=='Network':
            key='StaName';bools=np.array([d[key].split('.')[0] for d in catalog.Deployment]);filts=np.unique(bools)
            hexclrs=[ColorStandard.network[k] for k in filts]
        if color_by=='Instrument':
            key='Instrument_Design';bools=np.array([d[key].split('.')[0] for d in catalog.Deployment]);filts=np.unique(bools)
            hexclrs=[ColorStandard.instrument[k] for k in filts]
        if color_by=='Seismometer':
            key='Seismometer';bools=np.array([d[key].split('.')[0] for d in catalog.Deployment]);filts=np.unique(bools)
            markers = [ColorStandard.seismometer_marker[k] for k in filts]
        if color_by=='Pressure':key='Pressure_Gauge';bools=np.array([d[key].split('.')[0] for d in catalog.Deployment]);filts=np.unique(bools)
        sta_set=[];[sta_set.append(catalog[bools==b].copy()) for b in filts]
        lobnds = [catalog.Longitude.min(),catalog.Longitude.max()];labnds = [catalog.Latitude.min(),catalog.Latitude.max()]
        sta_set_orig = sta_set.copy()
        # Tiles--------------------------------------------------------
        tilekey = dict() # good map tiles: http://leaflet-extras.github.io/leaflet-providers/preview/
        tilekey[0] = ('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}','Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic')
        tilekey[1] = ('https://server.arcgisonline.com/ArcGIS/rest/services/World_Terrain_Base/MapServer/tile/{z}/{y}/{x}','Tiles &copy; Esri &mdash; Source: USGS, Esri, TANA, DeLorme, and NPS')
        tilekey[2] = ('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}','Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community')
        if hexclrs is None:
            colors = distinctipy.get_colors(len(sta_set),pastel_factor=0.0,exclude_colors=[(1.0, 0.0, 1.0),(1,1,1),(0,0,0),(0,1,0),(0,1,1)],colorblind_type='Deuteranomaly',n_attempts=1000)
            hexclrs = [matplotlib.colors.to_hex(colors[i]) for i in range(len(colors))][0:len(colors)]
            hexclrs = hexclrs[0:len(sta_set)]
        # Base Map--------------------------------------------------------
        if figsize is None:figsize = [2000,800]
        fig = Figure(width=figsize[0], height=figsize[1])
        m = folium.Map(tiles=tilekey[mapkey][0],attr=tilekey[mapkey][1],zoom_start=1,min_zoom=1,minlon = lobnds[0],maxlon=lobnds[1],control=False,max_bounds=True,zoom_control=False,scrollWheelZoom=True) #,width=figsize[0],height=figsize[1]
        m.add_child(folium.LatLngPopup())
        # Build Map--------------------------------------------------------
        if sta_set is not None:
                legend_dict = dict()
                for clri,(label,color,cur_sta_set) in enumerate(zip(filts,hexclrs,sta_set)):
                        legend_dict[label]  = color
                        locationlist = (np.array([cur_sta_set['Latitude'],cur_sta_set['Longitude']]).T.tolist())
                        g1 = folium.FeatureGroup(name=('<span style=>{txt}</span>').format( txt='All Stations'))
                        popups=list(cur_sta_set['Network'] + '.' + list(cur_sta_set['Station']) + ' (' + list(cur_sta_set['Experiment']) + ')')
                        for point in range(0, len(locationlist)):
                                marker = folium.Marker(locationlist[point],
                                popup=popups[point],
                                icon=folium.Icon(color='white',
                                icon_color=color,prefix='fa',
                                icon='caret-up',
                                shadow_size=0,
                                icon_size=(0,0),
                                icon_anchor=(0,0)),
                                legend_name='hello')
                                L = folium.FeatureGroup().add_child(marker)
                                marker.add_to(g1)
                        m.add_child(g1)
                if len(sta_set)==1:
                        labels = [color_by]
                else:
                        labels = list(legend_dict.keys())
                colors = list(legend_dict.values())
                macro = getLegendTemplate(colors,labels,title_str=legend_title,legx='"auto"',legy='"auto"')
                m.add_child(macro)
                m.fit_bounds(lobnds)
        fig.add_child(m)
        return m




def mapby(df,color_by=None,height=900,width=1000,map_set='NOAA/NGDC/ETOPO1',elevation=['black','white', 'black','white'],layercontrol=False,center=None,file=None):
        if color_by is not None:
                items = df[color_by].unique()
                sta_set = [df.iloc[np.in1d(df[color_by],item)] for item in items]
        else:
                sta_set = [df]
        MapObj = geemap.Map(height=height,width=width)
        landcover = ee.Image(map_set).select('bedrock')
        elevationVis = {'min': -7000.0,'max': 3000.0,'palette': elevation}
        MapObj.addLayer(landcover, elevationVis, 'NOAA')
        if layercontrol:
                # MapObj.addLayerControl()
                pass
        if len(sta_set)==1:
                hexclrs = ['#e41a1c'] #red
        else:
                colors = distinctipy.get_colors(len(sta_set),pastel_factor=0.0,exclude_colors=[(1.0, 0.0, 1.0),(1,1,1),(0,0,0),(0,1,0),(0,1,1)],colorblind_type='Deuteranomaly',n_attempts=1000)
                hexclrs = [matplotlib.colors.to_hex(colors[i]) for i in range(len(colors))]
        legend_dict = dict()
        for clri,cur_sta_set in enumerate(sta_set):
                color = hexclrs[clri]
                if color_by is not None:
                        label = cur_sta_set[color_by].tolist()[0]
                        legend_dict[label]  = color
                        MapObj.add_circle_markers_from_xy(cur_sta_set, x="Longitude", y="Latitude", radius=1, color=color, fill_color="black",label=label)
                else:
                        MapObj.add_circle_markers_from_xy(cur_sta_set, x="Longitude", y="Latitude", radius=1, color=color, fill_color="black")
                        # MapObj.add_legend(legend_title="Test", legend_dict=legend_dict,position="bottomright")
        if color_by is not None:
                MapObj.add_legend(title=color_by, legend_dict=legend_dict,position="bottomright")
        if center is not None:
                MapObj.set_center(lat=center[0],lon=center[1],zoom=center[2])
        else:
                MapObj.zoom_to_bounds((df.Longitude.min(),df.Latitude.min(),df.Longitude.max(),df.Latitude.max()))
        MapObj.set_control_visibility(layerControl=False,fullscreenControl=False,latLngPopup=False)
        if file is not None:
                MapObj
                MapObj.to_image(file)
        return MapObj






# -----MACRO-----
def getLegendTemplate(colors,labels,title_str='Legend',opacity = 0.7,legx='50px',legy='50px'):
    template = """
    {% macro html(this, kwargs) %}

    <!doctype html>
    <html lang="en">
    <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>jQuery UI Draggable - Default functionality</title>
    <link rel="stylesheet" href="//code.jquery.com/ui/1.12.1/themes/base/jquery-ui.css">

    <script src="https://code.jquery.com/jquery-1.12.4.js"></script>
    <script src="https://code.jquery.com/ui/1.12.1/jquery-ui.js"></script>
    
    <script>
    $( function() {
        $( "#maplegend" ).draggable({
                        start: function (event, ui) {
                            $(this).css({
                                right: "auto",
                                top: "auto",
                                bottom: "auto"
                            });
                        }
                    });
    });

    </script>
    </head>
    <body>

    
    <div id='maplegend' class='maplegend' 
        style='position: absolute; bottom: ##LEGY##; left: ##LEGX##; z-index:9999; border:2px solid grey; background-color:rgba(255, 255, 255, 0.8);
        border-radius:6px; padding: 10px; font-size:14px;'>
        
    <div class='legend-title'>##TITLE##</div>
    <div class='legend-scale'>
    <ul class='legend-labels'>
        ##ITEM##
        ##BOTTOM##

    </ul>
    </div>
    </div>
    
    </body>
    </html>

    <style type='text/css'>
    .maplegend .legend-title {
        text-align: left;
        margin-bottom: 5px;
        font-weight: bold;
        font-size: 90%;
        }
    .maplegend .legend-scale ul {
        margin: 0;
        margin-bottom: 5px;
        padding: 0;
        float: left;
        list-style: none;
        }
    .maplegend .legend-scale ul li {
        font-size: 80%;
        list-style: none;
        margin-left: 0;
        line-height: 18px;
        margin-bottom: 2px;
        }
    .maplegend ul.legend-labels li span {
        display: block;
        float: left;
        height: 16px;
        width: 30px;
        margin-right: 5px;
        margin-left: 0;
        border: 1px solid #999;
        }
    .maplegend .legend-source {
        font-size: 80%;
        color: #777;
        clear: both;
        }
    .maplegend a {
        color: #777;
        }
    </style>
    {% endmacro %}"""

    template = template.replace('##TITLE##',title_str)
    template = template.replace('##LEGX##',legx)
    template = template.replace('##LEGY##',legy)
    for c,l in zip(colors,labels):
        item = "<li><span style='background:" + str(c) + ";opacity:" + str(opacity) + ";'></span>" + str(l) + "</li>"
        template = template.replace('##ITEM##',item)
        template = template.replace('##BOTTOM##',"##ITEM##\n    ##BOTTOM##")
    template = template.replace('    ##ITEM##','')
    template = template.replace('    ##BOTTOM##','')
    
    macro = MacroElement()
    macro._template = Template(template)
    return macro




# ======



m=map()
k=1

