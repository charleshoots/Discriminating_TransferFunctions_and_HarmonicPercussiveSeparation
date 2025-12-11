### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *
from source.imports import *
from source.modules import *
import pygmt
# -*- coding: utf-8 -*-
import numpy as np, pandas as pd, re, pygmt
from pathlib import Path

# Options
MapsToMake = ['global','network']
# mapfolder value
mapfolder = dirs.Plots/'_Papers'/'ImageOutputs'/'_main_figures'/'Figure1_MapPlots' #Folder maps are saved to.

# Data and variables
# --- inputs from your codebase ---


# save_format='png'
save_format='pdf' #plot save save_format

cat = catalog.copy()
# TRN value
TRN = PyGMT_PLT_Scatter_Translator
# plot subarray center value
plot_subarray_center = True
# plot stations value
plot_stations = False
# plot_stations = True
# -------- event table (one row per event) --------
icat = cat.sr.copy()
# evcat value
evcat = pd.DataFrame([icat.aloc[evi].iloc[0] for evi in np.unique(icat.Name)])
# evcat value
evcat = evcat[["Name","Magnitude","Stations","EvDepth","LaLo"]].copy()
evcat["Longitude"] = np.array(list(evcat.LaLo))[:,1]
evcat["Latitude"]  = np.array(list(evcat.LaLo))[:,0]
evcat["Depth"]     = np.array(list(evcat.EvDepth)).astype(float)
# evcat value
evcat = evcat[["Name","Longitude","Latitude","Magnitude","Stations","Depth"]]

# -------- helpers --------
def wrap_lon_360(lon):lon=np.asarray(lon, float); return np.mod(lon, 360.0)
# function map width cm
def map_width_cm(proj, default=12.0):
    # m value
    m = re.search(r"/([-+]?\d+(?:\.\d+)?)([cip])$", proj)
    if not m: return default
    w,u = float(m.group(1)), m.group(2)
    return w if u=="c" else (w*2.54 if u=="i" else w*0.0352778)
# function size pt from projection
def size_pt_from_projection(proj, base_pt_at_12cm=8.0):
    return base_pt_at_12cm * (12.0 / map_width_cm(proj))
# function mag to points
def mag_to_points(mag, pmin=3.0, pmax=12.0, mmin=5.0, mmax=8.0):
    # mag value
    mag = np.asarray(mag, float)
    # t value
    t = np.clip((mag - mmin) / (mmax - mmin), 0, 1)
    return pmin + t*(pmax - pmin)
# function wrap lon 360
def wrap_lon_360(lon):lon = np.asarray(lon, float);return np.mod(lon, 360.0)
# function wrap lon 180
def wrap_lon_180(lon):lon = np.asarray(lon, float);return (lon + 180.0) % 360.0 - 180.0
# function lon span
def lon_span(lons):lmin, lmax = np.nanmin(lons), np.nanmax(lons);return (lmax - lmin, lmin, lmax)
# function choose region
def choose_region(lons_deg, lats_deg, pad_deg=1.2, pad_frac=0.07, min_pad=0.35):
    """
    Build a tight rectangular region around points with less extra map area.
    - pad_deg: absolute lat/lon padding in degrees (was 2.0)
    - pad_frac: fraction of the lon span added as extra pad (was 0.10)
    - min_pad: minimum extra pad (keeps room for markers/outline)
    Returns: (region [w,e,s,n], wrap_mode '360' or '180')
    """
    # lons value
    lons = np.asarray(lons_deg, float)
    # lats value
    lats = np.asarray(lats_deg, float)
    # Compare wraps and pick the one with smaller span
    lons180 = wrap_lon_180(lons)
    span180, min180, max180 = lon_span(lons180)

    # lons360 value
    lons360 = wrap_lon_360(lons)
    span360, min360, max360 = lon_span(lons360)

    if span360 < span180:
        # wrap mode value
        wrap_mode = "360"; lon_min, lon_max = min360, max360; span = span360
    else:
        # wrap mode value
        wrap_mode = "180"; lon_min, lon_max = min180, max180; span = span180
    # Tighter longitudinal padding: smaller fraction + floor
    lon_pad = max(min_pad, pad_deg, pad_frac * max(span, 1e-6))
    # Tighter latitudinal padding: use pad_deg; also ensure a min margin
    lat_min = np.nanmin(lats) - max(min_pad, pad_deg)
    # lat max value
    lat_max = np.nanmax(lats) + max(min_pad, pad_deg)
    # Clamp latitude for Mercator
    lat_min = max(lat_min, -80.0)
    # lat max value
    lat_max = min(lat_max,  80.0)
    # region value
    region = [lon_min - lon_pad, lon_max + lon_pad, lat_min, lat_max]
    return region, wrap_mode
# function subset events to region
def subset_events_to_region(evdf, region, wrap_mode):
    w,e,s,n = region
    # lon value
    lon = wrap_lon_360(evdf["Longitude"].values) if wrap_mode=="360" else wrap_lon_180(evdf["Longitude"].values)
    # lat value
    lat = evdf["Latitude"].values
    # msk value
    msk = (lon>=w)&(lon<=e)&(lat>=s)&(lat<=n)
    # out value
    out = evdf.loc[msk].copy()
    out["__lon_wrapped__"] = lon[msk]
    return out
# function wrap series for plot
def wrap_series_for_plot(lon_series, wrap_mode):
    # arr value
    arr = lon_series.values if hasattr(lon_series, "values") else np.asarray(lon_series)
    return wrap_lon_360(arr) if wrap_mode=="360" else wrap_lon_180(arr)
# function station keyframe
def station_keyframe(df):
    """Return a Series key for de-dup: prefer Station/Name/ID; fallback to rounded lon/lat."""
    # cols value
    cols = df.columns
    if "Station" in cols: return df["Station"].astype(str)
    if "Name" in cols:    return df["Name"].astype(str)
    if "ID" in cols:      return df["ID"].astype(str)
    # fallback: lon/lat rounding to ~100 m (3 decimals ~ 100 m)
    lon = np.asarray(df["Longitude"], float)
    # lat value
    lat = np.asarray(df["Latitude"],  float)
    return pd.Series(np.char.add(np.round(lon,3).astype(str), np.round(lat,3).astype(str)), index=df.index)

if 'global' in MapsToMake:
    print('Starting global map..')
    # -------- style / projection --------
    projection = "Ks-160/12c"     # Eckert VI, 160°W center
    # region value
    region = "g"                   # global 0–360
    # bathy cmap value
    bathy_cmap = "lapaz"
    # depth cmap value
    depth_cmap = "romaO"
    # depth cmap value
    depth_cmap = "seis"

    # before any plotting (right after creating `fig` is fine):
    pygmt.config(
    # MAP LINE STEP value
    MAP_LINE_STEP="0.20p",    # much finer resampling of curves
    # PS LINE JOIN value
    PS_LINE_JOIN="round",     # rounded joins
    # PS LINE CAP value
    PS_LINE_CAP="round",      # rounded caps
    COLOR_NAN="dimgrey",
    FONT_ANNOT="7p"
    )
    fig = pygmt.Figure()

    # -------- basemap + bathymetry --------
    grid = pygmt.datasets.load_earth_relief(data_source="synbath", resolution="03m", region=region)
    grid = grid.where(grid <= 0)     # keeps ocean (<=0), makes land NaN
    grid.data = grid.data / 1000.0  # km



    fig.basemap(region=region, projection=projection, frame=True)
    # fig.coast(land="dimgrey", water="dimgrey", shorelines=False,resolution="c", area_thresh=0)
    # fig.plot(x=[0, 360, 360, 0], y=[-90, -90, 90, 90],fill="dimgrey") #, pen="0.1p,none", close=True)


    fig.grdimage(grid=grid, cmap=bathy_cmap)

    fig.coast(resolution="c",
    # shorelines="1/0.25p,black",# thin, crisp shoreline
    land="dimgrey",
    area_thresh=2000)# keep tiny islands; raise to de-speckle

    fig.coast(resolution="h",
    shorelines="1/0.25p,black",# thin, crisp shoreline
    land="dimgrey",
    area_thresh=1000)# keep tiny islands; raise to de-speckle


    # ===================== EARTHQUAKES =====================
    elons = wrap_lon_360(evcat.Longitude.values)
    elats = evcat.Latitude.values
    edeps = evcat.Depth.values
    emags = evcat.Magnitude.values

    pygmt.makecpt(cmap=depth_cmap, series=[np.nanmin(edeps), np.nanmax(edeps)], continuous=True)

    # Size in cm; 0.4 factor makes them ~60% smaller
    pt_to_cm = 0.0352778
    eq_sizes_cm = mag_to_points(emags) * pt_to_cm * 1.2

    fig.plot(
    x=elons, y=elats,
    style="cc", size=eq_sizes_cm,
    fill=edeps, cmap=True,pen="0.2p,black",transparency=40)
    fig.colorbar(frame="xaf+lDepth (km)", position="JMR+o0.6c/0c+w7c/0.35c+n+mc")


    if plot_stations:
        ## ===================== STATIONS (refined X styling) =====================
        slons = wrap_lon_360(np.array(list(cat.r.Longitude)))
        slats = np.array(list(cat.r.Latitude))
        fills = [ColorStandard.instrument[s.Deployment.Instrument_Design] for s in cat.r.iloc]
        markers = [ColorStandard.seismometer_marker[s.Deployment.Seismometer] for s in cat.r.iloc]

        # Base station size (points), already 60% smaller upstream via *0.4
        st_pt_base = size_pt_from_projection(projection, base_pt_at_12cm=8.0)
        st_pt = st_pt_base * 0.4

        for lon, lat, mkr, fill in zip(slons, slats, markers, fills):
            code = TRN[mkr]
            if code == "x":
                # Very thin strokes + small size offset to avoid "square" look
                s_red_pt   = st_pt * 0.96   # slightly smaller red on top
                s_blck_pt  = st_pt * 1.06   # tiny larger black "halo" underneath
                fig.plot(x=lon, y=lat, style=f"x{ s_blck_pt:.3g}p", pen="0.45p,black")
                fig.plot(x=lon, y=lat, style=f"x{  s_red_pt:.3g}p", pen="0.28p,red")
            else:
                fig.plot(x=lon, y=lat, style=f"{code}{st_pt:.3g}p", fill=fill, pen="0.45p,black")


    if plot_subarray_center:
        font="1p,Helvetica-Bold,black" #keep it small so i can put a tiny marker over it in inDesign
        fill = 'black' #color
        code = 's' #shape used for center marker
        subarray_center_latlon = {n:[cat.r.aloc[n].Latitude.mean(),cat.r.aloc[n].Longitude.mean()] for n in cat.r.Network.unique()}
        for n in subarray_center_latlon.keys():

            fig.text(x=subarray_center_latlon[n][1], #lon
                y=subarray_center_latlon[n][0], #lat
                text=n,font=font,)


    # -------- save/show --------
    outdir=mapfolder
    outdir.mkdir(parents=True,exist_ok=True)

    outfile = outdir / f"'S03.01_GlobalMap{'.with.Stations' if plot_stations else '.without.Stations'}.{save_format}"
    if save_format=='png':fig.savefig(str(outfile), dpi=900)
    else:fig.savefig(str(outfile))
    fig.show()
    print(f"Global map done - Saved: {outfile}")
    plt.close();del fig


if 'network' in MapsToMake:
    print('Starting inset maps for each network..')
    cat  = catalog.copy()
    TRN  = PyGMT_PLT_Scatter_Translator
    jdf  = catalog.Janiszewski23.copy()  # extra station catalog
    # ---- events table (one row per event) ----
    icat  = cat.sr.copy()
    evcat = pd.DataFrame([icat.aloc[evi].iloc[0] for evi in np.unique(icat.Name)])
    evcat = evcat[["Name","Magnitude","Stations","EvDepth","LaLo"]].copy()
    evcat["Longitude"] = np.array(list(evcat.LaLo))[:,0]
    evcat["Latitude"]  = np.array(list(evcat.LaLo))[:,1]
    evcat["Depth"]     = np.array(list(evcat.EvDepth)).astype(float)
    evcat = evcat[["Name","Longitude","Latitude","Magnitude","Stations","Depth"]]
    print(f"Starting inset maps for each network..")
    # --------------- styles / constants ----------------
    bathy_cmap  = "lapaz"
    depth_cmap  = "romaO"
    grid_res    = "01m"
    projection  = "M12c"     # rectangular Mercator
    pt_to_cm    = 0.0352778

    # crisper lines
    pygmt.config(
        MAP_LINE_STEP="0.20p",
        PS_LINE_JOIN="round",
        PS_LINE_CAP="round",
        FONT_ANNOT="7p",
    )


    networks = np.unique(cat.r["Network"])
    for net in networks:
        sub = cat.r[cat.r["Network"]==net].copy()
        if sub.empty: continue

        # region/wrap per-network
        slons_raw = np.array(list(sub["Longitude"]), float)
        slats     = np.array(list(sub["Latitude"]),  float)
        region, wrap_mode = choose_region(slons_raw, slats, pad_deg=.5)

        # fig & background
        fig = pygmt.Figure()
        grid = pygmt.datasets.load_earth_relief(data_source="synbath", resolution=grid_res, region=region)
        grid.data[grid.data>0] = 0.0; grid.data = grid.data/1000.0
        fig.basemap(region=region, projection=projection, frame=True)
        fig.grdimage(grid=grid, cmap=bathy_cmap)
        fig.coast(land="dimgrey", shorelines="1/0.25p,black", resolution="full", area_thresh=500)

        # events inside region
        ev_sub = subset_events_to_region(evcat, region, wrap_mode)
        if not ev_sub.empty:
            elons = ev_sub["__lon_wrapped__"].values
            elats = ev_sub["Latitude"].values
            edeps = ev_sub["Depth"].values
            emags = ev_sub["Magnitude"].values
            pygmt.makecpt(cmap=depth_cmap, series=[np.nanmin(edeps), np.nanmax(edeps)], continuous=True)
            eq_sizes_cm = mag_to_points(emags) * pt_to_cm
            fig.plot(x=elons, y=elats, style="cc", size=eq_sizes_cm, fill=edeps, cmap=True, pen="0.2p,black", transparency=50)
            fig.colorbar(frame="xaf+lDepth (km)", position="JMB+o0.2c/0.3c+w5c/0.25c+n+mc")

        st_pt = size_pt_from_projection(projection, base_pt_at_12cm=7.0) * 0.8 * 3.5  # 3× larger
        hjan_st_pt=st_pt*.45
        hjan_transparency = 45

        # -------- supplemental stations (Janiszewski23) --------
        jnet = jdf[jdf["Network"]==net].copy()
        if not jnet.empty:
            # remove any stations that are already in your list (by Station/Name/ID if present; else by ~100 m coord key)
            key_primary = station_keyframe(sub)
            key_supp    = station_keyframe(jnet)
            keep_mask   = ~key_supp.isin(set(key_primary.astype(str)))
            jsub        = jnet.loc[keep_mask].copy()

            if not jsub.empty:
                j_lons = wrap_series_for_plot(jsub["Longitude"], wrap_mode)
                j_lats = np.array(list(jsub["Latitude"]), float)
                fills_j   = [ColorStandard.instrument[s.Instrument_Design] for s in jsub.iloc]
                markers_j = [ColorStandard.seismometer_marker[s.Seismometer] for s in jsub.iloc]

                # same size as primary (still 3×), but half the alpha → more transparent
                for lon, lat, mkr, fill in zip(j_lons, j_lats, markers_j, fills_j):
                    code = TRN[mkr]
                    if code == "x":
                        s_red_pt = hjan_st_pt * 1.46
                        s_blk_pt = hjan_st_pt * 1.64
                        fig.plot(x=lon, y=lat, style=f"x{s_blk_pt:.3g}p", pen="3.0p,black", transparency=hjan_transparency)
                        fig.plot(x=lon, y=lat, style=f"x{s_red_pt:.3g}p", pen="1.5p,red",   transparency=hjan_transparency)
                    else:
                        fig.plot(x=lon, y=lat, style=f"{code}{hjan_st_pt:.3g}p", fill=fill, pen="0.45p,black", transparency=hjan_transparency)

        # -------- primary stations (your list) --------
        edgec='black'
        edgec='white'
        fills_primary   = [ColorStandard.instrument[s.Deployment.Instrument_Design] for s in sub.iloc]
        markers_primary = [ColorStandard.seismometer_marker[s.Deployment.Seismometer] for s in sub.iloc]
        slons = wrap_series_for_plot(sub["Longitude"], wrap_mode)
        for lon, lat, mkr, fill in zip(slons, slats, markers_primary, fills_primary):
            code = TRN[mkr]
            # --- PRIMARY stations (X style) ---
            if code == "x":
                # keep the micro size offset to prevent the arms merging visually
                s_red_pt = st_pt * 1.13
                s_blk_pt = st_pt * 1.15
                # thicker strokes than before
                fig.plot(x=lon, y=lat, style=f"x{s_blk_pt:.3g}p", pen=f"3.30p,{edgec}")
                fig.plot(x=lon, y=lat, style=f"x{s_red_pt:.3g}p", pen="2.7p,red")
            else:
                fig.plot(x=lon, y=lat, style=f"{code}{st_pt:.3g}p", fill=fill, pen=f"1.0p,{edgec}")



        # title/save
        fig.basemap(frame=f"+t{net}")
        outdir=mapfolder
        outdir.mkdir(parents=True, exist_ok=True)
        outfile = outdir / f"S03.01_network_{net}.{save_format}"


        if save_format=='png':fig.savefig(str(outfile), dpi=900)
        else:fig.savefig(str(outfile))


        # pdffile = outdir / f"network_{net}.pdf"
        # fig.savefig(str(pdffile))
        fig.show()
        del fig
    print(f"Inset maps done - Saved to: {outdir}")

print("Done.")
