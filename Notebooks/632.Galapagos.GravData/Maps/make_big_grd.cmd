#! /bin/bash
mbgrid -I datalist_SeisArray.mb-1 \
    -A2 -F5 -N -C4 -E200/200 \
    -O SeismoArray_bathymetry_200m -V
mbgrdviz -I SeismoArray_bathymetry_200m.grd &
