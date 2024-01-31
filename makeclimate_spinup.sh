#!/bin/bash

basefile=/home/terraces/datasets/climate/NorthAm_hybrid/NAclimate.nc

outfile=NAclimate_spinup.nc

numyrs=1020

let tlen=$numyrs*12

cat << EOF > sub.sed
s|time = 12 ;|time = $tlen ;|g
s|_ChunkSizes = 1, 511, 513 ;|_ChunkSizes = 360, 25, 25 ;|g
s|365_day|proleptic_gregorian|g
s|days since 0000-01-01 00:00:00|days since 0001-01-01 00:00:00|g
EOF

jobfile=${outfile%%.*}.namelist

sed -e "s|BASEFILE|$basefile|g" -e "s|NUMYRS|$numyrs|g" template.namelist > $jobfile

ncdump -s -v x,y,lon,lat $basefile | sed -f sub.sed | ncgen -k nc4 -o $outfile

./makeclimate $jobfile $outfile
