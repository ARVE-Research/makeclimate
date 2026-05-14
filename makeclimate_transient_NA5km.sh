#!/bin/bash

module load Autotools pkg-config netCDF netCDF-Fortran

jobfile=namelist/NA5km/1980-2024_MERRA-2.namelist

tmp=${jobfile##*/}

outfile=${tmp%%.*}.nc

echo "generating output file:" $outfile

sed -e "s|ANOMTYPE|Transient|g" NAclimate_template.cdl | ncgen -k nc4 -o $outfile 

./makeclimate $jobfile $outfile
