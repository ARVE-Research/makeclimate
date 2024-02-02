#!/bin/bash

module purge ; module load Autotools pkg-config netCDF-Fortran

basefile=/home/terraces/datasets/climate/NorthAm_hybrid/NAclimate.nc

outfile=NAclimate_transient.nc

jobfile=${outfile%%.*}.namelist

echo "generating output file:" $outfile

sed -e "s|ANOMTYPE|Transient|g" NAclimate_template.cdl | ncgen -k nc4 -o $outfile 

./makeclimate $jobfile $outfile
