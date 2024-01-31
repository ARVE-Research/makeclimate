#!/bin/bash

module purge ; module load Autotools pkg-config netCDF-Fortran

basefile=/home/terraces/datasets/climate/NorthAm_hybrid/NAclimate.nc

outfile=NAclimate_transient.nc

jobfile=${outfile%%.*}.namelist

ncgen -k nc4 -o $outfile NAclimate_template.cdl

./makeclimate $jobfile $outfile
