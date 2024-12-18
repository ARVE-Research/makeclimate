#!/bin/bash

module purge ; module load Autotools pkg-config netCDF-Fortran

outfile=ModelE-LM_30m_spinup.nc

jobfile=${outfile%%.*}.namelist

echo "generating output file:" $outfile

sed -e "s|ANOMTYPE|Transient|g" global_30m_template.cdl | ncgen -k 4 -o $outfile 

./makeclimate $jobfile $outfile
