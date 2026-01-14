#!/usr/bin/env bash

module purge ; module load Autotools pkg-config netCDF-Fortran

jobfile=${1}

outfile=${jobfile%%.*}.nc

echo "generating output file:" $outfile

sed -e "s|ANOMTYPE|Spinup|g" global_30m_template.cdl | ncgen -k 4 -o $outfile 

./makeclimate $jobfile $outfile
