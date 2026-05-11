#!/usr/bin/env bash

module purge ; module load Autotools pkg-config netCDF-Fortran

for jobfile in namelist/NA5km/midHolocene_IPSL-CM5A-LR.namelist
do

  tmp=climate/NA5km/${jobfile##*/}
  
  outfile=${tmp%%.*}.nc
  
  echo "generating output file:" $outfile
  
  sed -e "s|ANOMTYPE|Spinup|g" NAclimate_template.cdl | ncgen -k 4 -o $outfile 
  
  ./makeclimate $jobfile $outfile

done

exit

for jobfile in namelist/NA5km/1836-1865_mean_20CR.namelist
do

  tmp=climate/NA5km/${jobfile##*/}
  
  outfile=${tmp%%.*}.nc
  
  echo "generating output file:" $outfile
  
  sed -e "s|ANOMTYPE|Spinup|g" NAclimate_template.cdl | ncgen -k 4 -o $outfile 
  
  ./makeclimate $jobfile $outfile

done

exit

for jobfile in `ls namelist/NA5km/midHolocene*.namelist`
do

  tmp=climate/NA5km/${jobfile##*/}
  
  outfile=${tmp%%.*}.nc
  
  echo "generating output file:" $outfile
  
  sed -e "s|ANOMTYPE|Spinup|g" NAclimate_template.cdl | ncgen -k 4 -o $outfile 
  
  ./makeclimate $jobfile $outfile

done
