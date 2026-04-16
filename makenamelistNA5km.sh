#!/usr/bin/env bash

while read model
do

  echo $model
  
  sed -e "s/MODEL/$model/g" namelist/NA5km/template.namelist > namelist/NA5km/midHolocene_$model.namelist

done < ../makepaleoclimate/modellist.txt
