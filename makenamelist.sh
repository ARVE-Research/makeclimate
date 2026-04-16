#!/usr/bin/env bash

while read model
do

  echo $model
  
  sed -e "s/MODEL/$model/g" namelist/template.namelist > namelist/midHolocene_$model.namelist

done < ../makepaleoclimate/modellist.txt
