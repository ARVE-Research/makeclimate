# makeclimate

`makeclimate` is a utility program to generate climate input files for [LPJ-LMfire](https://github.com/ARVE-Research/LPJ-LMfire "The LPJ-LMfire Dynamic Global Vegetation Model").

Required software: `pkg-config Autotools netCDF-Fortran`

After checkout, perform initial setup with

`autoreconf -if`

then

`./configure` (set flags manually as needed)

then

`make`
