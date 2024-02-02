# makeclimate

`makeclimate` is a utility program to generate climate input files for [LPJ-LMfire](https://github.com/ARVE-Research/LPJ-LMfire "The LPJ-LMfire Dynamic Global Vegetation Model").

Required software: `pkg-config Autotools netCDF-Fortran`

After checkout, perform initial setup with

`autoreconf -if`

then

`./configure` (set flags manually as needed)

then

`make`

Job options are specified in the `.namelist` files. When generating a climatology, the user should specify the calendar year in CE (negative for BCE ages) for the first year of the transient run (or the year following the last year of the spinup). In generating a spinup climatology, `makeclimate` will count years down to the final year of the spinup. Climate generated for the spinup will not be a true reflection of the climate over this period, but should be roughly similar to the mean climate state over the first part of the transient run (if any). For the transient period in contrast, climate should reflect the actual climate corresponding to the calendar year (insofar as the reanalysis is a good representation of past climate).

In principle, when generating spinup climates, the climatological mean climate used as the baseline should be pre-adjusted to reflect conditions around the end of the spinup period.