ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = makeclimate

makeclimate_SOURCES =  \
  parametersmod.f90    \
  netcdfmod.f90        \
  calcvarmod.f90       \
  calcwetmod.f90       \
  randomdistmod.f90    \
  calendarmod.f90      \
  makeclimate.f90      

clean-local:
	rm -r -f *.mod *.dSYM
