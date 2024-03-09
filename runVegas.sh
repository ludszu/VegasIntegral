#!/bin/sh

gfortran -c functions/Bessel.f90
gfortran -c functions/TightBinding.f90
gfortran -c functions/MonteCarlo.f90


gfortran -o prog programVegas.f90 bessel.o tightbinding.o montecarlo.o -I [library path]/linux/bin/intel64/ifort/fgsl -lfgsl
./prog
