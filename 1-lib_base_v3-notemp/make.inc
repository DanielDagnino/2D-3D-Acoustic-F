
LIB=../lib/libdag.a

F90  = ifort -warn all -fpp \
             -O3 -ipo -unroll-aggressive -fp-model fast=2 \
             -m64 -xSSE4.2 \
             -heap-arrays -mcmodel=large -shared-intel

# F90 = ifort -g -traceback -check all -fp-stack-check -gen-interfaces -warn interfaces

# F90  = ifort -g -traceback -check bounds -check format -check uninit -check output_conversion -fp-stack-check -gen-interfaces -warn interfaces
       
#       -Ofast -ipo -m64 -xSSE4.2 -unroll-aggressive -fp-model fast=2 -opt-prefetch \
#       -heap-arrays -mcmodel=large -shared-intel

#F90 := gfortran -ffree-line-length-none #-O3 -ffast-math -funroll-loops -m64 

MAKE := make
AR := xiar -qoffload-build

#AR := ar

