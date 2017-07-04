
LIB=../lib/libdfftpack.a

F90  = ifort -g -traceback -fpp \
             -r8 \
             -O3 -ipo -unroll-aggressive -fp-model fast=2 \
         	 -m64 -xSSE4.2 \
         	 -heap-arrays -mcmodel=large -shared-intel

#F90 := gfortran -ffree-line-length-none #-O3 -ffast-math -funroll-loops -m64 

MAKE := make
AR := xiar -qoffload-build

#AR := ar

