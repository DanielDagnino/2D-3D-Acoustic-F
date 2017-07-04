include ../make.inc

#SRC := $(wildcard *.f90)

SRC  = 000-type.f90 \
	 010-det.f90 \
	 010-mat.f90 \
	 010-split_work.f90 \
	 020-init_random_seed.f90 \
	 020-lin_decay.f90 \
	 020-near_small_prime_factoriz.f90 \
	 020-random_dist.f90 \
	 020-noise.f90 \
	 020-transfer_funct.f90 \
	 021-interpol.f90 \
	 022-interpol.f90 \
	 023-interpol.f90 \
	 024-interpol.f90 \
	 025-interpol.f90 \
	 026-interpol.f90 \
	 020-interpol.f90 \
	 100-filter_1D.f90 \
	 100-filter_1D_win.f90 \
	 100-fourier_trans_1D.f90 \
	 100-integrate.f90 \
	 100-space_2_filter_1D_win.f90 \
	 100-space_filter_1D_win.f90 \
	 200-filter_2D.f90 \
	 200-space_2_filter_2D.f90 \
	 200-space_filter_2D.f90 \
	 201-space_c_filter_2D.f90 \
	 200-space_filter_2D.f90 \
	 300-filter_1D_zp.f90 \
	 300-space_2_filter_1D.f90 \
	 300-space_2_filter_1D_win.f90 \
	 300-space_2_filter_2D.f90 \
	 300-time.f90 \
	 400-quick_sort.f90 \
	 400-filter_1D_zp_extra.f90 \
	 400-space_2_filter_1D_win_extra.f90
	 
OBJ=$(subst .f90,.o,$(SRC))
OBJS=$(addprefix ../objs/,$(OBJ))
MOD=$(wildcard ./*.mod)

#.SUFFIXES: .f90 .o
#OBJ=$(SRC:.f90=.o)

#.f90.o:
#	$(FC) $(FFLAGS) $(HEAD2) -c $<

$(LIB) : $(OBJS) $(MOD)
	$(AR) -rv $@ $? 

../objs/%.o : %.f90
	$(F90) -c $< -o $@

clean:
	rm -f $(LIB) $(OBJS) $(MOD)


