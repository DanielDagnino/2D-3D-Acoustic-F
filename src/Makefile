#!/bin/bash
include $(CWPROOT)/src/Makefile.config
D = $L/libcwp.a $L/libpar.a $L/libsu.a

#FC = gfortran
FC = mpif90

#FFLAGS = -Wall -cpp -ffree-line-length-none -fmax-errors=0 #-O3 -ffast-math -funroll-loops -m64 

#test_no_RK_loop
#-g -traceback -check all -fpp -Dusempi
#FFLAGS = -Wall -cpp -Dvalgrind -Dusempi -ffree-line-length-none -fmax-errors=0 #-O3 -ffast-math -funroll-loops -m64 
#FFLAGS = -g -Wall -cpp -Dvalgrind -Dusempi -ffree-line-length-none -fmax-errors=0 #-O3 -ffast-math -funroll-loops -m64 
#FFLAGS = -Wall -cpp -Dvalgrind -ffree-line-length-none -fmax-errors=0 #-O3 -ffast-math -funroll-loops -m64 
#FFLAGS = -Wall -cpp -ffree-line-length-none -fmax-errors=0 #-O3 -ffast-math -funroll-loops -m64 
#LIBS = -L./1-lib_base_v3_gnu/lib -ldag \
#       -L../../../fftpack5.1_gnu/lib -lfftpack \
#       -L../../../dfftpack5.1_gnu/lib -ldfftpack \
#	   -L../../../0-lapack-3.5.0_gnu -llapack -lrefblas -ltmglib \
#       sunmo_dag.o sudipfilt_dag.o -lstdc++ \
#       $(PRELFLAGS) -L$L  -lcwp -lsu -lpar -lm $(POSTLFLAGS)
#HEAD = -I./1-lib_base_v3_gnu/src

#FFLAGS = -g -traceback -check all -fpp -Dusempi

#FFLAGS = -g -traceback -fpp -Dusempi -Dtest_no_RK_loop -check bounds -check format -check arg_temp_created -check uninit -check output_conversion
#FFLAGS = -g -traceback -fpp -Dusempi -Dtest_no_RK_loop -check bounds -check format -check uninit -check output_conversion
#FFLAGS = -g -traceback -fpp -Dusempi -check bounds -check format -check uninit -check output_conversion -fp-stack-check -gen-interfaces -warn interfaces
#FFLAGS = -g -traceback -fpp -Dusempi -check bounds -check format -check arg_temp_created -check uninit -check output_conversion -fp-stack-check -gen-interfaces -warn interfaces
#FFLAGS = -g -traceback -fpp -Dusempi -check all -fp-stack-check -gen-interfaces -warn interfaces

#FFLAGS = -g -fpp -Dusempi \

#FFLAGS = -g -fpp -Dusempi \
#         -O3 -ipo-separate -unroll-aggressive -fp-model fast=2 -opt-prefetch=4 \
#         -m64 -xSSE4.2 -align \
#         -heap-arrays -mcmodel=large -shared-intel

FFLAGS = -fpp -Dusempi \
         -O3 -ipo -unroll-aggressive -fp-model fast=2 \
         -m64 -xSSE4.2 -align \
         -heap-arrays -mcmodel=large -shared-intel

#LIBS = -L./1-lib_base_v3-notemp/lib -ldag \
#       -L../../../fftpack5.1_gnu/lib -lfftpack \
#       -L../../../dfftpack5.1_gnu/lib -ldfftpack \
#	   -L../../../0-lapack-3.5.0_gnu -llapack -lrefblas -ltmglib \
#       sudipfilter_dag.o sunmo_dag.o -lstdc++ \
#       $(PRELFLAGS) -L$L  -lcwp -lsu -lpar -lm $(POSTLFLAGS)
         
LIBS = -L./1-lib_base_v3-notemp/lib -ldag \
       -L../../../fftpack5.1/lib -lfftpack \
       -L../../../dfftpack5.1/lib -ldfftpack \
	   -L../../../0-lapack-3.5.0 -llapack -lrefblas -ltmglib \
       sudipfilter_dag.o sunmo_dag.o -lstdc++ \
       $(PRELFLAGS) -L$L  -lcwp -lsu -lpar -lm $(POSTLFLAGS)
HEAD = -I./1-lib_base_v3-notemp/src

EXE = fwi

SRC = 00-m_data_kind.f90 \
	 m_mpi.f90 \
	 matrix_inversion.f90 \
	 00-module.f90 \
	 m_solv3Dahvc_data_kind.f90 \
	 01-correctness_criterion.f90 \
	 01-m_mat.f90 \
	 01-mod_definit.f90 \
	 m_splines.f90 \
	 m_support_funct.f90 \
	 m_sg_path_file.f90 \
	 m_sg_data.f90 \
	 m_geo.f90 \
	 m_sg_type.f90 \
	 m_grid_bath.f90 \
	 m_sg_func.f90 \
	 m_sg_get_pos.f90 \
	 m_sg_get_sou.f90 \
	 m_sg_get_adj_sou.f90 \
	 m_sg_get_tr.f90 \
	 m_sg_get_tr_picks.f90 \
	 m_sg.f90 \
	 m_sg_filt.f90 \
	 02-module_iH.f90 \
	 m_mygeo.f90 \
	 m_get_freq_inv.f90 \
	 m_solver_source.f90 \
	 m_interpol_sg.f90 \
	 m_solver.f90 \
	 m_deriv_b.f90 \
	 m_deriv.f90 \
	 m_PML_definitions_param_face.f90 \
	 m_PML_building.f90 \
	 m_coll_allo.f90 \
	 m_read_field.f90 \
	 m_write_field.f90 \
	 alpha_cota.f90 \
	 coef_norm_grad.f90 \
	 control_grad.f90 \
	 m_fix_edge.f90 \
	 m_calib.f90 \
	 m_sg_pr.f90 \
	 02-module_iH_red.f90 \
	 m_adj_sou_L1.f90 \
	 m_adj_sou_L2.f90 \
	 m_adj_sou_L2_array.f90 \
	 m_data_proc.f90 \
	 m_adj_sou_L2_water.f90 \
	 m_adj_sou_L2_water_Wd.f90 \
	 m_adj_sou_L2_water_Wd_phase.f90 \
	 m_adj_sou_WCC_EA.f90 \
	 m_adj_sou_L2_EA.f90 \
	 m_adj_sou_L2_dip_EA.f90 \
	 m_adj_sou_L2_dip_EA_weight.f90 \
	 m_adj_sou_L2_EA_cc.f90 \
	 m_adj_sou_p.f90 \
	 m_adj_sou_e.f90 \
	 m_adj_sou_CC_W.f90 \
	 m_apply_weight.f90 \
	 m_apply_DRO.f90 \
	 m_change_var.f90 \
	 m_clear_div_sou.f90 \
	 m_copy_mat_2d.f90 \
	 m_copy_mat_3D.f90 \
	 m_copy_mat_4.f90 \
	 m_get_id_shots.f90 \
	 m_get_initial_data.f90 \
	 m_get_strategy.f90 \
	 param_control.f90 \
	 m_get_model_and_QC.f90 \
	 m_prec_mod.f90 \
	 m_poly_fit_min.f90 \
	 m_model_modif.f90 \
	 m_source_ricker.f90 \
	 non_lin_conj_grad.f90 \
	 norm_grad.f90 \
	 relat_change.f90 \
	 RK_der_edge_BD.f90 \
	 RK_der_edge_BU.f90 \
	 RK_der_edge_FD.f90 \
	 RK_der_edge_FU.f90 \
	 RK_der_edge_LB.f90 \
	 RK_der_edge_LD.f90 \
	 RK_der_edge_LF.f90 \
	 RK_der_edge_LU.f90 \
	 RK_der_edge_RB.f90 \
	 RK_der_edge_RD.f90 \
	 RK_der_edge_RF.f90 \
	 RK_der_edge_RU.f90 \
	 RK_der.f90 \
	 RK_der_face_B.f90 \
	 RK_der_face_D.f90 \
	 RK_der_face_F.f90 \
	 RK_der_face_L.f90 \
	 RK_der_face_R.f90 \
	 RK_der_face_U.f90 \
	 RK_der_vertex_LBD.f90 \
	 RK_der_vertex_LBU.f90 \
	 RK_der_vertex_LFD.f90 \
	 RK_der_vertex_LFU.f90 \
	 RK_der_vertex_RBD.f90 \
	 RK_der_vertex_RBU.f90 \
	 RK_der_vertex_RFD.f90 \
	 RK_der_vertex_RFU.f90 \
	 RK_dif_op_FS.f90 \
	 RK_dif_op_FS_back.f90 \
	 RK_dif_op_split_back.f90 \
	 RK_dif_op_split_bd.f90 \
	 RK_dif_op_split.f90 \
	 RK_new_step_edge.f90 \
	 RK_new_step.f90 \
	 RK_new_step_face.f90 \
	 RK_new_step_FS.f90 \
	 RK_new_step_vertex.f90 \
	 RK_new_step_split.f90 \
	 RK_save_boundary.f90 \
	 RK_sou_adj.f90 \
	 RK_sou.f90 \
	 RK_step_edge.f90 \
	 RK_step.f90 \
	 RK_step_face.f90 \
	 RK_step_FS.f90 \
	 RK_step_vertex.f90 \
	 RK_step_split_bd.f90 \
	 RK_step_split.f90 \
	 RK_update_boundary_2.f90 \
	 RK_update_boundary.f90 \
	 RK_update_boundary_FS.f90 \
	 RK_loop_adj.f90 \
	 RK_loop_back_2.f90 \
	 RK_loop_back_aux.f90 \
	 RK_loop_back.f90 \
	 RK_loop.f90 \
	 time_shift.f90 \
	 m_smooth_op.f90 \
	 solver_acoustic_3d_vnh_RK4_CPML_adj.f90 \
	 solver_acoustic_3d_vnh_RK4_CPML_deallo.f90 \
	 solver_acoustic_3d_vnh_RK4_CPML_dir.f90 \
	 solver_acoustic_3d_vnh_RK4_CPML_test.f90 \
	 solv3Dahv_test_all_nmo.f90 \
	 m_comp_J.f90 \
	 m_comp_J_NMO.f90 \
	 solv3Dahv_dir_ad.f90 \
	 solv3Dahv_test_all.f90 \
	 source_inversion.f90 \
	 search_min.f90 \
	 correct_search.f90 \
	 m_minimize.f90 \
	 trial_source.f90 \
	 zz-fwi.f90 \
	 zz-main.f90
	 
# No need to edit below this line
.SUFFIXES: .f90 .o

OBJ = $(SRC:.f90=.o)

.f90.o:
	$(FC) $(FFLAGS) $(HEAD) -c $<

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) $(OBJ) $(HEAD) $(LIBS) -o $@ 

$(OBJ): $(MF)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) core *.mod *.MOD





