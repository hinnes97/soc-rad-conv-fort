# Makefile created by mkmf $Id: mkmf,v 18.0 2010/03/02 23:26:08 fms Exp $ 

SRCROOT = /home/i/innes/matrixrt/fortran/

CPPDEFS = -DSOC


include /home/i/innes/matrixrt/fortran/mkmf.template.oxford_ubuntu_1804


.DEFAULT:
	-echo $@ does not exist.
all: main.exe
adjust.o: $(SRCROOT)src/adjust.f90 phys_new.o params.o condense.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/adjust.f90
adjust_ir_radiance.o: $(SRCROOT)src/socrates/src/radiance_core/adjust_ir_radiance.F90 realtype_rd.o def_control.o def_dimen.o def_atm.o def_out.o def_planck.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/adjust_ir_radiance.F90
aggregate_cloud.o: $(SRCROOT)src/socrates/src/radiance_core/aggregate_cloud.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/aggregate_cloud.F90
augment_radiance.o: $(SRCROOT)src/socrates/src/radiance_core/augment_radiance.F90 realtype_rd.o def_control.o def_out.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/augment_radiance.F90
augment_tiled_radiance.o: $(SRCROOT)src/socrates/src/radiance_core/augment_tiled_radiance.F90 realtype_rd.o def_control.o def_out.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/augment_tiled_radiance.F90
band_grey_rad.o: $(SRCROOT)src/band_grey_rad.f90 io.o utils.o params.o twostr_module.o nc_wrap.o phys.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/band_grey_rad.f90
band_solver.o: $(SRCROOT)src/socrates/src/radiance_core/band_solver.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/band_solver.F90
build_sph_matrix.o: $(SRCROOT)src/socrates/src/radiance_core/build_sph_matrix.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/build_sph_matrix.F90
calc_brdf.o: $(SRCROOT)src/socrates/src/radiance_core/calc_brdf.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_brdf.F90
calc_cg_coeff.o: $(SRCROOT)src/socrates/src/radiance_core/calc_cg_coeff.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_cg_coeff.F90
calc_contrib_func.o: $(SRCROOT)src/socrates/src/radiance_core/calc_contrib_func.F90 rad_ccf.o diffusivity_factor.o realtype_rd.o def_ss_prop.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_contrib_func.F90
calc_flux_ipa.o: $(SRCROOT)src/socrates/src/radiance_core/calc_flux_ipa.F90 realtype_rd.o def_control.o def_bound.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_flux_ipa.F90
calc_gauss_weight_90.o: $(SRCROOT)src/socrates/src/radiance_core/calc_gauss_weight_90.F90 realtype_rd.o rad_pcf.o errormessagelength_mod.o ereport_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_gauss_weight_90.F90
calc_radiance_ipa.o: $(SRCROOT)src/socrates/src/radiance_core/calc_radiance_ipa.F90 realtype_rd.o def_ss_prop.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_radiance_ipa.F90
calc_surf_rad.o: $(SRCROOT)src/socrates/src/radiance_core/calc_surf_rad.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_surf_rad.F90
calc_top_rad.o: $(SRCROOT)src/socrates/src/radiance_core/calc_top_rad.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_top_rad.F90
calc_uplm_sol.o: $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_sol.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_uplm_sol.F90
calc_uplm_zero.o: $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_zero.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/calc_uplm_zero.F90
cg_kappa_ms.o: $(SRCROOT)src/socrates/src/radiance_core/cg_kappa_ms.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/cg_kappa_ms.F90
check_phf_term.o: $(SRCROOT)src/socrates/src/radiance_core/check_phf_term.F90 rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/check_phf_term.F90
circumsolar_fraction.o: $(SRCROOT)src/socrates/src/radiance_core/circumsolar_fraction.F90 realtype_rd.o legendre_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/circumsolar_fraction.F90
cloud_maxcs_split.o: $(SRCROOT)src/socrates/src/radiance_core/cloud_maxcs_split.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/cloud_maxcs_split.F90
column_solver.o: $(SRCROOT)src/socrates/src/radiance_core/column_solver.F90 realtype_rd.o def_control.o def_bound.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o spherical_solar_source.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/column_solver.F90
compress_spectrum.o: $(SRCROOT)src/socrates/interface/compress_spectrum.F90 def_control.o def_spectrum.o gas_list_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/compress_spectrum.F90
condense.o: $(SRCROOT)src/condense.f90 params.o phys_new.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/condense.f90
convection.o: $(SRCROOT)src/convection.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/convection.f90
copy_clr_full.o: $(SRCROOT)src/socrates/src/radiance_core/copy_clr_full.F90 realtype_rd.o yomhook.o parkind1.o def_control.o rad_pcf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/copy_clr_full.F90
copy_clr_sol.o: $(SRCROOT)src/socrates/src/radiance_core/copy_clr_sol.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/copy_clr_sol.F90
def_aer.o: $(SRCROOT)src/socrates/src/radiance_core/def_aer.F90 realtype_rd.o def_dimen.o def_spectrum.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_aer.F90
def_atm.o: $(SRCROOT)src/socrates/src/radiance_core/def_atm.F90 realtype_rd.o def_dimen.o def_spectrum.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_atm.F90
def_bound.o: $(SRCROOT)src/socrates/src/radiance_core/def_bound.F90 realtype_rd.o def_dimen.o def_spectrum.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_bound.F90
def_cld.o: $(SRCROOT)src/socrates/src/radiance_core/def_cld.F90 realtype_rd.o def_dimen.o def_spectrum.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_cld.F90
def_control.o: $(SRCROOT)src/socrates/src/radiance_core/def_control.F90 filenamelength_mod.o missing_data_mod.o realtype_rd.o def_spectrum.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_control.F90
def_dimen.o: $(SRCROOT)src/socrates/src/radiance_core/def_dimen.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_dimen.F90
def_mcica.o: $(SRCROOT)src/socrates/src/radiance_core/def_mcica.F90 realtype_rd.o def_spectrum.o yomhook.o parkind1.o ereport_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_mcica.F90
def_out.o: $(SRCROOT)src/socrates/src/radiance_core/def_out.F90 realtype_rd.o def_control.o def_dimen.o def_spectrum.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_out.F90
def_planck.o: $(SRCROOT)src/socrates/src/radiance_core/def_planck.F90 realtype_rd.o def_dimen.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_planck.F90
def_spectrum.o: $(SRCROOT)src/socrates/src/radiance_core/def_spectrum.F90 realtype_rd.o missing_data_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_spectrum.F90
def_spherical_geometry.o: $(SRCROOT)src/socrates/src/radiance_core/def_spherical_geometry.F90 realtype_rd.o def_dimen.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_spherical_geometry.F90
def_ss_prop.o: $(SRCROOT)src/socrates/src/radiance_core/def_ss_prop.F90 realtype_rd.o def_dimen.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/def_ss_prop.F90
diff_albedo_basis.o: $(SRCROOT)src/socrates/src/radiance_core/diff_albedo_basis.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/diff_albedo_basis.F90
diff_planck_source_mod.o: $(SRCROOT)src/socrates/src/radiance_core/diff_planck_source_mod.F90 def_control.o def_dimen.o def_spectrum.o def_atm.o def_bound.o def_planck.o realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o planck_flux_band_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/diff_planck_source_mod.F90
diffusivity_factor.o: $(SRCROOT)src/socrates/src/radiance_core/diffusivity_factor.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/diffusivity_factor.F90
dimensions_spec_ucf.o: $(SRCROOT)src/socrates/src/modules_core/dimensions_spec_ucf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/dimensions_spec_ucf.F90
eig_sys.o: $(SRCROOT)src/socrates/src/radiance_core/eig_sys.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/eig_sys.F90
eigenvalue_tri.o: $(SRCROOT)src/socrates/src/radiance_core/eigenvalue_tri.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/eigenvalue_tri.F90
ereport_mod.o: $(SRCROOT)src/socrates/src/modules_core/ereport_mod.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/ereport_mod.F90
errormessagelength_mod.o: $(SRCROOT)src/socrates/src/modules_core/errormessagelength_mod.F90 filenamelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/errormessagelength_mod.F90
eval_uplm.o: $(SRCROOT)src/socrates/src/radiance_core/eval_uplm.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/eval_uplm.F90
file_manager.o: $(SRCROOT)src/socrates/src/modules_core/file_manager.F90 ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/file_manager.F90
filenamelength_mod.o: $(SRCROOT)src/socrates/src/modules_core/filenamelength_mod.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/filenamelength_mod.F90
flux.o: $(SRCROOT)src/flux.F90 params.o condense.o new_socrates_interface.o tsrad_short_char.o toon.o short_char_split_layer.o ross.o short_char_ross.o rad_ross.o band_grey_rad.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/flux.F90
fsolve.o: $(SRCROOT)src/fsolve.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/fsolve.f90
gas_list_pcf.o: $(SRCROOT)src/socrates/src/radiance_core/gas_list_pcf.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/gas_list_pcf.F90
gas_optical_properties.o: $(SRCROOT)src/socrates/src/radiance_core/gas_optical_properties.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/gas_optical_properties.F90
gauss_angle.o: $(SRCROOT)src/socrates/src/radiance_core/gauss_angle.F90 realtype_rd.o rad_pcf.o rad_ccf.o gaussian_weight_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/gauss_angle.F90
gaussian_weight_pcf.o: $(SRCROOT)src/socrates/src/radiance_core/gaussian_weight_pcf.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/gaussian_weight_pcf.F90
grey_opt_prop.o: $(SRCROOT)src/socrates/src/radiance_core/grey_opt_prop.F90 realtype_rd.o def_control.o def_out.o def_ss_prop.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/grey_opt_prop.F90
hemi_sph_integ.o: $(SRCROOT)src/socrates/src/radiance_core/hemi_sph_integ.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/hemi_sph_integ.F90
increment_rad_cf.o: $(SRCROOT)src/socrates/src/radiance_core/increment_rad_cf.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/increment_rad_cf.F90
inter_k.o: $(SRCROOT)src/socrates/src/radiance_core/inter_k.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/inter_k.F90
inter_pt.o: $(SRCROOT)src/socrates/src/radiance_core/inter_pt.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/inter_pt.F90
inter_pt_lookup.o: $(SRCROOT)src/socrates/src/radiance_core/inter_pt_lookup.F90 realtype_rd.o rad_ccf.o gas_list_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/inter_pt_lookup.F90
inter_t_lookup.o: $(SRCROOT)src/socrates/src/radiance_core/inter_t_lookup.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/inter_t_lookup.F90
interp1d.o: $(SRCROOT)src/socrates/src/radiance_core/interp1d.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/interp1d.F90
io.o: $(SRCROOT)src/io.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/io.f90
ir_source.o: $(SRCROOT)src/socrates/src/radiance_core/ir_source.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/ir_source.F90
layer_part_integ.o: $(SRCROOT)src/socrates/src/radiance_core/layer_part_integ.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/layer_part_integ.F90
legendre_mod.o: $(SRCROOT)src/socrates/src/radiance_core/legendre_mod.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/legendre_mod.F90
legendre_weight.o: $(SRCROOT)src/socrates/src/radiance_core/legendre_weight.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/legendre_weight.F90
main.o: $(SRCROOT)src/main.F90 io.o utils.o params.o timestep.o matrix.o condense.o new_socrates_interface.o rad_ross.o band_grey_rad.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/main.F90
matrix.o: $(SRCROOT)src/matrix.f90 params.o flux.o utils.o condense.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/matrix.f90
mcica_column.o: $(SRCROOT)src/socrates/src/radiance_core/mcica_column.F90 realtype_rd.o def_control.o def_cld.o def_bound.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o set_n_source_coeff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/mcica_column.F90
mcica_sample.o: $(SRCROOT)src/socrates/src/radiance_core/mcica_sample.F90 realtype_rd.o def_control.o def_dimen.o def_atm.o def_cld.o def_bound.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/mcica_sample.F90
missing_data_mod.o: $(SRCROOT)src/socrates/src/modules_core/missing_data_mod.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/missing_data_mod.F90
mix_app_scat.o: $(SRCROOT)src/socrates/src/radiance_core/mix_app_scat.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/mix_app_scat.F90
mix_column.o: $(SRCROOT)src/socrates/src/radiance_core/mix_column.F90 realtype_rd.o def_control.o def_cld.o def_bound.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o set_n_source_coeff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/mix_column.F90
mixed_solar_source.o: $(SRCROOT)src/socrates/src/radiance_core/mixed_solar_source.F90 realtype_rd.o def_control.o def_bound.o rad_pcf.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/mixed_solar_source.F90
moisture_mod.o: $(SRCROOT)src/moisture_mod.f90 params.o condense.o phys_new.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/moisture_mod.f90
monochromatic_gas_flux.o: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_gas_flux.F90 realtype_rd.o rad_pcf.o vectlib_mod.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/monochromatic_gas_flux.F90
monochromatic_ir_radiance.o: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_ir_radiance.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/monochromatic_ir_radiance.F90
monochromatic_radiance.o: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance.F90 realtype_rd.o def_control.o def_atm.o def_cld.o def_bound.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o spherical_trans_coeff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance.F90
monochromatic_radiance_sph.o: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_sph.F90 realtype_rd.o def_control.o def_ss_prop.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_sph.F90
monochromatic_radiance_tseq.o: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_tseq.F90 realtype_rd.o def_control.o def_cld.o def_bound.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_tseq.F90
nc_wrap.o: $(SRCROOT)src/nc_wrap.f90 io.o params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/nc_wrap.f90
new_socrates_calc.o: $(SRCROOT)src/socrates/interface/new_socrates_calc.f90 rad_pcf.o rad_ccf.o def_control.o def_spectrum.o def_dimen.o def_atm.o def_bound.o def_cld.o def_aer.o def_out.o set_control.o set_dimen.o set_atm.o set_bound.o set_cld.o set_aer.o soc_constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/new_socrates_calc.f90
new_socrates_interface.o: $(SRCROOT)src/socrates/interface/new_socrates_interface.f90 read_control.o def_control.o def_spectrum.o soc_constants.o socrates_config_mod.o params.o rad_ccf.o realtype_rd.o new_socrates_calc.o compress_spectrum.o def_dimen.o def_atm.o def_cld.o def_aer.o def_bound.o def_out.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/new_socrates_interface.f90
opt_prop_aerosol.o: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_aerosol.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/opt_prop_aerosol.F90
opt_prop_ice_cloud.o: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ice_cloud.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/opt_prop_ice_cloud.F90
opt_prop_inhom_corr_cairns.o: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_inhom_corr_cairns.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/opt_prop_inhom_corr_cairns.F90
opt_prop_ukca_aerosol.o: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ukca_aerosol.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/opt_prop_ukca_aerosol.F90
opt_prop_water_cloud.o: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_water_cloud.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/opt_prop_water_cloud.F90
overlap_coupled.o: $(SRCROOT)src/socrates/src/radiance_core/overlap_coupled.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/overlap_coupled.F90
params.o: $(SRCROOT)src/params.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/params.f90
parkind1.o: $(SRCROOT)src/socrates/src/modules_core/parkind1.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/parkind1.F90
phys.o: $(SRCROOT)src/phys.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/phys.f90
phys_new.o: $(SRCROOT)src/phys_new.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/phys_new.f90
planck_flux_band_mod.o: $(SRCROOT)src/socrates/src/radiance_core/planck_flux_band_mod.F90 realtype_rd.o def_spectrum.o interp1d.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/planck_flux_band_mod.F90
prsc_gather_spline.o: $(SRCROOT)src/socrates/src/radiance_core/prsc_gather_spline.F90 realtype_rd.o rad_pcf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/prsc_gather_spline.F90
prsc_opt_prop.o: $(SRCROOT)src/socrates/src/radiance_core/prsc_opt_prop.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/prsc_opt_prop.F90
quicksort.o: $(SRCROOT)src/socrates/src/radiance_core/quicksort.F90 realtype_rd.o rad_pcf.o errormessagelength_mod.o ereport_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/quicksort.F90
rad_ccf.o: $(SRCROOT)src/socrates/src/modules_core/rad_ccf.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/rad_ccf.F90
rad_pcf.o: $(SRCROOT)src/socrates/src/radiance_core/rad_pcf.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/rad_pcf.F90
rad_ross.o: $(SRCROOT)src/rad_ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/rad_ross.f90
radiance_calc.o: $(SRCROOT)src/socrates/src/radiance_core/radiance_calc.F90 realtype_rd.o rad_pcf.o rad_ccf.o def_spectrum.o def_dimen.o def_control.o def_atm.o def_cld.o def_aer.o def_bound.o def_out.o def_ss_prop.o gas_list_pcf.o errormessagelength_mod.o yomhook.o parkind1.o ereport_mod.o def_planck.o def_spherical_geometry.o spherical_path.o diff_planck_source_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/radiance_calc.F90
read_control.o: $(SRCROOT)src/socrates/interface/read_control.F90 rad_pcf.o def_control.o def_spectrum.o socrates_config_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/read_control.F90
read_spectrum.o: $(SRCROOT)src/socrates/src/radiance_core/read_spectrum.F90 realtype_rd.o def_spectrum.o rad_pcf.o dimensions_spec_ucf.o gas_list_pcf.o file_manager.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/read_spectrum.F90
realtype_rd.o: $(SRCROOT)src/socrates/src/modules_core/realtype_rd.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/realtype_rd.f90
rebin_esft_terms.o: $(SRCROOT)src/socrates/src/radiance_core/rebin_esft_terms.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/rebin_esft_terms.F90
rescale_continuum.o: $(SRCROOT)src/socrates/src/radiance_core/rescale_continuum.F90 realtype_rd.o rad_pcf.o rad_ccf.o def_control.o vectlib_mod.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/rescale_continuum.F90
rescale_phase_fnc.o: $(SRCROOT)src/socrates/src/radiance_core/rescale_phase_fnc.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/rescale_phase_fnc.F90
rescale_tau_csr.o: $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_csr.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/rescale_tau_csr.F90
rescale_tau_omega.o: $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_omega.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/rescale_tau_omega.F90
ross.o: $(SRCROOT)src/ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/ross.f90
scale_absorb.o: $(SRCROOT)src/socrates/src/radiance_core/scale_absorb.F90 realtype_rd.o rad_pcf.o vectlib_mod.o scale_wenyi.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/scale_absorb.F90
scale_wenyi.o: $(SRCROOT)src/socrates/src/radiance_core/scale_wenyi.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/scale_wenyi.F90
ses_rescale_contm.o: $(SRCROOT)src/socrates/src/radiance_core/ses_rescale_contm.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/ses_rescale_contm.F90
set_aer.o: $(SRCROOT)src/socrates/interface/set_aer.F90 rad_pcf.o def_spectrum.o def_dimen.o def_control.o def_aer.o soc_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/set_aer.F90
set_atm.o: $(SRCROOT)src/socrates/interface/set_atm.F90 rad_pcf.o def_control.o def_dimen.o def_spectrum.o def_atm.o socrates_config_mod.o gas_list_pcf.o soc_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/set_atm.F90
set_bound.o: $(SRCROOT)src/socrates/interface/set_bound.F90 rad_pcf.o def_spectrum.o def_dimen.o def_control.o def_bound.o soc_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/set_bound.F90
set_cld.o: $(SRCROOT)src/socrates/interface/set_cld.F90 rad_pcf.o def_spectrum.o def_dimen.o def_control.o def_cld.o soc_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/set_cld.F90
set_cloud_geometry.o: $(SRCROOT)src/socrates/src/radiance_core/set_cloud_geometry.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_cloud_geometry.F90
set_cloud_pointer.o: $(SRCROOT)src/socrates/src/radiance_core/set_cloud_pointer.F90 rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_cloud_pointer.F90
set_control.o: $(SRCROOT)src/socrates/interface/set_control.F90 def_control.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/set_control.F90
set_dimen.o: $(SRCROOT)src/socrates/interface/set_dimen.F90 rad_pcf.o def_control.o def_dimen.o def_spectrum.o soc_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/set_dimen.F90
set_dirn_weights.o: $(SRCROOT)src/socrates/src/radiance_core/set_dirn_weights.F90 realtype_rd.o rad_pcf.o rad_ccf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_dirn_weights.F90
set_level_weights.o: $(SRCROOT)src/socrates/src/radiance_core/set_level_weights.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_level_weights.F90
set_matrix_pentadiagonal.o: $(SRCROOT)src/socrates/src/radiance_core/set_matrix_pentadiagonal.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_matrix_pentadiagonal.F90
set_n_cloud_parameter.o: $(SRCROOT)src/socrates/src/radiance_core/set_n_cloud_parameter.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_n_cloud_parameter.F90
set_n_source_coeff.o: $(SRCROOT)src/socrates/src/radiance_core/set_n_source_coeff.F90 rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_n_source_coeff.F90
set_rad_layer.o: $(SRCROOT)src/socrates/src/radiance_core/set_rad_layer.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_rad_layer.F90
set_truncation.o: $(SRCROOT)src/socrates/src/radiance_core/set_truncation.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/set_truncation.F90
shell_sort.o: $(SRCROOT)src/socrates/src/radiance_core/shell_sort.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/shell_sort.F90
short_char_ross.o: $(SRCROOT)src/short_char_ross.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/short_char_ross.f90
short_char_split_layer.o: $(SRCROOT)src/short_char_split_layer.f90 params.o tau_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/short_char_split_layer.f90
single_scat_sol.o: $(SRCROOT)src/socrates/src/radiance_core/single_scat_sol.F90 realtype_rd.o rad_ccf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/single_scat_sol.F90
single_scattering.o: $(SRCROOT)src/socrates/src/radiance_core/single_scattering.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/single_scattering.F90
single_scattering_all.o: $(SRCROOT)src/socrates/src/radiance_core/single_scattering_all.F90 realtype_rd.o def_ss_prop.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/single_scattering_all.F90
soc_constants.o: $(SRCROOT)src/socrates/interface/soc_constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/soc_constants.F90
socrates_config_mod.o: $(SRCROOT)src/socrates/interface/socrates_config_mod.f90 def_control.o def_spectrum.o soc_constants.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/interface/socrates_config_mod.f90
sol_scat_cos.o: $(SRCROOT)src/socrates/src/radiance_core/sol_scat_cos.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/sol_scat_cos.F90
solar_coefficient_basic.o: $(SRCROOT)src/socrates/src/radiance_core/solar_coefficient_basic.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o def_control.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solar_coefficient_basic.F90
solar_source.o: $(SRCROOT)src/socrates/src/radiance_core/solar_source.F90 realtype_rd.o def_control.o def_bound.o rad_pcf.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solar_source.F90
solve_band_k_eqv.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv.F90 realtype_rd.o def_control.o def_dimen.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o diffusivity_factor.o vectlib_mod.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv.F90
solve_band_k_eqv_scl.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv_scl.F90 realtype_rd.o def_control.o def_dimen.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o diffusivity_factor.o vectlib_mod.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv_scl.F90
solve_band_one_gas.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_one_gas.F90 realtype_rd.o def_control.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_one_gas.F90
solve_band_random_overlap.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap.F90 realtype_rd.o def_control.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap.F90
solve_band_random_overlap_resort_rebin.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap_resort_rebin.F90 realtype_rd.o def_control.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap_resort_rebin.F90
solve_band_ses.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_ses.F90 realtype_rd.o def_control.o def_dimen.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_ses.F90
solve_band_without_gas.o: $(SRCROOT)src/socrates/src/radiance_core/solve_band_without_gas.F90 realtype_rd.o def_control.o def_atm.o def_cld.o def_bound.o def_out.o def_planck.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solve_band_without_gas.F90
solver_homogen_direct.o: $(SRCROOT)src/socrates/src/radiance_core/solver_homogen_direct.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_homogen_direct.F90
solver_mix_direct.o: $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct.F90
solver_mix_direct_hogan.o: $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct_hogan.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct_hogan.F90
solver_no_scat.o: $(SRCROOT)src/socrates/src/radiance_core/solver_no_scat.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_no_scat.F90
solver_triple.o: $(SRCROOT)src/socrates/src/radiance_core/solver_triple.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_triple.F90
solver_triple_app_scat.o: $(SRCROOT)src/socrates/src/radiance_core/solver_triple_app_scat.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_triple_app_scat.F90
solver_triple_hogan.o: $(SRCROOT)src/socrates/src/radiance_core/solver_triple_hogan.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/solver_triple_hogan.F90
sph_matrix_solver.o: $(SRCROOT)src/socrates/src/radiance_core/sph_matrix_solver.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/sph_matrix_solver.F90
sph_solver.o: $(SRCROOT)src/socrates/src/radiance_core/sph_solver.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/sph_solver.F90
spherical_path.o: $(SRCROOT)src/socrates/src/radiance_core/spherical_path.F90 realtype_rd.o def_dimen.o def_atm.o def_bound.o def_spherical_geometry.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/spherical_path.F90
spherical_solar_source.o: $(SRCROOT)src/socrates/src/radiance_core/spherical_solar_source.F90 realtype_rd.o def_control.o def_bound.o def_spherical_geometry.o rad_pcf.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/spherical_solar_source.F90
spherical_trans_coeff.o: $(SRCROOT)src/socrates/src/radiance_core/spherical_trans_coeff.F90 realtype_rd.o def_ss_prop.o def_spherical_geometry.o vectlib_mod.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/spherical_trans_coeff.F90
spline_evaluate.o: $(SRCROOT)src/socrates/src/radiance_core/spline_evaluate.F90 realtype_rd.o rad_pcf.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/spline_evaluate.F90
spline_fit.o: $(SRCROOT)src/socrates/src/radiance_core/spline_fit.F90 realtype_rd.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/spline_fit.F90
sum_k.o: $(SRCROOT)src/socrates/src/radiance_core/sum_k.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/sum_k.F90
tau_mod.o: $(SRCROOT)src/tau_mod.f90 params.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/tau_mod.f90
timestep.o: $(SRCROOT)src/timestep.f90 params.o flux.o utils.o convection.o condense.o moisture_mod.o adjust.o io.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/timestep.f90
toon.o: $(SRCROOT)src/toon.F90 params.o tau_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/toon.F90
trans_source_coeff.o: $(SRCROOT)src/socrates/src/radiance_core/trans_source_coeff.F90 realtype_rd.o rad_pcf.o vectlib_mod.o yomhook.o parkind1.o def_control.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/trans_source_coeff.F90
triple_column.o: $(SRCROOT)src/socrates/src/radiance_core/triple_column.F90 realtype_rd.o def_control.o def_cld.o def_bound.o def_ss_prop.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o set_n_source_coeff.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/triple_column.F90
triple_solar_source.o: $(SRCROOT)src/socrates/src/radiance_core/triple_solar_source.F90 realtype_rd.o def_control.o def_bound.o rad_pcf.o parkind1.o yomhook.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/triple_solar_source.F90
tsrad_short_char.o: $(SRCROOT)src/tsrad_short_char.f90 params.o tau_mod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/tsrad_short_char.f90
two_coeff.o: $(SRCROOT)src/socrates/src/radiance_core/two_coeff.F90 realtype_rd.o rad_pcf.o vectlib_mod.o yomhook.o parkind1.o def_control.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_coeff.F90
two_coeff_basic.o: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_basic.F90 realtype_rd.o rad_pcf.o diffusivity_factor.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_coeff_basic.F90
two_coeff_cloud.o: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_cloud.F90 realtype_rd.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o def_control.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_coeff_cloud.F90
two_coeff_fast_lw.o: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_fast_lw.F90 realtype_rd.o vectlib_mod.o rad_pcf.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_coeff_fast_lw.F90
two_coeff_region.o: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region.F90 realtype_rd.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o def_control.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_coeff_region.F90
two_coeff_region_fast_lw.o: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region_fast_lw.F90 realtype_rd.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o errormessagelength_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_coeff_region_fast_lw.F90
two_stream.o: $(SRCROOT)src/socrates/src/radiance_core/two_stream.F90 realtype_rd.o def_control.o def_bound.o def_spherical_geometry.o rad_pcf.o yomhook.o parkind1.o ereport_mod.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/radiance_core/two_stream.F90
twostr_module.o: $(SRCROOT)src/twostr_module.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/twostr_module.f90
utils.o: $(SRCROOT)src/utils.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/utils.f90
vectlib_mod.o: $(SRCROOT)src/socrates/src/modules_core/vectlib_mod.F90 realtype_rd.o yomhook.o parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/vectlib_mod.F90
yomhook.o: $(SRCROOT)src/socrates/src/modules_core/yomhook.F90 parkind1.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)src/socrates/src/modules_core/yomhook.F90
./spherical_path.F90: $(SRCROOT)src/socrates/src/radiance_core/spherical_path.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/spherical_path.F90 .
./spherical_trans_coeff.F90: $(SRCROOT)src/socrates/src/radiance_core/spherical_trans_coeff.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/spherical_trans_coeff.F90 .
./mcica_column.F90: $(SRCROOT)src/socrates/src/radiance_core/mcica_column.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/mcica_column.F90 .
./monochromatic_radiance_sph.F90: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_sph.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_sph.F90 .
./single_scattering.F90: $(SRCROOT)src/socrates/src/radiance_core/single_scattering.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/single_scattering.F90 .
./rescale_phase_fnc.F90: $(SRCROOT)src/socrates/src/radiance_core/rescale_phase_fnc.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/rescale_phase_fnc.F90 .
./matrix.f90: $(SRCROOT)src/matrix.f90
	cp $(SRCROOT)src/matrix.f90 .
./eigenvalue_tri.F90: $(SRCROOT)src/socrates/src/radiance_core/eigenvalue_tri.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/eigenvalue_tri.F90 .
./hemi_sph_integ.F90: $(SRCROOT)src/socrates/src/radiance_core/hemi_sph_integ.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/hemi_sph_integ.F90 .
./calc_contrib_func.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_contrib_func.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_contrib_func.F90 .
./opt_prop_ice_cloud.F90: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ice_cloud.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ice_cloud.F90 .
./opt_prop_ukca_aerosol.F90: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ukca_aerosol.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ukca_aerosol.F90 .
./calc_gauss_weight_90.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_gauss_weight_90.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_gauss_weight_90.F90 .
./solve_band_one_gas.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_one_gas.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_one_gas.F90 .
./condense.f90: $(SRCROOT)src/condense.f90
	cp $(SRCROOT)src/condense.f90 .
./prsc_gather_spline.F90: $(SRCROOT)src/socrates/src/radiance_core/prsc_gather_spline.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/prsc_gather_spline.F90 .
./single_scattering_all.F90: $(SRCROOT)src/socrates/src/radiance_core/single_scattering_all.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/single_scattering_all.F90 .
./opt_prop_water_cloud.F90: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_water_cloud.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/opt_prop_water_cloud.F90 .
./diff_albedo_basis.F90: $(SRCROOT)src/socrates/src/radiance_core/diff_albedo_basis.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/diff_albedo_basis.F90 .
./solver_homogen_direct.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_homogen_direct.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_homogen_direct.F90 .
./short_char_split_layer.f90: $(SRCROOT)src/short_char_split_layer.f90
	cp $(SRCROOT)src/short_char_split_layer.f90 .
./eval_uplm.F90: $(SRCROOT)src/socrates/src/radiance_core/eval_uplm.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/eval_uplm.F90 .
./solar_coefficient_basic.F90: $(SRCROOT)src/socrates/src/radiance_core/solar_coefficient_basic.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solar_coefficient_basic.F90 .
./two_coeff_region_fast_lw.F90: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region_fast_lw.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region_fast_lw.F90 .
./rad_ross.f90: $(SRCROOT)src/rad_ross.f90
	cp $(SRCROOT)src/rad_ross.f90 .
./set_level_weights.F90: $(SRCROOT)src/socrates/src/radiance_core/set_level_weights.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_level_weights.F90 .
./aggregate_cloud.F90: $(SRCROOT)src/socrates/src/radiance_core/aggregate_cloud.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/aggregate_cloud.F90 .
./phys_new.f90: $(SRCROOT)src/phys_new.f90
	cp $(SRCROOT)src/phys_new.f90 .
./compress_spectrum.F90: $(SRCROOT)src/socrates/interface/compress_spectrum.F90
	cp $(SRCROOT)src/socrates/interface/compress_spectrum.F90 .
./inter_pt.F90: $(SRCROOT)src/socrates/src/radiance_core/inter_pt.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/inter_pt.F90 .
./solver_no_scat.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_no_scat.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_no_scat.F90 .
./calc_surf_rad.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_surf_rad.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_surf_rad.F90 .
./ross.f90: $(SRCROOT)src/ross.f90
	cp $(SRCROOT)src/ross.f90 .
./set_cld.F90: $(SRCROOT)src/socrates/interface/set_cld.F90
	cp $(SRCROOT)src/socrates/interface/set_cld.F90 .
./flux.F90: $(SRCROOT)src/flux.F90
	cp $(SRCROOT)src/flux.F90 .
./radiance_calc.F90: $(SRCROOT)src/socrates/src/radiance_core/radiance_calc.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/radiance_calc.F90 .
./rescale_tau_csr.F90: $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_csr.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_csr.F90 .
./trans_source_coeff.F90: $(SRCROOT)src/socrates/src/radiance_core/trans_source_coeff.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/trans_source_coeff.F90 .
./mcica_sample.F90: $(SRCROOT)src/socrates/src/radiance_core/mcica_sample.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/mcica_sample.F90 .
./gauss_angle.F90: $(SRCROOT)src/socrates/src/radiance_core/gauss_angle.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/gauss_angle.F90 .
./solve_band_without_gas.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_without_gas.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_without_gas.F90 .
./read_control.F90: $(SRCROOT)src/socrates/interface/read_control.F90
	cp $(SRCROOT)src/socrates/interface/read_control.F90 .
./check_phf_term.F90: $(SRCROOT)src/socrates/src/radiance_core/check_phf_term.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/check_phf_term.F90 .
./sph_solver.F90: $(SRCROOT)src/socrates/src/radiance_core/sph_solver.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/sph_solver.F90 .
./set_n_cloud_parameter.F90: $(SRCROOT)src/socrates/src/radiance_core/set_n_cloud_parameter.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_n_cloud_parameter.F90 .
./opt_prop_inhom_corr_cairns.F90: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_inhom_corr_cairns.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/opt_prop_inhom_corr_cairns.F90 .
./mixed_solar_source.F90: $(SRCROOT)src/socrates/src/radiance_core/mixed_solar_source.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/mixed_solar_source.F90 .
./legendre_mod.F90: $(SRCROOT)src/socrates/src/radiance_core/legendre_mod.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/legendre_mod.F90 .
./band_solver.F90: $(SRCROOT)src/socrates/src/radiance_core/band_solver.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/band_solver.F90 .
./tsrad_short_char.f90: $(SRCROOT)src/tsrad_short_char.f90
	cp $(SRCROOT)src/tsrad_short_char.f90 .
./rad_pcf.F90: $(SRCROOT)src/socrates/src/radiance_core/rad_pcf.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/rad_pcf.F90 .
./solve_band_random_overlap_resort_rebin.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap_resort_rebin.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap_resort_rebin.F90 .
./def_dimen.F90: $(SRCROOT)src/socrates/src/radiance_core/def_dimen.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_dimen.F90 .
./solve_band_k_eqv_scl.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv_scl.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv_scl.F90 .
./monochromatic_radiance.F90: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance.F90 .
./fsolve.f90: $(SRCROOT)src/fsolve.f90
	cp $(SRCROOT)src/fsolve.f90 .
./mix_column.F90: $(SRCROOT)src/socrates/src/radiance_core/mix_column.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/mix_column.F90 .
./calc_cg_coeff.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_cg_coeff.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_cg_coeff.F90 .
./two_stream.F90: $(SRCROOT)src/socrates/src/radiance_core/two_stream.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_stream.F90 .
./adjust.f90: $(SRCROOT)src/adjust.f90
	cp $(SRCROOT)src/adjust.f90 .
./two_coeff_fast_lw.F90: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_fast_lw.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_coeff_fast_lw.F90 .
./solver_mix_direct_hogan.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct_hogan.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct_hogan.F90 .
./short_char_ross.f90: $(SRCROOT)src/short_char_ross.f90
	cp $(SRCROOT)src/short_char_ross.f90 .
./moisture_mod.f90: $(SRCROOT)src/moisture_mod.f90
	cp $(SRCROOT)src/moisture_mod.f90 .
./adjust_ir_radiance.F90: $(SRCROOT)src/socrates/src/radiance_core/adjust_ir_radiance.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/adjust_ir_radiance.F90 .
./two_coeff_basic.F90: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_basic.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_coeff_basic.F90 .
./set_control.F90: $(SRCROOT)src/socrates/interface/set_control.F90
	cp $(SRCROOT)src/socrates/interface/set_control.F90 .
./def_bound.F90: $(SRCROOT)src/socrates/src/radiance_core/def_bound.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_bound.F90 .
./solver_triple_app_scat.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_triple_app_scat.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_triple_app_scat.F90 .
./cloud_maxcs_split.F90: $(SRCROOT)src/socrates/src/radiance_core/cloud_maxcs_split.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/cloud_maxcs_split.F90 .
./monochromatic_gas_flux.F90: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_gas_flux.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/monochromatic_gas_flux.F90 .
./inter_t_lookup.F90: $(SRCROOT)src/socrates/src/radiance_core/inter_t_lookup.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/inter_t_lookup.F90 .
./monochromatic_radiance_tseq.F90: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_tseq.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_tseq.F90 .
./scale_absorb.F90: $(SRCROOT)src/socrates/src/radiance_core/scale_absorb.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/scale_absorb.F90 .
./file_manager.F90: $(SRCROOT)src/socrates/src/modules_core/file_manager.F90
	cp $(SRCROOT)src/socrates/src/modules_core/file_manager.F90 .
./legendre_weight.F90: $(SRCROOT)src/socrates/src/radiance_core/legendre_weight.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/legendre_weight.F90 .
./spline_fit.F90: $(SRCROOT)src/socrates/src/radiance_core/spline_fit.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/spline_fit.F90 .
./solver_triple_hogan.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_triple_hogan.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_triple_hogan.F90 .
./yomhook.F90: $(SRCROOT)src/socrates/src/modules_core/yomhook.F90
	cp $(SRCROOT)src/socrates/src/modules_core/yomhook.F90 .
./def_cld.F90: $(SRCROOT)src/socrates/src/radiance_core/def_cld.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_cld.F90 .
./dimensions_spec_ucf.F90: $(SRCROOT)src/socrates/src/modules_core/dimensions_spec_ucf.F90
	cp $(SRCROOT)src/socrates/src/modules_core/dimensions_spec_ucf.F90 .
./set_rad_layer.F90: $(SRCROOT)src/socrates/src/radiance_core/set_rad_layer.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_rad_layer.F90 .
./set_dimen.F90: $(SRCROOT)src/socrates/interface/set_dimen.F90
	cp $(SRCROOT)src/socrates/interface/set_dimen.F90 .
./two_coeff.F90: $(SRCROOT)src/socrates/src/radiance_core/two_coeff.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_coeff.F90 .
./toon.F90: $(SRCROOT)src/toon.F90
	cp $(SRCROOT)src/toon.F90 .
./socrates_config_mod.f90: $(SRCROOT)src/socrates/interface/socrates_config_mod.f90
	cp $(SRCROOT)src/socrates/interface/socrates_config_mod.f90 .
./io.f90: $(SRCROOT)src/io.f90
	cp $(SRCROOT)src/io.f90 .
./opt_prop_aerosol.F90: $(SRCROOT)src/socrates/src/radiance_core/opt_prop_aerosol.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/opt_prop_aerosol.F90 .
./def_aer.F90: $(SRCROOT)src/socrates/src/radiance_core/def_aer.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_aer.F90 .
./triple_column.F90: $(SRCROOT)src/socrates/src/radiance_core/triple_column.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/triple_column.F90 .
./gaussian_weight_pcf.F90: $(SRCROOT)src/socrates/src/radiance_core/gaussian_weight_pcf.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/gaussian_weight_pcf.F90 .
./shell_sort.F90: $(SRCROOT)src/socrates/src/radiance_core/shell_sort.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/shell_sort.F90 .
./diffusivity_factor.F90: $(SRCROOT)src/socrates/src/radiance_core/diffusivity_factor.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/diffusivity_factor.F90 .
./copy_clr_sol.F90: $(SRCROOT)src/socrates/src/radiance_core/copy_clr_sol.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/copy_clr_sol.F90 .
./sph_matrix_solver.F90: $(SRCROOT)src/socrates/src/radiance_core/sph_matrix_solver.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/sph_matrix_solver.F90 .
./circumsolar_fraction.F90: $(SRCROOT)src/socrates/src/radiance_core/circumsolar_fraction.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/circumsolar_fraction.F90 .
./gas_list_pcf.F90: $(SRCROOT)src/socrates/src/radiance_core/gas_list_pcf.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/gas_list_pcf.F90 .
./copy_clr_full.F90: $(SRCROOT)src/socrates/src/radiance_core/copy_clr_full.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/copy_clr_full.F90 .
./ir_source.F90: $(SRCROOT)src/socrates/src/radiance_core/ir_source.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/ir_source.F90 .
./set_n_source_coeff.F90: $(SRCROOT)src/socrates/src/radiance_core/set_n_source_coeff.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_n_source_coeff.F90 .
./two_coeff_cloud.F90: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_cloud.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_coeff_cloud.F90 .
./tau_mod.f90: $(SRCROOT)src/tau_mod.f90
	cp $(SRCROOT)src/tau_mod.f90 .
./def_out.F90: $(SRCROOT)src/socrates/src/radiance_core/def_out.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_out.F90 .
./inter_pt_lookup.F90: $(SRCROOT)src/socrates/src/radiance_core/inter_pt_lookup.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/inter_pt_lookup.F90 .
./layer_part_integ.F90: $(SRCROOT)src/socrates/src/radiance_core/layer_part_integ.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/layer_part_integ.F90 .
./set_truncation.F90: $(SRCROOT)src/socrates/src/radiance_core/set_truncation.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_truncation.F90 .
./filenamelength_mod.F90: $(SRCROOT)src/socrates/src/modules_core/filenamelength_mod.F90
	cp $(SRCROOT)src/socrates/src/modules_core/filenamelength_mod.F90 .
./eig_sys.F90: $(SRCROOT)src/socrates/src/radiance_core/eig_sys.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/eig_sys.F90 .
./set_aer.F90: $(SRCROOT)src/socrates/interface/set_aer.F90
	cp $(SRCROOT)src/socrates/interface/set_aer.F90 .
./ses_rescale_contm.F90: $(SRCROOT)src/socrates/src/radiance_core/ses_rescale_contm.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/ses_rescale_contm.F90 .
./set_cloud_pointer.F90: $(SRCROOT)src/socrates/src/radiance_core/set_cloud_pointer.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_cloud_pointer.F90 .
./def_spectrum.F90: $(SRCROOT)src/socrates/src/radiance_core/def_spectrum.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_spectrum.F90 .
./soc_constants.F90: $(SRCROOT)src/socrates/interface/soc_constants.F90
	cp $(SRCROOT)src/socrates/interface/soc_constants.F90 .
./new_socrates_interface.f90: $(SRCROOT)src/socrates/interface/new_socrates_interface.f90
	cp $(SRCROOT)src/socrates/interface/new_socrates_interface.f90 .
./calc_top_rad.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_top_rad.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_top_rad.F90 .
./augment_tiled_radiance.F90: $(SRCROOT)src/socrates/src/radiance_core/augment_tiled_radiance.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/augment_tiled_radiance.F90 .
./set_dirn_weights.F90: $(SRCROOT)src/socrates/src/radiance_core/set_dirn_weights.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_dirn_weights.F90 .
./monochromatic_ir_radiance.F90: $(SRCROOT)src/socrates/src/radiance_core/monochromatic_ir_radiance.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/monochromatic_ir_radiance.F90 .
./set_cloud_geometry.F90: $(SRCROOT)src/socrates/src/radiance_core/set_cloud_geometry.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_cloud_geometry.F90 .
./calc_brdf.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_brdf.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_brdf.F90 .
./params.f90: $(SRCROOT)src/params.f90
	cp $(SRCROOT)src/params.f90 .
./def_atm.F90: $(SRCROOT)src/socrates/src/radiance_core/def_atm.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_atm.F90 .
./parkind1.F90: $(SRCROOT)src/socrates/src/modules_core/parkind1.F90
	cp $(SRCROOT)src/socrates/src/modules_core/parkind1.F90 .
./rescale_continuum.F90: $(SRCROOT)src/socrates/src/radiance_core/rescale_continuum.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/rescale_continuum.F90 .
./phys.f90: $(SRCROOT)src/phys.f90
	cp $(SRCROOT)src/phys.f90 .
./def_ss_prop.F90: $(SRCROOT)src/socrates/src/radiance_core/def_ss_prop.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_ss_prop.F90 .
./solar_source.F90: $(SRCROOT)src/socrates/src/radiance_core/solar_source.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solar_source.F90 .
./two_coeff_region.F90: $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region.F90 .
./def_planck.F90: $(SRCROOT)src/socrates/src/radiance_core/def_planck.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_planck.F90 .
./increment_rad_cf.F90: $(SRCROOT)src/socrates/src/radiance_core/increment_rad_cf.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/increment_rad_cf.F90 .
./missing_data_mod.F90: $(SRCROOT)src/socrates/src/modules_core/missing_data_mod.F90
	cp $(SRCROOT)src/socrates/src/modules_core/missing_data_mod.F90 .
./set_bound.F90: $(SRCROOT)src/socrates/interface/set_bound.F90
	cp $(SRCROOT)src/socrates/interface/set_bound.F90 .
./scale_wenyi.F90: $(SRCROOT)src/socrates/src/radiance_core/scale_wenyi.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/scale_wenyi.F90 .
./spline_evaluate.F90: $(SRCROOT)src/socrates/src/radiance_core/spline_evaluate.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/spline_evaluate.F90 .
./errormessagelength_mod.F90: $(SRCROOT)src/socrates/src/modules_core/errormessagelength_mod.F90
	cp $(SRCROOT)src/socrates/src/modules_core/errormessagelength_mod.F90 .
./convection.f90: $(SRCROOT)src/convection.f90
	cp $(SRCROOT)src/convection.f90 .
./interp1d.F90: $(SRCROOT)src/socrates/src/radiance_core/interp1d.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/interp1d.F90 .
./twostr_module.f90: $(SRCROOT)src/twostr_module.f90
	cp $(SRCROOT)src/twostr_module.f90 .
./band_grey_rad.f90: $(SRCROOT)src/band_grey_rad.f90
	cp $(SRCROOT)src/band_grey_rad.f90 .
./realtype_rd.f90: $(SRCROOT)src/socrates/src/modules_core/realtype_rd.f90
	cp $(SRCROOT)src/socrates/src/modules_core/realtype_rd.f90 .
./calc_uplm_sol.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_sol.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_sol.F90 .
./ereport_mod.F90: $(SRCROOT)src/socrates/src/modules_core/ereport_mod.F90
	cp $(SRCROOT)src/socrates/src/modules_core/ereport_mod.F90 .
./vectlib_mod.F90: $(SRCROOT)src/socrates/src/modules_core/vectlib_mod.F90
	cp $(SRCROOT)src/socrates/src/modules_core/vectlib_mod.F90 .
./calc_radiance_ipa.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_radiance_ipa.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_radiance_ipa.F90 .
./solver_triple.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_triple.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_triple.F90 .
./single_scat_sol.F90: $(SRCROOT)src/socrates/src/radiance_core/single_scat_sol.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/single_scat_sol.F90 .
./gas_optical_properties.F90: $(SRCROOT)src/socrates/src/radiance_core/gas_optical_properties.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/gas_optical_properties.F90 .
./spherical_solar_source.F90: $(SRCROOT)src/socrates/src/radiance_core/spherical_solar_source.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/spherical_solar_source.F90 .
./sol_scat_cos.F90: $(SRCROOT)src/socrates/src/radiance_core/sol_scat_cos.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/sol_scat_cos.F90 .
./calc_uplm_zero.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_zero.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_zero.F90 .
./new_socrates_calc.f90: $(SRCROOT)src/socrates/interface/new_socrates_calc.f90
	cp $(SRCROOT)src/socrates/interface/new_socrates_calc.f90 .
./planck_flux_band_mod.F90: $(SRCROOT)src/socrates/src/radiance_core/planck_flux_band_mod.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/planck_flux_band_mod.F90 .
./rebin_esft_terms.F90: $(SRCROOT)src/socrates/src/radiance_core/rebin_esft_terms.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/rebin_esft_terms.F90 .
./solve_band_k_eqv.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv.F90 .
./solve_band_ses.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_ses.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_ses.F90 .
./nc_wrap.f90: $(SRCROOT)src/nc_wrap.f90
	cp $(SRCROOT)src/nc_wrap.f90 .
./rescale_tau_omega.F90: $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_omega.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_omega.F90 .
./timestep.f90: $(SRCROOT)src/timestep.f90
	cp $(SRCROOT)src/timestep.f90 .
./column_solver.F90: $(SRCROOT)src/socrates/src/radiance_core/column_solver.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/column_solver.F90 .
./diff_planck_source_mod.F90: $(SRCROOT)src/socrates/src/radiance_core/diff_planck_source_mod.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/diff_planck_source_mod.F90 .
./def_mcica.F90: $(SRCROOT)src/socrates/src/radiance_core/def_mcica.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_mcica.F90 .
./main.F90: $(SRCROOT)src/main.F90
	cp $(SRCROOT)src/main.F90 .
./cg_kappa_ms.F90: $(SRCROOT)src/socrates/src/radiance_core/cg_kappa_ms.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/cg_kappa_ms.F90 .
./utils.f90: $(SRCROOT)src/utils.f90
	cp $(SRCROOT)src/utils.f90 .
./rad_ccf.F90: $(SRCROOT)src/socrates/src/modules_core/rad_ccf.F90
	cp $(SRCROOT)src/socrates/src/modules_core/rad_ccf.F90 .
./grey_opt_prop.F90: $(SRCROOT)src/socrates/src/radiance_core/grey_opt_prop.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/grey_opt_prop.F90 .
./sum_k.F90: $(SRCROOT)src/socrates/src/radiance_core/sum_k.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/sum_k.F90 .
./mix_app_scat.F90: $(SRCROOT)src/socrates/src/radiance_core/mix_app_scat.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/mix_app_scat.F90 .
./set_atm.F90: $(SRCROOT)src/socrates/interface/set_atm.F90
	cp $(SRCROOT)src/socrates/interface/set_atm.F90 .
./build_sph_matrix.F90: $(SRCROOT)src/socrates/src/radiance_core/build_sph_matrix.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/build_sph_matrix.F90 .
./triple_solar_source.F90: $(SRCROOT)src/socrates/src/radiance_core/triple_solar_source.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/triple_solar_source.F90 .
./calc_flux_ipa.F90: $(SRCROOT)src/socrates/src/radiance_core/calc_flux_ipa.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/calc_flux_ipa.F90 .
./prsc_opt_prop.F90: $(SRCROOT)src/socrates/src/radiance_core/prsc_opt_prop.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/prsc_opt_prop.F90 .
./set_matrix_pentadiagonal.F90: $(SRCROOT)src/socrates/src/radiance_core/set_matrix_pentadiagonal.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/set_matrix_pentadiagonal.F90 .
./read_spectrum.F90: $(SRCROOT)src/socrates/src/radiance_core/read_spectrum.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/read_spectrum.F90 .
./quicksort.F90: $(SRCROOT)src/socrates/src/radiance_core/quicksort.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/quicksort.F90 .
./overlap_coupled.F90: $(SRCROOT)src/socrates/src/radiance_core/overlap_coupled.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/overlap_coupled.F90 .
./solver_mix_direct.F90: $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct.F90 .
./def_spherical_geometry.F90: $(SRCROOT)src/socrates/src/radiance_core/def_spherical_geometry.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_spherical_geometry.F90 .
./augment_radiance.F90: $(SRCROOT)src/socrates/src/radiance_core/augment_radiance.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/augment_radiance.F90 .
./inter_k.F90: $(SRCROOT)src/socrates/src/radiance_core/inter_k.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/inter_k.F90 .
./solve_band_random_overlap.F90: $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap.F90 .
./def_control.F90: $(SRCROOT)src/socrates/src/radiance_core/def_control.F90
	cp $(SRCROOT)src/socrates/src/radiance_core/def_control.F90 .
SRC = $(SRCROOT)src/socrates/src/radiance_core/calc_surf_rad.F90 $(SRCROOT)src/params.f90 $(SRCROOT)src/socrates/src/radiance_core/scale_wenyi.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_basic.F90 $(SRCROOT)src/socrates/src/radiance_core/def_out.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_radiance_ipa.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_ses.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_phase_fnc.F90 $(SRCROOT)src/socrates/src/radiance_core/sph_solver.F90 $(SRCROOT)src/rad_ross.f90 $(SRCROOT)src/socrates/src/radiance_core/calc_contrib_func.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_triple.F90 $(SRCROOT)src/socrates/src/radiance_core/def_ss_prop.F90 $(SRCROOT)src/socrates/src/radiance_core/single_scattering.F90 $(SRCROOT)src/socrates/src/modules_core/missing_data_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/def_spherical_geometry.F90 $(SRCROOT)src/socrates/src/radiance_core/sol_scat_cos.F90 $(SRCROOT)src/socrates/src/radiance_core/layer_part_integ.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_omega.F90 $(SRCROOT)src/utils.f90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff.F90 $(SRCROOT)src/twostr_module.f90 $(SRCROOT)src/socrates/src/radiance_core/set_dirn_weights.F90 $(SRCROOT)src/socrates/src/radiance_core/legendre_mod.F90 $(SRCROOT)src/socrates/src/modules_core/yomhook.F90 $(SRCROOT)src/socrates/src/modules_core/rad_ccf.F90 $(SRCROOT)src/socrates/src/radiance_core/gas_list_pcf.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_fast_lw.F90 $(SRCROOT)src/socrates/src/radiance_core/read_spectrum.F90 $(SRCROOT)src/socrates/src/radiance_core/diffusivity_factor.F90 $(SRCROOT)src/short_char_split_layer.f90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_water_cloud.F90 $(SRCROOT)src/socrates/src/radiance_core/spherical_solar_source.F90 $(SRCROOT)src/socrates/src/radiance_core/solar_coefficient_basic.F90 $(SRCROOT)src/socrates/src/radiance_core/set_matrix_pentadiagonal.F90 $(SRCROOT)src/socrates/src/radiance_core/triple_solar_source.F90 $(SRCROOT)src/phys.f90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_cloud.F90 $(SRCROOT)src/socrates/src/radiance_core/copy_clr_full.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_t_lookup.F90 $(SRCROOT)src/socrates/interface/new_socrates_interface.f90 $(SRCROOT)src/socrates/src/radiance_core/sph_matrix_solver.F90 $(SRCROOT)src/socrates/src/radiance_core/mcica_sample.F90 $(SRCROOT)src/socrates/src/radiance_core/rebin_esft_terms.F90 $(SRCROOT)src/flux.F90 $(SRCROOT)src/socrates/src/radiance_core/cloud_maxcs_split.F90 $(SRCROOT)src/socrates/src/modules_core/dimensions_spec_ucf.F90 $(SRCROOT)src/socrates/src/radiance_core/column_solver.F90 $(SRCROOT)src/socrates/interface/soc_constants.F90 $(SRCROOT)src/socrates/src/radiance_core/set_cloud_geometry.F90 $(SRCROOT)src/matrix.f90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance.F90 $(SRCROOT)src/io.f90 $(SRCROOT)src/socrates/src/radiance_core/diff_planck_source_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/augment_radiance.F90 $(SRCROOT)src/socrates/interface/compress_spectrum.F90 $(SRCROOT)src/socrates/interface/set_control.F90 $(SRCROOT)src/nc_wrap.f90 $(SRCROOT)src/main.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_continuum.F90 $(SRCROOT)src/tau_mod.f90 $(SRCROOT)src/socrates/src/radiance_core/def_spectrum.F90 $(SRCROOT)src/socrates/src/radiance_core/eval_uplm.F90 $(SRCROOT)src/socrates/src/radiance_core/scale_absorb.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_gas_flux.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_triple_hogan.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_brdf.F90 $(SRCROOT)src/socrates/src/radiance_core/prsc_gather_spline.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_ir_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/gaussian_weight_pcf.F90 $(SRCROOT)src/socrates/src/modules_core/realtype_rd.f90 $(SRCROOT)src/socrates/src/radiance_core/shell_sort.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_gauss_weight_90.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_sph.F90 $(SRCROOT)src/socrates/interface/new_socrates_calc.f90 $(SRCROOT)src/socrates/src/radiance_core/radiance_calc.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_cg_coeff.F90 $(SRCROOT)src/short_char_ross.f90 $(SRCROOT)src/socrates/src/radiance_core/triple_column.F90 $(SRCROOT)src/socrates/src/radiance_core/augment_tiled_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/sum_k.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_k.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_pt.F90 $(SRCROOT)src/socrates/src/modules_core/filenamelength_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/build_sph_matrix.F90 $(SRCROOT)src/socrates/src/radiance_core/def_planck.F90 $(SRCROOT)src/socrates/src/radiance_core/legendre_weight.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region.F90 $(SRCROOT)src/moisture_mod.f90 $(SRCROOT)src/socrates/src/radiance_core/calc_flux_ipa.F90 $(SRCROOT)src/socrates/interface/read_control.F90 $(SRCROOT)src/socrates/src/radiance_core/gas_optical_properties.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_aerosol.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_without_gas.F90 $(SRCROOT)src/socrates/src/radiance_core/trans_source_coeff.F90 $(SRCROOT)src/socrates/interface/set_bound.F90 $(SRCROOT)src/socrates/src/radiance_core/def_control.F90 $(SRCROOT)src/socrates/src/radiance_core/mixed_solar_source.F90 $(SRCROOT)src/socrates/interface/set_aer.F90 $(SRCROOT)src/socrates/src/modules_core/errormessagelength_mod.F90 $(SRCROOT)src/socrates/src/modules_core/file_manager.F90 $(SRCROOT)src/socrates/src/radiance_core/increment_rad_cf.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_tseq.F90 $(SRCROOT)src/socrates/src/radiance_core/set_truncation.F90 $(SRCROOT)src/socrates/src/radiance_core/two_stream.F90 $(SRCROOT)src/socrates/src/radiance_core/spherical_trans_coeff.F90 $(SRCROOT)src/toon.F90 $(SRCROOT)src/socrates/src/radiance_core/spline_fit.F90 $(SRCROOT)src/socrates/interface/set_cld.F90 $(SRCROOT)src/socrates/src/radiance_core/planck_flux_band_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/def_cld.F90 $(SRCROOT)src/socrates/src/radiance_core/interp1d.F90 $(SRCROOT)src/socrates/interface/set_atm.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv_scl.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_csr.F90 $(SRCROOT)src/socrates/src/radiance_core/aggregate_cloud.F90 $(SRCROOT)src/socrates/src/radiance_core/grey_opt_prop.F90 $(SRCROOT)src/socrates/interface/set_dimen.F90 $(SRCROOT)src/condense.f90 $(SRCROOT)src/socrates/src/radiance_core/set_rad_layer.F90 $(SRCROOT)src/socrates/src/radiance_core/ir_source.F90 $(SRCROOT)src/socrates/src/modules_core/ereport_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/eigenvalue_tri.F90 $(SRCROOT)src/socrates/src/radiance_core/def_aer.F90 $(SRCROOT)src/socrates/src/radiance_core/single_scattering_all.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_one_gas.F90 $(SRCROOT)src/socrates/src/radiance_core/quicksort.F90 $(SRCROOT)src/socrates/src/radiance_core/solar_source.F90 $(SRCROOT)src/socrates/src/radiance_core/diff_albedo_basis.F90 $(SRCROOT)src/socrates/src/radiance_core/band_solver.F90 $(SRCROOT)src/socrates/src/radiance_core/set_level_weights.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_pt_lookup.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_sol.F90 $(SRCROOT)src/fsolve.f90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv.F90 $(SRCROOT)src/socrates/src/radiance_core/single_scat_sol.F90 $(SRCROOT)src/socrates/src/radiance_core/check_phf_term.F90 $(SRCROOT)src/socrates/src/radiance_core/def_bound.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct_hogan.F90 $(SRCROOT)src/socrates/src/radiance_core/def_atm.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ice_cloud.F90 $(SRCROOT)src/socrates/src/radiance_core/rad_pcf.F90 $(SRCROOT)src/socrates/src/radiance_core/adjust_ir_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/cg_kappa_ms.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct.F90 $(SRCROOT)src/socrates/src/radiance_core/spline_evaluate.F90 $(SRCROOT)src/socrates/src/radiance_core/spherical_path.F90 $(SRCROOT)src/socrates/src/modules_core/vectlib_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/copy_clr_sol.F90 $(SRCROOT)src/socrates/src/radiance_core/set_cloud_pointer.F90 $(SRCROOT)src/socrates/src/modules_core/parkind1.F90 $(SRCROOT)src/socrates/src/radiance_core/overlap_coupled.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_inhom_corr_cairns.F90 $(SRCROOT)src/socrates/src/radiance_core/def_dimen.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap_resort_rebin.F90 $(SRCROOT)src/convection.f90 $(SRCROOT)src/socrates/src/radiance_core/set_n_cloud_parameter.F90 $(SRCROOT)src/socrates/src/radiance_core/def_mcica.F90 $(SRCROOT)src/socrates/src/radiance_core/circumsolar_fraction.F90 $(SRCROOT)src/band_grey_rad.f90 $(SRCROOT)src/timestep.f90 $(SRCROOT)src/socrates/src/radiance_core/solver_no_scat.F90 $(SRCROOT)src/socrates/src/radiance_core/gauss_angle.F90 $(SRCROOT)src/ross.f90 $(SRCROOT)src/socrates/src/radiance_core/mcica_column.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ukca_aerosol.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region_fast_lw.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_top_rad.F90 $(SRCROOT)src/tsrad_short_char.f90 $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_zero.F90 $(SRCROOT)src/socrates/src/radiance_core/mix_column.F90 $(SRCROOT)src/phys_new.f90 $(SRCROOT)src/socrates/src/radiance_core/eig_sys.F90 $(SRCROOT)src/socrates/src/radiance_core/prsc_opt_prop.F90 $(SRCROOT)src/socrates/src/radiance_core/set_n_source_coeff.F90 $(SRCROOT)src/socrates/src/radiance_core/ses_rescale_contm.F90 $(SRCROOT)src/socrates/src/radiance_core/mix_app_scat.F90 $(SRCROOT)src/socrates/src/radiance_core/hemi_sph_integ.F90 $(SRCROOT)src/adjust.f90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_homogen_direct.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_triple_app_scat.F90 $(SRCROOT)src/socrates/interface/socrates_config_mod.f90
OBJ = calc_surf_rad.o params.o scale_wenyi.o two_coeff_basic.o def_out.o calc_radiance_ipa.o solve_band_ses.o rescale_phase_fnc.o sph_solver.o rad_ross.o calc_contrib_func.o solver_triple.o def_ss_prop.o single_scattering.o missing_data_mod.o def_spherical_geometry.o sol_scat_cos.o layer_part_integ.o rescale_tau_omega.o utils.o two_coeff.o twostr_module.o set_dirn_weights.o legendre_mod.o yomhook.o rad_ccf.o gas_list_pcf.o two_coeff_fast_lw.o read_spectrum.o diffusivity_factor.o short_char_split_layer.o opt_prop_water_cloud.o spherical_solar_source.o solar_coefficient_basic.o set_matrix_pentadiagonal.o triple_solar_source.o phys.o two_coeff_cloud.o copy_clr_full.o inter_t_lookup.o new_socrates_interface.o sph_matrix_solver.o mcica_sample.o rebin_esft_terms.o flux.o cloud_maxcs_split.o dimensions_spec_ucf.o column_solver.o soc_constants.o set_cloud_geometry.o matrix.o monochromatic_radiance.o io.o diff_planck_source_mod.o augment_radiance.o compress_spectrum.o set_control.o nc_wrap.o main.o rescale_continuum.o tau_mod.o def_spectrum.o eval_uplm.o scale_absorb.o monochromatic_gas_flux.o solver_triple_hogan.o calc_brdf.o prsc_gather_spline.o monochromatic_ir_radiance.o gaussian_weight_pcf.o realtype_rd.o shell_sort.o calc_gauss_weight_90.o monochromatic_radiance_sph.o new_socrates_calc.o radiance_calc.o calc_cg_coeff.o short_char_ross.o triple_column.o augment_tiled_radiance.o sum_k.o inter_k.o inter_pt.o filenamelength_mod.o build_sph_matrix.o def_planck.o legendre_weight.o two_coeff_region.o moisture_mod.o calc_flux_ipa.o read_control.o gas_optical_properties.o opt_prop_aerosol.o solve_band_without_gas.o trans_source_coeff.o set_bound.o def_control.o mixed_solar_source.o set_aer.o errormessagelength_mod.o file_manager.o increment_rad_cf.o monochromatic_radiance_tseq.o set_truncation.o two_stream.o spherical_trans_coeff.o toon.o spline_fit.o set_cld.o planck_flux_band_mod.o def_cld.o interp1d.o set_atm.o solve_band_k_eqv_scl.o rescale_tau_csr.o aggregate_cloud.o grey_opt_prop.o set_dimen.o condense.o set_rad_layer.o ir_source.o ereport_mod.o eigenvalue_tri.o def_aer.o single_scattering_all.o solve_band_one_gas.o quicksort.o solar_source.o diff_albedo_basis.o band_solver.o set_level_weights.o inter_pt_lookup.o calc_uplm_sol.o fsolve.o solve_band_k_eqv.o single_scat_sol.o check_phf_term.o def_bound.o solver_mix_direct_hogan.o def_atm.o opt_prop_ice_cloud.o rad_pcf.o adjust_ir_radiance.o cg_kappa_ms.o solver_mix_direct.o spline_evaluate.o spherical_path.o vectlib_mod.o copy_clr_sol.o set_cloud_pointer.o parkind1.o overlap_coupled.o opt_prop_inhom_corr_cairns.o def_dimen.o solve_band_random_overlap_resort_rebin.o convection.o set_n_cloud_parameter.o def_mcica.o circumsolar_fraction.o band_grey_rad.o timestep.o solver_no_scat.o gauss_angle.o ross.o mcica_column.o opt_prop_ukca_aerosol.o two_coeff_region_fast_lw.o calc_top_rad.o tsrad_short_char.o calc_uplm_zero.o mix_column.o phys_new.o eig_sys.o prsc_opt_prop.o set_n_source_coeff.o ses_rescale_contm.o mix_app_scat.o hemi_sph_integ.o adjust.o solve_band_random_overlap.o solver_homogen_direct.o solver_triple_app_scat.o socrates_config_mod.o
OFF = $(SRCROOT)src/socrates/src/radiance_core/spherical_path.F90 $(SRCROOT)src/socrates/src/radiance_core/spherical_trans_coeff.F90 $(SRCROOT)src/socrates/src/radiance_core/mcica_column.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_sph.F90 $(SRCROOT)src/socrates/src/radiance_core/single_scattering.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_phase_fnc.F90 $(SRCROOT)src/matrix.f90 $(SRCROOT)src/socrates/src/radiance_core/eigenvalue_tri.F90 $(SRCROOT)src/socrates/src/radiance_core/hemi_sph_integ.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_contrib_func.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ice_cloud.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_ukca_aerosol.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_gauss_weight_90.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_one_gas.F90 $(SRCROOT)src/condense.f90 $(SRCROOT)src/socrates/src/radiance_core/prsc_gather_spline.F90 $(SRCROOT)src/socrates/src/radiance_core/single_scattering_all.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_water_cloud.F90 $(SRCROOT)src/socrates/src/radiance_core/diff_albedo_basis.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_homogen_direct.F90 $(SRCROOT)src/short_char_split_layer.f90 $(SRCROOT)src/socrates/src/radiance_core/eval_uplm.F90 $(SRCROOT)src/socrates/src/radiance_core/solar_coefficient_basic.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region_fast_lw.F90 $(SRCROOT)src/rad_ross.f90 $(SRCROOT)src/socrates/src/radiance_core/set_level_weights.F90 $(SRCROOT)src/socrates/src/radiance_core/aggregate_cloud.F90 $(SRCROOT)src/phys_new.f90 $(SRCROOT)src/socrates/interface/compress_spectrum.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_pt.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_no_scat.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_surf_rad.F90 $(SRCROOT)src/ross.f90 $(SRCROOT)src/socrates/interface/set_cld.F90 $(SRCROOT)src/flux.F90 $(SRCROOT)src/socrates/src/radiance_core/radiance_calc.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_csr.F90 $(SRCROOT)src/socrates/src/radiance_core/trans_source_coeff.F90 $(SRCROOT)src/socrates/src/radiance_core/mcica_sample.F90 $(SRCROOT)src/socrates/src/radiance_core/gauss_angle.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_without_gas.F90 $(SRCROOT)src/socrates/interface/read_control.F90 $(SRCROOT)src/socrates/src/radiance_core/check_phf_term.F90 $(SRCROOT)src/socrates/src/radiance_core/sph_solver.F90 $(SRCROOT)src/socrates/src/radiance_core/set_n_cloud_parameter.F90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_inhom_corr_cairns.F90 $(SRCROOT)src/socrates/src/radiance_core/mixed_solar_source.F90 $(SRCROOT)src/socrates/src/radiance_core/legendre_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/band_solver.F90 $(SRCROOT)src/tsrad_short_char.f90 $(SRCROOT)src/socrates/src/radiance_core/rad_pcf.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap_resort_rebin.F90 $(SRCROOT)src/socrates/src/radiance_core/def_dimen.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv_scl.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance.F90 $(SRCROOT)src/fsolve.f90 $(SRCROOT)src/socrates/src/radiance_core/mix_column.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_cg_coeff.F90 $(SRCROOT)src/socrates/src/radiance_core/two_stream.F90 $(SRCROOT)src/adjust.f90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_fast_lw.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct_hogan.F90 $(SRCROOT)src/short_char_ross.f90 $(SRCROOT)src/moisture_mod.f90 $(SRCROOT)src/socrates/src/radiance_core/adjust_ir_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_basic.F90 $(SRCROOT)src/socrates/interface/set_control.F90 $(SRCROOT)src/socrates/src/radiance_core/def_bound.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_triple_app_scat.F90 $(SRCROOT)src/socrates/src/radiance_core/cloud_maxcs_split.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_gas_flux.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_t_lookup.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_radiance_tseq.F90 $(SRCROOT)src/socrates/src/radiance_core/scale_absorb.F90 $(SRCROOT)src/socrates/src/modules_core/file_manager.F90 $(SRCROOT)src/socrates/src/radiance_core/legendre_weight.F90 $(SRCROOT)src/socrates/src/radiance_core/spline_fit.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_triple_hogan.F90 $(SRCROOT)src/socrates/src/modules_core/yomhook.F90 $(SRCROOT)src/socrates/src/radiance_core/def_cld.F90 $(SRCROOT)src/socrates/src/modules_core/dimensions_spec_ucf.F90 $(SRCROOT)src/socrates/src/radiance_core/set_rad_layer.F90 $(SRCROOT)src/socrates/interface/set_dimen.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff.F90 $(SRCROOT)src/toon.F90 $(SRCROOT)src/socrates/interface/socrates_config_mod.f90 $(SRCROOT)src/io.f90 $(SRCROOT)src/socrates/src/radiance_core/opt_prop_aerosol.F90 $(SRCROOT)src/socrates/src/radiance_core/def_aer.F90 $(SRCROOT)src/socrates/src/radiance_core/triple_column.F90 $(SRCROOT)src/socrates/src/radiance_core/gaussian_weight_pcf.F90 $(SRCROOT)src/socrates/src/radiance_core/shell_sort.F90 $(SRCROOT)src/socrates/src/radiance_core/diffusivity_factor.F90 $(SRCROOT)src/socrates/src/radiance_core/copy_clr_sol.F90 $(SRCROOT)src/socrates/src/radiance_core/sph_matrix_solver.F90 $(SRCROOT)src/socrates/src/radiance_core/circumsolar_fraction.F90 $(SRCROOT)src/socrates/src/radiance_core/gas_list_pcf.F90 $(SRCROOT)src/socrates/src/radiance_core/copy_clr_full.F90 $(SRCROOT)src/socrates/src/radiance_core/ir_source.F90 $(SRCROOT)src/socrates/src/radiance_core/set_n_source_coeff.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_cloud.F90 $(SRCROOT)src/tau_mod.f90 $(SRCROOT)src/socrates/src/radiance_core/def_out.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_pt_lookup.F90 $(SRCROOT)src/socrates/src/radiance_core/layer_part_integ.F90 $(SRCROOT)src/socrates/src/radiance_core/set_truncation.F90 $(SRCROOT)src/socrates/src/modules_core/filenamelength_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/eig_sys.F90 $(SRCROOT)src/socrates/interface/set_aer.F90 $(SRCROOT)src/socrates/src/radiance_core/ses_rescale_contm.F90 $(SRCROOT)src/socrates/src/radiance_core/set_cloud_pointer.F90 $(SRCROOT)src/socrates/src/radiance_core/def_spectrum.F90 $(SRCROOT)src/socrates/interface/soc_constants.F90 $(SRCROOT)src/socrates/interface/new_socrates_interface.f90 $(SRCROOT)src/socrates/src/radiance_core/calc_top_rad.F90 $(SRCROOT)src/socrates/src/radiance_core/augment_tiled_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/set_dirn_weights.F90 $(SRCROOT)src/socrates/src/radiance_core/monochromatic_ir_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/set_cloud_geometry.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_brdf.F90 $(SRCROOT)src/params.f90 $(SRCROOT)src/socrates/src/radiance_core/def_atm.F90 $(SRCROOT)src/socrates/src/modules_core/parkind1.F90 $(SRCROOT)src/socrates/src/radiance_core/rescale_continuum.F90 $(SRCROOT)src/phys.f90 $(SRCROOT)src/socrates/src/radiance_core/def_ss_prop.F90 $(SRCROOT)src/socrates/src/radiance_core/solar_source.F90 $(SRCROOT)src/socrates/src/radiance_core/two_coeff_region.F90 $(SRCROOT)src/socrates/src/radiance_core/def_planck.F90 $(SRCROOT)src/socrates/src/radiance_core/increment_rad_cf.F90 $(SRCROOT)src/socrates/src/modules_core/missing_data_mod.F90 $(SRCROOT)src/socrates/interface/set_bound.F90 $(SRCROOT)src/socrates/src/radiance_core/scale_wenyi.F90 $(SRCROOT)src/socrates/src/radiance_core/spline_evaluate.F90 $(SRCROOT)src/socrates/src/modules_core/errormessagelength_mod.F90 $(SRCROOT)src/convection.f90 $(SRCROOT)src/socrates/src/radiance_core/interp1d.F90 $(SRCROOT)src/twostr_module.f90 $(SRCROOT)src/band_grey_rad.f90 $(SRCROOT)src/socrates/src/modules_core/realtype_rd.f90 $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_sol.F90 $(SRCROOT)src/socrates/src/modules_core/ereport_mod.F90 $(SRCROOT)src/socrates/src/modules_core/vectlib_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_radiance_ipa.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_triple.F90 $(SRCROOT)src/socrates/src/radiance_core/single_scat_sol.F90 $(SRCROOT)src/socrates/src/radiance_core/gas_optical_properties.F90 $(SRCROOT)src/socrates/src/radiance_core/spherical_solar_source.F90 $(SRCROOT)src/socrates/src/radiance_core/sol_scat_cos.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_uplm_zero.F90 $(SRCROOT)src/socrates/interface/new_socrates_calc.f90 $(SRCROOT)src/socrates/src/radiance_core/planck_flux_band_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/rebin_esft_terms.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_k_eqv.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_ses.F90 $(SRCROOT)src/nc_wrap.f90 $(SRCROOT)src/socrates/src/radiance_core/rescale_tau_omega.F90 $(SRCROOT)src/timestep.f90 $(SRCROOT)src/socrates/src/radiance_core/column_solver.F90 $(SRCROOT)src/socrates/src/radiance_core/diff_planck_source_mod.F90 $(SRCROOT)src/socrates/src/radiance_core/def_mcica.F90 $(SRCROOT)src/main.F90 $(SRCROOT)src/socrates/src/radiance_core/cg_kappa_ms.F90 $(SRCROOT)src/utils.f90 $(SRCROOT)src/socrates/src/modules_core/rad_ccf.F90 $(SRCROOT)src/socrates/src/radiance_core/grey_opt_prop.F90 $(SRCROOT)src/socrates/src/radiance_core/sum_k.F90 $(SRCROOT)src/socrates/src/radiance_core/mix_app_scat.F90 $(SRCROOT)src/socrates/interface/set_atm.F90 $(SRCROOT)src/socrates/src/radiance_core/build_sph_matrix.F90 $(SRCROOT)src/socrates/src/radiance_core/triple_solar_source.F90 $(SRCROOT)src/socrates/src/radiance_core/calc_flux_ipa.F90 $(SRCROOT)src/socrates/src/radiance_core/prsc_opt_prop.F90 $(SRCROOT)src/socrates/src/radiance_core/set_matrix_pentadiagonal.F90 $(SRCROOT)src/socrates/src/radiance_core/read_spectrum.F90 $(SRCROOT)src/socrates/src/radiance_core/quicksort.F90 $(SRCROOT)src/socrates/src/radiance_core/overlap_coupled.F90 $(SRCROOT)src/socrates/src/radiance_core/solver_mix_direct.F90 $(SRCROOT)src/socrates/src/radiance_core/def_spherical_geometry.F90 $(SRCROOT)src/socrates/src/radiance_core/augment_radiance.F90 $(SRCROOT)src/socrates/src/radiance_core/inter_k.F90 $(SRCROOT)src/socrates/src/radiance_core/solve_band_random_overlap.F90 $(SRCROOT)src/socrates/src/radiance_core/def_control.F90
clean: neat
	-rm -f .main.exe.cppdefs $(OBJ) main.exe
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
main.exe: $(OBJ) 
	$(LD) $(OBJ) -o main.exe  $(LDFLAGS)
