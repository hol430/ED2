#Makefile objects.mk

# Define main source.

MAIN = $(CORE)/rammain.F90
MAINOBJ = rammain.o


# Define objects.

OBJ_MODEL =                         \
	adap_init.o                 \
	altera_dia.o                \
	an_header.o                 \
	aobj.o                      \
	asgen.o                     \
	asnc.o                      \
	asti.o                      \
	asti2.o                     \
	astp.o                      \
	avarf.o                     \
	catt_start.o                \
	charutils.o                 \
	cond_read.o                 \
	cond_update.o               \
	conv_coms.o                 \
	coriolis.o                  \
	cu_read.o                   \
	cup_dn.o                    \
	cup_env.o                   \
	cup_grell2.o                \
	cup_grell2_shcu.o           \
	cup_up.o                    \
	cyclic_mod.o                \
	dateutils.o                 \
	dealloc.o                   \
	diffsclr.o                  \
	diffuse.o                   \
	domain_decomp.o             \
	dry_dep.o                   \
	dted.o                      \
	eenviron.o                  \
	emission_source_map.o       \
	error_mess.o                \
	extra.o                     \
	file_inv.o                  \
	filelist.o                  \
	first_rams.o                \
	gaspart.o                   \
	geodat.o                    \
	getvar.o                    \
	great_circle.o              \
	grell_coms.o                \
	grell_cupar_aux.o           \
	grell_cupar_downdraft.o     \
	grell_cupar_driver.o        \
	grell_cupar_dynamic.o       \
	grell_cupar_ensemble.o      \
	grell_cupar_environment.o   \
	grell_cupar_feedback.o      \
	grell_cupar_static.o        \
	grell_cupar_updraft.o       \
	grell_extras_catt.o         \
	grid_dims.o                 \
	grid_struct.o               \
	gridset.o                   \
	harr_coms.o                 \
	harr_rad.o                  \
	harr_raddriv.o              \
	harr_radinit.o              \
	hdf5_coms.o                 \
	hdf5_utils.o                \
	hemi2.o                     \
	htint-opt.o                 \
	inithis.o                   \
	interp_lib.o                \
	io_params.o                 \
	isan_coms.o                 \
	isan_io.o                   \
	ke_coms.o                   \
	kuo_cupar_driver.o          \
	landuse_input.o             \
	leaf_coms.o                 \
	leaf3.o                     \
	leaf3_bc.o                  \
	leaf3_can.o                 \
	leaf3_hyd.o                 \
	leaf3_init.o                \
	leaf3_ocean.o               \
	leaf3_teb.o                 \
	leaf3_tw.o                  \
	leaf3_utils.o               \
	local_proc.o                \
	machine_arq.o               \
	map_proj.o                  \
	mem_aerad.o                 \
	mem_all.o                   \
	mem_basic.o                 \
	mem_carma.o                 \
	mem_cuparm.o                \
	mem_emiss.o                 \
	mem_ensemble.o              \
	mem_gaspart.o               \
	mem_globaer.o               \
	mem_globrad.o               \
	mem_grell_param2.o          \
	mem_grid.o                  \
	mem_grid_dim_defs.o         \
	mem_harr.o                  \
	mem_leaf.o                  \
	mem_mass.o                  \
	mem_mclat.o                 \
	mem_micro.o                 \
	mem_mnt_advec.o             \
	mem_mksfc.o                 \
	mem_nestb.o                 \
	mem_oda.o                   \
	mem_opt_scratch.o           \
	mem_radiate.o               \
	mem_scalar.o                \
	mem_scratch.o               \
	mem_scratch_grell.o         \
	mem_scratch1_brams.o        \
	mem_scratch1_grell.o        \
	mem_scratch2_grell.o        \
	mem_scratch2_grell_sh.o     \
	mem_scratch3_grell.o        \
	mem_scratch3_grell_sh.o     \
	mem_soil_moisture.o         \
	mem_tconv.o                 \
	mem_teb.o                   \
	mem_teb_common.o            \
	mem_teb_vars_const.o        \
	mem_tend.o                  \
	mem_turb.o                  \
	mem_turb_scalar.o           \
	mem_varinit.o               \
	mic_coll.o                  \
	mic_driv.o                  \
	mic_gamma.o                 \
	mic_init.o                  \
	mic_misc.o                  \
	mic_nuc.o                   \
	mic_tabs.o                  \
	mic_vap.o                   \
	micphys.o                   \
	micro_coms.o                \
	mksfc_driver.o              \
	mksfc_fuso.o                \
	mksfc_ndvi.o                \
	mksfc_sfc.o                 \
	mksfc_sst.o                 \
	mksfc_top.o                 \
	mnt_advec_aux.o             \
	mnt_advec_main.o            \
	mod_advect_kit.o            \
	mod_GhostBlock.o            \
	mod_GhostBlockPartition.o   \
	mod_ozone.o                 \
	model.o                     \
	modsched.o                  \
	mpass_advec.o               \
	mpass_cyclic.o              \
	mpass_dtl.o                 \
	mpass_feed.o                \
	mpass_full.o                \
	mpass_init.o                \
	mpass_lbc.o                 \
	mpass_nest.o                \
	mpass_oda.o                 \
	mpass_st.o                  \
	ncarg_dummy.o               \
	ndvi_read.o                 \
	nest_drivers.o              \
	nest_feed.o                 \
	nest_filldens.o             \
	nest_geosst.o               \
	nest_init_aux.o             \
	nest_intrp.o                \
	nest_move.o                 \
	node_mod.o                  \
	nud_analysis.o              \
	nud_read.o                  \
	nud_update.o                \
	numutils.o                  \
	obs_input.o                 \
	oda_krig.o                  \
	oda_nudge.o                 \
	oda_proc_obs.o              \
	oda_read.o                  \
	oda_sta_count.o             \
	oda_sta_input.o             \
	old_grell_cupar_driver.o    \
	opspec.o                    \
	ozone.o                     \
	par_decomp.o                \
	para_init.o                 \
	paral.o                     \
	polarst.o                   \
	plumerise_vector.o          \
	raco.o                      \
	raco_adap.o                 \
	rad_carma.o                 \
	rad_ccmp.o                  \
	rad_driv.o                  \
	rad_mclat.o                 \
	rad_stable.o                \
	radvc.o                     \
	radvc_adap.o                \
	radvc_new.o                 \
	rams_grid.o                 \
	rams_master.o               \
	rams_mem_alloc.o            \
	rams_read_header.o          \
	ranlavg.o                   \
	rbnd.o                      \
	rbnd_adap.o                 \
	rcio.o                      \
	rconstants.o                \
	rconv_driver.o              \
	rdint.o                     \
	read_ralph.o                \
	recycle.o                   \
	ref_sounding.o              \
	refstate.o                  \
	rexev.o                     \
	rgrad.o                     \
	rhhi.o                      \
	rhdf5.o                     \
	rinit.o                     \
	rio.o                       \
	rmass.o                     \
	rname.o                     \
	rnest_par.o                 \
	rnode.o                     \
	rpara.o                     \
	rprnt.o                     \
	rthrm.o                     \
	rtimh.o                     \
	rtimi.o                     \
	rsys.o                      \
	ruser.o                     \
	shcu_vars_const.o           \
	soil_moisture_init.o        \
	souza_cupar_driver.o        \
	sst_read.o                  \
	teb_spm_start.o             \
	therm_lib.o                 \
	therm_lib8.o                \
	tkenn.o                     \
	tmpname.o                   \
	turb_coms.o                 \
	turb_derivs.o               \
	turb_diff.o                 \
	turb_diff_adap.o            \
	turb_k.o                    \
	turb_k_adap.o               \
	turb_ke.o                   \
	urban.o                     \
	urban_canopy.o              \
	utils_c.o                   \
	utils_f.o                   \
	v_interps.o                 \
	var_tables.o                \
	varf_read.o                 \
	varf_update.o               \
	varutils.o                  \
	vformat.o                   \
	vtab_fill.o                 \
	edcp_driver.o               \
	edcp_init.o                 \
	edcp_lake_driver.o          \
	edcp_lake_misc.o            \
	edcp_lake_stepper.o         \
	edcp_load_namelist.o        \
	edcp_met.o                  \
	edcp_met_init.o             \
	edcp_model.o                \
	edcp_mpiutils.o             \
	edcp_para_init.o            \
	mem_edcp.o                  \
	lake_coms.o                 \
	allometry.o                 \
	average_utils.o             \
	budget_utils.o              \
	c34constants.o              \
	canopy_air_coms.o           \
	canopy_layer_coms.o         \
	canopy_radiation_coms.o     \
	canopy_struct_dynamics.o    \
	consts_coms.o               \
	decomp_coms.o               \
	disturb_coms.o              \
	disturbance.o               \
	ed_bigleaf_init.o           \
	ed_filelist.o               \
	ed_grid.o                   \
	ed_init_full_history.o      \
	ed_init.o                   \
	ed_init_atm.o               \
	ed_max_dims.o               \
	ed_mem_grid_dim_defs.o      \
	ed_misc_coms.o              \
	ed_nbg_init.o               \
	ed_node_coms.o              \
	ed_opspec.o                 \
	ed_para_coms.o              \
	ed_params.o                 \
	ed_print.o                 \
	ed_read_ed10_20_history.o   \
	ed_read_ed21_history.o      \
	ed_state_vars.o             \
	ed_therm_lib.o              \
	ed_type_init.o              \
	ed_var_tables.o             \
	ed_work_vars.o              \
	ed_xml_config.o             \
	edio.o                      \
	ename_coms.o                \
	euler_driver.o              \
	events.o                    \
	farq_leuning.o              \
	fatal_error.o               \
	fire.o                      \
	forestry.o                  \
	fuse_fiss_utils.o           \
	fusion_fission_coms.o       \
	grid_coms.o                 \
	growth_balive.o             \
	h5_output.o                 \
	heun_driver.o               \
	hydrology_coms.o            \
	hydrology_constants.o       \
	init_hydro_sites.o          \
	invmondays.o                \
	landuse_init.o              \
	lapse.o                     \
	leaf_database.o             \
	libxml2f90.f90_pp.o         \
	lsm_hyd.o                   \
	mem_polygons.o              \
	met_driver_coms.o           \
	mortality.o                 \
	multiple_scatter.o          \
	old_twostream_rad.o         \
	optimiz_coms.o              \
	phenology_aux.o             \
	phenology_coms.o            \
	phenology_driv.o            \
	phenology_startup.o         \
	photosyn_driv.o             \
	physiology_coms.o           \
	pft_coms.o                  \
	radiate_driver.o            \
	radiate_utils.o             \
	random_utils.o              \
	reproduction.o              \
	rk4_coms.o                  \
	rk4_derivs.o                \
	rk4_driver.o                \
	rk4_integ_utils.o           \
	rk4_misc.o                  \
	rk4_stepper.o               \
	soil_coms.o                 \
	soil_respiration.o          \
	stable_cohorts.o            \
	structural_growth.o         \
	twostream_rad.o             \
	update_derived_props.o      \
	vegetation_dynamics.o       \
	detailed_coms.o             \
	bdf2_solver.o               \
	hybrid_driver.o

