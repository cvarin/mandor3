//configuration files reader

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"params.h"
#include"utils.h"
#include"io_utils.h"
#include"types.h"
#include"data_loader.h"
#include"report.h"
#include"cfgreader.h"

void read_cfg_common(char *bin_data_dir, int_array *checkpoints_list, int_array *snapshots_list, int *hor_pic_size, subdomain *subdmn, int *coarsen_data){
	//reads the parameters from the CFG_COMMON file and some data stored in the Mandor binData folder
	int start_chkp, n_skipped_chkp, start_snsht, n_skipped_snsht, n_saved_chkp, n_saved_snsht, i, n = 0, l;
	int enable_subdomain, imin, imax, jmin, jmax, kmin, kmax;

	cfg_read_str("bin_data_dir", bin_data_dir, CFG_COMMON);

	start_chkp = cfg_read_int("start_chkp", CFG_COMMON);
	n_skipped_chkp = cfg_read_int("n_skipped_chkp", CFG_COMMON);
	start_snsht = cfg_read_int("start_snsht", CFG_COMMON);
	n_skipped_snsht = cfg_read_int("n_skipped_snsht", CFG_COMMON);

	*hor_pic_size = cfg_read_int("hor_pic_size", CFG_COMMON);
	*coarsen_data = cfg_read_int("coarsen_data", CFG_COMMON);
	*coarsen_data = (*coarsen_data > 0)?1:0;

	enable_subdomain = cfg_read_int("enable_subdomain", CFG_COMMON);
	imin = cfg_read_int("imin", CFG_COMMON);
	imax = cfg_read_int("imax", CFG_COMMON);
	jmin = cfg_read_int("jmin", CFG_COMMON);
	jmax = cfg_read_int("jmax", CFG_COMMON);
	kmin = cfg_read_int("kmin", CFG_COMMON);
	kmax = cfg_read_int("kmax", CFG_COMMON);

	l = strlen(bin_data_dir);
	if(bin_data_dir[l-1] != '/'){
		bin_data_dir[l] = '/';
		bin_data_dir[l+1] = '\0';
	}
	if(*hor_pic_size < 2) error(FUNCTION, "hor_pic_size = %d < 2", hor_pic_size);
	if(start_chkp < 0) error(FUNCTION, "start_chkp = %d < 0. This is forbidden.", start_chkp);
	if(n_skipped_chkp < 0) error(FUNCTION, "n_skipped_chkp = %d < 0. This is forbidden.", n_skipped_chkp);
	if(start_snsht < 0) error(FUNCTION, "start_snsht = %d < 0. This is forbidden.", start_snsht);
	if(n_skipped_snsht < 0) error(FUNCTION, "n_skipped_snsht = %d < 0. This is forbidden.", n_skipped_snsht);

	enable_subdomain = (enable_subdomain > 0)?1:0;
	if(imin < 0) error(FUNCTION, "(imin = %d) < 0. This is forbidden.", imin);
	if(jmin < 0) error(FUNCTION, "(jmin = %d) < 0. This is forbidden.", jmin);
	if(kmin < 0) error(FUNCTION, "(kmin = %d) < 0. This is forbidden.", kmin);
	if(imax <= imin) error(FUNCTION, "(imax = %d) <= (imin = %d). This is forbidden.", imax, imin);
	if(jmax <= jmin) error(FUNCTION, "(jmax = %d) <= (jmin = %d). This is forbidden.", jmax, jmin);
	if(kmax <= kmin) error(FUNCTION, "(kmax = %d) <= (kmin = %d). This is forbidden.", kmax, kmin);
	subdmn->enabled = enable_subdomain;
	subdmn->imin = imin;
	subdmn->jmin = jmin;
	subdmn->kmin = kmin;
	subdmn->imax = imax;
	subdmn->jmax = jmax;
	subdmn->kmax = kmax;
	
	//loading number of saved checkpoints and snapshot to build the appropriate lists
	//accounting for the parameters just read
	n_saved_chkp = get_checkpoints_num(bin_data_dir, "sys");
	n_saved_snsht = get_checkpoints_num(bin_data_dir, "tecplot");
	printf("\to  total %d saved checkpoints and %d saved plotting snapshots detected\n", n_saved_chkp, n_saved_snsht);
	for(i = start_chkp; i < n_saved_chkp; i += (n_skipped_chkp+1)) n++;
	checkpoints_list->N = n;
	checkpoints_list->data = calloc(n, sizeof(int));
	n = 0;
	for(i = start_chkp; i < n_saved_chkp; i += (n_skipped_chkp+1)){
		checkpoints_list->data[n] = i;
		n++;
	}
	n = 0;
	for(i = start_snsht; i < n_saved_snsht; i += (n_skipped_snsht+1)) n++;
	snapshots_list->N = n;
	snapshots_list->data = calloc(n, sizeof(int));
	n = 0;
	for(i = start_snsht; i < n_saved_snsht; i += (n_skipped_snsht+1)){
		snapshots_list->data[n] = i;
		n++;
	}

	report_what_to_be_processed(checkpoints_list, snapshots_list);
}

void read_cfg_mdata(int *fs, fields_task *f_task, charges_task *c_task){
	//reads meshed data postprocessing parameters
	*fs = cfg_read_int("field_source", CFG_MDATA);
	f_task->ex = cfg_read_int("process_ex", CFG_MDATA);
	f_task->ey = cfg_read_int("process_ey", CFG_MDATA);
	f_task->ez = cfg_read_int("process_ez", CFG_MDATA);
	f_task->hx = cfg_read_int("process_hx", CFG_MDATA);
	f_task->hy = cfg_read_int("process_hy", CFG_MDATA);
	f_task->hz = cfg_read_int("process_hz", CFG_MDATA);
	f_task->I = cfg_read_int("process_I", CFG_MDATA);
	c_task->rho = cfg_read_int("process_rho", CFG_MDATA);
	c_task->jx = cfg_read_int("process_jx", CFG_MDATA);
	c_task->jy = cfg_read_int("process_jy", CFG_MDATA);
	c_task->jz = cfg_read_int("process_jz", CFG_MDATA);
	if((*fs != 1)&&(*fs != 2)) error(FUNCTION, "unvalid field source: field_source = %d\n", *fs);
	f_task->ex = (f_task->ex > 0)?1:0;
	f_task->ey = (f_task->ey > 0)?1:0;
	f_task->ez = (f_task->ez > 0)?1:0;
	f_task->hx = (f_task->hx > 0)?1:0;
	f_task->hy = (f_task->hy > 0)?1:0;
	f_task->hz = (f_task->hz > 0)?1:0;
	f_task->I = (f_task->I > 0)?1:0;
	c_task->rho = (c_task->rho > 0)?1:0;
	c_task->jx = (c_task->jx > 0)?1:0;
	c_task->jy = (c_task->jy > 0)?1:0;
	c_task->jz = (c_task->jz > 0)?1:0;
	f_task->any = ( f_task->ex + f_task->ey + f_task->ez + f_task->hx + f_task->hy + f_task->hz + f_task->I > 0)?1:0;
	c_task->any = ( c_task->jx + c_task->jy + c_task->jz + c_task->rho > 0)?1:0;
}

void read_cfg_mslicing(int *mp, int *cm, int_array *cuts){
	//reads and preprocesses the cutting info
	//for automatic cutting just leaves cuts.data filled with zeros
	int cn, i, l, ncommas = 0;
	char mcuts[STRLEN], mcuts_backup[STRLEN];
	char delim[] = ",";
	char *tmp = NULL;

	*mp = cfg_read_int("main_plane", CFG_MCUTS);
	cn = cfg_read_int("cuts_number", CFG_MCUTS);
	cuts->N = cn;
	cuts->data = calloc(cn, sizeof(int));
	*cm = cfg_read_int("cutting_method", CFG_MCUTS);
	if(*cm == 2){
		cfg_read_str("manual_cuts", mcuts, CFG_MCUTS);
		strcpy(mcuts_backup, mcuts);
		l = strlen(mcuts);
		for(i = 0; i < l; i++) if(mcuts[i] == ',') ncommas++;
		if((ncommas+1) != cn){
			cn = ncommas+1;
			cuts->N = cn;
			cuts->data = realloc(cuts->data, cn*sizeof(int));
		}

		tmp = strtok(mcuts, delim);
		for(i = 0; i < cn; i++){
			if(tmp == NULL) error(FUNCTION, "inconsistency in string manual_cuts = \"%s\" in file %s", mcuts_backup, CFG_MCUTS);
			cuts->data[i] = atoi(tmp);
			tmp = strtok(NULL, delim);
		}
	}
}

void read_cfg_markers(phase_space_task *ps_task){
	//reades the parameters of the phase space visualization
	
	int sum;

	ps_task->q_div_m = cfg_read_double("q_div_m", CFG_PS);

	ps_task->energy_filter = cfg_read_int("energy_filter", CFG_PS);
	ps_task->energy_filter = (ps_task->energy_filter > 0)?1:0;
	if(ps_task->energy_filter){
		ps_task->ef.min = cfg_read_double("ef_min_gamma", CFG_PS);
		ps_task->ef.max = cfg_read_double("ef_max_gamma", CFG_PS);
	}

	ps_task->x = cfg_read_int("ps_x", CFG_PS);
	ps_task->y = cfg_read_int("ps_y", CFG_PS);
	ps_task->z = cfg_read_int("ps_z", CFG_PS);
	ps_task->px = cfg_read_int("ps_px", CFG_PS);
	ps_task->py = cfg_read_int("ps_py", CFG_PS);
	ps_task->pz = cfg_read_int("ps_pz", CFG_PS);
	//ps_task->coarsen = cfg_read_int("s_coarsening_enable", CFG_PS);
	ps_task->coarsen = cfg_read_int("coarsen_data", CFG_COMMON);

	ps_task->x = (ps_task->x > 0)?1:0;
	ps_task->y = (ps_task->y > 0)?1:0;
	ps_task->z = (ps_task->z > 0)?1:0;
	ps_task->px = (ps_task->px > 0)?1:0;
	ps_task->py = (ps_task->py > 0)?1:0;
	ps_task->pz = (ps_task->pz > 0)?1:0;
	ps_task->coarsen = (ps_task->coarsen > 0)?1:0;

	sum = ps_task->x + ps_task->y + ps_task->z + ps_task->px + ps_task->py + ps_task->pz;
	if(sum == 0){
		ps_task->any = 0;
		return;
	}
	else ps_task->any = 1;
	if(sum < 2) error(FUNCTION, "less than 2 dimensions in phase space is requested for visualization. The number of dimensions should be 0 or 2. Right now it's %d.", sum);
	if(sum > 2) error(FUNCTION, "more than 2 dimensions in phase space is requested for visualization. The number of dimensions should be 0 or 2. Right now it's %d.", sum);

	if(ps_task->px){
		ps_task->dpx.min = cfg_read_double("px_min", CFG_PS);
		ps_task->dpx.max = cfg_read_double("px_max", CFG_PS);
	}
	if(ps_task->py){
		ps_task->dpy.min = cfg_read_double("py_min", CFG_PS);
		ps_task->dpy.max = cfg_read_double("py_max", CFG_PS);
	}
	if(ps_task->pz){
		ps_task->dpz.min = cfg_read_double("pz_min", CFG_PS);
		ps_task->dpz.max = cfg_read_double("pz_max", CFG_PS);
	}

	report_phase_space(ps_task);
}

void read_cfg_distrib(distrib_function_task *df_task){
	//reads info for creating distribution function
	int sum = 0;
	df_task->any = 0;

	df_task->q_div_m = cfg_read_double("q_div_m", CFG_PS);
	df_task->energy_filter = cfg_read_int("energy_filter", CFG_PS);
	df_task->energy_filter = (df_task->energy_filter > 0)?1:0;
	if(df_task->energy_filter){
		df_task->ef.min = cfg_read_double("ef_min_gamma", CFG_PS);
		df_task->ef.max = cfg_read_double("ef_max_gamma", CFG_PS);
	}

	df_task->x = ( cfg_read_int("distrib_x", CFG_PS) > 0 ) ? 1 : 0;
	df_task->y = ( cfg_read_int("distrib_y", CFG_PS) > 0 ) ? 1 : 0;
	df_task->z = ( cfg_read_int("distrib_z", CFG_PS) > 0 ) ? 1 : 0;
	df_task->px = ( cfg_read_int("distrib_px", CFG_PS) > 0 ) ? 1 : 0;
	df_task->py = ( cfg_read_int("distrib_py", CFG_PS) > 0 ) ? 1 : 0;
	df_task->pz = ( cfg_read_int("distrib_pz", CFG_PS) > 0 ) ? 1 : 0;
	df_task->gamma = ( cfg_read_int("distrib_gamma", CFG_PS) > 0 ) ? 1 : 0;

	df_task->theta_pos_x = ( cfg_read_int("distrib_theta_pos_x", CFG_PS) > 0 ) ? 1 : 0;
	df_task->theta_pos_y = ( cfg_read_int("distrib_theta_pos_y", CFG_PS) > 0 ) ? 1 : 0;
	df_task->theta_pos_z = ( cfg_read_int("distrib_theta_pos_z", CFG_PS) > 0 ) ? 1 : 0;
	df_task->theta_mom_x = ( cfg_read_int("distrib_theta_mom_x", CFG_PS) > 0 ) ? 1 : 0;
	df_task->theta_mom_y = ( cfg_read_int("distrib_theta_mom_y", CFG_PS) > 0 ) ? 1 : 0;
	df_task->theta_mom_z = ( cfg_read_int("distrib_theta_mom_z", CFG_PS) > 0 ) ? 1 : 0;
	df_task->phi_pos_xy = ( cfg_read_int("distrib_phi_pos_xy", CFG_PS) > 0 ) ? 1 : 0;
	df_task->phi_pos_yz = ( cfg_read_int("distrib_phi_pos_yz", CFG_PS) > 0 ) ? 1 : 0;
	df_task->phi_pos_xz = ( cfg_read_int("distrib_phi_pos_xz", CFG_PS) > 0 ) ? 1 : 0;
	df_task->phi_mom_xy = ( cfg_read_int("distrib_phi_mom_xy", CFG_PS) > 0 ) ? 1 : 0;
	df_task->phi_mom_yz = ( cfg_read_int("distrib_phi_mom_yz", CFG_PS) > 0 ) ? 1 : 0;
	df_task->phi_mom_xz = ( cfg_read_int("distrib_phi_mom_xz", CFG_PS) > 0 ) ? 1 : 0;

	sum = (df_task->x) + (df_task->y) + (df_task->z) + (df_task->px) + (df_task->py) + (df_task->pz) + (df_task->gamma) + (df_task->theta_pos_x) + 
		(df_task->theta_pos_y) + (df_task->theta_pos_z) + (df_task->theta_mom_x) + (df_task->theta_mom_y) + (df_task->theta_mom_z) + 
		(df_task->phi_pos_xy) + (df_task->phi_pos_yz) + (df_task->phi_pos_xz) + (df_task->phi_mom_xy) + (df_task->phi_mom_yz) + (df_task->phi_mom_xz);
	if(sum > 2) error(FUNCTION, "Histogram creation over %d parameters is requested. Only 1 or 2 are allowed", sum);

	if(sum > 0){
		df_task->nhist1 = cfg_read_int("nhist_parameter_1", CFG_PS);
		df_task->nhist2 = cfg_read_int("nhist_parameter_2", CFG_PS);

		if(df_task->nhist1 < 2) error(FUNCTION, "nhist_parameter_1 = %d is specified. Minumum allowed value is 2", df_task->nhist1);
		if(sum == 2){if(df_task->nhist2 < 2) error(FUNCTION, "nhist_parameter_2 = %d is specified. Minumum allowed value is 2", df_task->nhist2);}
		else df_task->nhist2 = 1;

		df_task->i1.min = cfg_read_double("parameter_1_min", CFG_PS);
		df_task->i1.max = cfg_read_double("parameter_1_max", CFG_PS);
		if(df_task->i1.min >= df_task->i1.max) error(FUNCTION, "(parameter_1_min = %e) >= (parameter_1_max = %e)", df_task->i1.min, df_task->i1.max);
		df_task->i2.min = cfg_read_double("parameter_2_min", CFG_PS);
		df_task->i2.max = cfg_read_double("parameter_2_max", CFG_PS);
		if(df_task->i2.min >= df_task->i2.max) error(FUNCTION, "(parameter_2_min = %e) >= (parameter_2_max = %e)", df_task->i2.min, df_task->i2.max);

		df_task->log_distr = ( cfg_read_int("logarithmic_distribution", CFG_PS) > 0 ) ? 1 : 0;

		df_task->any = 1;
	}

	df_task->sum = sum;

	df_report(df_task, sum);
}

void read_cfg_special_task(char *bin_data_dir, special_task *st, domain dmn){
	//reads the special tasks parameters from the config file
	
	char probes_nums[STRLEN], probes_nums_bak[STRLEN], *tmp = NULL;
	int i, l, ncommas = 0;
	char delim[] = ",";
	FILE *fp, *fp1;

	st->WE = ( cfg_read_int("show_WE", CFG_ST) > 0 ) ? 1 : 0;
	st->WM = ( cfg_read_int("show_WM", CFG_ST) > 0 ) ? 1 : 0;
	st->WTx = ( cfg_read_int("show_WTx", CFG_ST) > 0 ) ? 1 : 0;
	st->WTy = ( cfg_read_int("show_WTy", CFG_ST) > 0 ) ? 1 : 0;
	st->WTz = ( cfg_read_int("show_WTz", CFG_ST) > 0 ) ? 1 : 0;
	st->WT_perp = ( cfg_read_int("show_WT_perp", CFG_ST) > 0 ) ? 1 : 0;
	st->W_kin = ( cfg_read_int("show_W_kin", CFG_ST) > 0 ) ? 1 : 0;
	st->WS = ( cfg_read_int("show_WS", CFG_ST) > 0 ) ? 1 : 0;
	st->e_diag = (st->WE + st->WM + st->WTx + st->WTy + st->WTz + st->WT_perp + st->W_kin + st->WS);

	st->vis_probes = ( cfg_read_int("probes_processor_on", CFG_ST) > 0 ) ? 1 : 0;
	if(st->vis_probes){
		//reading, what probes are supposed to be processed
		cfg_read_str("probes_nums", probes_nums, CFG_ST);
		strcpy(probes_nums_bak, probes_nums);
		l = strlen(probes_nums);
		for(i = 0; i < l; i++) if(probes_nums[i] == ',') ncommas++;

		st->probes.N = ncommas+1;
		if(!(st->probes.data = calloc(st->probes.N, sizeof(int)))) error(FUNCTION, "cannot allocate %d kbytes for storing probes numbers", st->probes.N*sizeof(int)/1024);

		//parsing probes_num
		tmp = strtok(probes_nums, delim);
		for(i = 0; i < st->probes.N; i++){
			if(tmp == NULL) error(FUNCTION, "inconsistency in string manual_cuts = \"%s\" in file %s", probes_nums_bak, CFG_ST);
			st->probes.data[i] = atoi(tmp);
			tmp = strtok(NULL, delim);
		}

		sprintf(probes_nums, "%s../output/probe_locations.dat", bin_data_dir);
		OPEN_FILE(fp, probes_nums, "rt");
		
		for(i = 0; i < 11; i++) fgets(probes_nums, STRLEN, fp);		//skipping unnecessary data

		l = 0;
		while(fgets(probes_nums, STRLEN, fp)) if(probe_is_inside_domain(probes_nums, dmn)) l++;
		st->probes_num_domain = l;

		//check consistency
		for(i = 0; i < st->probes.N; i++){
			if(st->probes.data[i] > st->probes_num_domain) error(FUNCTION, "probe number %d requested which is larger than the total number of probes = %d", st->probes.data[i], st->probes_num_domain);
			if(st->probes.data[i] < 0) error(FUNCTION, "probe number #d is requested. This number should be greater than or equal to zero", st->probes.data[i]);
		}

		if((st->probes.N == 1)&&(st->probes.data[0] == 0)){
			st->probes.N = l;
			if(!(st->probes.data = calloc(st->probes.N, sizeof(int)))) error(FUNCTION, "cannot allocate %d bytes for storing probes numbers", st->probes.N*sizeof(int));
			for(i = 0; i < st->probes.N; i++) st->probes.data[i] = i+1;
		}

		fseek(fp, 0, SEEK_SET);
		for(i = 0; i < 11; i++) fgets(probes_nums, STRLEN, fp);	//skipping unnecessary data

		if(!(st->probe_locations = calloc(l, sizeof(float_vector)))) error(FUNCTION, "cannot allocate %d bytes for storing probes numbers", l*sizeof(float_vector));
		OPEN_FILE(fp1, "dat/probes.log", "wt");
		fprintf(fp1, "the following probese were detected:\n");
		i = 0;
		while(fgets(probes_nums, STRLEN, fp)){
			if(probe_is_inside_domain(probes_nums, dmn)){
				sscanf(probes_nums, "%f %f %f", &(st->probe_locations[i].x), &(st->probe_locations[i].y), &(st->probe_locations[i].z));
				fprintf(fp1, "(%d) [%.3f %.3f %.3f]\n", i+1, (st->probe_locations[i].x), (st->probe_locations[i].y), (st->probe_locations[i].z));
				i++;
			}
		}

		fclose(fp1);
		fclose(fp);

		st->vis_probe_locations = ( cfg_read_int("probes_locations_visualizer", CFG_ST) > 0 ) ? 1 : 0;
		st->probe_Ex = ( cfg_read_int("probe_Ex", CFG_ST) > 0 ) ? 1 : 0;
		st->probe_Ey = ( cfg_read_int("probe_Ey", CFG_ST) > 0 ) ? 1 : 0;
		st->probe_Ez = ( cfg_read_int("probe_Ez", CFG_ST) > 0 ) ? 1 : 0;
		st->probe_Hx = ( cfg_read_int("probe_Hx", CFG_ST) > 0 ) ? 1 : 0;
		st->probe_Hy = ( cfg_read_int("probe_Hy", CFG_ST) > 0 ) ? 1 : 0;
		st->probe_Hz = ( cfg_read_int("probe_Hz", CFG_ST) > 0 ) ? 1 : 0;
	}

	st->vis_fourier = ( cfg_read_int("fourier_processor_on", CFG_ST) > 0 ) ? 1 : 0;
	if(st->vis_fourier){
		st->enable_intervals = ( cfg_read_int("enable_k_intervals", CFG_ST) > 0 ) ? 1 : 0;
		st->ix.min = cfg_read_double("kx_min", CFG_ST);
		st->ix.max = cfg_read_double("kx_max", CFG_ST);
		st->iy.min = cfg_read_double("ky_min", CFG_ST);
		st->iy.max = cfg_read_double("ky_max", CFG_ST);
		st->iz.min = cfg_read_double("kz_min", CFG_ST);
		st->iz.max = cfg_read_double("kz_max", CFG_ST);

		if(st->ix.max < st->ix.min) error(FUNCTION, "kx_max = %f < kx_min = %f", st->ix.max, st->ix.min);
		if(st->iy.max < st->iy.min) error(FUNCTION, "ky_max = %f < ky_min = %f", st->iy.max, st->iy.min);
		if(st->iz.max < st->iz.min) error(FUNCTION, "kz_max = %f < kz_min = %f", st->iz.max, st->iz.min);
	}
	st->vis_gamma = 1; // TEMP

	st_report(st);
}

int probe_is_inside_domain(char *probe, domain dmn){
	//checks if the specified probe is inside the domain
	double x, y, z;
	sscanf(probe, "%le %le %le", &x, &y, &z);
	if( (x<dmn.Lx)&&(y<dmn.Ly)&&(z<dmn.Lz) ) return 1;
	else return 0;
}



