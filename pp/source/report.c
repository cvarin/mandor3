//functions reporting anything

#include<stdio.h>
#include<string.h>
#include"types.h"
#include"utils.h"
#include"params.h"

void report_what_to_be_processed(int_array *checkpoints_list, int_array *snapshots_list){
	//reports the program's view of what binary data is going to be processed
	int i;
	if(checkpoints_list->N > 0){
		printf("\to  the following checkpoints numbers will be processed: ");
		for(i = 0; i < (checkpoints_list->N - 1); i++) printf("%d, ", checkpoints_list->data[i]);
		printf("%d\n", checkpoints_list->data[checkpoints_list->N - 1]);
	}
	if(snapshots_list->N > 0){
		printf("\to  the following plotting snapshots numbers will be processed: ");
		for(i = 0; i < (snapshots_list->N - 1); i++) printf("%d, ", snapshots_list->data[i]);
		printf("%d\n", snapshots_list->data[snapshots_list->N - 1]);
	}

}

void report_domain(domain dmn){
	//reports the loaded domain information
	printf("Domain information of the Mandor data being loaded:\n");
	printf("\to  h1 = %f, h2 = %f, h3 = %f, Lx = %f, Ly = %f, Lz = %f, tau = %f\n", dmn.h1, dmn.h2, dmn.h3, dmn.Lx, dmn.Ly, dmn.Lz, dmn.tau);
	printf("\to  imin = %d, jmin = %d, kmin = %d, imax = %d, jmax = %d, kmax = %d\n", dmn.imin, dmn.jmin, dmn.kmin, dmn.imax, dmn.jmax, dmn.kmax);
	printf("\to  the domain dimensionality is %dD", dmn.dim);
	if(dmn.dim == 2){
		if(dmn.mp_2D == 1) printf(", main plane is XY");
		else if(dmn.mp_2D == 2) printf(", main plane is YZ");
		else printf(", main plane is XZ");
	}
	if(dmn.dim == 1){
		if(dmn.aa_1D == 1) printf(", active axis is X");
		else if(dmn.aa_1D == 2) printf(", active axis is Y");
		else printf(", active axis is Z");
	}
	printf("\n");
}

void report_mesh_slicing(int_array mcuts, int main_plane, domain dmn){
	//reports slicing for 3D
	char plane[3], axis[2];
	int i;
	double h = 0;

	if(main_plane == 1){
		strcpy(plane, "XY");
		strcpy(axis, "z");
		h = dmn.h3;
	}
	else if(main_plane == 2){
		strcpy(plane, "YZ");
		strcpy(axis, "x");
		h = dmn.h1;
	}
	else if(main_plane == 3){
		strcpy(plane, "XZ");
		strcpy(axis, "y");
		h = dmn.h2;
	}
	else error(FUNCTION, "unknown main_plane = %d", main_plane);

	printf("I will draw pictures for %d plane(s) parallel to the %s plane and having coordinates %s = ", mcuts.N, plane, axis);
	for(i = 0; i < (mcuts.N-1); i++) printf("%f, ", h*mcuts.data[i]);
	printf("%f lambdas\n", h*mcuts.data[mcuts.N-1]);
}

void report_subdomain(domain dmn, subdomain subdmn){
	//reports the just read subdomain information
	if(subdmn.enabled){
		printf("The following subdomain will be attempted to cut from the domain:\n");
		printf("\to  imin = %d, jmin = %d, kmin = %d <==> xmin = %f, ymin = %f, zmin = %f lambdas\n", subdmn.imin, subdmn.jmin, subdmn.kmin, dmn.h1*subdmn.imin, dmn.h2*subdmn.jmin, dmn.h3*subdmn.kmin);
		printf("\to  imax = %d, jmax = %d, kmax = %d <==> xmax = %f, ymax = %f, zmax = %f lambdas\n", subdmn.imax, subdmn.jmax, subdmn.kmax, dmn.h1*subdmn.imax, dmn.h2*subdmn.jmax, dmn.h3*subdmn.kmax);
	}
}

void report_mcoarsening(int_vector *cf){
	//reports the coarsening state of meshed data if necessary
	if( (cf->i > 1)||(cf->j > 1)||(cf->k > 1) ){
		printf("The meshed data will be coarsened by the following factors: cf.i = %d, cf.j = %d, cf.k = %d\n", cf->i, cf->j, cf->k);
	}
}

void report_phase_space(phase_space_task *pt){
	//reports parameters of the phase space visualization
	char a1[3], a2[3];
	int tmp = 0;
	if(!pt->any) return;

	printf("The plasma in phase space will be visualized with the following parameters:\n");
	if(pt->q_div_m > MAX_Q_DIV_M) printf("\to  all the species will be shown\n");
	else printf("\to  only specie with q/m = %f will be shown\n", pt->q_div_m);

	if(pt->energy_filter) printf("\to  energy filter is on: only particles having gamma-factors from interval [%f:%f] are accounted for\n", pt->ef.min, pt->ef.max);
	else printf("\to  energy filter is off\n");

	if(pt->x){
		sprintf(a1, "X");
		tmp = 1;
	}
	if(pt->y){
		if(!tmp) sprintf(a1, "Y");
		else sprintf(a2, "Y");
		tmp = 1;
	}
	if(pt->z){
		if(!tmp) sprintf(a1, "Z");
		else sprintf(a2, "Z");
		tmp = 1;
	}
	if(pt->px){
		if(!tmp) sprintf(a1, "Px");
		else sprintf(a2, "Px");
		tmp = 1;
	}
	if(pt->py){
		if(!tmp) sprintf(a1, "Py");
		else sprintf(a2, "Py");
		tmp = 1;
	}
	if(pt->pz) sprintf(a2, "Pz");

	printf("\to  the projection onto the %s-%s plane of the phase space will be visualized\n", a1, a2);
	if(pt->px) printf("\to  the visualization domain is limited to Px in [%f:%f]\n", pt->dpx.min, pt->dpx.max);
	if(pt->px) printf("\to  the visualization domain is limited to Py in [%f:%f]\n", pt->dpy.min, pt->dpy.max);
	if(pt->px) printf("\to  the visualization domain is limited to Pz in [%f:%f]\n", pt->dpz.min, pt->dpz.max);

	if(pt->coarsen) printf("\to  an automatic coarsening algorithm will be applied to the data\n");
}

void df_report(distrib_function_task *dt, int max){
	//reports the distribution function creation task
	int sum = 0;
	char andsymb[5];

	if(!dt->any) return;
	printf("The distribution function will be shown for ");
	if(dt->q_div_m > MAX_Q_DIV_M) printf("all particles in (sub)domain ");
	else printf("particles in (sub)domain having q/m = %e ", dt->q_div_m);
	if(dt->energy_filter) printf("and gamma-factors between %e and %e\n", dt->ef.min, dt->ef.max);
	else printf("\n");
	printf("\to  histogram over ");
#define REPORT_AXIS(symb_) {							\
	if(sum == (max-2)) sprintf(andsymb, "and "); else andsymb[0] = '\0';	\
	printf("%s %s", symb_, andsymb);					\
	sum++;									\
}
	if(dt->x) REPORT_AXIS("particles x/lambda");
	if(dt->y) REPORT_AXIS("particles y/lambda");
	if(dt->z) REPORT_AXIS("particles z/lambda");
	if(dt->px) REPORT_AXIS("particles px/mc");
	if(dt->py) REPORT_AXIS("particles py/mc");
	if(dt->pz) REPORT_AXIS("particles pz/mc");
	if(dt->gamma) REPORT_AXIS("particles gamma-factor");
	if(dt->theta_pos_x) REPORT_AXIS("particles theta_pos_x");
	if(dt->theta_pos_y) REPORT_AXIS("particles theta_pos_y");
	if(dt->theta_pos_z) REPORT_AXIS("particles theta_pos_z");
	if(dt->theta_mom_x) REPORT_AXIS("particles theta_mom_x");
	if(dt->theta_mom_y) REPORT_AXIS("particles theta_mom_y");
	if(dt->theta_mom_z) REPORT_AXIS("particles theta_mom_z");
	if(dt->phi_pos_xy) REPORT_AXIS("particles phi_pos_xy");
	if(dt->phi_pos_yz) REPORT_AXIS("particles phi_pos_yz)");
	if(dt->phi_pos_xz) REPORT_AXIS("particles phi_pos_xz)");
	if(dt->phi_mom_xy) REPORT_AXIS("particles phi_mom_xy)");
	if(dt->phi_mom_yz) REPORT_AXIS("particles phi_mom_yz)");
	if(dt->phi_mom_xz) REPORT_AXIS("particles phi_mom_xz)");

	printf("\n");
	printf("\to  number of histogram mesh points is %d", dt->nhist1);
	if(sum == 2) printf("x%d", dt->nhist2);
	printf("\n");

	sum = 0;
#define REPORT_AXIS_INTERVAL(symb_) {													\
	if(sum == 0) printf("\to  histogram over %s will be created in the interval [%.3e, %.3e]\n", symb_, dt->i1.min, dt->i1.max);	\
	else printf("\to  histogram over %s will be created in the interval [%.3e, %.3e]\n", symb_, dt->i2.min, dt->i2.max);		\
	sum++;																\
}
	if(dt->x) REPORT_AXIS_INTERVAL("x/lambda");
	if(dt->y) REPORT_AXIS_INTERVAL("y/lambda");
	if(dt->z) REPORT_AXIS_INTERVAL("z/lambda");
	if(dt->px) REPORT_AXIS_INTERVAL("px/mc");
	if(dt->py) REPORT_AXIS_INTERVAL("py/mc");
	if(dt->pz) REPORT_AXIS_INTERVAL("pz/mc");
	if(dt->gamma) REPORT_AXIS_INTERVAL("gamma-factor");
	if(dt->theta_pos_x) REPORT_AXIS_INTERVAL("theta_pos_x");
	if(dt->theta_pos_y) REPORT_AXIS_INTERVAL("theta_pos_y");
	if(dt->theta_pos_z) REPORT_AXIS_INTERVAL("theta_pos_z");
	if(dt->theta_mom_x) REPORT_AXIS_INTERVAL("theta_mom_x");
	if(dt->theta_mom_y) REPORT_AXIS_INTERVAL("theta_mom_y");
	if(dt->theta_mom_z) REPORT_AXIS_INTERVAL("theta_mom_z");
	if(dt->phi_pos_xy) REPORT_AXIS_INTERVAL("phi_pos_xy");
	if(dt->phi_pos_yz) REPORT_AXIS_INTERVAL("phi_pos_yz)");
	if(dt->phi_pos_xz) REPORT_AXIS_INTERVAL("phi_pos_xz)");
	if(dt->phi_mom_xy) REPORT_AXIS_INTERVAL("phi_mom_xy)");
	if(dt->phi_mom_yz) REPORT_AXIS_INTERVAL("phi_mom_yz)");
	if(dt->phi_mom_xz) REPORT_AXIS_INTERVAL("phi_mom_xz)");

	if(dt->log_distr) printf("\to  the logarithmic histogram will be created\n");
}

void st_report(special_task *st){
	//reports the special tasks to be performed
	int i;

	printf("The following special visualization tasks have been requested:\n");
	if(st->WE) printf("\to  electric energy diagnostic\n");
	if(st->WM) printf("\to  magnetic energy diagnostic\n");
	if(st->WTx) printf("\to  particles x-component of temperature\n");
	if(st->WTy) printf("\to  particles y-component of temperature\n");
	if(st->WTz) printf("\to  particles z-component of temperature\n");
	if(st->WT_perp) printf("\to  particles perpendicular temperature\n");
	if(st->W_kin) printf("\to  particles total kinetic energy\n");
	if(st->WS) printf("\to  total energy in the domain (fields + particles)\n");

	if(st->vis_probes){
		if(st->vis_probe_locations) printf("\to  probes locations visualization\n");
		printf("\to  data visualization from the following probes:");
		for(i = 0; i < st->probes.N; i++){
			printf(" (%d) [%.2f, %.2f, %.2f]", st->probes.data[i], st->probe_locations[st->probes.data[i]-1].x, st->probe_locations[st->probes.data[i]-1].y, st->probe_locations[st->probes.data[i]-1].z);
			if(i != (st->probes.N - 1)) printf(",");
			else printf(" (for the list of the available see dat/probes.log)\n");
		}
		if(st->probes_num_domain == 0) printf(" no probes were found inside the domain\n");
	}

	if(st->vis_fourier) printf("\to  fourier transformation in space\n");
}






