//important routines used for data processing

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<fftw3.h>
#include"utils.h"
#include"types.h"
#include"params.h"
#include"data_loader.h"
#include"processor.h"
#include"matrix.h"
#include"report.h"
#include"plot.h"

void auto_slice(subdomain subdmn, int_array *mcuts, int mmain_plane){
	//performs auto slicing
	int i, n = 0, m = 0;
	if(mmain_plane == 1){
		n = subdmn.kmax - subdmn.kmin + 1;
		m = subdmn.kmin;
	}
	else if(mmain_plane == 2){
		n = subdmn.imax - subdmn.imin + 1;
		m = subdmn.imin;
	}
	else if(mmain_plane == 3){
		n = subdmn.jmax - subdmn.jmin + 1;
		m = subdmn.jmin;
	}
	else error(FUNCTION, "cannot recognize mmain_plane code (mmain_plane = %d)", mmain_plane);
	for(i = 0; i < mcuts->N; i++) mcuts->data[i] = m + n*( (i+0.5)/mcuts->N );
}

void check_slicing(int_array mcuts, int mmain_plane, subdomain subdmn){
	//checks if all the slices are inside the subdomain
	int i, nmax = 0, nmin = 0;
	if(mmain_plane == 1){
		nmax = subdmn.kmax;
		nmin = subdmn.kmin;
	}
	if(mmain_plane == 2){
		nmax = subdmn.imax;
		nmin = subdmn.imin;
	}
	if(mmain_plane == 3){
		nmax = subdmn.jmax;
		nmin = subdmn.jmin;
	}
	for(i = 0; i < mcuts.N; i++){
		if(mcuts.data[i] > nmax) error(FUNCTION, "slice is outside subdomain: mcuts.data[%d] = %d; subdomain upper boundary is n = %d", i, mcuts.data[i], nmax);
		if(mcuts.data[i] < nmin) error(FUNCTION, "slice is outside subdomain: mcuts.data[%d] = %d; subdomain lower boundary is n = %d", i, mcuts.data[i], nmin);
	}
}

void correct_subdomain(domain dmn, subdomain *subdmn){
	//corrects subdomain to make sure it properly describes the task to be done
	if(!subdmn->enabled){
		subdmn->imin = dmn.imin;
		subdmn->jmin = dmn.jmin;
		subdmn->kmin = dmn.kmin;
		subdmn->imax = dmn.imax;
		subdmn->jmax = dmn.jmax;
		subdmn->kmax = dmn.kmax;
		subdmn->enabled = 1;
	}
	else{
#define CHECK(min_, max1_, max2_, spec_){if(min_ >= max2_) error(FUNCTION, "(subdmn.%smin = %d) >= (dmn.%smax = %d)", spec_, min_, spec_, max2_); if(max2_ < max1_) max1_ = max2_;}
		if(dmn.dim == 3){
			CHECK( subdmn->imin, subdmn->imax, dmn.imax, "i" );
			CHECK( subdmn->jmin, subdmn->jmax, dmn.jmax, "i" );
			CHECK( subdmn->kmin, subdmn->kmax, dmn.kmax, "i" );
		}
		else if(dmn.dim == 2){
			if(dmn.mp_2D == 1){
				subdmn->kmin = dmn.kmin;
				subdmn->kmax = dmn.kmax;
				CHECK( subdmn->imin, subdmn->imax, dmn.imax, "i" );
				CHECK( subdmn->jmin, subdmn->jmax, dmn.jmax, "j" );
			}
			else if(dmn.mp_2D == 2){
				subdmn->imin = dmn.imin;
				subdmn->imax = dmn.imax;
				CHECK( subdmn->jmin, subdmn->jmax, dmn.jmax, "j" );
				CHECK( subdmn->kmin, subdmn->kmax, dmn.kmax, "k" );
			}
			else{
				subdmn->jmin = dmn.jmin;
				subdmn->jmax = dmn.jmax;
				CHECK( subdmn->imin, subdmn->imax, dmn.imax, "i" );
				CHECK( subdmn->kmin, subdmn->kmax, dmn.kmax, "k" );
			}
		}
		else{
			if(dmn.aa_1D == 1){
				subdmn->jmin = dmn.jmin;
				subdmn->kmin = dmn.kmin;
				subdmn->jmax = dmn.jmax;
				subdmn->kmax = dmn.kmax;
				CHECK( subdmn->imin, subdmn->imax, dmn.imax, "i" );
			}
			else if(dmn.aa_1D == 2){
				subdmn->imin = dmn.imin;
				subdmn->kmin = dmn.kmin;
				subdmn->imax = dmn.imax;
				subdmn->kmax = dmn.kmax;
				CHECK( subdmn->jmin, subdmn->jmax, dmn.jmax, "j" );
			}
			else{
				subdmn->imin = dmn.imin;
				subdmn->jmin = dmn.jmin;
				subdmn->imax = dmn.imax;
				subdmn->jmax = dmn.jmax;
				CHECK( subdmn->kmin, subdmn->kmax, dmn.kmax, "k" );
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// MESHED DATA VISUALIZATION //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void get_coarse_factors(int_vector *cf, int enable, domain dmn, int  mmain_plane, subdomain subdmn, int hor_pic_size, int vert_pic_size){
	//adjusts coarse factor
	if(!enable){
		cf->i = 1;
		cf->j = 1;
		cf->k = 1;
	}
	else{
#define COARSEN(f1_, f2_, f3_) {cf->i = f1_+1; cf->j = f2_+1; cf->k = f3_+1;}
		hor_pic_size *= WEAKEN_COARSENING;
		vert_pic_size *= WEAKEN_COARSENING;
		if(dmn.dim == 3){
			if(mmain_plane == 1) COARSEN( (subdmn.imax - subdmn.imin)/hor_pic_size, (subdmn.jmax - subdmn.jmin)/vert_pic_size, 0 );
			if(mmain_plane == 2) COARSEN( (subdmn.jmax - subdmn.jmin)/hor_pic_size, (subdmn.kmax - subdmn.kmin)/vert_pic_size, 0 );
			if(mmain_plane == 3) COARSEN( (subdmn.imax - subdmn.imin)/hor_pic_size, (subdmn.kmax - subdmn.kmin)/vert_pic_size, 0 );
		}
		else if(dmn.dim == 2){
			if(dmn.mp_2D == 1) COARSEN( (subdmn.imax - subdmn.imin)/hor_pic_size, (subdmn.jmax - subdmn.jmin)/vert_pic_size, 0 );
			if(dmn.mp_2D == 2) COARSEN( (subdmn.jmax - subdmn.jmin)/hor_pic_size, (subdmn.kmax - subdmn.kmin)/vert_pic_size, 0 );
			if(dmn.mp_2D == 3) COARSEN( (subdmn.imax - subdmn.imin)/hor_pic_size, (subdmn.kmax - subdmn.kmin)/vert_pic_size, 0 );
		}
		else{
			if(dmn.aa_1D == 1) COARSEN( (subdmn.imax - subdmn.imin)/hor_pic_size, 0, 0 );
			if(dmn.aa_1D == 2) COARSEN( (subdmn.jmax - subdmn.jmin)/hor_pic_size, 0, 0 );
			if(dmn.aa_1D == 3) COARSEN( (subdmn.kmax - subdmn.kmin)/hor_pic_size, 0, 0 );
		}
	}
	report_mcoarsening(cf);
}


void process_fields(char *bin_data_dir, int num, domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, char *source_spec, FILE *plt, int hor_pic_size, int vert_pic_size, fields_task f_task, int_vector coarse_factor, special_task st){
	//does all the necessary operations to postprocess fields
	partition p;
	int i, j, k, cpu, process_E = 0, process_H = 0, point_readable = 0, nx, ny;
	vector E, H;
	double time;
	arr2d ex[mcuts.N], ey[mcuts.N], ez[mcuts.N], hx[mcuts.N], hy[mcuts.N], hz[mcuts.N], I[mcuts.N], edensity[mcuts.N];
	arr2d spectral_edensity;

	p = load_partition_info(bin_data_dir, num, source_spec, &time);
	FILE *fs[2*p.N];
	open_mesh_binary_files(bin_data_dir, source_spec, "E", num, fs, p.N);
	open_mesh_binary_files(bin_data_dir, source_spec, "H", num, fs+p.N, p.N);

	calc_cut_size(dmn, mmain_plane, mcuts, subdmn, coarse_factor, &nx, &ny);

	if(f_task.ex||f_task.ey||f_task.ez||f_task.I||st.vis_fourier) process_E = 1;
	if(f_task.hx||f_task.hy||f_task.hz||f_task.I||st.vis_fourier) process_H = 1;
	for(i = 0; i < mcuts.N; i++){
		if(f_task.ex) ex[i] = arr2d_create(nx, ny);
		if(f_task.ey) ey[i] = arr2d_create(nx, ny);
		if(f_task.ez) ez[i] = arr2d_create(nx, ny);
		if(f_task.hx) hx[i] = arr2d_create(nx, ny);
		if(f_task.hy) hy[i] = arr2d_create(nx, ny);
		if(f_task.hz) hz[i] = arr2d_create(nx, ny);
		if(f_task.I) I[i] = arr2d_create(nx, ny);
		if(st.vis_fourier) edensity[i] = arr2d_create(nx, ny);
	}

	for(i = dmn.imin; i <= dmn.imax; i++){
		for(j = dmn.jmin; j <= dmn.jmax; j++){
			for(k = dmn.kmin; k <= dmn.kmax; k++){
				if(process_E) E = read_vector_point(i, j, k, p, fs);
				if(process_H) H = read_vector_point(i, j, k, p, fs+p.N);
				point_readable = point_kind(i, j, k, dmn, mmain_plane, mcuts, subdmn, coarse_factor);
				if(point_readable){
					if(f_task.ex) savepoint(E.x/(2*PI), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, ex);
					if(f_task.ey) savepoint(E.y/(2*PI), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, ey);
					if(f_task.ez) savepoint(E.z/(2*PI), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, ez);
					if(f_task.hx) savepoint(H.x/(2*PI), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, hx);
					if(f_task.hy) savepoint(H.y/(2*PI), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, hy);
					if(f_task.hz) savepoint(H.z/(2*PI), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, hz);
					if(f_task.I) savepoint(intensity(E, H), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, I);
					if(st.vis_fourier) savepoint(energy_density(E, H), i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, edensity);
				}
			}
		}
	}
	//save data and update plt file
	if(f_task.ex) save_dat_files(dmn, subdmn, mmain_plane, mcuts, ex, num, "Ex", plt, time);
	if(f_task.ey) save_dat_files(dmn, subdmn, mmain_plane, mcuts, ey, num, "Ey", plt, time);
	if(f_task.ez) save_dat_files(dmn, subdmn, mmain_plane, mcuts, ez, num, "Ez", plt, time);
	if(f_task.hx) save_dat_files(dmn, subdmn, mmain_plane, mcuts, hx, num, "Hx", plt, time);
	if(f_task.hy) save_dat_files(dmn, subdmn, mmain_plane, mcuts, hy, num, "Hy", plt, time);
	if(f_task.hz) save_dat_files(dmn, subdmn, mmain_plane, mcuts, hz, num, "Hz", plt, time);
	if(f_task.I) save_dat_files(dmn, subdmn, mmain_plane, mcuts, I, num, "I", plt, time);

	for(cpu = 0; cpu < p.N; cpu++) fclose(fs[cpu]);
	for(i = 0; i < mcuts.N; i++){
		if(f_task.ex) arr2d_destroy(ex[i]);
		if(f_task.ey) arr2d_destroy(ey[i]);
		if(f_task.ez) arr2d_destroy(ez[i]);
		if(f_task.hx) arr2d_destroy(hx[i]);
		if(f_task.hy) arr2d_destroy(hy[i]);
		if(f_task.hz) arr2d_destroy(hz[i]);
		if(f_task.I) arr2d_destroy(I[i]);
	}
	printf("\t  - fields loaded\n");

	//performing fourier transformation if neccessary
	if(st.vis_fourier){
		for(i = 0; i < mcuts.N; i++){
			if(ny == 1) spectral_edensity = arr2d_create(nx, 1);
			else spectral_edensity = arr2d_create(nx, ny);
			perform_fourier_transform(edensity[i], &spectral_edensity);
			plot_fourier_transform(dmn, st, mmain_plane, mcuts, spectral_edensity, i, num, plt, hor_pic_size, vert_pic_size, time);
			arr2d_destroy(spectral_edensity);
			arr2d_destroy(edensity[i]);
		}
		printf("\t  - fourier transformation done\n");
	}

	partiton_erase(p);
}

void calc_cut_size(domain dmn, int mmain_plane, int_array mcuts, subdomain subdmn, int_vector coarse_factor, int *nx, int *ny){
	//calculates size of a single cut
	if(dmn.dim == 3){
		if(mmain_plane == 1){*nx = (subdmn.imax - subdmn.imin)/coarse_factor.i + 1; *ny = (subdmn.jmax - subdmn.jmin)/coarse_factor.j + 1;}
		if(mmain_plane == 2){*nx = (subdmn.jmax - subdmn.jmin)/coarse_factor.j + 1; *ny = (subdmn.kmax - subdmn.kmin)/coarse_factor.k + 1;}
		if(mmain_plane == 3){*nx = (subdmn.imax - subdmn.imin)/coarse_factor.i + 1; *ny = (subdmn.kmax - subdmn.kmin)/coarse_factor.k + 1;}
	}
	else if(dmn.dim == 2){
		if(dmn.mp_2D == 1){*nx = (subdmn.imax - subdmn.imin)/coarse_factor.i + 1; *ny = (subdmn.jmax - subdmn.jmin)/coarse_factor.j + 1;}
		if(dmn.mp_2D == 2){*nx = (subdmn.jmax - subdmn.jmin)/coarse_factor.j + 1; *ny = (subdmn.kmax - subdmn.kmin)/coarse_factor.k + 1;}
		if(dmn.mp_2D == 3){*nx = (subdmn.imax - subdmn.imin)/coarse_factor.i + 1; *ny = (subdmn.kmax - subdmn.kmin)/coarse_factor.k + 1;}
	}
	else{
		*ny = 1;
		if(dmn.aa_1D == 1) *nx = (subdmn.imax - subdmn.imin)/coarse_factor.i + 1;
		if(dmn.aa_1D == 2) *nx = (subdmn.jmax - subdmn.jmin)/coarse_factor.j + 1;
		if(dmn.aa_1D == 3) *nx = (subdmn.kmax - subdmn.kmin)/coarse_factor.k + 1;
	}

}

int point_kind(int i, int j, int k, domain dmn, int mmain_plane, int_array mcuts, subdomain subdmn, int_vector coarse_factor){
	//decides whether the point is to be read from the binary files or just can be skipped
	int cn;
	if(inside_subdomain(i, j, k, subdmn)){
		if(dmn.dim == 3){
			for(cn = 0; cn< mcuts.N; cn++){
				if((mmain_plane == 1)&&(k == mcuts.data[cn])&&((i-subdmn.imin)%coarse_factor.i == 0)&&((j-subdmn.jmin)%coarse_factor.j == 0)) return 1;
				if((mmain_plane == 2)&&(i == mcuts.data[cn])&&((j-subdmn.jmin)%coarse_factor.j == 0)&&((k-subdmn.kmin)%coarse_factor.k == 0)) return 1;
				if((mmain_plane == 3)&&(j == mcuts.data[cn])&&((i-subdmn.imin)%coarse_factor.i == 0)&&((k-subdmn.kmin)%coarse_factor.k == 0)) return 1;
			}
		}
		else if(dmn.dim == 2){
			if((dmn.mp_2D == 1)&&((i-subdmn.imin)%coarse_factor.i == 0)&&((j-subdmn.jmin)%coarse_factor.j == 0)) return 1;
			if((dmn.mp_2D == 2)&&((j-subdmn.jmin)%coarse_factor.j == 0)&&((k-subdmn.kmin)%coarse_factor.k == 0)) return 1;
			if((dmn.mp_2D == 3)&&((i-subdmn.imin)%coarse_factor.i == 0)&&((k-subdmn.kmin)%coarse_factor.k == 0)) return 1;
		}
		else{
			if((dmn.aa_1D == 1)&&((i-subdmn.imin)%coarse_factor.i == 0)) return 1;
			if((dmn.aa_1D == 2)&&((j-subdmn.jmin)%coarse_factor.j == 0)) return 1;
			if((dmn.aa_1D == 3)&&((k-subdmn.kmin)%coarse_factor.k == 0)) return 1;
		}
	}
	return 0;
}

int inside_subdomain(int i, int j, int k, subdomain subdmn){
	//decides whether the point is inside the domain or not
	if((i >= subdmn.imin)&&(i <= subdmn.imax)&&(j >= subdmn.jmin)&&(j <= subdmn.jmax)&&(k >= subdmn.kmin)&&(k <= subdmn.kmax)) return 1;
	else return 0;
}

int get_cut_num(int i, int j, int k, int mmain_plane, int_array mcuts){
	//figures out the number of the cut the point belongs to
	int n;
	for(n = 0; n < mcuts.N; i++){
		if((mmain_plane == 1)&&(i == mcuts.data[n])) return n;
		if((mmain_plane == 2)&&(j == mcuts.data[n])) return n;
		if((mmain_plane == 3)&&(k == mcuts.data[n])) return n;
	}
	error(FUNCTION, "point i,j,k = (%d,%d,%d) cannot be found in the mcuts", i, j, k);
	return 0;
}

void savepoint(double val, int i, int j, int k, domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, int_vector cf, arr2d *A){
	//determines the number of the cut (if necessary) and saves the data point into the array
#define SAVE_3D(i_, j_, k_) for(cnt = 0; cnt < mcuts.N; cnt++) if(k_ == mcuts.data[cnt]) arr2d_set_ij(val, A+cnt, i_, j_);
#define SAVE_2D(i_, j_) arr2d_set_ij(val, A, i_, j_);
#define SAVE_1D(i_) arr2d_set_ij(val, A, i_, 0);
	int cnt;
	if(dmn.dim == 3){
		if(mmain_plane == 1) SAVE_3D((i-subdmn.imin)/cf.i, (j-subdmn.jmin)/cf.j, k);
		if(mmain_plane == 2) SAVE_3D((j-subdmn.jmin)/cf.j, (k-subdmn.kmin)/cf.k, i);
		if(mmain_plane == 3) SAVE_3D((i-subdmn.imin)/cf.i, (k-subdmn.kmin)/cf.k, j);
	}
	else if(dmn.dim == 2){
		if(dmn.mp_2D == 1) SAVE_2D((i-subdmn.imin)/cf.i, (j-subdmn.jmin)/cf.j);
		if(dmn.mp_2D == 2) SAVE_2D((j-subdmn.jmin)/cf.j, (k-subdmn.kmin)/cf.k);
		if(dmn.mp_2D == 3) SAVE_2D((i-subdmn.imin)/cf.i, (k-subdmn.kmin)/cf.k);
	}
	else{
		if(dmn.aa_1D == 1) SAVE_1D((i-subdmn.imin)/cf.i);
		if(dmn.aa_1D == 2) SAVE_1D((j-subdmn.jmin)/cf.j);
		if(dmn.aa_1D == 3) SAVE_1D((k-subdmn.kmin)/cf.k);
	}
}

double intensity(vector E, vector H){
	//calculates intensity given the vectors E, H
	vector S;
	S.x = E.y*H.z - E.z*H.y;
	S.y = E.z*H.x - E.x*H.z;
	S.z = E.x*H.y - E.y*H.x;
	return(sqrt(S.x*S.x + S.y*S.y + S.z*S.z)/(16*PI*PI*PI));
}

double energy_density(vector E, vector H){
	//calculates energy density given the vectors E, H
	return((E.x*E.x + E.y*E.y + E.z*E.z + H.x*H.x + H.y*H.y + H.z*H.z)/(32*PI*PI*PI));
}


void save_dat_files(domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, arr2d *A, int num, char *spec, FILE *plt, double time){
	//saves a number of cuts into dat-files using procedures save_dat_file
	int i, mp;
	double xmin = 0, ymin = 0, xmax = 0, ymax = 0;
	double z[mcuts.N];
	char str[STRLEN];
#define BOUNDS(hx_, hy_, hz_, imin_, imax_, jmin_, jmax_){				\
	xmin = hx_*imin_; xmax = hx_*imax_; ymin = hy_*jmin_; ymax = hy_*jmax_;		\
	for(i = 0; i < mcuts.N; i++) z[i] = hz_*mcuts.data[i];				\
}
	if(dmn.dim == 3){
		if(mmain_plane == 1) BOUNDS(dmn.h1, dmn.h2, dmn.h2, subdmn.imin, subdmn.imax, subdmn.jmin, subdmn.jmax);
		if(mmain_plane == 2) BOUNDS(dmn.h2, dmn.h3, dmn.h1, subdmn.jmin, subdmn.jmax, subdmn.kmin, subdmn.kmax);
		if(mmain_plane == 3) BOUNDS(dmn.h1, dmn.h3, dmn.h2, subdmn.imin, subdmn.imax, subdmn.kmin, subdmn.kmax);
		mp = mmain_plane;
	}
	else if(dmn.dim == 2){
		if(dmn.mp_2D == 1) BOUNDS(dmn.h1, dmn.h2, 0, subdmn.imin, subdmn.imax, subdmn.jmin, subdmn.jmax);
		if(dmn.mp_2D == 2) BOUNDS(dmn.h2, dmn.h3, 0, subdmn.jmin, subdmn.jmax, subdmn.kmin, subdmn.kmax);
		if(dmn.mp_2D == 3) BOUNDS(dmn.h1, dmn.h3, 0, subdmn.imin, subdmn.imax, subdmn.kmin, subdmn.kmax);
		mp = dmn.mp_2D;
	}
	else{
		if(dmn.aa_1D == 1) BOUNDS(dmn.h1, 0, 0, subdmn.imin, subdmn.imax, 0, 0);
		if(dmn.aa_1D == 2) BOUNDS(dmn.h2, 0, 0, subdmn.jmin, subdmn.jmax, 0, 0);
		if(dmn.aa_1D == 3) BOUNDS(dmn.h3, 0, 0, subdmn.kmin, subdmn.kmax, 0, 0);
		mp = dmn.aa_1D;
	}
	for(i = 0; i < mcuts.N; i++){
		sprintf(str, "dat/%s_%06d(%03d).dat", spec, num, i);
		save_dat_file(A[i], str, dmn.dim, xmin, ymin, xmax, ymax, mp, z[i], spec, num, i, plt, time);
	}
}

void save_dat_file(arr2d A, char *fn, int dim, double xmin, double ymin, double xmax, double ymax, int mp, double z, char *spec, int num, int cut, FILE *plt, double time){
	//saves the array A into a dat file having name fn and updates the plt-file
	int i, j;
	double x, y;
	FILE *fp;
	OPEN_FILE(fp, fn, "wt");
	char xplane[2], yplane[2], zplane[2];
#define SET_PLANES(xp_, yp_, zp_) {sprintf(xplane, "%s", xp_); sprintf(yplane, "%s", yp_); sprintf(zplane, "%s", zp_);}
	if(mp == 1) SET_PLANES("x", "y", "z");
	if(mp == 2) SET_PLANES("y", "z", "x");
	if(mp == 3) SET_PLANES("x", "z", "y");

	fprintf(plt, "set output 'pics/%s_%06d(%03d).png'\n", spec, num, cut);

	if(dim == 1){
		for(i = 0; i < A.nx; i++){
			x = xmin + (i+0.5)*(xmax-xmin)/A.nx;
			fprintf(fp, "%e %e\n", x, arr2d_ij(A, i, 0));
		}
		fprintf(plt, "set title '%s: time = %.2f waveperiods'\n", spec, time);
		fprintf(plt, "set xlabel '%s/lambda'\n", (mp == 3)?"z":xplane);
		fprintf(plt, "set ylabel '%s'\n", spec);
		fprintf(plt, "set grid\n");
		fprintf(plt, "plot '%s' w l notitle\n", fn);
		fprintf(plt, "unset grid\n");
	}
	else{
		arr2d_set_ij(arr2d_ij(A, 0, 0) + 1e-100, &A, 0, 0);
		for(i = 0; i < A.nx; i++){
			for(j = 0; j < A.ny; j++){
				x = xmin + (i+0.5)*(xmax-xmin)/A.nx;
				y = ymin + (j+0.5)*(ymax-ymin)/A.ny;
				fprintf(fp, "%e %e %e\n", x, y, arr2d_ij(A, i, j));
			}
			fprintf(fp, "\n");
		}
		if(dim == 3) fprintf(plt, "set title '%s: cut at %s = %.2f lambdas, time = %.2f waveperiods'\n", spec, zplane, z, time);
		else fprintf(plt, "set title '%s: time = %.2f waveperiods'\n", spec, time);
		fprintf(plt, "set xlabel '%s/lambda'\n", xplane);
		fprintf(plt, "set ylabel '%s/lambda'\n", yplane);
		fprintf(plt, "set pm3d map\n");
		fprintf(plt, "splot '%s' notitle\n", fn);
		fprintf(plt, "unset pm3d\n");
	}
	fclose(fp);


}

void perform_fourier_transform(arr2d in, arr2d *out){
	//performs fourier transfromation of the field energy density stored in the array in
	//the resulting spectral energy density is saved in the array out
	fftw_complex *tmp_in, *tmp_out;
	fftw_plan p;
	int i, j, cnt;

	if(!in.allocated) error(FUNCTION, "array in is not allocated");

	if( !(tmp_in = (fftw_complex*) fftw_malloc(in.nx*in.ny*sizeof(fftw_complex))) ) error(FUNCTION, "cannot allocate %d kbytes to store fourier transformation", in.nx*in.ny*sizeof(fftw_complex)/1024);
	if( !(tmp_out = (fftw_complex*) fftw_malloc(in.nx*in.ny*sizeof(fftw_complex))) ) error(FUNCTION, "cannot allocate %d kbytes to store fourier transformation", in.nx*in.ny*sizeof(fftw_complex)/1024);

	p = fftw_plan_dft_2d(in.nx, in.ny, tmp_in, tmp_out, FFTW_FORWARD, FFTW_ESTIMATE);

	for(i = 0; i < in.nx; i++){
		for(j = 0; j < in.ny; j++){
			cnt = i*in.ny + j;
			tmp_in[cnt][0] = arr2d_ij(in, i, j);
			tmp_in[cnt][1] = 0;
		}
	}


	fftw_execute(p);

	for(i = 0; i < in.nx; i++){
		for(j = 0; j < in.ny; j++){
			cnt = i*in.ny + j;
			if(i < in.nx/2){
				if(j < in.ny/2) arr2d_set_ij(sqrt(tmp_out[cnt][0]*tmp_out[cnt][0] + tmp_out[cnt][1]*tmp_out[cnt][1]), out, in.nx/2-1-i, in.ny/2-1-j);
				else arr2d_set_ij(sqrt(tmp_out[cnt][0]*tmp_out[cnt][0] + tmp_out[cnt][1]*tmp_out[cnt][1]), out, in.nx/2-1-i, 3*in.ny/2-1-j);
			}
			else{
				if(j < in.ny/2) arr2d_set_ij(sqrt(tmp_out[cnt][0]*tmp_out[cnt][0] + tmp_out[cnt][1]*tmp_out[cnt][1]), out, 3*in.nx/2-1-i, in.ny/2-1-j);
				else arr2d_set_ij(sqrt(tmp_out[cnt][0]*tmp_out[cnt][0] + tmp_out[cnt][1]*tmp_out[cnt][1]), out, 3*in.nx/2-1-i, 3*in.ny/2-1-j);
			}

		}
	}

	fftw_destroy_plan(p);
	fftw_free(tmp_in);
	fftw_free(tmp_out);
}


void process_mplasma(char *bin_data_dir, int num, domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, char *source_spec, FILE *plt, charges_task c_task, int_vector coarse_factor){
	//does all the necessary operations to postprocess meshed plasma information (rho and J)
	partition p;
	int i, j, k, cpu, process_J = 0, point_readable = 0, nx, ny;
	vector J;
	double r = 0, time;
	arr2d rho[mcuts.N], jx[mcuts.N], jy[mcuts.N], jz[mcuts.N];

	p = load_partition_info(bin_data_dir, num, source_spec, &time);
	FILE *fs[2*p.N];
	open_mesh_binary_files(bin_data_dir, source_spec, "rho", num, fs, p.N);
	open_mesh_binary_files(bin_data_dir, source_spec, "J", num, fs+p.N, p.N);

	calc_cut_size(dmn, mmain_plane, mcuts, subdmn, coarse_factor, &nx, &ny);

	if(c_task.jx||c_task.jy||c_task.jz) process_J = 1;
	for(i = 0; i < mcuts.N; i++){
		if(c_task.rho) rho[i] = arr2d_create(nx, ny);
		if(c_task.jx) jx[i] = arr2d_create(nx, ny);
		if(c_task.jy) jy[i] = arr2d_create(nx, ny);
		if(c_task.jz) jz[i] = arr2d_create(nx, ny);
	}

	for(i = dmn.imin; i <= dmn.imax; i++){
		for(j = dmn.jmin; j <= dmn.jmax; j++){
			for(k = dmn.kmin; k <= dmn.kmax; k++){
				if(c_task.rho) r = read_double_point(i, j, k, p, fs);
				if(process_J) J = read_vector_point(i, j, k, p, fs+p.N);
				point_readable = point_kind(i, j, k, dmn, mmain_plane, mcuts, subdmn, coarse_factor);
				if(point_readable){
					if(c_task.rho) savepoint(r, i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, rho);
					if(c_task.jx) savepoint(J.x, i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, jx);
					if(c_task.jy) savepoint(J.y, i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, jy);
					if(c_task.jz) savepoint(J.z, i, j, k, dmn, subdmn, mmain_plane, mcuts, coarse_factor, jz);
				}
			}
		}
	}
	//save data and update plt file
	if(c_task.rho) save_dat_files(dmn, subdmn, mmain_plane, mcuts, rho, num, "rho", plt, time);
	if(c_task.jx) save_dat_files(dmn, subdmn, mmain_plane, mcuts, jx, num, "Jx", plt, time);
	if(c_task.jy) save_dat_files(dmn, subdmn, mmain_plane, mcuts, jy, num, "Jy", plt, time);
	if(c_task.jz) save_dat_files(dmn, subdmn, mmain_plane, mcuts, jz, num, "Jz", plt, time);

	for(cpu = 0; cpu < p.N; cpu++) fclose(fs[cpu]);
	for(i = 0; i < mcuts.N; i++){
		if(c_task.rho) arr2d_destroy(rho[i]);
		if(c_task.jx) arr2d_destroy(jx[i]);
		if(c_task.jy) arr2d_destroy(jy[i]);
		if(c_task.jz) arr2d_destroy(jz[i]);
	}

	partiton_erase(p);
	printf("\t  - plasma loaded\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// SCATTERED DATA VISUALIZATION ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void process_markers(char *bin_data_dir, int num, domain dmn, subdomain subdmn, phase_space_task ps_task, distrib_function_task df_task, int hor_pic_size, int vert_pic_size, long int *n, int *n_species, int *df_success, FILE *plt){
	//creates the phase space plots and processes distribution functions
	partition p;
	marker *m;
	long int N;
	int cpu, i, vps;
	double time;
	double_array q_div_ms_global, q_div_ms_local;
	q_div_ms_global.N = 0;
	q_div_ms_local.N = 0;
	arr2d_char *A;
	arr2d df;

	*df_success = 0;
	*n = 0;
	vps = markers_get_vert_pic_size(ps_task, dmn, subdmn, hor_pic_size);

	p = load_partition_info(bin_data_dir, num, "sys", &time);
	if(df_task.any) df = arr2d_create(df_task.nhist1, df_task.nhist2);

	for(cpu = 0; cpu < p.N; cpu++){
		N = get_markers_num(bin_data_dir, num, cpu);
		if(N > 0){
			if(!(m = calloc((size_t)N, sizeof(marker)))) error(FUNCTION, "cannot allocate %ld kbytes of memory for loading markers from cpu #%d", N*sizeof(marker)/1024, cpu);
			load_plasma(m, bin_data_dir, num, cpu);
			if(ps_task.q_div_m > MAX_Q_DIV_M){
				find_local_qm(m, N, &q_div_ms_local, dmn, subdmn, ps_task);
				if(q_div_ms_local.N > 0){
					update_global_qm(q_div_ms_local, &q_div_ms_global);
					free(q_div_ms_local.data);
				}
				q_div_ms_local.N = 0;
			}
			else{
				q_div_ms_global.N = 1;
				q_div_ms_global.data = calloc(1, sizeof(double));
				q_div_ms_global.data[0] = ps_task.q_div_m;
			}
			if(ps_task.any && (q_div_ms_global.N > 0)){
				if(!ps_task.coarsen) save_markers_direct(m, N, num, dmn, subdmn, n, q_div_ms_global, ps_task);
				else{
					A = calloc(q_div_ms_global.N, sizeof(arr2d_char));
					for(i = 0; i < q_div_ms_global.N; i++){
						A[i] = arr2d_char_create(hor_pic_size, vps);
						carray_fill(A+i, dmn, subdmn, ps_task, q_div_ms_global.data[i], m, N, n);
						save_markers_coarsened(A[i], num, i, dmn, subdmn, ps_task);
						arr2d_char_destroy(A[i]);
					}
					free(A);
				}
			}
			if(df_task.any) update_df(m, N, df, df_task, dmn, subdmn);
			free(m);
		}
	}

	if(*n > 0){
		markers_update_plt(ps_task, dmn, subdmn, hor_pic_size, vert_pic_size, num, q_div_ms_global, time, plt);
		*n_species = q_div_ms_global.N;
	}
	else *n_species = 0;

	if(q_div_ms_global.N > 0) free(q_div_ms_global.data);
	if(df_task.any){
		df_plot(df_task, df, hor_pic_size, vert_pic_size, num, time, plt);
		arr2d_destroy(df);
		*df_success = 1;
	}
}

long int get_markers_num(char *bin_data_dir, int num, int cpu){
	//finds out the number of the saved markers in the checkpoint num for this cpu
	long int N;
	char fn[STRLEN];
	FILE *fp;
	sprintf(fn, "%splasma_%06d_%03d.bin", bin_data_dir, num, cpu);
	OPEN_FILE(fp, fn, "rb");
	fread(&N, sizeof(long int), 1, fp);
	fclose(fp);
	return N;
}

void find_local_qm(marker *m, long int N, double_array *a, domain dmn, subdomain subdmn, phase_space_task ps_task){
	//creates q/m array for the markers loaded from a single cpu. N is supposed to be > 0
	long int i, i0 = 0;
	int j, flag;
	double qm;

	if(N < 1) error(FUNCTION, "N = %d. It should be at least 1", N);
	for(i = 1; i < N; i++){
		if( phase_point_inside_bounds(m[i], dmn, subdmn, ps_task) && marker_is_filtered_by_ef(m[i], ps_task.energy_filter, ps_task.ef) ){
			a->N = 1;
			a->data = calloc(1, sizeof(double));
			a->data[0] = m[i].q_div_m;
			i0 = i;
			break;
		}
	}

	for(i = i0+1; i < N; i++){
		qm = m[i].q_div_m;
		if( phase_point_inside_bounds(m[i], dmn, subdmn, ps_task) && marker_is_filtered_by_ef(m[i], ps_task.energy_filter, ps_task.ef) ){
			flag = 0;
			for(j = 0; j < a->N; j++) if(qm == a->data[j]) flag = 1;
			if(!flag){
				a->N++;
				a->data = realloc(a->data, a->N*sizeof(double));
				assert(a->data);
				a->data[a->N - 1] = qm;
			}
		}
	}
}

void update_global_qm(double_array local, double_array *global){
	//creates or updates global q/m list
	int i, j, flag, i0 = 0;
	if(local.N == 0) error(FUNCTION, "local.N = 0");
	if(global->N == 0){
		global->N = 1;
		global->data = calloc(1, sizeof(double));
		global->data[0] = local.data[0];
		i0 = 1;
	}
	for(i = i0; i < local.N; i++){
		flag = 0;
		for(j = 0; j < global->N; j++) if(local.data[i] == global->data[j]) flag = 1;
		if(!flag){
			global->N++;
			global->data = realloc(global->data, global->N*sizeof(double));
			assert(global->data);
			global->data[global->N - 1] = local.data[i];
		}
	}
}

void save_markers_direct(marker *m, long int N, int num, domain dmn, subdomain subdmn, long int *n, double_array qm, phase_space_task ps_task){
	//directly saves markers into the .dat file and returns the number of the q/m saved
	char fn[STRLEN];
	FILE *fps[qm.N];
	int i;
	long int j;

	for(i = 0; i < qm.N; i++){
		sprintf(fn, "dat/markers_%06d_%03d.dat", num, i);
		OPEN_FILE(fps[i], fn, "at");
		for(j = 0; j < N; j++) if(( (m[j].q_div_m <= qm.data[i]+0.05*fabs(qm.data[i]))&&(m[j].q_div_m >= qm.data[i]-0.05*fabs(qm.data[i])) )&&phase_point_inside_bounds(m[j], dmn, subdmn, ps_task)){
			if(marker_is_filtered_by_ef(m[j], ps_task.energy_filter, ps_task.ef)){
				fprintf(fps[i], "%e %e\n", extract_val(m[j], ps_task, 1), extract_val(m[j], ps_task, 2));
				(*n)++;
			}
		}
		fclose(fps[i]);
	}
}

int phase_point_inside_bounds(marker m, domain dmn, subdomain subdmn, phase_space_task pt){
	//figures out whether the marker is inside the subdomain and specified bounds
	int inside = 0, inside_subdmn = 0;
	if( (pt.x)&&(m.x >= subdmn.imin*dmn.h1)&&(m.x <= subdmn.imax*dmn.h1) ) inside++;
	if( (pt.y)&&(m.y >= subdmn.jmin*dmn.h2)&&(m.y <= subdmn.jmax*dmn.h2) ) inside++;
	if( (pt.z)&&(m.z >= subdmn.kmin*dmn.h3)&&(m.z <= subdmn.kmax*dmn.h3) ) inside++;
	if(dmn.dim == 3){
		if( (m.x >= subdmn.imin*dmn.h1)&&(m.x <= subdmn.imax*dmn.h1) ) inside_subdmn++;
		if( (m.y >= subdmn.jmin*dmn.h2)&&(m.y <= subdmn.jmax*dmn.h2) ) inside_subdmn++;
		if( (m.z >= subdmn.kmin*dmn.h3)&&(m.z <= subdmn.kmax*dmn.h3) ) inside_subdmn++;
	}
	else if(dmn.dim == 2){
		if( ((dmn.mp_2D == 1)||(dmn.mp_2D == 3))&&(m.x >= subdmn.imin*dmn.h1)&&(m.x <= subdmn.imax*dmn.h1) ) inside_subdmn++;
		if( ((dmn.mp_2D == 1)||(dmn.mp_2D == 2))&&(m.y >= subdmn.jmin*dmn.h2)&&(m.y <= subdmn.jmax*dmn.h2) ) inside_subdmn++;
		if( ((dmn.mp_2D == 2)||(dmn.mp_2D == 3))&&(m.z >= subdmn.kmin*dmn.h3)&&(m.z <= subdmn.kmax*dmn.h3) ) inside_subdmn++;
	}
	else{
		if( (dmn.aa_1D == 1)&&(m.x >= subdmn.imin*dmn.h1)&&(m.x <= subdmn.imax*dmn.h1) ) inside_subdmn++;
		if( (dmn.aa_1D == 2)&&(m.y >= subdmn.jmin*dmn.h2)&&(m.y <= subdmn.jmax*dmn.h2) ) inside_subdmn++;
		if( (dmn.aa_1D == 3)&&(m.z >= subdmn.kmin*dmn.h3)&&(m.z <= subdmn.kmax*dmn.h3) ) inside_subdmn++;
	}
	if( (pt.px)&&(m.px >= pt.dpx.min)&&(m.px <= pt.dpx.max) ) inside++;
	if( (pt.py)&&(m.py >= pt.dpy.min)&&(m.py <= pt.dpy.max) ) inside++;
	if( (pt.pz)&&(m.pz >= pt.dpz.min)&&(m.pz <= pt.dpz.max) ) inside++;
	if(inside > 2) error(FUNCTION, "variable inside is %d, must be 2 or less", inside);
	if((inside == 2)&&(inside_subdmn >= dmn.dim)) return 1;
	else return 0;
}

double extract_val(marker m, phase_space_task pt, int index){
	//extracts the appropriate value from the marker
	int i = 0;
	if(pt.x){
		i++;
		if(i == index) return m.x;
	}
	if(pt.y){
		i++;
		if(i == index) return m.y;
	}
	if(pt.z){
		i++;
		if(i == index) return m.z;
	}
	if(pt.px){
		i++;
		if(i == index) return m.px;
	}
	if(pt.py){
		i++;
		if(i == index) return m.py;
	}
	if(pt.pz){
		i++;
		if(i == index) return m.pz;
	}
	return(-1e100);			//dummy return
}

void carray_fill(arr2d_char *A, domain dmn, subdomain subdmn, phase_space_task pt, double qm, marker *m, long int N, long int *n){
	//fills the markers coarsening array 
	long int i;
	interval i1, i2;
	int j = 0, k = 0;
	markers_get_interval(dmn, subdmn, pt, 1, &i1);
	markers_get_interval(dmn, subdmn, pt, 2, &i2);
	for(i = 0; i < N; i++){
		if(( (m[i].q_div_m <= qm+0.05*fabs(qm))&&(m[i].q_div_m >= qm-0.05*fabs(qm)) )&&phase_point_inside_bounds(m[i], dmn, subdmn, pt)){
			if(marker_is_filtered_by_ef(m[i], pt.energy_filter, pt.ef)){
				j = (A->nx)*( (extract_val(m[i], pt, 1) - i1.min)/(i1.max - i1.min) );
				k = (A->ny)*( (extract_val(m[i], pt, 2) - i2.min)/(i2.max - i2.min) );
				(*n)++;
				arr2d_char_set_ij(1, A, j, k);
			}
		}
	}
}

void markers_get_interval(domain dmn, subdomain subdmn, phase_space_task pt, int ind, interval *i){
	//calculates the bounds of the phase space coordinate i
	int j = 0;
	if(pt.x){
		j++;
		if(j == ind){i->min = subdmn.imin*dmn.h1; i->max = subdmn.imax*dmn.h1;}
	}
	if(pt.y){
		j++;
		if(j == ind){i->min = subdmn.jmin*dmn.h2; i->max = subdmn.jmax*dmn.h2;}
	}
	if(pt.z){
		j++;
		if(j == ind){i->min = subdmn.kmin*dmn.h3; i->max = subdmn.kmax*dmn.h3;}
	}
	if(pt.px){
		j++;
		if(j == ind){i->min = pt.dpx.min; i->max = pt.dpx.max;}
	}
	if(pt.py){
		j++;
		if(j == ind){i->min = pt.dpy.min; i->max = pt.dpy.max;}
	}
	if(pt.pz){
		j++;
		if(j == ind){i->min = pt.dpz.min; i->max = pt.dpz.max;}
	}
	i->min -= VERY_SMALL_VALUE;	//to process the degenerate axes correctly
	i->max += VERY_SMALL_VALUE;
}

void save_markers_coarsened(arr2d_char A, int num, int specie, domain dmn, subdomain subdmn, phase_space_task pt){
	//saves the markers from the coarsening array
	char fn[STRLEN];
	FILE *fp;
	interval i1, i2;
	int i, j;
	double dx, dy;

	markers_get_interval(dmn, subdmn, pt, 1, &i1);
	markers_get_interval(dmn, subdmn, pt, 2, &i2);
	dx = (i1.max - i1.min)/A.nx;
	dy = (i2.max - i2.min)/A.ny;

	sprintf(fn, "dat/markers_%06d_%03d.dat", num, specie);
	OPEN_FILE(fp, fn, "at");
	for(i = 0; i < A.nx; i++){
		for(j = 0; j < A.ny; j++){
			if(arr2d_char_ij(A, i, j)) fprintf(fp, "%e %e\n", i1.min + (i+0.5)*dx, i2.min + (j+0.5)*dy);
		}
	}
	fclose(fp);
}

void update_df(marker *m, long int N, arr2d df, distrib_function_task dt, domain dmn, subdomain subdmn){
	//updates distribution function histogram df with the markers information m read from a cpu
	long int i;
	double qm = dt.q_div_m, param1, param2;
	int ind1, ind2;
	for(i = 0; i < N; i++){
		if( ((qm > MAX_Q_DIV_M)||( (m[i].q_div_m <= qm+0.05*fabs(qm))&&(m[i].q_div_m >= qm-0.05*fabs(qm)) )) && (marker_is_filtered_by_ef(m[i], dt.energy_filter, dt.ef)) && marker_inside_subdomain(m[i], dmn, subdmn) ){
			param1 = df_extract_param(m[i], dt, 1);
			if((param1 >= dt.i1.min)&&(param1 < dt.i1.max)){
				ind1 = (int)(dt.nhist1*(param1-dt.i1.min)/(dt.i1.max-dt.i1.min));
				if(dt.sum == 2){
					param2 = df_extract_param(m[i], dt, 2);
					if((param2 >= dt.i2.min)&&(param2 < dt.i2.max)) ind2 = (int)(dt.nhist2*(param2-dt.i2.min)/(dt.i2.max-dt.i2.min));
					else ind2 = -1;
				}
				else ind2 = 0;
			}
			else ind1 = -1;
			if((ind1 >= 0)&&(ind2 >= 0)) arr2d_set_ij(arr2d_ij(df, ind1, ind2)+1, &df, ind1, ind2);
			//if((ind1 >= 0)&&(ind2 >= 0)) arr2d_set_ij(arr2d_ij(df, ind1, ind2)+(m[i].rho), &df, ind1, ind2);
		}
	}
}

double df_extract_param(marker m, distrib_function_task dt, int n){
	//extracts the n-th parameter requested by dt
	int sum = 0;
	if(dt.sum < n) error(FUNCTION, "(dt.sum = %d) < (requested parameter number = %d)", dt.sum, n);
	if(n < 1) error(FUNCTION, "requested parameter number = %d < 1!", n);
	if(dt.x) if(++sum == n) return m.x;
	if(dt.y) if(++sum == n) return m.y;
	if(dt.z) if(++sum == n) return m.z;
	if(dt.px) if(++sum == n) return m.px;
	if(dt.py) if(++sum == n) return m.py;
	if(dt.pz) if(++sum == n) return m.pz;
	if(dt.gamma) if(++sum == n) return(sqrt(1 + m.px*m.px + m.py*m.py + m.pz*m.pz));
//theta_pos_axis is the angle between the axis and the radius-vector of the particles
	if(dt.theta_pos_x) if(++sum == n) return(calc_theta(sqrt(m.y*m.y + m.z*m.z), m.x));
	if(dt.theta_pos_y) if(++sum == n) return(calc_theta(sqrt(m.x*m.x + m.z*m.z), m.y));
	if(dt.theta_pos_z) if(++sum == n) return(calc_theta(sqrt(m.y*m.y + m.x*m.x), m.z));
//theta_mom_axis is the angle between the axis and the momentum of the particles
	if(dt.theta_mom_x) if(++sum == n) return(calc_theta(sqrt(m.py*m.py + m.pz*m.pz), m.px));
	if(dt.theta_mom_y) if(++sum == n) return(calc_theta(sqrt(m.px*m.px + m.pz*m.pz), m.py));
	if(dt.theta_mom_z) if(++sum == n) return(calc_theta(sqrt(m.py*m.py + m.px*m.px), m.pz));
//phi_pos_plane is the angle between the radius-vector of the particles and its projection onto the plane
	if(dt.phi_pos_xy) if(++sum == n) return(calc_theta(m.z, sqrt(m.x*m.x + m.y*m.y)));
	if(dt.phi_pos_yz) if(++sum == n) return(calc_theta(m.x, sqrt(m.y*m.y + m.z*m.z)));
	if(dt.phi_pos_xz) if(++sum == n) return(calc_theta(m.y, sqrt(m.x*m.x + m.z*m.z)));
//phi_mom_plane is the angle between the momentum of the particles and its projection onto the plane
	if(dt.phi_mom_xy) if(++sum == n) return(calc_theta(m.pz, sqrt(m.px*m.px + m.py*m.py)));
	if(dt.phi_mom_yz) if(++sum == n) return(calc_theta(m.px, sqrt(m.py*m.py + m.pz*m.pz)));
	if(dt.phi_mom_xz) if(++sum == n) return(calc_theta(m.py, sqrt(m.px*m.px + m.pz*m.pz)));
	return(-1E100);	//dummy return
}

double calc_theta(double r, double x){
	//returns theta: tg(theta) = r/x;
	if(x == 0){
		if(r == 0) return 0;
		else return(PI/2);
	}
	else{
		if(r == 0) return((x>0)?0:PI);
		else{
			if(x > 0) return atan(r/x);
			else return(PI - atan(fabs(r/x)));
		}
	}
	return 0;	//dummy return
}

int marker_is_filtered_by_ef(marker m, int ef, interval i){
	//returns 1 if energy filter permits this marker and 0 otherwise
	if(!ef) return 1;
	double gamma = sqrt(1 + m.px*m.px + m.py*m.py + m.pz*m.pz);
	if((gamma >= i.min)&&(gamma <= i.max)) return 1;
	return 0;
}

int marker_inside_subdomain(marker m, domain dmn, subdomain subdmn){
	//checks if the marker is inside subdomain
	int inside_subdmn = 0;
	if( (m.x >= subdmn.imin*dmn.h1)&&(m.x <= subdmn.imax*dmn.h1) ) inside_subdmn++;
	if( (m.y >= subdmn.jmin*dmn.h2)&&(m.y <= subdmn.jmax*dmn.h2) ) inside_subdmn++;
	if( (m.z >= subdmn.kmin*dmn.h3)&&(m.z <= subdmn.kmax*dmn.h3) ) inside_subdmn++;
	if(inside_subdmn == 3) return 1;
	else return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// SPECIAL TASKS VISUALIZATION ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void process_energy_diagnostic(char *bin_data_dir, special_task st, FILE *plt, int hor_pic_size){
	//visualizes energy diagnostic
	FILE *fp, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8;
	char str[STRLEN];
	float t, x1, x2, x3, x4, x5, x6, x7, x8;
	
	sprintf(str, "%s/diag.dat", bin_data_dir);
	OPEN_FILE(fp, str, "rt");
	
	if(st.WE) OPEN_FILE(fp1, "dat/we.dat", "w");
	if(st.WM) OPEN_FILE(fp2, "dat/wm.dat", "w");
	if(st.WTx) OPEN_FILE(fp3, "dat/wtx.dat", "w");
	if(st.WTy) OPEN_FILE(fp4, "dat/wty.dat", "w");
	if(st.WTz) OPEN_FILE(fp5, "dat/wtz.dat", "w");
	if(st.WT_perp) OPEN_FILE(fp6, "dat/wtperp.dat", "w");
	if(st.W_kin) OPEN_FILE(fp7, "dat/wkin.dat", "w");
	if(st.WS) OPEN_FILE(fp8, "dat/ws.dat", "w");

	fgets(str, STRLEN, fp);		//skip first two strings
	fgets(str, STRLEN, fp);
	while(fgets(str, STRLEN, fp)){
		float x9;
		sscanf(str, "%f %f %f %f %f %f %f %f %f %f", &t, &x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9);
		if(st.WE) fprintf(fp1, "%e %e\n", t, x1);
		if(st.WM) fprintf(fp2, "%e %e\n", t, x2);
		if(st.WTx) fprintf(fp3, "%e %e\n", t, x3);
		if(st.WTy) fprintf(fp4, "%e %e\n", t, x4);
		if(st.WTz) fprintf(fp5, "%e %e\n", t, x5);
		if(st.WT_perp) fprintf(fp6, "%e %e\n", t, x6);
		if(st.W_kin) fprintf(fp7, "%e %e\n", t, x7);
		if(st.WS) fprintf(fp8, "%e %e\n", t, x9);
	}

	if(st.WE) fclose(fp1);
	if(st.WM) fclose(fp2);
	if(st.WTx) fclose(fp3);
	if(st.WTy) fclose(fp4);
	if(st.WTz) fclose(fp5);
	if(st.WT_perp) fclose(fp6);
	if(st.W_kin) fclose(fp7);
	if(st.WS) fclose(fp8);
	fclose(fp);

	plot_energy_diag(st, hor_pic_size, plt);
}

void vis_probe_locations(domain dmn, special_task st, FILE *plt, int hor_pic_size, int vert_pic_size){
	//visualizes probes locations
	FILE *fp;
	int i;

	OPEN_FILE(fp, "dat/probes_locations.dat", "w");
	for(i = 0; i < st.probes_num_domain; i++){
		if(dmn.dim == 1){
			if(dmn.aa_1D == 1) fprintf(fp, "%e 0\n", st.probe_locations[i].x);
			if(dmn.aa_1D == 2) fprintf(fp, "%e 0\n", st.probe_locations[i].y);
			if(dmn.aa_1D == 3) fprintf(fp, "%e 0\n", st.probe_locations[i].z);
		}
		else if(dmn.dim == 2){
			if(dmn.mp_2D == 1) fprintf(fp, "%e %e\n", st.probe_locations[i].x, st.probe_locations[i].y);
			if(dmn.mp_2D == 2) fprintf(fp, "%e %e\n", st.probe_locations[i].y, st.probe_locations[i].z);
			if(dmn.mp_2D == 3) fprintf(fp, "%e %e\n", st.probe_locations[i].x, st.probe_locations[i].z);
		}
		else{
			fprintf(fp, "%e %e %e\n", st.probe_locations[i].x, st.probe_locations[i].y, st.probe_locations[i].z);
		}
	}
	fclose(fp);
	
	plot_probe_locations(dmn, plt, hor_pic_size, vert_pic_size);
}

void vis_probes(domain dmn, special_task st, int hor_pic_size, FILE *plt, char *bin_data_dir){
	//visualizes information from probes
	int_vector probe_coord;
	int i, j, n_steps = 0, plan_created = 0;
	FILE *fp, *fp_Ex, *fp_Ey, *fp_Ez, *fp_Hx, *fp_Hy, *fp_Hz;
	char str[STRLEN];
	double_array pd;
	fftw_complex *in, *out;
	fftw_plan p;
	float omega;

	if(st.probes_num_domain == 0){
		printf("\t-  no probes found in the domain\n");
		return;
	}

	pd.N = 7;
	if(!(pd.data = calloc(7, sizeof(double)))) error(FUNCTION, "cannot allocate %d bytes of memory for probes data holder", 7*sizeof(double));

#define PROBE_IND(i_, j_, k_) {			\
		probe_coord.i = (int)(i_);	\
		probe_coord.j = (int)(j_);	\
		probe_coord.k = (int)(k_);	\
}
#define OPEN_FILE_FIELD(fp_, num_, spec_){			\
	sprintf(str, "dat/probe_%03d_%s.dat", num_, spec_);	\
	OPEN_FILE(fp_, str, "w");				\
}

	for(i = 0; i < st.probes.N; i++){
		if(st.probe_Ex) OPEN_FILE_FIELD(fp_Ex, st.probes.data[i], "Ex");
		if(st.probe_Ey) OPEN_FILE_FIELD(fp_Ey, st.probes.data[i], "Ey");
		if(st.probe_Ez) OPEN_FILE_FIELD(fp_Ez, st.probes.data[i], "Ez");
		if(st.probe_Hx) OPEN_FILE_FIELD(fp_Hx, st.probes.data[i], "Hx");
		if(st.probe_Hy) OPEN_FILE_FIELD(fp_Hy, st.probes.data[i], "Hy");
		if(st.probe_Hz) OPEN_FILE_FIELD(fp_Hz, st.probes.data[i], "Hz");
		if(dmn.dim == 1){
			if(dmn.aa_1D == 1) PROBE_IND( (st.probe_locations[st.probes.data[i]-1].x)/dmn.h1 + SMALL_VALUE, 0, 0 );
			if(dmn.aa_1D == 2) PROBE_IND( 0, (st.probe_locations[st.probes.data[i]-1].y)/dmn.h2 + SMALL_VALUE, 0 );
			if(dmn.aa_1D == 3) PROBE_IND( 0, 0, (st.probe_locations[st.probes.data[i]-1].z)/dmn.h3 + SMALL_VALUE );

		}
		else if(dmn.dim == 2){
			if(dmn.mp_2D == 1) PROBE_IND( (st.probe_locations[st.probes.data[i]-1].x)/dmn.h1 + SMALL_VALUE, (st.probe_locations[st.probes.data[i]-1].y)/dmn.h2 + SMALL_VALUE, 0 );
			if(dmn.mp_2D == 2) PROBE_IND( 0, (st.probe_locations[st.probes.data[i]-1].y)/dmn.h2 + SMALL_VALUE, (st.probe_locations[st.probes.data[i]-1].z)/dmn.h3 + SMALL_VALUE );
			if(dmn.mp_2D == 3) PROBE_IND( (st.probe_locations[st.probes.data[i]-1].x)/dmn.h1 + SMALL_VALUE, 0, (st.probe_locations[st.probes.data[i]-1].z)/dmn.h3 + SMALL_VALUE );
		}
		else{
			PROBE_IND((st.probe_locations[st.probes.data[i]-1].x)/dmn.h1 + SMALL_VALUE, (st.probe_locations[st.probes.data[i]-1].y)/dmn.h2 + SMALL_VALUE, (st.probe_locations[st.probes.data[i]-1].z)/dmn.h3 + SMALL_VALUE);
		}
		sprintf(str, "%sprobe_%d_%d_%d.bin", bin_data_dir, probe_coord.i, probe_coord.j, probe_coord.k);
		OPEN_FILE(fp, str, "rb");
		fseek(fp, (long)(8*sizeof(char) + 3*sizeof(int)), SEEK_SET);

		//preparing data structures for Fourier transformation in time
		if(!plan_created){
			//first, lets get the number of timesteps
			n_steps = 0;
			while(!feof(fp)){
				fread(pd.data, sizeof(double), 7, fp);
				n_steps++;
			}
			fseek(fp, (long)(8*sizeof(char) + 3*sizeof(int)), SEEK_SET);

			if(!(in = (fftw_complex*) fftw_malloc(n_steps*sizeof(fftw_complex)))) error(FUNCTION, "cannot allocate memory for fourier transformation");
			if(!(out = (fftw_complex*) fftw_malloc(n_steps*sizeof(fftw_complex)))) error(FUNCTION, "cannot allocate memory for fourier transformation");
			p = fftw_plan_dft_1d(n_steps, in, out, FFTW_FORWARD, FFTW_MEASURE);
			plan_created = 1;
		}

		//now reading the data
		j = 0;
		while(!feof(fp)){
			fread(pd.data, sizeof(double), 7, fp);
			if(st.probe_Ex) fprintf(fp_Ex, "%e %e\n", pd.data[0], pd.data[1]/(2*PI));
			if(st.probe_Ey) fprintf(fp_Ey, "%e %e\n", pd.data[0], pd.data[2]/(2*PI));
			if(st.probe_Ez) fprintf(fp_Ez, "%e %e\n", pd.data[0], pd.data[3]/(2*PI));
			if(st.probe_Hx) fprintf(fp_Hx, "%e %e\n", pd.data[0], pd.data[4]/(2*PI));
			if(st.probe_Hy) fprintf(fp_Hy, "%e %e\n", pd.data[0], pd.data[5]/(2*PI));
			if(st.probe_Hz) fprintf(fp_Hz, "%e %e\n", pd.data[0], pd.data[6]/(2*PI));
			in[j][0] = sqrt( pd.data[1]*pd.data[1] + pd.data[2]*pd.data[2] + pd.data[3]*pd.data[3] + pd.data[4]*pd.data[4] + pd.data[5]*pd.data[5] + pd.data[6]*pd.data[6] ) / (16*PI*PI*PI);
			in[j][1] = 0;
			j++;
		}
		fclose(fp);
		if(st.probe_Ex) fclose(fp_Ex);
		if(st.probe_Ey) fclose(fp_Ey);
		if(st.probe_Ez) fclose(fp_Ez);
		if(st.probe_Hx) fclose(fp_Hx);
		if(st.probe_Hy) fclose(fp_Hy);
		if(st.probe_Hz) fclose(fp_Hz);

		fftw_execute(p);
		OPEN_FILE_FIELD(fp, st.probes.data[i], "spectr");
		for(j = 0; j < n_steps; j++){
			omega = 1./dmn.tau;
			omega *= ((float)(j-n_steps/2))/n_steps;
			if(j <= n_steps/2) fprintf(fp, "%e %e\n", omega, sqrt( out[n_steps/2 - j][0]*out[n_steps/2 - j][0] + out[n_steps/2 - j][1]*out[n_steps/2 - j][1] ));
			else fprintf(fp, "%e %e\n", omega, sqrt( out[3*n_steps/2 - j][0]*out[3*n_steps/2 - j][0] + out[3*n_steps/2 - j][1]*out[3*n_steps/2 - j][1] ));
		}
		fclose(fp);
	}
	pd.N = 0;
	free(pd.data);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);

	plot_probes_data(st, hor_pic_size, plt);
}





