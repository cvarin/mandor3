//routines for final plotting

#include<stdio.h>
#include<string.h>
#include<math.h>
#include"utils.h"
#include"types.h"
#include"params.h"
#include"plot.h"

int plot_init(FILE **plt, int hor_pic_size, domain dmn, subdomain subdmn, int mmain_plane){
	//initializes .plt file
	int vert_pic_size = 0;

	OPEN_FILE(*plt, "plot.plt", "wt");
	fprintf(*plt, "#gnuplot plot file\n");

	if(dmn.dim == 1) vert_pic_size = hor_pic_size*0.75;
	if(dmn.dim == 2){
		if(dmn.mp_2D == 1) vert_pic_size = (hor_pic_size*(subdmn.jmax - subdmn.jmin)*dmn.h2)/((subdmn.imax - subdmn.imin)*dmn.h1) + VERT_EXTRA;
		if(dmn.mp_2D == 2) vert_pic_size = (hor_pic_size*(subdmn.kmax - subdmn.kmin)*dmn.h3)/((subdmn.jmax - subdmn.jmin)*dmn.h2) + VERT_EXTRA;
		if(dmn.mp_2D == 3) vert_pic_size = (hor_pic_size*(subdmn.kmax - subdmn.kmin)*dmn.h3)/((subdmn.imax - subdmn.imin)*dmn.h1) + VERT_EXTRA;
	}
	if(dmn.dim == 3){
		if(mmain_plane == 1) vert_pic_size = (hor_pic_size*(subdmn.jmax - subdmn.jmin)*dmn.h2)/((subdmn.imax - subdmn.imin)*dmn.h1) + VERT_EXTRA;
		if(mmain_plane == 2) vert_pic_size = (hor_pic_size*(subdmn.kmax - subdmn.kmin)*dmn.h3)/((subdmn.jmax - subdmn.jmin)*dmn.h2) + VERT_EXTRA;
		if(mmain_plane == 3) vert_pic_size = (hor_pic_size*(subdmn.kmax - subdmn.kmin)*dmn.h3)/((subdmn.imax - subdmn.imin)*dmn.h1) + VERT_EXTRA;
	}
	fprintf(*plt, "set term png large size %d,%d\n", hor_pic_size, vert_pic_size);
	return vert_pic_size;
}

void markers_update_plt(phase_space_task pt, domain dmn, subdomain subdmn, int hor_pic_size, int vert_pic_size, int num, double_array qm, double time, FILE *plt){
	//plots markers
	int vps, ind = 0, i;
	double minx = 0, maxx = 0, miny = 0, maxy = 0;
	char xlabel[9], ylabel[9];
	int indices[qm.N];

	get_sorted_indices(indices, qm.N, qm.data);

	vps = markers_get_vert_pic_size(pt, dmn, subdmn, hor_pic_size);

#define GET_MAX_MINS(min_, max_, label_){										\
	if(!ind){minx = min_ - VERY_SMALL_VALUE; maxx = max_ + VERY_SMALL_VALUE; ind++; strcpy(xlabel, label_);}	\
	else{miny = min_ - VERY_SMALL_VALUE; maxy = max_ + VERY_SMALL_VALUE; ind++; strcpy(ylabel, label_);}		\
}
	if(pt.x) GET_MAX_MINS(subdmn.imin*dmn.h1, subdmn.imax*dmn.h1, "x/lambda");
	if(pt.y) GET_MAX_MINS(subdmn.jmin*dmn.h2, subdmn.jmax*dmn.h2, "y/lambda");
	if(pt.z) GET_MAX_MINS(subdmn.kmin*dmn.h3, subdmn.kmax*dmn.h3, "z/lambda");
	if(pt.px) GET_MAX_MINS(pt.dpx.min, pt.dpx.max, "p_x/(mc)");
	if(pt.py) GET_MAX_MINS(pt.dpy.min, pt.dpy.max, "p_y/(mc)");
	if(pt.pz) GET_MAX_MINS(pt.dpz.min, pt.dpz.max, "p_z/(mc)");

	if(vps != vert_pic_size) fprintf(plt, "set term png large size %d, %d\n", hor_pic_size, vps);
	fprintf(plt, "set output 'pics/markers_%06d.png'\n", num);
	fprintf(plt, "set title 'Markers distribution in the phase subspace: t = %.2f waveperiods'\n", time);
	fprintf(plt, "set xlabel '%s'\n", xlabel);
	fprintf(plt, "set ylabel '%s'\n", ylabel);
	fprintf(plt, "set grid\n");
	fprintf(plt, "plot [%e:%e][%e:%e] ", minx, maxx, miny, maxy);
	for(i = 0; i < qm.N; i++){
		fprintf(plt, "'dat/markers_%06d_%03d.dat' w d title 'q/m = %.3e'", num, indices[i], qm.data[indices[i]]);
		if(i != (qm.N-1)) fprintf(plt, ", ");
	}
	fprintf(plt, "\n");
	fprintf(plt, "unset grid\n");

	if(vps != vert_pic_size) fprintf(plt, "set term png large size %d, %d\n", hor_pic_size, vert_pic_size);
}

void get_sorted_indices(int *indices, int n, double *vals){
	//fills indices array with indices of a sorted array vals -- the procedure is used in the markers plotter
	double arr[n], a;
	int i, j;
	for(i = 0; i < n; i++) arr[i] = vals[i];
	for(j = 1; j < n; j++){
		a = arr[j];
		i = j-1;
		while((i >= 0)&&(arr[i] > a)){
			arr[i+1] = arr[i];
			i--;
		}
		arr[i+1] = a;
	}
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(vals[i] == arr[j]){
				indices[i] = j;
				break;
			}
		}
	}
}

int markers_get_vert_pic_size(phase_space_task pt, domain dmn, subdomain subdmn, int hor_pic_size){
	//calculates and returns vertical picture size for markers
	int vps = 0.67*hor_pic_size;
	if((pt.x + pt.y + pt.z) == 2){
		if(dmn.dim == 3){
			if(pt.x){
				if(pt.y) vps = hor_pic_size*( dmn.h2*(subdmn.jmax-subdmn.jmin)/(dmn.h1*(subdmn.imax-subdmn.imin)) );
				if(pt.z) vps = hor_pic_size*( dmn.h3*(subdmn.kmax-subdmn.kmin)/(dmn.h1*(subdmn.imax-subdmn.imin)) );
			}
			else vps = hor_pic_size*( dmn.h3*(subdmn.kmax-subdmn.kmin)/(dmn.h2*(subdmn.jmax-subdmn.jmin)) );
		}
		else if(dmn.dim == 2){
			if((dmn.mp_2D == 1)&&pt.x&&pt.y) vps = hor_pic_size*( dmn.h2*(subdmn.jmax-subdmn.jmin)/(dmn.h1*(subdmn.imax-subdmn.imin)) );
			if((dmn.mp_2D == 2)&&pt.y&&pt.z) vps = hor_pic_size*( dmn.h3*(subdmn.kmax-subdmn.kmin)/(dmn.h2*(subdmn.jmax-subdmn.jmin)) );
			if((dmn.mp_2D == 3)&&pt.x&&pt.z) vps = hor_pic_size*( dmn.h3*(subdmn.kmax-subdmn.kmin)/(dmn.h1*(subdmn.imax-subdmn.imin)) );
		}
	}
	return (vps + VERT_EXTRA);
}

void df_plot(distrib_function_task dt, arr2d df, int hor_pic_size, int vert_pic_size, int num, double time, FILE *plt){
	//saves distribution function data and updates plot file
	FILE *fp;
	char str[STRLEN], xlabel[32], ylabel[32];
	int i, j, sum = 0;
	double x, y;

	sprintf(str, "dat/df_%06d.dat", num);
	OPEN_FILE(fp, str, "w");
	if((df.nx > 1)&&(df.ny > 1)) arr2d_set_ij(arr2d_ij(df, 0, 0) + 1, &df, 0, 0);
	for(i = 0; i < dt.nhist1; i++){
		x = dt.i1.min + (i+0.5)*(dt.i1.max-dt.i1.min)/dt.nhist1;
		if(dt.sum == 2){
			for(j = 0; j < dt.nhist2; j++){
				y = dt.i2.min + (j+0.5)*(dt.i2.max-dt.i2.min)/dt.nhist2;
				if(!dt.log_distr) fprintf(fp, "%e %e %.0f\n", x, y, arr2d_ij(df, i, j));
				else fprintf(fp, "%e %e %e\n", x, y, log(1 + arr2d_ij(df, i, j)));
			}
			fprintf(fp, "\n");
		}
		else{
			if(!dt.log_distr) fprintf(fp, "%e %.0f\n", x, arr2d_ij(df, i, 0));
			else fprintf(fp, "%e %e\n", x, log(1 + arr2d_ij(df, i, 0)));
		}
	}
	fclose(fp);

	fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, (int)(hor_pic_size*0.75));
#define MKLABEL(label_) {				\
	if(sum == 0) sprintf(xlabel, "%s", label_);	\
	else sprintf(ylabel, "%s", label_);		\
	sum++;						\
}
	if(dt.x) MKLABEL("x/lambda");
	if(dt.y) MKLABEL("y/lambda");
	if(dt.z) MKLABEL("z/lambda");
	if(dt.px) MKLABEL("p_x/mc");
	if(dt.py) MKLABEL("p_y/mc");
	if(dt.pz) MKLABEL("p_z/mc");
	if(dt.gamma) MKLABEL("particles gamma factor");
	if(dt.theta_pos_x) MKLABEL("theta_pos_x");
	if(dt.theta_pos_y) MKLABEL("theta_pos_y");
	if(dt.theta_pos_z) MKLABEL("theta_pos_z");
	if(dt.theta_mom_x) MKLABEL("theta_mom_x");
	if(dt.theta_mom_y) MKLABEL("theta_mom_y");
	if(dt.theta_mom_z) MKLABEL("theta_mom_z");
	if(dt.phi_pos_xy) MKLABEL("phi_pos_xy");
	if(dt.phi_pos_yz) MKLABEL("phi_pos_yz");
	if(dt.phi_pos_xz) MKLABEL("phi_pos_xz");
	if(dt.phi_mom_xy) MKLABEL("phi_mom_xy");
	if(dt.phi_mom_yz) MKLABEL("phi_mom_yz");
	if(dt.phi_mom_xz) MKLABEL("phi_mom_xz");
	if(dt.sum == 1){
		if(!dt.log_distr) sprintf(ylabel, "delta N");
		else sprintf(ylabel, "ln ( 1 + delta N )");
	}

	fprintf(plt, "set output 'pics/df_%06d.png'\n", num);
	sprintf(str, "q/m = %f", dt.q_div_m);
	fprintf(plt, "set title 'Particle %shistogram", dt.log_distr?"logarithmic ":"");
	if((dt.q_div_m <= MAX_Q_DIV_M)||(dt.energy_filter)) fprintf(plt, " for particles having ");
	fprintf(plt, "%s", (dt.q_div_m > MAX_Q_DIV_M)?"":str);
	sprintf(str, "gamma in [%.3e, %.3e]", dt.ef.min, dt.ef.max);
	fprintf(plt, "%s%s", ((dt.q_div_m <= MAX_Q_DIV_M)?" and ":""), (dt.energy_filter?str:""));
	fprintf(plt, ": t = %.2f'\n", time);
	fprintf(plt, "set xlabel '%s'\n", xlabel);
	fprintf(plt, "set ylabel '%s'\n", ylabel);
	if(dt.sum == 2){
		fprintf(plt, "set pm3d map\n");
		fprintf(plt, "splot 'dat/df_%06d.dat' using 1:2:3 notitle\n", num);
		fprintf(plt, "unset pm3d\n");
	}
	else{
		fprintf(plt, "set grid\n");
		fprintf(plt, "plot 'dat/df_%06d.dat' w l notitle\n", num);
		fprintf(plt, "unset grid\n");
	}
	fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, vert_pic_size);
}

void plot_energy_diag(special_task st, int hor_pic_size, FILE *plt){
	//plots energy diagnostic
	int nprinted = 0;
	
	fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, (int)(0.75*hor_pic_size));
	fprintf(plt, "set output 'pics/energy_daig.png'\n");
	fprintf(plt, "set title 'energy diagnostic'\n");
	fprintf(plt, "set xlabel 'time, waveperiods'\n");
	fprintf(plt, "set ylabel 'energy/w0, w0 = m_e^2 c^2 omega_0^2 / ( 16 pi^3 e^2 )'\n");
	fprintf(plt, "set grid\n");
	fprintf(plt, "plot ");
#define APPEND_PLT(datfile_, title_) {				\
	fprintf(plt, "'%s' w l t '%s'", datfile_, title_);	\
	nprinted++;						\
	if(nprinted < st.e_diag) fprintf(plt, ", ");		\
	else fprintf(plt, "\n");				\
}
	if(st.WE) APPEND_PLT("dat/we.dat", "Electric energy");
	if(st.WM) APPEND_PLT("dat/wm.dat", "Magnetic energy");
	if(st.WTx) APPEND_PLT("dat/wtx.dat", "Particles T_x");
	if(st.WTy) APPEND_PLT("dat/wty.dat", "Particles T_y");
	if(st.WTz) APPEND_PLT("dat/wtz.dat", "Particles T_z");
	if(st.WT_perp) APPEND_PLT("dat/wtperp.dat", "Particles T_perp");
	if(st.W_kin) APPEND_PLT("dat/wkin.dat", "Particles kinetic energy");
	if(st.WS) APPEND_PLT("dat/ws.dat", "Total energy");
	fprintf(plt, "unset grid\n");
}

void plot_probe_locations(domain dmn, FILE *plt, int hor_pic_size, int vert_pic_size){
	//plots the probes locations
#define PL_APPEND_PLT_1D(label_, L_) {fprintf(plt, "set xlabel '%s'\n", label_); fprintf(plt, "plot [0:%e] 'dat/probes_locations.dat' w p pt 18 notitle\n", L_);}
#define PL_APPEND_PLT_2D(labelx_, labely_, Lx_, Ly_) {							\
	fprintf(plt, "set xlabel '%s'\n", labelx_);							\
	fprintf(plt, "set ylabel '%s'\n", labely_);							\
	fprintf(plt, "plot [0:%e][0:%e] 'dat/probes_locations.dat' w p pt 18 notitle\n", Lx_, Ly_);	\
}

	if(dmn.dim != 2) fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, (int)(0.75*hor_pic_size));
	else  fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, vert_pic_size);
	fprintf(plt, "set output 'pics/probes_locations.png'\n");
	fprintf(plt, "set title 'Probes locations'\n");
	fprintf(plt, "set grid\n");
	fprintf(plt, "unset pm3d\n");
	fprintf(plt, "set view 60,30,1,1\n");
	if(dmn.dim == 1){
		if(dmn.aa_1D == 1) PL_APPEND_PLT_1D("x/lambda", dmn.Lx);
		if(dmn.aa_1D == 2) PL_APPEND_PLT_1D("y/lambda", dmn.Ly);
		if(dmn.aa_1D == 3) PL_APPEND_PLT_1D("z/lambda", dmn.Lz);
	}
	else if(dmn.dim == 2){
		if(dmn.mp_2D == 1) PL_APPEND_PLT_2D("x/lambda", "y/lambda", dmn.Lx, dmn.Ly);
		if(dmn.mp_2D == 2) PL_APPEND_PLT_2D("y/lambda", "z/lambda", dmn.Ly, dmn.Lz);
		if(dmn.mp_2D == 3) PL_APPEND_PLT_2D("x/lambda", "z/lambda", dmn.Lx, dmn.Lz);
	}
	else{
		fprintf(plt, "set xlabel 'x/lambda'\n");
		fprintf(plt, "set ylabel 'y/lambda'\n");
		fprintf(plt, "set zlabel 'z/lambda'\n");
		fprintf(plt, "splot [0:%e][0:%e][0:%e] 'dat/probes_locations.dat' w p pt 18 notitle\n", dmn.Lx, dmn.Ly, dmn.Lz);
	}
	fprintf(plt, "unset grid\n");
}

void plot_probes_data(special_task st, int hor_pic_size, FILE *plt){
	//plots probes data
	int i;

	fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, (int)(0.75*hor_pic_size + ( (st.probes.N > 4) ? st.probes.N*3 : 0 )));
	fprintf(plt, "set grid\n");
	if(st.probes.N > 4) fprintf(plt, "set key below\n");
	fprintf(plt, "set xlabel 'time (waveperiods)'\n");

#define PLOT_COMPONENT(E_, title_) {				\
	fprintf(plt, "set output 'pics/probes_%s.png'\n", E_);	\
	fprintf(plt, "set title '%s from probes'\n", title_);	\
	fprintf(plt, "set ylabel '%s'\n", E_);			\
	fprintf(plt, "plot");					\
	for(i = 0; i < st.probes.N; i++){			\
		fprintf(plt, " 'dat/probe_%03d_%s.dat' w l t 'probe #%d at [%.2f, %.2f, %.2f]'", st.probes.data[i], E_, st.probes.data[i], st.probe_locations[st.probes.data[i]-1].x, st.probe_locations[st.probes.data[i]-1].y, st.probe_locations[st.probes.data[i]-1].z);\
		if(i < (st.probes.N - 1)) fprintf(plt, ", ");	\
		else fprintf(plt, "\n");			\
	}							\
}

	if(st.probe_Ex) PLOT_COMPONENT("Ex", "Ex");
	if(st.probe_Ey) PLOT_COMPONENT("Ey", "Ey");
	if(st.probe_Ez) PLOT_COMPONENT("Ez", "Ez");
	if(st.probe_Hx) PLOT_COMPONENT("Hx", "Hx");
	if(st.probe_Hy) PLOT_COMPONENT("Hy", "Hy");
	if(st.probe_Hz) PLOT_COMPONENT("Hz", "Hz");
	fprintf(plt, "set xlabel 'omega/omega_0'\n");
	PLOT_COMPONENT("spectr", "Spectral energy density");

	fprintf(plt, "unset grid\n");
	fprintf(plt, "set key default\n");

}

void plot_fourier_transform(domain dmn, special_task st, int mmain_plane, int_array mcuts, arr2d A, int i_slice, int num, FILE *plt, int hor_pic_size, int vert_pic_size, double time){
	//saves a number of cuts for fourier components into dat-files using procedures save_dat_file
	int i, j, flag;
	double kxmin = 0, kymin = 0, kxmax = 0, kymax = 0, kx, ky, kxmin_bound = 0, kxmax_bound = 0, kymin_bound = 0, kymax_bound = 0;
	double z = 0;
	char str[STRLEN];
	char xlabel[16], ylabel[16], zlabel[2];
	FILE *fp;

#define BOUNDS_F(kxmin_, kxmax_, kymin_, kymax_, hx_, hy_, hz_, xlabel_, ylabel_, zlabel_){	\
	kxmin = -1./(2*hx_);									\
	kxmax = 1./(2*hx_);									\
	kymin = -1./(2*hy_);									\
	kymax = 1./(2*hy_);									\
	if(st.enable_intervals){								\
		kxmin_bound = kxmin_;								\
		kxmax_bound = kxmax_;								\
		kymin_bound = kymin_;								\
		kymax_bound = kymax_;								\
	}											\
	else{											\
		kxmin_bound = kxmin;								\
		kxmax_bound = kxmax;								\
		kymin_bound = kymin;								\
		kymax_bound = kymax;								\
	}											\
	z = hz_*mcuts.data[i_slice];								\
	sprintf(xlabel, "%s", xlabel_);								\
	sprintf(ylabel, "%s", ylabel_);								\
	sprintf(zlabel, "%s", zlabel_);								\
}
	if(dmn.dim == 3){
		if(mmain_plane == 1) BOUNDS_F(st.ix.min, st.ix.max, st.iy.min, st.iy.max, dmn.h1, dmn.h2, dmn.h3, "lambda k_x", "lambda k_y", "z");
		if(mmain_plane == 2) BOUNDS_F(st.iy.min, st.iy.max, st.iz.min, st.iz.max, dmn.h2, dmn.h3, dmn.h1, "lambda k_y", "lambda k_z", "x");
		if(mmain_plane == 3) BOUNDS_F(st.ix.min, st.ix.max, st.iz.min, st.iz.max, dmn.h1, dmn.h3, dmn.h2, "lambda k_x", "lambda k_z", "y");
	}
	else if(dmn.dim == 2){
		if(dmn.mp_2D == 1) BOUNDS_F(st.ix.min, st.ix.max, st.iy.min, st.iy.max, dmn.h1, dmn.h2, 0, "lambda k_x", "lambda k_y", "");
		if(dmn.mp_2D == 2) BOUNDS_F(st.iy.min, st.iy.max, st.iz.min, st.iz.max, dmn.h2, dmn.h3, 0, "lambda k_y", "lambda k_z", "");
		if(dmn.mp_2D == 3) BOUNDS_F(st.ix.min, st.ix.max, st.iz.min, st.iz.max, dmn.h1, dmn.h3, 0, "lambda k_x", "lambda k_z", "");
	}
	else{
		if(dmn.aa_1D == 1) BOUNDS_F(st.ix.min, st.ix.max, 0, 1, dmn.h1, 1, 0, "lambda k_x", "", "");
		if(dmn.aa_1D == 2) BOUNDS_F(st.iy.min, st.iy.max, 0, 1, dmn.h2, 1, 0, "lambda k_y", "", "");
		if(dmn.aa_1D == 3) BOUNDS_F(st.iz.min, st.iz.max, 0, 1, dmn.h3, 1, 0, "lambda k_z", "", "");
	}
	sprintf(str, "dat/spectr_%06d(%03d).dat", num, i_slice);
	OPEN_FILE(fp, str, "w");
	if(dmn. dim > 1){
		flag = 0;
		for(i = 0; i < A.nx; i++){
			for(j = 0; j < A.ny; j++){
				kx = kxmax*(i+0.5-A.nx/2)*2/A.nx;
				ky = kymax*(j+0.5-A.ny/2)*2/A.ny;
				if(in(kx, kxmin_bound, kxmax_bound) && in(ky, kymin_bound, kymax_bound)){
					if(flag == 0){
						flag = 1;
						arr2d_set_ij(arr2d_ij(A, i, j) + VERY_SMALL_VALUE, &A, i, j);	//for gnuplot
					}
					fprintf(fp, "%e %e %e\n", kx, ky, arr2d_ij(A, i, j));
				}
			}
			fprintf(fp, "\n");
		}
	}
	else{
		for(i = 0; i < A.nx; i++){
			kx = kxmax*(i+0.5-A.nx/2)*2/A.nx;
			if(in(kx, kxmin_bound, kxmax_bound)) fprintf(fp, "%e %e\n", kx, arr2d_ij(A, i, 0));
		}
	}
	fclose(fp);

	fprintf(plt, "set output 'pics/spectr_%06d(%03d).png'\n", num, i_slice);
	if(dmn.dim == 3) fprintf(plt, "set title 'energy spectral density: cut at %s = %.2f lambdas, time = %.2f waveperiods'\n", zlabel, z, time);
	else fprintf(plt, "set title 'energy spectral density: time = %.2f waveperiods'\n", time);
	fprintf(plt, "set xlabel '%s'\n", xlabel);
	fprintf(plt, "set ylabel '%s'\n", ylabel);
	if(dmn.dim > 1){
		fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, (int)((hor_pic_size*kymax_bound)/kxmax_bound + VERT_EXTRA));
		fprintf(plt, "set pm3d map\n");
		fprintf(plt, "splot 'dat/spectr_%06d(%03d).dat' notitle\n", num, i_slice);
	}
	else{
		fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, (int)(hor_pic_size*0.75));
		fprintf(plt, "set grid\n");
		fprintf(plt, "plot 'dat/spectr_%06d(%03d).dat' w l notitle\n", num, i_slice);
		fprintf(plt, "unset grid\n");
	}
	fprintf(plt, "set term png large size %d,%d\n", hor_pic_size, vert_pic_size);
}










