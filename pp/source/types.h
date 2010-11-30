// types used

#ifndef TYPES_LOADED
	#define TYPES_LOADED 1
	
	typedef struct{
		int N;
		int *data;
	} int_array;

	typedef struct{
		int N;
		double *data;
	} double_array;

	typedef struct{
		int i;
		int j;
		int k;
	} int_vector;

	typedef struct{
		float x;
		float y;
		float z;
	} float_vector;


	typedef struct{
		double h1;
		double h2;
		double h3;
		double Lx;
		double Ly;
		double Lz;
		double tau;
		int imin;
		int jmin;
		int kmin;
		int imax;
		int jmax;
		int kmax;
		int dim;
		int mp_2D;	//main plane for 2D: 1 -- XY, 2 -- YZ, 3 -- XZ
		int aa_1D;	//active axis for 1D: 1 -- X, 2 -- Y, 3 -- Z
	} domain;

	typedef struct{
		int enabled;
		int imin;
		int jmin;
		int kmin;
		int imax;
		int jmax;
		int kmax;
	} subdomain;

	typedef struct{
		int ex;
		int ey;
		int ez;
		int hx;
		int hy;
		int hz;
		int I;
		int any;	//shows if there is at least one nonzero
	} fields_task;

	typedef struct{
		int rho;
		int jx;
		int jy;
		int jz;
		int any;	//shows if there is at least one nonzero
	} charges_task;

	typedef struct{
		int imin;
		int jmin;
		int kmin;
		int imax;
		int jmax;
		int kmax;
	} cpudmn;

	typedef struct{
		int N;
		cpudmn *nodes;
	} partition;

	typedef struct{
		double x;
		double y;
		double z;
	} vector;

	typedef struct{
		double min;
		double max;
	} interval;

	typedef struct{
		double q_div_m;
		int energy_filter;
		interval ef;
		int x;
		int y;
		int z;
		int px;
		int py;
		int pz;
		interval dpx;
		interval dpy;
		interval dpz;
		int coarsen;
		int any;
	} phase_space_task;

	typedef struct{
		double q_div_m;
		int energy_filter;
		interval ef;
		int x;
		int y;
		int z;
		int px;
		int py;
		int pz;
		int gamma;
		int theta_pos_x;		//theta_pos_axis is the angle between the axis and the radius-vector of the particles
		int theta_pos_y;
		int theta_pos_z;
		int theta_mom_x;		//theta_mom_axis is the angle between the axis and the momentum of the particles
		int theta_mom_y;
		int theta_mom_z;
		int phi_pos_xy;			//phi_pos_plane is the angle between the radius-vector of the particles and its projection onto the plane
		int phi_pos_yz;
		int phi_pos_xz;
		int phi_mom_xy;			//phi_mom_plane is the angle between the momentum of the particles and its projection onto the plane
		int phi_mom_yz;
		int phi_mom_xz;
		int nhist1;
		int nhist2;
		interval i1;
		interval i2;
		int log_distr;
		int any;
		int sum;			//dimensionality of the histogram to be created
	} distrib_function_task;

	typedef struct{
	        double x;
	        double y;
	        double z;
	        double px;
	        double py;
	        double pz;
	        float rho;
	        float q_div_m;
	} marker;

	typedef struct{
		int WE;
		int WM;
		int WTx;
		int WTy;
		int WTz;
		int WT_perp;
		int W_kin;
		int WS;
		int e_diag;

		int vis_probes;			//flag turning on the diagnostic
		int_array probes;		//array containing probes numbers for visualization
		int probes_num_domain;		//total number of probes in the domain
		float_vector *probe_locations;	//locations of the probes inside the entire domain
		int vis_probe_locations;	//whether visualize locations or not
		int probe_Ex;			//components of the field to visualize
		int probe_Ey;
		int probe_Ez;
		int probe_Hx;
		int probe_Hy;
		int probe_Hz;

		int vis_fourier;		//flag turning the fourier diagnostic on
		int enable_intervals;		//whether the fourier components inside some explicitly specified interval is to be visualized or just inside the natural interval = [0; 1/h]
		interval ix;			//the intervals
		interval iy;
		interval iz;

		int vis_gamma; // TEMP
	} special_task;


#endif


