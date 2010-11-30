//the Mandor postprocessor reads Mandor core.out generated binary files (checkpoints and plotting snapshots)
//and uses the information to create fields/currents/markers distribution/rho/energy graphs and so on

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"io_utils.h"
#include"utils.h"
#include"params.h"
#include"types.h"
#include"cfgreader.h"
#include"data_loader.h"
#include"report.h"
#include"processor.h"
#include"plot.h"

int main(){
	
	char bin_data_dir[STRLEN];
	int hor_pic_size, vert_pic_size, coarsen_data;
	int_array checkpoints_list, snapshots_list;
	checkpoints_list.N = 0, snapshots_list.N = 0;

	int field_source;
	int mmain_plane, mslicing_method;
	int_array mcuts;
	mcuts.N = 1;				//it's necessary for correct work with dimensionality < 3 (although mcuts.data[0] contains garbage)

	fields_task f_task;
	charges_task c_task;
	phase_space_task ps_task;
	distrib_function_task df_task;
	special_task s_task;

	domain dmn;
	subdomain subdmn;
	int i, df_success = 0, n_species;
	long int n;
	int_vector coarse_factor;

	FILE *plt;
	
	printf("Reading postprocessor parameters and Mandor binary data configuration...\n");

	//preparing folders
	preparefolders();

	//reading parameters from cfg-file
	read_cfg_common(bin_data_dir, &checkpoints_list, &snapshots_list, &hor_pic_size, &subdmn, &coarsen_data);
	read_cfg_mdata(&field_source, &f_task, &c_task);
	read_cfg_markers(&ps_task);
	read_cfg_distrib(&df_task);

	//gathering domain information and performing automatic slicing if necessary
	load_domain_info(&dmn, bin_data_dir);
	read_cfg_special_task(bin_data_dir, &s_task, dmn);
	report_domain(dmn);
	report_subdomain(dmn, subdmn);
	correct_subdomain(dmn, &subdmn);		//corrects subdomain if it overlaps with the domain
	if(dmn.dim == 3){
		read_cfg_mslicing(&mmain_plane, &mslicing_method, &mcuts);
		if(mslicing_method == 1) auto_slice(subdmn, &mcuts, mmain_plane);
		report_mesh_slicing(mcuts, mmain_plane, dmn);
		check_slicing(mcuts, mmain_plane, subdmn);
	}
	if(f_task.any){
		printf("Fields will be loaded from ");
		(field_source == 2)?printf("plot snapshots\n"):printf("checkpoints\n");
	}

	vert_pic_size = plot_init(&plt, hor_pic_size, dmn, subdmn, mmain_plane);
	get_coarse_factors(&coarse_factor, coarsen_data, dmn, mmain_plane, subdmn, hor_pic_size, vert_pic_size - VERT_EXTRA);

	printf("Data processing:\n");
	for(i = 0; i < checkpoints_list.N; i++){
		printf("\tprocessing checkpoint #%d (%d of %d)...\n", checkpoints_list.data[i], i+1, checkpoints_list.N);
		if((field_source == 1)&&f_task.any){
			process_fields(bin_data_dir, checkpoints_list.data[i], dmn, subdmn, mmain_plane, mcuts, "sys", plt, hor_pic_size, vert_pic_size, f_task, coarse_factor, s_task);
		}
		if(ps_task.any||df_task.any){
			process_markers(bin_data_dir, checkpoints_list.data[i], dmn, subdmn, ps_task, df_task, hor_pic_size, vert_pic_size, &n, &n_species, &df_success, plt);
			if(ps_task.any){
				if(n > 0) printf("\t  - scattered data loaded: %ld particles of totally %d species will be visualized\n", n, n_species);
				else printf("\t  - scattered data loaded: no particles for visualization in the specified subdomain\n");
			}
			if(df_success) printf("\t  - distribution function successfully created\n");
			else if(df_task.any) printf("\t  - distribution function creation FAILED\n");
		}
	}

	printf("\n");
	for(i = 0; i < snapshots_list.N; i++){
		printf("\tprocessing plot snapshot #%d (%d of %d)...\n", snapshots_list.data[i], i+1, snapshots_list.N);
		if((field_source == 2)&&f_task.any){
			process_fields(bin_data_dir, snapshots_list.data[i], dmn, subdmn, mmain_plane, mcuts, "tec", plt, hor_pic_size, vert_pic_size, f_task, coarse_factor, s_task);
		}
		if(c_task.any){
			process_mplasma(bin_data_dir, snapshots_list.data[i], dmn, subdmn, mmain_plane, mcuts, "tec", plt, c_task, coarse_factor);
		}
	}

	printf("Visualizing special tasks...\n");
	if(s_task.e_diag) process_energy_diagnostic(bin_data_dir, s_task, plt, hor_pic_size);
	if(s_task.vis_probes){
		vis_probes(dmn, s_task, hor_pic_size, plt, bin_data_dir);
		if((s_task.vis_probe_locations)&&(s_task.probes_num_domain > 0)) vis_probe_locations(dmn, s_task, plt, hor_pic_size, vert_pic_size);
	}

#ifdef HUI
	if (s_task.vis_gamma) { // TEMP
		double time;
		double deltaN[100];
		double deltaN_plus[100];
		double gamma[100];
		double gamma_plus[100];
		double total[100];
		double total_plus[100];
		char name[300];
		double tmp, tmp2;
		FILE *fp;

		for (int i = 0; i < 100; i++) total_plus[i] = total[i] = 0.0;

		for (int i = 0; i < checkpoints_list.N; i++) {
			int num = checkpoints_list.data[i];

			if (num == 0) continue;
			
			partition p = load_partition_info(bin_data_dir, num, "sys", &time);
			
			for (int j = 0; j < 100; j++) deltaN_plus[j] = deltaN[j] = 0.0;

			for (int j = 0; j < p.N; j++) {
				sprintf(name, "%sgamma_%06d_%03d.bin", bin_data_dir, num, j);
				fp = fopen(name, "rt");
				
				for (int k = 0; k < 100; k++) {
					fscanf(fp, "%lg %*lg %lg %lg", &gamma[k], &tmp, &tmp2);
					deltaN[k] += tmp;
					total[k] += tmp;
					deltaN_plus[k] += tmp2;
					total_plus[k] += tmp2;
				}

				fclose(fp);
			}

			sprintf(name, "dat/gamma_%06d.dat", num);
			fp = fopen(name, "wt");

			for (int j = 0; j < 100; j++) fprintf(fp, "%e %e\n", gamma[j], log(1 + deltaN[j]));

			fclose(fp);
/*
			sprintf(name, "dat/gamma_plus_%06d.dat", num);
			fp = fopen(name, "wt");

			for (int j = 0; j < 100; j++) fprintf(fp, "%e %e\n", gamma[j], log(1 + deltaN_plus[j]));

			fclose(fp);
			*/
			
			fprintf(plt, "set output 'pics/gamma_%06d.png'\n", num);
			fprintf(plt, "set title 'gamma: time = %.2f waveperiods'\n", time);
			fprintf(plt, "set ylabel 'ln (1 + delta N)'\n");
			fprintf(plt, "set xlabel 'gamma'\n");
			fprintf(plt, "set grid\n");
			fprintf(plt, "set size square\n");
			fprintf(plt, "plot '%s' w l notitle\n", name);
			fprintf(plt, "unset grid\n");

		}

		fp = fopen("dat/gamma_total.dat", "wt");
		
		for (int i = 0; i < 100; i++) fprintf(fp, "%e %e\n", gamma[i], log(1 + total[i]));

		fclose(fp);
		
		fprintf(plt, "set output 'pics/gamma_total.png'\n");
		fprintf(plt, "set title 'total gamma'\n");
		fprintf(plt, "set ylabel 'ln (1 + delta N)'\n");
		fprintf(plt, "set xlabel 'gamma'\n");
		fprintf(plt, "set grid\n");
		fprintf(plt, "plot 'dat/gamma_total.dat' w l notitle\n");
		fprintf(plt, "unset grid\n");

		fp = fopen("dat/gamma_plus_total.dat", "wt");
		
		for (int i = 0; i < 100; i++) fprintf(fp, "%e %e\n", gamma[i], log(1 + total_plus[i]));

		fclose(fp);
	} // TEMP
#endif

	free_int_array(&checkpoints_list);
	free_int_array(&snapshots_list);
	if(dmn.dim == 3) free_int_array(&mcuts);
	if(s_task.vis_probes){
		free_int_array(&(s_task.probes));
		free(s_task.probe_locations);
	}

	fclose(plt);
	printf("Plotting pictures...\n");
	if(system("/usr/bin/gnuplot plot.plt 2>dat/gnuplot_message.log")) printf("Warning: gnuplot has exited with non-zero return code. Probably, the plotting failed. See file dat/gnuplot_message.log for details\n");
	else printf("Done\n");

	return(0);
}












