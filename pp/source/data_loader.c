//contains utils dealing with loading data

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"utils.h"
#include"params.h"
#include"data_loader.h"

int get_checkpoints_num(char *bin_data_dir, char *suffix){
	//reads the number of saved checkpoints
	FILE *fp;
	char fn[STRLEN];
	int num;

	sprintf(fn, "%s%s.N", bin_data_dir, suffix);
	OPEN_FILE(fp, fn, "r");
	fscanf(fp, "@ %d", &num);
	fclose(fp);
	return num;
}

void load_domain_info(domain *dmn, char *bin_data_dir){
	//loads domain information
	FILE *fp;
	char fn[STRLEN];
	sprintf(fn, "%s/domain.bin", bin_data_dir);
	OPEN_FILE(fp, fn, "rb");
	fread(&(dmn->h1), sizeof(double), 1, fp);
	fread(&(dmn->h2), sizeof(double), 1, fp);
	fread(&(dmn->h3), sizeof(double), 1, fp);
	fread(&(dmn->Lx), sizeof(double), 1, fp);
	fread(&(dmn->Ly), sizeof(double), 1, fp);
	fread(&(dmn->Lz), sizeof(double), 1, fp);
	fread(&(dmn->tau), sizeof(double), 1, fp);
	fread(&(dmn->imin), sizeof(int), 1, fp);
	fread(&(dmn->jmin), sizeof(int), 1, fp);
	fread(&(dmn->kmin), sizeof(int), 1, fp);
	fread(&(dmn->imax), sizeof(int), 1, fp);
	fread(&(dmn->jmax), sizeof(int), 1, fp);
	fread(&(dmn->kmax), sizeof(int), 1, fp);
	fclose(fp);
	if( (dmn->imin != dmn->imax)&&(dmn->jmin != dmn->jmax)&&(dmn->kmin != dmn->kmax) ){
		dmn->dim = 3;
		dmn->mp_2D = 0;
		dmn->aa_1D = 0;
	}
	else if( (dmn->imin == dmn->imax)&&(dmn->jmin != dmn->jmax)&&(dmn->kmin != dmn->kmax) ){
		dmn->dim = 2;
		dmn->mp_2D = 2;
		dmn->aa_1D = 0;
	}
	else if( (dmn->imin != dmn->imax)&&(dmn->jmin == dmn->jmax)&&(dmn->kmin != dmn->kmax) ){
		dmn->dim = 2;
		dmn->mp_2D = 3;
		dmn->aa_1D = 0;
	}
	else if( (dmn->imin != dmn->imax)&&(dmn->jmin != dmn->jmax)&&(dmn->kmin == dmn->kmax) ){
		dmn->dim = 2;
		dmn->mp_2D = 1;
		dmn->aa_1D = 0;
	}
	else if( (dmn->imin == dmn->imax)&&(dmn->jmin == dmn->jmax)&&(dmn->kmin != dmn->kmax) ){
		dmn->dim = 1;
		dmn->mp_2D = 0;
		dmn->aa_1D = 3;
	}
	else if( (dmn->imin == dmn->imax)&&(dmn->jmin != dmn->jmax)&&(dmn->kmin == dmn->kmax) ){
		dmn->dim = 1;
		dmn->mp_2D = 0;
		dmn->aa_1D = 2;
	}
	else{
		dmn->dim = 1;
		dmn->mp_2D = 0;
		dmn->aa_1D = 1;
	}
}

partition load_partition_info(char *bin_data_dir, int num, char *spec, double *time){
	//loads partitioning information
	partition p;
	int i, n, id;
	FILE *fp;
	char fn[STRLEN], tmp[MAXSTRLEN];
	float tm;

	sprintf(fn, "%s%s_%06d.txt", bin_data_dir, spec, num);
	if(!strcmp(spec, "sys")){
		OPEN_FILE(fp, fn, "rt");
		fscanf(fp, "@ %e", &tm);
		*time = (double)tm;
		fgets(tmp, MAXSTRLEN, fp);	//skip the remainder of the line
		fscanf(fp, "@ %d", &n);
		fgets(tmp, MAXSTRLEN, fp);	//skip the remainder of the line
		fscanf(fp, "@ %d", &id);
		fclose(fp);
	}
	else if(!strcmp(spec, "tec")){
		OPEN_FILE(fp, fn, "rb");
		fread(time, sizeof(double), 1, fp);
		fread(&n, sizeof(int), 1, fp);
		fread(&id, sizeof(int), 1, fp);
		fclose(fp);
	}
	else error(FUNCTION, "unknown spec = \"%s\"", spec);

	p.N = n;
	p.nodes = calloc(n, sizeof(cpudmn));

	sprintf(fn, "%sfileMap_%06d.txt", bin_data_dir, id);
	OPEN_FILE(fp, fn, "rt");
	fgets(tmp, MAXSTRLEN, fp);	//skip first line
	for(i = 0; i < p.N; i++){
		fgets(tmp, MAXSTRLEN, fp);
		if(strlen(tmp) < 2) error(FUNCTION, "inconsistent structure of file %s", fn);
		sscanf(tmp, "%d: (%d %d %d) - (%d %d %d)", &n, &(p.nodes[i].imin), &(p.nodes[i].jmin), &(p.nodes[i].kmin), &(p.nodes[i].imax), &(p.nodes[i].jmax), &(p.nodes[i].kmax));
		//printf("cpu %d: (%d %d %d) - (%d %d %d)\n", n, p.nodes[i].imin, p.nodes[i].jmin, p.nodes[i].kmin, p.nodes[i].imax, p.nodes[i].jmax, p.nodes[i].kmax);
	}
	fclose(fp);

	return p;
}

void open_mesh_binary_files(char *bin_data_dir, char *source_spec, char *target_spec, int num, FILE **fs, int N){
	//opens the binary files containing meshed data and positions heads at the beginning of useful data
	char fn[STRLEN];
	int cpu;
#if CPU_64_BIT == 1
	#define HEADER_SIZE 112
#else
	#define HEADER_SIZE 108
#endif
	for(cpu = 0; cpu < N; cpu++){
		sprintf(fn, "%s%s_%s_%06d(%03d).bin", bin_data_dir, source_spec, target_spec, num, cpu);
		OPEN_FILE(fs[cpu], fn, "rb");
		fseek(fs[cpu], HEADER_SIZE, SEEK_SET);
	}

}

vector read_vector_point(int i, int j, int k, partition p, FILE **fs){
	//reads the data element from the binary files
	int cpu;
	vector vec;
#if SMART_READING == 1
	//smart reading would keep part(s) of the binary file(s) in fuffer or so
	error(FUNCTION, "the SMART_READING option is currently not implemented");
#else
	for(cpu = 0; cpu < p.N; cpu++) if(inside_cpu_subdomain(i, j, k, p, cpu)) fread(&vec, sizeof(vector), 1, fs[cpu]);
#endif
	return vec;
}

double read_double_point(int i, int j, int k, partition p, FILE **fs){
	//reads the data element from the binary files
	int cpu;
	double val;
#if SMART_READING == 1
	error(FUNCTION, "the SMART_READING option is currently not implemented");
#else
	for(cpu = 0; cpu < p.N; cpu++) if(inside_cpu_subdomain(i, j, k, p, cpu)) fread(&val, sizeof(double), 1, fs[cpu]);
#endif
	return val;
}

int inside_cpu_subdomain(int i, int j, int k, partition p, int cpu){
	//decides whether the point belongs this cpu or not
	if((i >= p.nodes[cpu].imin)&&(i <= p.nodes[cpu].imax)&&(j >= p.nodes[cpu].jmin)&&(j <= p.nodes[cpu].jmax)&&(k >= p.nodes[cpu].kmin)&&(k <= p.nodes[cpu].kmax)) return 1;
	else return 0;
}

void load_plasma(marker *m, char *bin_data_dir, int num, int cpu){
	//loads the plasma from checpoint for the current cpu
	long int N;
	char fn[STRLEN];
	FILE *fp;
	sprintf(fn, "%splasma_%06d_%03d.bin", bin_data_dir, num, cpu);
	OPEN_FILE(fp, fn, "rb");
	fread(&N, sizeof(long int), 1, fp);
	fread(m, sizeof(marker), (size_t)N, fp);
	fclose(fp);
}


