//header file for processor.c

#include"types.h"
#include"matrix.h"

void auto_slice(subdomain subdmn, int_array *mcuts, int mmain_plane);
void check_slicing(int_array mcuts, int mmain_plane, subdomain subdmn);
void process_fields(char *bin_data_dir, int num, domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, char *source_spec, FILE *plt, int hor_pic_size, int vert_pic_size, fields_task f_task, int_vector coarse_factor, special_task st);
void process_mplasma(char *bin_data_dir, int num, domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, char *source_spec, FILE *plt, charges_task c_task, int_vector coarse_factor);
int point_kind(int i, int j, int k, domain dmn, int mmain_plane, int_array mcuts, subdomain subdmn, int_vector coarse_factor);
int inside_subdomain(int i, int j, int k, subdomain subdmn);
void correct_subdomain(domain dmn, subdomain *subdmn);
void calc_cut_size(domain dmn, int mmain_plane, int_array mcuts, subdomain subdmn, int_vector coarse_factor, int *nx, int *ny);
int get_cut_num(int i, int j, int k, int mmain_plane, int_array mcuts);
void get_coarse_factors(int_vector *cf, int enable, domain dmn, int  mmain_plane, subdomain subdmn, int hor_pic_size, int vert_pic_size);
void savepoint(double val, int i, int j, int k, domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, int_vector cf, arr2d *A);
double intensity(vector E, vector H);
double energy_density(vector E, vector H);
void save_dat_files(domain dmn, subdomain subdmn, int mmain_plane, int_array mcuts, arr2d *A, int num, char *spec, FILE *plt, double time);
void save_dat_file(arr2d A, char *fn, int dim, double xmin, double ymin, double xmax, double ymax, int mp, double z, char *spec, int num, int cut, FILE *plt, double time);
void perform_fourier_transform(arr2d in, arr2d *out);
void process_markers(char *bin_data_dir, int num, domain dmn, subdomain subdmn, phase_space_task ps_task, distrib_function_task df_task, int hor_pic_size, int vert_pic_size, long int *n, int *n_species, int *df_success, FILE *plt);
long int get_markers_num(char *bin_data_dir, int num, int cpu);
void find_local_qm(marker *m, long int N, double_array *a, domain dmn, subdomain subdmn, phase_space_task ps_task);
void update_global_qm(double_array local, double_array *global);
void save_markers_direct(marker *m, long int N, int num, domain dmn, subdomain subdmn, long int *n, double_array qm, phase_space_task ps_task);
int phase_point_inside_bounds(marker m, domain dmn, subdomain subdmn, phase_space_task pt);
double extract_val(marker m, phase_space_task ps_task, int index);
void markers_get_interval(domain dmn, subdomain subdmn, phase_space_task pt, int ind, interval *i);
void carray_fill(arr2d_char *A, domain dmn, subdomain subdmn, phase_space_task pt, double qm, marker *m, long int N, long int *n);
void save_markers_coarsened(arr2d_char A, int num, int specie, domain dmn, subdomain subdmn, phase_space_task ps_task);
void update_df(marker *m, long int N, arr2d df, distrib_function_task df_task, domain dmn, subdomain subdmn);
double df_extract_param(marker m, distrib_function_task dt, int n);
double calc_theta(double r, double x);
int marker_is_filtered_by_ef(marker m, int ef, interval i);
int marker_inside_subdomain(marker m, domain dmn, subdomain subdmn);
void process_energy_diagnostic(char *bin_data_dir, special_task st, FILE *plt, int hor_pic_size);
void vis_probe_locations(domain dmn, special_task st, FILE *plt, int hor_pic_size, int vert_pic_size);
void vis_probes(domain dmn, special_task st, int hor_pic_size, FILE *plt, char *bin_data_dir);




