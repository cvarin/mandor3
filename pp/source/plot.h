//header file plot.c

#include"types.h"
#include"matrix.h"

int plot_init(FILE **plt, int hor_pic_size, domain dmn, subdomain subdmn, int mmain_plane);
void markers_update_plt(phase_space_task pt, domain dmn, subdomain subdmn, int hor_pic_size, int vert_pic_size, int num, double_array qm, double time, FILE *plt);
int markers_get_vert_pic_size(phase_space_task pt, domain dmn, subdomain subdmn, int hor_pic_size);
void get_sorted_indices(int *indices, int n, double *vals);
void df_plot(distrib_function_task dt, arr2d df, int hor_pic_size, int vert_pic_size, int num, double time, FILE *plt);
void plot_energy_diag(special_task st, int hor_pic_size, FILE *plt);
void plot_probe_locations(domain dmn, FILE *plt, int hor_pic_size, int vert_pic_size);
void plot_probes_data(special_task st, int hor_pic_size, FILE *plt);
void plot_fourier_transform(domain dmn, special_task st, int mmain_plane, int_array mcuts, arr2d A, int i, int num, FILE *plt, int hor_pic_size, int vert_pic_size, double time);

