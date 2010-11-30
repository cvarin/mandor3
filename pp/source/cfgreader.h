//header file for cfgreader.c
#include"types.h"

void read_cfg_common(char *bin_data_dir, int_array *checkpoints_list, int_array *snapshots_list, int *hor_pic_size, subdomain *subdmn, int *coarsen_data);
void read_cfg_mdata(int *fs, fields_task *f_task, charges_task *c_task);
void read_cfg_mslicing(int *mp, int *cm, int_array *cuts);
void read_cfg_markers(phase_space_task *ps_task);
void read_cfg_distrib(distrib_function_task *df_task);
void read_cfg_special_task(char *bin_data_dir, special_task *s_task, domain dmn);
int probe_is_inside_domain(char *probe, domain dmn);

