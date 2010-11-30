//header file for reports.c

#include"types.h"

void report_what_to_be_processed(int_array *checkpoints_list, int_array *snapshots_list);
void report_domain(domain dmn);
void report_mesh_slicing(int_array mcuts, int mp, domain dmn);
void report_subdomain(domain dmn, subdomain subdmn);
void report_mcoarsening(int_vector *cf);
void report_phase_space(phase_space_task *pt);
void df_report(distrib_function_task *dt, int max);
void st_report(special_task *st);

