//header file for data_loader.c

int get_checkpoints_num(char *bin_data_dir, char *suffix);
void load_domain_info(domain *dmn, char *bin_data_dir);
partition load_partition_info(char *bin_data_dir, int num, char *spec, double *time);
void open_mesh_binary_files(char *bin_data_dir, char *source_spec, char *target_spec, int num, FILE **fs, int N);
vector read_vector_point(int i, int j, int k, partition p, FILE **fs);
double read_double_point(int i, int j, int k, partition p, FILE **fs);
int inside_cpu_subdomain(int i, int j, int k, partition p, int cpu);
void load_plasma(marker *m, char *bin_data_dir, int num, int cpu);



