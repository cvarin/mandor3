//header file for utils.c

#include"types.h"

void error(char *source, char *errmsg, ...);
void free_int_array(int_array *array);
void partiton_erase(partition prtn);
int in(double x, double xmin, double xmax);

#define FUNCTION (char*)__func__

#define OPEN_FILE(fp_, filename_, mode_) {						\
	fp_ = fopen(filename_, mode_);							\
	if(!fp_) error(FUNCTION, "cannot open file %s in %s mode", filename_, mode_);	\
}
















