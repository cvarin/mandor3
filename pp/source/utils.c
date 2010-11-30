//general utilities for Mandor postprocessor

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdarg.h>
#include"utils.h"
#include"params.h"
#include"types.h"

void error(char *source, char *errmsg, ...){
	//prints the error message and terminates the program
	va_list argptr;
	/*va_start(argptr, errmsg);
	char *source;
	source = va_arg(argptr, char*);
	va_end(argptr);*/
	va_start(argptr, errmsg);
	printf("\nerror in module %s: ", source);
	vprintf(errmsg, argptr);
	printf("\n");
	va_end(argptr);
	exit(1);
}

void free_int_array(int_array *array){
	//empties the array
	array->N = 0;
	free(array->data);
}

void partiton_erase(partition prtn){
	//erases partition
	prtn.N = 0;
	free(prtn.nodes);
}

int in(double x, double xmin, double xmax){
	//returns 1 if x \in [xmin, xmax] and 0 otherwise
	if((x >= xmin)&&(x <= xmax)) return 1;
	else return 0;
}










