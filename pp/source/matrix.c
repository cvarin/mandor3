//matrix library for further usage
//
//All the matrices start from zero

#include<stdio.h>
#include<stdlib.h>
#include"matrix.h"
#include"utils.h"

#define ARRAYS_OUB_CHECK 1	//defines whether the program has to always check the out of bounds criterion or not


arr2d arr2d_create(int Ni, int Nj){
	//creates and allocates 2D an Ni x Nj array
	arr2d result;
	arr2d_allocate(&result, Ni, Nj);
	return result;
}

void arr2d_allocate(arr2d *A, int Ni, int Nj){
	//allocates the Ni x Nj array
	A->nx = Ni;
	A->ny = Nj;
	A->allocated = 1;
	A->storage = calloc(Ni*Nj, sizeof(double));
	if(!A->storage) error(FUNCTION, "cannot allocate %d kbytes of memory to store data\n", (Ni*Nj)/1024);
}

double arr2d_ij(arr2d A, int i, int j){
	//extracts (i,j)-th element from the 2D array A
	int cnt = i*A.ny + j;
#if ARRAYS_OUB_CHECK == 1
	if(i > (A.nx - 1)) error(FUNCTION, "index i = %u exceeds maximum value of %u\n", i, A.nx - 1);
	if(j > (A.ny - 1)) error(FUNCTION, "index j = %u exceeds maximum value of %u\n", j, A.ny - 1);
#endif
	return(A.storage[cnt]);
}

void arr2d_set_ij(double A_ij, arr2d *A, int i, int j){
	//sets (i,j)-th element of the 2D array A to A_ij
	int cnt = i*(A->ny) + j;
#if ARRAYS_OUB_CHECK == 1
	if(i > (A->nx - 1)) error(FUNCTION, "index i = %u exceeds maximum value of %u\n", i, A->nx - 1);
	if(j > (A->ny - 1)) error(FUNCTION, "index j = %u exceeds maximum value of %u\n", j, A->ny - 1);
#endif
	A->storage[cnt] = A_ij;
}

void arr2d_print(arr2d A){
	//prints the contents of the array
	int i, j;
	for(j = 0; j < A.ny; j++){
		for(i = 0; i < A.nx; i++){
			printf("%.2f ", arr2d_ij(A, i, j));
		}
		printf("\n");
	}
}

void arr2d_cpy(arr2d A, arr2d B){
	//copies array A to B
	if(A.nx != B.nx) error(FUNCTION, "A.nx(%d) != B.nx(%d)\n", A.nx, B.nx);
	if(A.ny != B.ny) error(FUNCTION, "A.ny(%d) != B.ny(%d)\n", A.ny, B.ny);
	int i, j;
	for(i = 0; i < A.nx; i++){
		for(j = 0; j < A.ny; j++){
			arr2d_set_ij(arr2d_ij(A, i, j), &B, i, j);
		}
	}
}

void arr2d_destroy(arr2d A){
	//destroys 2D array
	if(A.allocated == 1){
		A.allocated = 0;
		free(A.storage);
	}
	else printf("arr2d_destroy warning: the array is not allocated\n");
}




// non-double vals:

arr2d_char arr2d_char_create(int Ni, int Nj){
	//creates and allocates 2D an Ni x Nj array
	arr2d_char result;
	arr2d_char_allocate(&result, Ni, Nj);
	return result;
}

void arr2d_char_allocate(arr2d_char *A, int Ni, int Nj){
	//allocates the Ni x Nj array
	A->nx = Ni;
	A->ny = Nj;
	A->allocated = 1;
	A->storage = calloc(Ni*Nj, sizeof(char));
	if(!A->storage) error(FUNCTION, "cannot allocate %d kbytes of memory to store data\n", (Ni*Nj)/1024);
}

double arr2d_char_ij(arr2d_char A, int i, int j){
	//extracts (i,j)-th element from the 2D array A
	int cnt = i*A.ny + j;
#if ARRAYS_OUB_CHECK == 1
	if(i > (A.nx - 1)) error(FUNCTION, "index i = %u exceeds maximum value of %u\n", i, A.nx - 1);
	if(j > (A.ny - 1)) error(FUNCTION, "index j = %u exceeds maximum value of %u\n", j, A.ny - 1);
#endif
	return(A.storage[cnt]);
}

void arr2d_char_set_ij(char A_ij, arr2d_char *A, int i, int j){
	//sets (i,j)-th element of the 2D array A to A_ij
	int cnt = i*(A->ny) + j;
#if ARRAYS_OUB_CHECK == 1
	if(i > (A->nx - 1)) error(FUNCTION, "index i = %u exceeds maximum value of %u\n", i, A->nx - 1);
	if(j > (A->ny - 1)) error(FUNCTION, "index j = %u exceeds maximum value of %u\n", j, A->ny - 1);
#endif
	A->storage[cnt] = A_ij;
}

void arr2d_char_destroy(arr2d_char A){
	//destroys 2D array
	if(A.allocated == 1){
		A.allocated = 0;
		free(A.storage);
	}
	else printf("arr2d_destroy warning: the array is not allocated\n");
}








