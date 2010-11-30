//header file for matrix.c

#ifndef MATRIX_H_DEFINED

#define MATRIX_H_DEFINED
	typedef struct{
		int nx;
		int ny;
		char allocated;
		double *storage;
	} arr2d;

	typedef struct{
		int nx;
		int ny;
		char allocated;
		char *storage;
	} arr2d_char;

	arr2d arr2d_create(int Ni, int Nj);
	void arr2d_allocate(arr2d *A, int Ni, int Nj);
	void arr2d_destroy(arr2d A);
	double arr2d_ij(arr2d A, int i, int j);
	void arr2d_set_ij(double A_ij, arr2d *A, int i, int j);
	void arr2d_print(arr2d A);
	void arr2d_cpy(arr2d A, arr2d B);

	arr2d_char arr2d_char_create(int Ni, int Nj);
	void arr2d_char_allocate(arr2d_char *A, int Ni, int Nj);
	double arr2d_char_ij(arr2d_char A, int i, int j);
	void arr2d_char_set_ij(char A_ij, arr2d_char *A, int i, int j);
	void arr2d_char_destroy(arr2d_char A);

#endif

