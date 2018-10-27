#include "Matrix.h"
#include <cblas.h>
#include <iostream>
using namespace std;


extern "C" {
	// LU decomoposition of a general matrix
	void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
	// void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);

	// generate inverse of a matrix given its LU decomposition
	void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
	// void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);
}

//Implementation of each function in the class Matrix

/*
		Accessibility
*/

Matrix::Matrix(int m, int n, double * a)
{
	row = m;
	column = n;
	array = new double[m*n];
	long int work = m * n;
	for (int i = 0; i < work; i++) array[i] = a[i];
}

Matrix::Matrix(int m, int n)
{
	row = m;
	column = n;
	array = new double[m*n];
	long int work = m * n;
	for (int i = 0; i < work; i++) array[i] = 0.0;
}

Matrix::Matrix(int n)
{
	column = row = n;
	array = new double[n*n];
	long int work = n * n;
	for (int i = 0, j = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j) array[i*n + j] = 1.0;
			else array[i*n + j] = 0.0;
		}
	}
}

void Matrix::setrow(int m)
{
	row = m;
}

void Matrix::setcol(int n)
{
	column = n;
}

void Matrix::setarray(double * a)
{
	long int work = column * row;
	for (int i = 0; i < work; i++) array[i] = a[i];
}



int Matrix::getrow()
{
	return row;
}

int Matrix::getcol()
{
	return column;
}

double * Matrix::getarray()
{
	return array;
}

Matrix::~Matrix()
{
	//delete []array;
}

void Matrix::show()
{
	long int work = row * column;
	for (int i = 0; i < work; i++) {
		cout << array[i] <<" ";
		if ( (i + 1) % column == 0) cout << "\n";
	}
	cout << endl;
}


/*
		Basic operation of the matrix
*/


void Matrix::MatrixAdd(int m, int n, double * a, double * c, double beta)
{	//C = A + C
	cblas_dgeadd(CblasRowMajor, m, n, 1, a, n, beta, c, n);
}

Matrix Matrix::operator+(const Matrix & B)
{
	Matrix C(this->row, B.column, B.array);//C = B -> C=A+B;
	this->MatrixAdd(C.row, C.column, this->array, C.array, 1);
	return C;
}

void Matrix::MatrixLess(int m, int n, double * a, double * c, double beta)
{	//C = A - C
	cblas_dgeadd(CblasRowMajor, m, n, 1, a, n, beta, c, n);
}

Matrix Matrix::operator-(const Matrix & B)
{
	Matrix C(this->row, B.column, B.array);//C = B -> C=A-B;
	this->MatrixLess(C.row, C.column, this->array, C.array, -1);
	return C;
}


void Matrix::MatrixMulMatrix(int m, int n, int k, double * a, double * b, double * c)
{
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, a, k, b, n, 0, c, n);
}

Matrix Matrix::operator*(const Matrix & B)
{
	Matrix C(this->row, B.column);
	this->MatrixMulMatrix(this->row, B.column, B.row, this->array, B.array, C.array);
	return C;
}


void Matrix::MatrixMulVector(int m, int n, double * a, double * x, double * y)
{
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1, a, n, x, 1, 0, y, 1);
}

void Matrix::MatrixMulNumber(int m, int n, double * a, double * c, double alpha)
{
	cblas_dgeadd(CblasRowMajor, m, n, alpha, a, n, 0, c, n);
}

Matrix Matrix::MatrixAlpha(double alpha)
{
	Matrix C(this->getrow(), this->getcol());
	MatrixMulNumber(this->getrow(), this->getcol(), this->getarray(), C.getarray(), alpha);
	return C;
}

void Matrix::MatrixInverse(double * A, int N)
{
	int *IPIV = new int[N];
	int LWORK = N * N;
	double *WORK = new double[LWORK];
	int INFO;

	dgetrf_(&N, &N, A, &N, IPIV, &INFO);
	dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO);

	this->array = A;
	delete[] IPIV;
	delete[] WORK;
}
//A \ B = inv(A) * B
void Matrix::MatrixLeftDiv(Matrix B)
{
	//function ld() is ok!
}

//A / B = A * inv(B)
void Matrix::MatrixRightDiv()
{
	// function operator/() is ok!
}

Matrix Matrix::operator/(Matrix & A)
{
	int m = this->row;
	int n = A.column;
	double *barr = new double[n*n];
	double *temp = A.getarray();
	int work = n * n;
	for (int i = 0; i < work; i++) {
		barr[i] = temp[i];
	}
	Matrix B(n, n, barr);
	B.MatrixInverse(barr, n);
	//B is changed! but A is not changed!
	Matrix D(m, n);
	D = *this * B;
	delete[] barr;
	return D;
}


void Matrix::MatrixTranspose(int m, int n, double * a, double *c)
{
	Matrix B(m);
	//B.show();
	//this->show();//lda is A's lda, but not op(A)'s lda
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, m, 1, a, n, B.array, m, 0, c, m);
	delete[]B.array;
}




/*
		Eigenvalue calculation of matrix
*/


/*
		Matrix decomposition
*/


/*
		Matrix decomposition
*/

/*
	Not Class member function
*/

Matrix inv(Matrix A)
{
	int n = A.getrow();
	double *barr = new double[n*n];
	double *temp = A.getarray();
	int work = n * n;
	for (int i = 0; i < work; i++) {
		barr[i] = temp[i];
	}
	Matrix B(n, n, barr);
	B.MatrixInverse(barr, n);
	////cout << "Inverse matrix :\n";
	//B.show();
	//delete[] barr;
	return B;
}

Matrix ld(Matrix A, Matrix B)
{
	int m = A.getrow();
	int n = B.getcol();
	Matrix D(m, n); //D = A \ B = inv(A) * B = C * B
	double *barr = new double[m*m];
	double *temp = A.getarray();
	int work = m * m;
	for (int i = 0; i < work; i++) {
		barr[i] = temp[i];
	}
	Matrix C(m, m, barr);
	C.MatrixInverse(barr, m);//C = A, C -> inv(C)
	D = C * B;
	delete[]barr;
	return D;
}

Matrix tran(Matrix A)
{
	int m = A.getrow();
	int n = A.getcol();
	Matrix C(n,m);
	A.MatrixTranspose(m, n, A.getarray(), C.getarray());
	return C;
}
