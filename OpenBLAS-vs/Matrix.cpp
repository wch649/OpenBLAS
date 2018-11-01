#include "Matrix.h"
#include <cblas.h>
#include <iostream>
#include <math.h>
#include <lapacke.h>
#include <string.h>
using namespace std;


#define min(a,b) ((a)>(b)?(b):(a))
double max(double *s) {
	double M = s[0];
	int len = sizeof(*s) / sizeof(double);
	for (int i = 1; i < len; i++) {
		cout << s[i] << endl;
		if (s[i] > M) M = s[i];
	}
	return M;
}


#ifdef __cplusplus
extern "C" {
#endif

	// LU decomoposition of a general matrix
	//void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
	// void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);

	// generate inverse of a matrix given its LU decomposition
	//void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
	// void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK, int* lwork, int* INFO);

#ifdef __cplusplus
}
#endif
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

Matrix::Matrix()
{
	row = 0;
	column = 0;
	array = NULL;
}

void Matrix::setsize(int m, int n)
{
	if (array) delete[]array;
	array = new double[m*n];
	row = m;
	column = n;
	for (long int i = m * n - 1; i >= 0; i--) array[i] = 0.0;
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
	array = new double[work];
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

void Matrix::modifyM(int m, int n, double * a)
{
	setrow(m);
	setcol(n);
	if (array) delete[]array;
	setarray(a);
}

void Matrix::free()
{
	row = 0;
	column = 0;
	delete[]array;
}


/*
		Basic operation of the matrix
*/


void Matrix::operator=(const Matrix & A)
{
	if (array) delete[]array;
	row = A.row;
	column = A.column;
	array = new double[row*column];
	long int work = row * column;
	for (int i = 0; i < work; i++) array[i] = A.array[i];
}

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

Vector Matrix::operator*(Vector &v)
{
	Vector x(this->row);
	this->MatrixMulVector(this->row, this->column, this->getarray(), v.getv(), x.getv());
	return x;
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
	//ld(*this, B);
}

//A / B = A * inv(B)
void Matrix::MatrixRightDiv(Matrix B)
{
	//*this / B;
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

bool Matrix::MatrixIsEqual(Matrix B)
{
	if (this->row == B.row && this->column == B.column) {
		long int work = B.row * B.column;
		//cout <<work<<" "<< memcmp(this->array, B.array, work * sizeof(double)) << endl;
		if (!memcmp(this->array, B.array, work * sizeof(double))) {
			return true;
		}
	}	
	return false;
}

bool Matrix::operator==(const Matrix &B)
{
	return this->MatrixIsEqual(B);
}


void Matrix::MatrixTranspose(int m, int n, double * a, double *c)
{
	Matrix B(m);
	//lda is A's lda, but not op(A)'s lda
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, m, 1, a, n, B.array, m, 0, c, m);
	delete[]B.array;
}

double Matrix::MatrixDeterminant(int n, double *a)
{
	double ret = 1.0;
	for (int i = 0; i < n; i++) {
		if (a[i*n + i] < 0) {
			ret = -ret;
			cblas_dscal(n - i, -1, a + i * n + i, 1);
		}
		for (int j = i + 1; j < n; j++) {
			while (fabs(a[j*n + i]-0) >= 1e-12) {
				if (a[j*n + i] < 0) {
					ret = -ret;
					cblas_dscal(n - i, -1, a + j * n + i, 1);
				}
				double t = a[i*n + i] / a[j*n + i];
				cblas_daxpy(n, -t, a + j * n + i, 1, a + i * n + i, 1);
				cblas_dswap(n - i, a + i * n + i, 1, a + j * n + i, 1);
				ret = -ret;
			}
		}
		if (fabs(a[i*n + i] -0)<1e-12) return 0;
		ret = ret * a[i*n + i];
	}
	/*for(int i=0;i<n*n;i++)
		cout<<a[i]<<endl;*/
	return ret;
}

int Matrix::MatrixRank(int m, int n, double *a)
{
	for (int i = 0; i < n; i++) {
		if (a[i*n + i] < 0) {
			cblas_dscal(n - i, -1, a + i * n + i, 1);
		}
		for (int j = i + 1; j < n; j++) {
			while (fabs(a[j*n + i] - 0) >= 1e-12) {
				if (a[j*n + i] < 0) {
					cblas_dscal(n - i, -1, a + j * n + i, 1);
				}
				double t = a[i*n + i] / a[j*n + i];
				cblas_daxpy(n, -t, a + j * n + i, 1, a + i * n + i, 1);
				cblas_dswap(n - i, a + i * n + i, 1, a + j * n + i, 1);
			}
		}
	}
	int cnt = 0;
	for (int i = 0; i < m; i++) {
		bool isHave = false;
		for (int j = 0; j < n; j++) {
			if (a[i*m + j] != 0) {
				isHave = true;
				break;
			}
		}
		if (!isHave) break;
		else cnt++;
	}
	return cnt;
}


/*
		Eigenvalue calculation of matrix
*/

int Matrix::MatrixIsPosDef()
{
	if (!this->MatrixIsSym()) return 0;
	Matrix v;
	double *er = new double[row];
	double *ei = new double[row];
	this->MatrixEig(er, ei, v);
	delete[]ei;
	for (int i = row - 1; i >= 0; i--)
		if (er[i] <= 0.0)
		{
			delete[]er;
			return 0;
		}
	delete[]er;
	return 1;
}

int Matrix::MatrixIsSym()
{
	if (row != column) return 0;
	double *p0 = array;
	for (int i = row - 1; i > 0; i--)
	{
		double *rp = p0 + row, *cp = p0 + 1;
		for (int j = 0; j < i; j++)
			if (*rp == *cp)
			{
				rp += row;
				cp++;
			}
			else return 0;
		p0 += row + 1;
	}
	return 1;
}

double Matrix::MatrixTrace()
{
	double t = 0.0;
	double *p = array;
	for (int i = row; i > 0; i--)
	{
		t += *p;
		p += row + 1;
	}
	return t;
}

int Matrix::MatrixEig(double * er, double * ei, Matrix & v) // Matrix will be overwrite
{
	if (row != column)
	{
		printf("Not square matrix.");
		return -1;
	}
	v.setsize(row, row);
	Matrix vl(row, row);
	lapack_int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', row, array, row, er, ei, vl.array, row, v.array, row);
	/*
	er and ei contain the real and imaginary parts, respectively, of the computed eigenvalues.
	Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.

	the right eigenvectors vr(j) are stored one after another in the columns of v, in the same order as their eigenvalues.
	If the j-th eigenvalue is real, then vr(j) = v(:,j), the j-th column of v. If the j-th and (j+1)-st eigenvalues form a complex
	conjugate pair, then vr(j) = (:,j) + i*v(:,j+1) and vr(j+1) = v(:,j) - i*v(:,j+1).
	*/
	return info;
}


/*
		Matrix decomposition
*/

int Matrix::MatrixQRfactorization(Matrix & q, Matrix & r) // Matrix won't be overwrite
{
	lapack_int info;
	if (row <= column)
	{
		r.setsize(row, column);
		memcpy(r.array, array, row*column * sizeof(double));
		double *tau = new double[row];
		info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, row, column, r.array, column, tau); // get R and tau
		q.setsize(row, row);
		double *qp = q.array, *rp = r.array;
		*qp = 1.0;
		for (int i = 1; i < row; i++)
		{
			qp += row;
			rp += column;
			memcpy(qp, rp, i * sizeof(double));
			for (int j = 0; j < i; j++) rp[j] = 0.0;
			qp[i] = 1.0;
		}
		info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, row, row, row, q.array, row, tau); // get Q
		delete[] tau;
		return info;
	}
	else
	{
		// plan A : q is row * column
		q.setsize(row, column);
		memcpy(q.array, array, row*column * sizeof(double));
		double *tau = new double[row];
		info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, row, column, q.array, column, tau); // get R and tau
		r.setsize(column, column);
		double *qp = q.array, *rp = r.array;
		for (int i = column; i > 1; i--)
		{
			memcpy(rp, qp, i * sizeof(double));
			qp += column + 1;
			rp += column + 1;
		}
		*rp = *qp;
		info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, row, column, column, q.array, column, tau); // get Q
		delete[] tau;

		// plan B : q is column * column
		/*

		*/
	}

	return info;
}

int Matrix::MatrixLUfactorization(Matrix & l, Matrix & u, int * ipiv) // Matrix won't be overwrite
{
	lapack_int info;
	//int *ipiv = new int[min(row, column)];
	if (row > column)
	{
		l.setsize(row, column);
		memcpy(l.array, array, row*column * sizeof(double));
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, row, column, l.array, column, ipiv);
		// ?? I have no idea what ipiv meant......
		u.setsize(column, column);
		double *lp = l.array, *up = u.array;
		for (int i = column; i > 1; i--)
		{
			memcpy(up, lp, i * sizeof(double));
			*lp = 1.0;
			for (int j = 1; j < i; j++) lp[j] = 0.0;
			lp += column + 1;
			up += column + 1;
		}
		*up = *lp;
		*lp = 1.0;
	}
	else
	{
		u.setsize(row, column);
		memcpy(u.array, array, row*column * sizeof(double));
		info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, row, column, u.array, column, ipiv);
		l.setsize(row, row);
		double *lp = l.array, *up = u.array;
		*lp = 1.0;
		for (int i = 1; i < row; i++)
		{
			lp += row;
			up += column;
			memcpy(lp, up, i * sizeof(double));
			for (int j = 0; j < i; j++) up[j] = 0.0;
			lp[i] = 1.0;
		}
	}
	/*
	info = 0:  successful exit
		 > 0 : if INFO = i, U(i, i) is exactly zero.The factorization has been completed, but the factor U is exactly singular,
			   and division by zero will occur if it is used to solve a system of equations.
	*/
	return info;
}

int Matrix::MatrixSVDfactorization(double * s, Matrix & u, Matrix & vt)// Matrix will be rewrite
{
	u.setsize(row, row);
	vt.setsize(column, column);
	double *superb = new double[min(row, column) - 1];
	lapack_int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', row, column, array, column, s, u.array, row, vt.array, column, superb);
	// ?? I have no idea what superb meant......
	/*
	info = 0:  successful exit.
         > 0:  The algorithm computing did not converge, info specifies how many superdiagonals of an intermediate bidiagonal form B did not converge to zero.
	*/

	delete[]superb;
	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm computing SVD failed to converge.\n");
		return info;
	}
	return info;
}


/*
		Norm of the matrix
*/

double Matrix::Matrix_1Norm(Matrix A)
{
	int m = A.row;
	int n = A.column;
	double *a = A.array;
	double sum = cblas_dasum(m, a, n);
	for (int i = 1; i < n; i++) {
		double tmp = cblas_dasum(m, a + i, n);
		//cout<<tmp<<endl;
		if (tmp - sum > 1e-12)
			sum = tmp;
	}
	return sum;
}

double Matrix::Matrix_2Norm(Matrix A)
{	//2-norm = max(svd(A))
	Matrix u(row, row);
	Matrix vt(column, column);
	double *superb = new double[min(row, column) - 1];
	double *s = new double[min(row, column)];
	lapack_int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'N', 'N', row, column, array, column, s, u.array, row, vt.array, column, superb);
	//求*s的最大值
	double ans = max(s);
	delete[]s;
	delete[]superb;
	return ans;
}


double Matrix::Matrix_InfNorm(Matrix A)
{
	int m = A.row;
	int n = A.column;
	double *a = A.array;
	double sum = cblas_dasum(n, a, 1);
	for (int i = 1; i < m; i++) {
		double tmp = cblas_dasum(n, a + i * n, 1);
		//cout<<tmp<<endl;
		if (tmp - sum > 1e-12)
			sum = tmp;
	}
	return sum;
}



/*
	Not Class member function
*/

Matrix operator*(const Vector &)
{
	return Matrix();
}

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

Vector ld(Matrix A, Vector x)
{
	int m = A.getrow();//A: m*m
	int n = x.getn();//x:n*1
	double *barr = new double[m*m];
	double *temp = A.getarray();
	int work = m * m;
	for (int i = 0; i < work; i++) {
		barr[i] = temp[i];
	}
	Matrix C(m, m, barr);
	C.MatrixInverse(barr, m);//C = A, C -> inv(C)
	Vector y = C * x;
	delete[]barr;
	return y;
}

Matrix tran(Matrix A)
{
	int m = A.getrow();
	int n = A.getcol();
	Matrix C(n,m);
	A.MatrixTranspose(m, n, A.getarray(), C.getarray());
	return C;
}

double norm(Matrix A, int count)
{
	Matrix B(A.getrow(),A.getcol(),A.getarray());
	double ans = 0.0;
	switch (count) {
	case 1: ans = B.Matrix_1Norm(B); break;
	case 2: ans = B.Matrix_2Norm(B); break;
	case INF: ans = B.Matrix_InfNorm(B); break;
	default:cout << "Parameter error!\n";
	}
	delete[]B.getarray();
	return ans;
}

double det(Matrix A)
{
	int n = A.getrow();
	double *barr = new double[n*n];
	double *temp = A.getarray();
	int work = n * n;
	for (int i = 0; i < work; i++) {
		barr[i] = temp[i];
	}
	Matrix B(n, n, barr);
	//B.show();
	double detB = B.MatrixDeterminant(n,B.getarray());
	//cout << detB << endl;
	delete[]barr;
	return detB;	
}

int ranka(Matrix A)
{
	int ran = 0;
	int m = A.getrow();
	int n = A.getcol();
	if (m > n) {
		Matrix B(n,m);
		B = tran(A);
		//B.show();
		ran = B.MatrixRank(n, m, B.getarray());	
	}
	else {
		double *barr = new double[m*n];
		double *temp = A.getarray();
		int work = m * n;
		for (int i = 0; i < work; i++) {
			barr[i] = temp[i];
		}
		Matrix B(m, n, barr);
		//B.show();
		ran = B.MatrixRank(m, n, B.getarray());
		delete[]barr;
	}
	
	return ran;
}
