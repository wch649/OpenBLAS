#pragma once
constexpr auto inf = 100;
#ifndef MATRIX
class Matrix {

protected:
	int column;
	int row;
	double *array;

public:
	/*
		Accessibility
	*/
	//create A={a[0],a[1]....}
	Matrix(int m, int n, double *a);
	//create A={0,0,0,0....}
	Matrix(int m, int n);
	/*create A={1,0,0,
				0,1,0,
				0,0,1}
	*/
	Matrix(int n);

	void setrow(int m);
	void setcol(int n);
	void setarray(double *a);
	int getrow();
	int getcol();
	double *getarray();
	void show();
	~Matrix();
	
	/*
		Basic operation of the matrix
	*/

	/*
		Description: A+B 
		Input:	m: rows of C、A
				n: columns of C、A
				a: A(m*n)
				beta: 1
				c: C(m*n)
		Output: C = a*A + b*C;
		Others:
	*/
	void MatrixAdd(int m, int n, double * a, double * c, double beta);
	Matrix operator+(const Matrix &A);
	/*
		Description: A-B equal A+(-1)B
		Input:	m: rows of C、A
				n: columns of C、A
				a: A(m*n)
				beta: -1
				c: C(m*n)
		Output: C = A - C;
		Others:
	*/
	void MatrixLess(int m, int n, double * a, double * c, double beta);
	Matrix operator-(const Matrix &A);
	/*
		Description: A*B
		C := a*op(A)*op(B) + b*C，a、b is scalar
		Input:	m: the row of A and C
				n: the column of B and C
				k: the column of B and A
				a: 	A(m*k)
				b:  B(k*n)
				c:	C(m*n)
		Output: C = A*B
	*/
	void MatrixMulMatrix(int m, int n, int k, double *a, double *b, double *c);
	Matrix operator*(const Matrix &A);
	/*
		Description: A*x
		Input:
			m: row of A
			n: column of A
			a: A(m*n)
			x: vector
			y: vector
		Output: y = a*op(A)*x + beta*y
				y = A*x
	*/	
	void MatrixMulVector(int m, int n, double *a, double *x, double *y);
	//friend Matrix operator*(const Vector &);
	/*
		Description: A*k
		Input:
		Output:
		Others:
	*/
	void MatrixMulNumber(int m, int n, double * a, double * c, double alpha);
	Matrix MatrixAlpha(double alpha);
	/*
		Description: inv(A)
		Input: A(n*n) 
			n: matrix's size
		Output: A = inv(A)
	*/
	void MatrixInverse(double *A, int N);
	//inv(A) to use this method
	/*
		Description: left division A\B = inv(A) * B
		Input: Matrix B
		Output: Matrix D = A\B
		Others: \ is not operator and not overloaded
	*/
	void MatrixLeftDiv(Matrix B);
	//Matrix operator\(); '\' is not overloaded
	/*
		Description: right division A/B = A * inv(B)
		Input: Matrix B
		Output: Matrix D = A/B
		Others: function operator/() is ok!
	*/
	void MatrixRightDiv();
	Matrix operator/(Matrix &A);
	/*
		Description: judge A = B
		Input:
		Output:
		Others:
	*/
	bool MatrixIsEqual();
	/*
		Description: C = A'
		Input: Matrix A's row col and array
		Output: C = A'
		Others: 
			C := a*op(A)*op(B) + b*C，a=1,b=0,op(A) = CblasTrans, B = E
	*/
	void MatrixTranspose(int m, int n, double * a, double *c);
	/*
		Description: |A|
		Input:
		Output:
		Others:
	*/
	int MatrixDeterminant();
	/*
		Description: rank(A)
		Input:
		Output:
		Others:
	*/
	int MatrixRank();
	

	/*
		Eigenvalue calculation of matrix
	*/
	
	/*
		Matrix decomposition
	*/
	
	/*
		Norm of the matrix
	*/
	double Matrix_1Norm(Matrix A);
	double Matrix_2Norm(Matrix A);
	double Matrix_InfNorm(Matrix A);
};

/*
	Description: inv(A)
	Input: Matrix A
	Output: Matrix B = inv(A);
*/
Matrix inv( Matrix A);
/*
	Description: A\B
	Input: Matrix A、B
	Output: Matrix C = A\B
*/
Matrix ld(Matrix A, Matrix B);
/*
	Description: A'
	Input: Matrix A
	Output: Matrix C = A'
*/
Matrix tran(Matrix A);

/*
	Description: ?-norm
	Input: Matrix A 
		?: 1、2、Inf
	Output: ?-norm
*/
double norm(Matrix A, int count);
#endif