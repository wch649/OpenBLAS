#include <stdio.h>
#include "Matrix.h"
extern "C" {
#include "cblas.h"
#include <lapacke.h>
}
#include <iostream>
#include <stdlib.h>
using namespace std;


#define M 6
#define N 5
#define LDA N
#define LDU M
#define LDVT N
extern void print_matrix(const char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda);
int main(void) {
	////�����޸�����������
	//double b[6] = { 1,0,1,0,1,0 };
	//Vector y(6, b);
	//y.printv();
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//y.modifyv(9, a);
	//y.printv();
	//y.free();
	//y.printv();
	////���Ծ���
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//double b[6] = { 1,0,1,0,1,0 };
	//Matrix A(3, 3, a);
	//A.show();
	//A.modifyM(2, 3, b);
	//A.show();
	//A.free();
	//A.show();

	////���Ծ��󿽱�
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//double b[6] = { 1,0,1,0,1,0 };
	//Matrix A(3, 3, a);
	//Matrix B(2, 3, b);
	//A.show();
	//B.show();
	//B = A;
	//A.show();
	//B.show();
	//��������
	//Vector x(9, a);
	//Vector y(6, b);
	//x.printv();
	//y.printv();
	//y = x;
	//x.printv();
	//y.printv();


	////���ԼӼ���
	//double a[9] = { 1,2,3,4,5,6,7,8,9};
	//double b[9] = { 1,0,1,0,1,0,0,1,0};
	//double c[9] = { 0,0,0,0,0,0,0,0,0};
	//Matrix A(3, 3,a);
	//Matrix B(3, 3,b);
	//Matrix C(3, 3);
	//C = A - B;
	//A.show(); 
	//B.show(); 
	//C.show();

	////���Ծ���*����
	//double a[9] = { 1,2,3,4,5,6,7,8,9};
	//double b[9] = { 1,0,1,0,1,0,0,1,0};
	//double c[9] = { 0,0,0,0,0,0,0,0,0};
	//Matrix A(3, 3,a);
	//Matrix B(3, 3,b);
	//Matrix C(3, 3);
	//C = A * B;
	//A.show(); 
	//B.show(); 
	//C.show(); 

	////���Ծ���*����
	//double a[9] = { 1,2,3,4,5,6,7,8,9};
	//double b[3] = { 1,2,3 };
	//Matrix A(3, 3, a);
	//Vector B(3, b);
	//Vector C(3);
	//C = A * B;
	//C.printv();

	////���Ծ���*��ֵ
	//double a[9] = { 1,2,3,4,5,6,7,8,9};
	//Matrix A(3, 3, a);
	//Matrix B = A.MatrixAlpha(3);
	//B.show();


	////��������
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//Matrix A(3, 3, a);
	//Matrix B(3, 3, a);
	//cout << "A:" << endl;
	////A.MatrixInverse(A.getarray(), A.getrow());
	//A.show();
	//cout << "B=inv(A):" << endl;
	//B = inv(A);
	//cout << "A:" << endl;
	//A.show();
	//B.show();

	////�����ҳ�
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//double b[12] = { 1,0,1,0,1,0,0,1,0,1,1,0 };
	//Matrix A(3, 3, a);
	//Matrix B(4, 3, b);
	//Matrix C(4, 3);
	//C = B/A;
	//cout << "A:\n";
	//A.show();
	//cout << "B:\n";
	//B.show();
	//cout << "C:\n";
	//C.show();


	////���Ծ������
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//double b[12] = { 1,0,1,0,1,0,0,1,0,1,1,0 };
	//Matrix A(3, 3, a);
	//Matrix B(3, 4, b);
	//Matrix C(3, 4);
	//C = ld(A,B); // C = A \ B = inv(A) * B
	//cout << "A:\n";
	//A.show();
	//cout << "B:\n";
	//B.show();
	//cout << "C=A\\B:\n";
	//C.show();

	////���Ծ����������
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//double b[3] = { 1,2,3 };
	//Matrix A(3, 3, a);
	//Matrix B(3, 3, a);
	//B = inv(A);
	//B.show();
	//Vector x(3,b);
	//Vector y(3);
	//y = ld(A, x);
	//x.printv();
	//A.show();
	//y.printv();

	//�����������

	////�����������
	//double a[9] = { 1,2,3,4,5,6,7,8,9};
	//double b[9] = { 1,0,1,0,1,0,0,1,0};
	//Vector x(9, a);
	//Vector y(9, b);
	//x.printv();
	//y.printv();
	//cout <<x * y<<endl;
	
	////������������
	//double a[9] = { 1,2,3,4,5,6,7,8,9};
	//Vector x(9, a);
	//x.printv();
	//cout << "1-norm: "<< norm(x, 1) << endl;
	//cout << "2-norm: " << norm(x, 2) << endl;
	//cout << "inf-norm: " << norm(x, INF)<<endl;
	//x.printv();

	////�������
	//double a[6] = { 1.0,2.0,3.0,4.0,5.0,6.0 };
	//double b[6] = { 1,2,3,4,5,6 };
	//Matrix A(2, 3, a);
	//Matrix B(2, 3, b);
	////cout << A.MatrixIsEqual(B) << endl;
	//if (A == B) cout << "equal\n";
	//else cout << "not equal" << endl;


	////����ת��
	//double a[12] = { 1.23456,2,3,4,5,6,7.465456132,8,9,8,7,6 };
	//Matrix A(3, 4, a);
	//Matrix C(4, 3);
	//C = tran(A);
	//cout << "A:\n";
	//A.show();
	//cout << "C=A':\n";
	//C.show();


	////���Ծ��������ʽ
	//double a[25] = { 0.4218  ,  0.0357   , 0.7431  ,  0.0318  ,  0.6948,
	//0.9157  ,  0.8491  ,  0.3922  ,  0.2769 ,   0.3171,
	//0.7922  ,  0.9340  ,  0.6555   , 0.0462  ,  0.9502,
	//0.9595  ,  0.6787  ,  0.1712  ,  0.0971  ,  0.0344,
	//0.6557  ,  0.7577  ,  0.7060   , 0.8235  ,  0.4387 };
	//Matrix A(5, 5, a);
	//cout << det(A)<<endl;
	////cout << A.MatrixDeterminant(A.getrow(), A.getarray());
	//A.show();

	////���Ծ������
	///*double a[25] = { 0.4218  ,  0.0357   , 0.7431  ,  0.0318  ,  0.6948,
	//0.9157  ,  0.8491  ,  0.3922  ,  0.2769 ,   0.3171,
	//0.7922  ,  0.9340  ,  0.6555   , 0.0462  ,  0.9502,
	//0.9595  ,  0.6787  ,  0.1712  ,  0.0971  ,  0.0344,
	//0.6557  ,  0.7577  ,  0.7060   , 0.8235  ,  0.4387 };*/
	//double a[12] = {
	//	 1   ,  2  ,   3   ,  4,
	//	 4   ,  5  ,   9   ,  6,
	//	 3   ,  5  ,   8   ,  6
	//};
	//Matrix A(4, 3, a);
	//cout << ranka(A)<<endl;
	//A.show();

	////���Ծ���ķ���
	//double a[12] = {
	//	 1   ,  2  ,   3   ,  4,
	//	 4   ,  5  ,   -9   ,  6,
	//	 3   ,  5  ,   8   ,  6
	//};
	//double b[12] = {
	//		 2   ,  2  ,   3   ,  4,
	//		 4   ,  5  ,   -9   ,  6,
	//		 3   ,  5  ,   8   ,  6
	//};	
	//Matrix A(3,4, a);
	////cout << A.Matrix_InfNorm(A)<<endl;
	//cout << "1-norm: "<< norm(A, 1) << endl;
	//cout << "2-norm: " << norm(A, 2) << endl;
	//cout << "inf-norm: " << norm(A, INF)<<endl;
	//A.show();	

	////�����������ֵ����������
	//const int N = 3;
	//double a[N*N] = { 1  ,   2   ,  3,
	//			 3  ,   1   ,  2,
	//			 2 ,    3  ,   1 },
	//er[N] = { 0 }, ei[N] = { 0 };

	//int ma = N;

	//Matrix A(ma, ma, a), v;
	//A.MatrixEig(er, ei, v);
	//printf("\nA:\n");
	//A.show();
	//print_matrix("er:", 1, N, er, ma);
	//print_matrix("ei:", 1, N, ei, ma);
	//printf("\nV:\n");
	//v.show();
	
	////���Ծ���ļ�
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//Matrix A(3, 3, a);
	//cout << A.MatrixTrace()<<endl;

	////���Ծ���LU�ֽ�
	///* Locals */
	//lapack_int m = M, n = N, lda = LDA, info;
	///* Local arrays */
	///*
	//double a[LDA*M] = {
	//8.79,  9.93,  9.83, 5.45,  3.16,
	//6.11,  6.91,  5.04, -0.27,  7.98,
	//-9.15, -7.93,  4.86, 4.85,  3.01,
	//9.57,  1.64,  8.83, 0.74,  5.80,
	//-3.49,  4.02,  9.80, 10.00,  4.27,
	//9.84,  0.15, -8.99, -6.02, -5.31
	//};*/

	//double a[LDA*M] = {
	//6.11,  6.91,  5.04, -0.27,  7.98,
	//-9.15, -7.93,  4.86, 4.85,  3.01,
	//9.57,  1.64,  8.83, 0.74,  5.80,
	//-3.49,  4.02,  9.80, 10.00,  4.27,
	//9.84,  0.15, -8.99, -6.02, -5.31
	//};

	//lapack_int ipiv[M + 1] = { 0,0,0,0,0,0,0 };
	//Matrix A(n, n, a), l, u;
	//printf("A\n");
	//A.show();
	//printf("A2\n");
	//info = A.MatrixLUfactorization(l, u, ipiv);
	//A.show();
	//printf("L\n");
	//l.show();
	//printf("U\n");
	//u.show();
	//printf("\n P\n");
	//for (lapack_int j = 0; j < M + 1; j++) printf("%6d", ipiv[j]);

	////���Ծ���QR�ֽ�
	////�𰸲�Ψһ
	//double a[5 * 3] = { 1,2,3,
	//				  4,5,6,
	//				  7,8,9,
	//				  10,11,12,
	//				  13,14,15 };

	//int ma = 5, na = 3;

	//Matrix A(ma, na, a), q, r;
	//A.MatrixQRfactorization(q, r);
	//printf("A:\n");
	//A.show();
	//printf("Q:\n");
	//q.show();
	//printf("R:\n");
	//r.show();

	////���Ծ���SVD�ֽ�
	//double s[N];
	//double a[LDA*M] = {
	//	8.79,  9.93,  9.83, 5.45,  3.16,
	//	6.11,  6.91,  5.04, -0.27,  7.98,
	//   -9.15, -7.93,  4.86, 4.85,  3.01,
	//	9.57,  1.64,  8.83, 0.74,  5.80,
	//   -3.49,  4.02,  9.80, 10.00,  4.27,
	//	9.84,  0.15, -8.99, -6.02, -5.31
	//};
	//Matrix A(M, N, a), U, Vt;
	//A.show();
	//A.MatrixSVDfactorization(s, U, Vt);
	////getchar();
	//cout << "A:\n";
	//A.show();
	//cout << "U:\n";
	//U.show();
	//cout << "V:\n";
	//Vt.show();
	//print_matrix("Singular values", 1, N, s, 1);

	return 0;
}

void print_matrix(const char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda) {
	lapack_int i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) printf(" %10.6f", a[i*lda + j]);
		printf("\n");
	}
}