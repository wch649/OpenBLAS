#include <stdio.h>
#include "Matrix.h"
extern "C" {
#include "cblas.h"
}
#include <iostream>
using namespace std;

int main(void) {

	////≤‚ ‘øΩ±¥
	//double a[9] = { 0,1,0, 1,0,1, -1,0,1 };
	//double b[9] = { 1,0,1,0,1,0,0,1,0 };
	//Matrix A(3, 3, a);
	//Matrix B(3, 3, b);
	//A.show();
	//B.show();
	//B = A;
	//B.show();
	//B = B + A;
	//B.show();
	//A.show();


	////≤‚ ‘º”ºı∑®
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

	////≤‚ ‘≥À∑®
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

	////≤‚ ‘«ÛƒÊ
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

	////≤‚ ‘”“≥˝
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


	////≤‚ ‘◊Û≥˝
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
	//cout << "C=A\B:\n";
	//C.show();

	
	////≤‚ ‘◊™÷√
	//double a[12] = { 1.23456,2,3,4,5,6,7.465456132,8,9,8,7,6 };
	//Matrix A(3, 4, a);
	//Matrix C(4, 3);
	//C = tran(A);
	//cout << "A:\n";
	//A.show();
	//cout << "C=A':\n";
	//C.show();

	////≤‚ ‘æÿ’Ûµƒ––¡– Ω
	//double a[25] = { 0.4218  ,  0.0357   , 0.7431  ,  0.0318  ,  0.6948,
	//0.9157  ,  0.8491  ,  0.3922  ,  0.2769 ,   0.3171,
	//0.7922  ,  0.9340  ,  0.6555   , 0.0462  ,  0.9502,
	//0.9595  ,  0.6787  ,  0.1712  ,  0.0971  ,  0.0344,
	//0.6557  ,  0.7577  ,  0.7060   , 0.8235  ,  0.4387 };
	//Matrix A(5, 5, a);
	//cout << det(A)<<endl;
	////cout << A.MatrixDeterminant(A.getrow(), A.getarray());
	//A.show();

	////≤‚ ‘æÿ’Ûµƒ÷»
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

	////≤‚ ‘æÿ’Ûµƒ∑∂ ˝
	//double a[12] = {
	//	 1   ,  2  ,   3   ,  4,
	//	 4   ,  5  ,   -9   ,  6,
	//	 3   ,  5  ,   8   ,  6
	//};
	//Matrix A(3,4, a);
	////Matrix A(4,3, a);
	////cout << A.Matrix_InfNorm(A)<<endl;
	//cout << "1-norm: "<< norm(A, 1) << endl;
	//cout << "2-norm: " << norm(A, 2) << endl;
	//cout << "inf-norm: " << norm(A, INF)<<endl;
	//A.show();	
	return 0;
}