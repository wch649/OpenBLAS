#include <stdio.h>
#include "Matrix.h"
#include <cblas.h>
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
	//A.MatrixInverse(A.getarray(), A.getrow());
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
	//double a[12] = { 1,2,3,4,5,6,7,8,9,8,7,6 };
	//Matrix A(3, 4, a);
	//Matrix C(4, 3);
	//C = tran(A);
	//cout << "A:\n";
	//A.show();
	//cout << "C=A':\n";
	//C.show();


	/*double a[4] = { 1,2,3,4 };
	Vector ()
	cout << norm();*/
	return 0;
}