#include "Vector.h"
#include <iostream>
#include <malloc.h>
#include <math.h>
using namespace std;


Vector::Vector(int n1)
{
	n = n1;
	v = (double *)malloc(n * sizeof(double));
}
Vector::Vector(int n1,double *a)
{
	n = n1;
	v = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) v[i] = a[i];
}
Vector::~Vector()
{
	//delete []v;
}


void Vector::setn(int n1)
{
	n = n1;
}

void Vector::setv(double *a)
{
	v = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) v[i] = a[i];
}
int Vector::getn()
{
	return n;
}

double* Vector::getv()
{
	return v;
}

void Vector::printv()
{
	int t = n;
	double *w = v;
	cout << "the numbers of this vector" << endl;
	while (t) {
		cout << *w << " ";
		w++;
		t--;
	}
	cout << endl;
}
void Vector::modifyv(int nl, double *a)
{
	setn(nl);
	if(v) delete[]v;
	setv(a);
}

void Vector::free()
{
	n = 0;
	delete[]v;
}


void Vector::operator = (const Vector &a)
{
	n = a.n;
	v = (double *)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) v[i] = a.v[i];
}

double operator *(const Vector &a, const Vector &b)
{
	return cblas_ddot(a.n, a.v, 1, b.v, 1);
}

double Vector::vector_1_norm()
{
	return cblas_dasum(this->n, this->v, 1);
}

double Vector::vector_2_norm()
{
	return cblas_dnrm2(this->n, this->v, 1);
}

double Vector::vector_inf_norm()
{
	return fabs(v[cblas_idamax(this->n, this->v, 1)]);
}

/*
	Not Class member function
*/

double norm(Vector x, int type)
{
	switch (type) {

	case 1: return x.vector_1_norm(); break;
	case 2: return x.vector_2_norm(); break;
	case INF:return x.vector_inf_norm(); break;
	default:cout << "Parameter error!\n";
	}
	return 0.0;
}
