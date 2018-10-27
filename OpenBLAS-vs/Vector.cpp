#include "Vector.h"
#include <iostream>
#include <malloc.h>
#include <math.h>
using namespace std;

#define INF -1

Vector::Vector(int n1)
{
	n = n1;
	v = (double *)malloc(n * sizeof(double));
}

void Vector::setn(int n1)
{
	n = n1;
	v = (double *)malloc(n * sizeof(double));
}

Vector::~Vector()
{
	//delete []v;
}
void Vector::setv()
{
	int t = n;
	double *w = v;
	cout << "input " << n << " numbers of this vector" << endl;
	while (t) {
		cin >> *w;
		w++;
		t--;
	}
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
void Vector::modifyv()
{
	delete[]v;
	setn(3);
	setv();
}


Vector Vector::operator = (Vector &a)
{
	cblas_dcopy(a.getn(), a.getv(), 1, this->getv(), 1);
	return *this;
}

double operator *(Vector &a, Vector &b)
{
	return cblas_ddot(a.getn(), a.getv(), 1, b.getv(), 1);
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
