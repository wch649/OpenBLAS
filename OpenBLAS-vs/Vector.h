#pragma once
#ifndef VECTOR_H
#define VECTOR_H
#include "cblas.h"
#define INF -1
class Vector {
private:
	int n;
	double *v;
public:
	Vector(int); //two methods to Initialize,choose one
	void setn(int);
	~Vector();
	void setv();
	int getn();
	double* getv(); // get the first node of the vector
	void printv();
	void modifyv(); //modify the length of the vector,and reset the value of it
	
	Vector operator = (Vector &);
	friend double operator *(Vector &, Vector &);

	double vector_1_norm();
	double vector_2_norm();
	double vector_inf_norm();
};

double norm(Vector x, int type);
#endif