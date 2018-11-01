#pragma once
#ifndef VECTOR_H
#define VECTOR_H
#include "cblas.h"
constexpr auto INF = 100;
class Vector {
private:
	int n;
	double *v;
public:
	//two methods to Initialize,choose one
	Vector(int); 
	Vector(int n1, double * a);
	~Vector();

	void setn(int);
	void setv(double *);

	int getn();
	double* getv(); // get the first node of the vector

	void printv();
	//modify the length of the vector,and reset the value of it
	void modifyv(int, double*);
	//destroy the vector
	void free();
	
	/*
		Description: y = x
		Input:	Vector x
		Output: y = x
	*/
	void operator = (const Vector &);
	/*
		Description: z = dot(x,y)
		Input:	Vector x、y
		Output: z = x*y
	*/
	friend double operator *(const Vector &, const Vector &);
	/*
		Description: norm
		Input:	null
		Output: norm-1、2、INF
	*/
	double vector_1_norm();
	double vector_2_norm();
	double vector_inf_norm();
};

/*
	Not Class member function
*/

/*
	Description: ?-norm
	Input: Vector x
		?: 1、2、Inf
	Output: ?-norm
*/
double norm(Vector x, int type);
#endif