#include "potentials.h"
#include <cmath>

double particleinabox(double x, double l, double r){
	if(x<l || x >r) return 9999;
	return 0.0;
}

double finitewell(double x, double l, double r, double d){
	if(x<l || x > r) return d;
	return 0.0;
}

double pbrectangular(double x, double l1, double r1, double l2, double r2, double d){
	if(x<l1||x>r1) return d*5;
	if(x>l2&&x<r2) return d;
	return 0.0;
}

double harmonic(double x, double k, double Rc){
	double diff=x-Rc;
	return 0.5*k*(diff*diff);
}

double morse(double x, double a, double d, double h,double Rc){
	return d*std::pow(1.0-std::exp(-a*(x-Rc)),2)+h;
}

