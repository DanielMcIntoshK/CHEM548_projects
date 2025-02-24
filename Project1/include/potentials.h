#ifndef POTENTIALS__H_
#define POTENTIALS__H_

double particleinabox(double x, double l, double r);
double finitewell(double x, double l, double r, double d);
double pbrectangular(double x, double l1, double r1, double l2, double r2,double d);
double harmonic(double x, double k, double Rc);
double morse(double x, double a, double d, double h,double Rc);



#endif

