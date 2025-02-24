#ifndef FD__H_
#define FD__H_
#include "mymath.h"

class FD{
public:
	FD();

	void compute(double m, double dx, dmath::vec V);

	dmath::mat eigenvals, eigenvecs;
};

#endif

