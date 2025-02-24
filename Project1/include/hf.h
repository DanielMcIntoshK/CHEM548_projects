#ifndef HF__H_
#define HF__H_
#include "mymath.h"
#include "mol.h"
#include "myints.h"

class HF{
public:
	struct result{
		double eelec;
		double enuc;

		dmath::mat C;
		dmath::mat D;
		dmath::mat G;
		dmath::mat F;
		dmath::mat E;

		int nelec;
		int norbs;
	};
public:
	HF(){}

	HF::result restrictedHF(molecule & mol, Integrals & ic);

private:
	dmath::mat computeGMat(dmath::mat & D, Integrals & ic);
};

#endif

