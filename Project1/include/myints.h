#ifndef INTEGRALS__H_
#define INTEGRALS__H_
#include <cmath>
#include "mymath.h"
#include "mol.h"
#include <map>

typedef std::vector<std::vector<std::vector<std::vector<double>>>> tblist;

class Integrals{
public:
	Integrals();

	void computeIntegrals(std::vector<atOrb> & orbitals);

	void computeOverlap();
	void computeKinetic();
	void computePotential();
	void computeERIs();

	dmath::mat S,T,V;
	tblist ERIs;

	std::vector<atOrb> orbs;
private:
	double computeOverlap_el(primitive a, primitive b);

	std::map<unsigned int,double> precalc;
	
	std::vector<double> primes;

	unsigned int hash(std::vector<primitive> orbs);
};


#endif

