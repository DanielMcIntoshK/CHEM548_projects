#ifndef INTEGRALS__H_
#define INTEGRALS__H_
#include <cmath>
#include <array>
#include "mymath.h"
#include "mol.h"
#include <map>

typedef std::vector<std::vector<std::vector<std::vector<double>>>> tblist;

class Integrals{
public:
	Integrals();

	void computeIntegrals(molecule mol);

	void computeOverlap();
	void computeKinetic();
	void computePotential();
	void computeERIs();

	dmath::mat S,T,V;
	tblist ERIs;

	std::vector<atOrb> orbs;
	molecule mol;

	double getERI(int a, int b, int c, int d);
private:
	double computeOverlap_el(primitive a, primitive b);
	double computeKinetic_el(primitive a, primitive b);
	double computePotential_el(primitive a, primitive b,int aux,Atom at);

	double computeERI_el(primitive a, primitive b, primitive c, primitive d, int aux);

	std::map<unsigned long,double> precalcOverlap;
	std::map<unsigned long,double> precalcKinetic;
	std::map<unsigned long,double> precalcPotential;
	std::map<unsigned long,double> precalcERI;
	
	std::vector<double> primes;

	std::array<int,4> getcord(int a, int b, int c, int d);

	unsigned long hash(std::vector<primitive> orbs,int aux=0);
	double Fa(int m, double a, double T);
	double K(double sa, double sb, dmath::vec v1, dmath::vec v2);
};


#endif

