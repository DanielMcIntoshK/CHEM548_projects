#ifndef MOL__H_
#define MOL__H_
#include "mymath.h"
#include <array>
#include <string>

struct Atom{
	dmath::vec center;
	int an;
};

struct primitive{
	int ang[3];
	int l(){return ang[0]+ang[1]+ang[2];}

	double norm;
	double exp;
	double co;

	Atom at;
};

struct atOrb{
	Atom at;
	
	std::vector<primitive> contraction;
};

struct molecule{
	std::vector<Atom> atoms;
	std::vector<atOrb> orbitals;
	int charge;
	//For now just assume spin=0
	int spin;

	void read(std::string filename);
};

#endif

