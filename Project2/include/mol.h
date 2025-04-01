#ifndef MOL__H_
#define MOL__H_
#include "mymath.h"
#include <array>
#include <string>
#include <map>

struct Atom{
	dmath::vec center;
	int an;

	static std::map<std::string,int> elements;
};

struct primitive{
	int ang[3];
	int l(){return ang[0]+ang[1]+ang[2];}

	primitive r(int k){
		primitive p;
		p.norm=norm;
		p.exp=exp;
		p.co=co;
		p.at=at;
		for(int i = 0; i < 3; i++){
			p.ang[i]=ang[i];
		}
		p.ang[k]-=1;
		return p;
	}
	primitive u(int k){
		primitive p;
		p.norm=norm;
		p.exp=exp;
		p.co=co;
		p.at=at;
		for(int i = 0; i < 3; i++){
			p.ang[i]=ang[i];
		}
		p.ang[k]+=1;
		return p;
	}

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
	int charge=0;
	//For now just assume spin=1
	int spin=1;

	int nelec;

	void read(std::string filename);
	void buildbasis(std::string filename);

	void swapbs(std::vector<int> swapvec);
};



#endif

