#include "fd.h"

FD::FD(){}

void FD::compute(double m, double dx, dmath::vec V){
	double dx2=dx*dx;

	dmath::mat H=dmath::mat::Zero(V.vals.size(),V.vals.size());
	
	for(int i = 0; i < V.vals.size(); i++){
		H(i,i)=1.0/(m*dx2)+V(i);
		if(i+1<V.vals.size()) H(i,i+1)=-1.0/(2.0*m*dx2);
		if(i-1>=0)H(i,i-1)=-1.0/(2.0*m*dx2);
	}

	//std::cout << "HMAT\n" << H <<std::endl;

	dmath::JacobiDiagonalizer jd;
	eigenvals=jd.diagonalize(H,0.000000001,eigenvecs);
}

