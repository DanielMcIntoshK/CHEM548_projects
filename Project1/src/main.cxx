#include <iostream>
#include "mymath.h"
#include "myints.h"
#include "mol.h"
#include "hf.h"

int main(int argc, char ** argv){

	if(argc==1) {
		std::cout << "NO INPUT FILE\n";
		return 1;
	}
	molecule mol;
	mol.read(argv[1]);

	Integrals ic;
	ic.computeIntegrals(mol);
	std::cout << "OVERLAP\n"<<ic.S<<std::endl<<std::endl;

	/*
	std::cout << "KINETIC\n"<<ic.T<<std::endl<<std::endl;
	std::cout << "POTENTIAL\n"<<ic.V<<std::endl<<std::endl;
	
	std::cout << "TWO BODY INTS\n";
	int n = mol.orbitals.size();
	for(int i = 0; i < n;i++){
	for(int j = 0; j < n;j++){
	for(int k = 0; k < n; k++){
	for(int l = 0; l < n; l++){
		std::cout << "("<<i<<j<<"|"<<k<<l<<") "<<ic.getERI(i,j,k,l)<<std::endl;
	}}}}
	*/

	dmath::mat U;
	dmath::JacobiDiagonalizer jd;
	dmath::mat d=jd.diagonalize(ic.S,0.0000001,U);
	for(int i = 0; i < d.r; i++){
	for(int j = 0; j < d.c; j++){
		if(std::abs(d(i,j))<0.0000001) d(i,j)=0.0;
	}
	}
	dmath::mat ident=dmath::mat::Identity(U.r,U.c);
	dmath::mat check = U*d*U.transpose();
	std::cout << d << std::endl << std::endl << U << std::endl << std::endl << check << std::endl;

	HF hartreesolve;
	HF::result rs=hartreesolve.restrictedHF(mol,ic);

	std::cout << "HF RESULT: ELECTRIC: " << rs.eelec << " TOTAL: " << rs.eelec+rs.enuc <<std::endl;

	return 0;
}
