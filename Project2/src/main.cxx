#include <iostream>
#include <fstream>
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
	std::cout << "MOLECULE FILE READ\n";

	
	//double sqrt3 =std::sqrt(3);
	//std::vector<int> dr{14,15,17};
	//for(int k = 0; k < dr.size(); k++){
	//for(int i = 0; i < mol.orbitals[dr[k]].contraction.size();i++){
	//	mol.orbitals[dr[k]].contraction[i].norm/=sqrt3;
	//}}
	
	//mol.swapbs(std::vector<int>{0,2,4,5,10,8,9,6,7,11,12,13,14,15,16,17,18,1,3});

	
	Integrals ic;
	
	ic.computeIntegrals(mol);

	std::cout << "OVERLAP\n"<<ic.S<<std::endl;
	std::cout << "KINETIC\n"<<ic.T<<std::endl;
	std::cout << "POTENTIAL\n"<<ic.V<<std::endl;

	HF hartreesolve;
	HF::result rs=hartreesolve.restrictedHF(mol,ic);

	std::cout << "HF RESULT: ELECTRIC: " << rs.eelec << " TOTAL: " << rs.eelec+rs.enuc <<std::endl;

	return 0;
}
