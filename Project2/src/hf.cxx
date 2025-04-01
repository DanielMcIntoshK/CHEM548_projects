#include "hf.h"
#include <iostream>
#include <iomanip>

HF::result HF::restrictedHF(molecule & mol, Integrals & ic){
	HF::result res;

	res.nelec=0;
	for(auto & atom: mol.atoms) res.nelec+=atom.an;
	res.nelec-=mol.charge;

	res.enuc=0.0;
	for(int i = 0; i < mol.atoms.size(); i++){
		for(int j = i+1; j< mol.atoms.size();j++){
			Atom & a1=mol.atoms[i],&a2=mol.atoms[j];
			dmath::vec diff=a1.center-a2.center;
			double dot=diff*diff;
			double r = std::sqrt(dot);

			res.enuc+=a1.an*a2.an/r;
		}
	}

	dmath::mat & S = ic.S;
	dmath::mat & T = ic.T;
	dmath::mat & V = ic.V;

	dmath::JacobiDiagonalizer jd;
	dmath::mat U;
	dmath::mat s=jd.diagonalize(S,0.00000001,U);


	for(int i = 0; i < s.r; i++) s(i,i)=1.0/std::sqrt(s(i,i));
	dmath::mat X = U*s*U.transpose();

	dmath::mat H = T+V;

	std::cout << "CORE HAMILTONIAN\n" << H << std::endl;
	std::cout << "XMat:\n"<<X<<std::endl;

	//s=dmath::mat(0,0);
	//U=dmath::mat(0,0);

	dmath::mat D = dmath::mat::Zero(H.r,H.c);

	int maxiter=200;
	double converge=0.000001;
	int iter=0;
	double rmsd=0.0;
	double ediff=0.0;
	double ehf=0.0;

	dmath::mat C,G,F,E;

	std::cout << "STARTING SCF\n";
	do{
		if(++iter>maxiter){
			std::cout << "SCF FAILED TO CONVERGE IN " << iter-1 << "STEPS" << std::endl;
			break;
		}

		double ehf_last=ehf;
		dmath::mat D_last = D;
		
		G=computeGMat(D,ic);
		F=H+G;

		dmath::mat Fx=X.transpose()*F*X;

		dmath::mat Cx;
		E=jd.diagonalize(Fx, 0.0000001,Cx,true);
		//Cx=-Cx;
		C=X*Cx;


		int ndocc=res.nelec/2;
		dmath::mat Ct=C.transpose();
		for(int i = 0; i < D.r; i++){
			for(int j = 0; j < D.c; j++){
				D(i,j)=0.0;
				for(int a = 0; a < ndocc; a++){
					D(i,j)+=2*C(i,a)*C(j,a);
				}
			}
		}

		ehf=0.0;
		for(int i = 0; i < mol.orbitals.size(); i++){
			for(int j = 0; j < mol.orbitals.size(); j++){
				ehf+=.5*D(i,j)*(H(i,j)+F(i,j));
			}
		}

		ediff=ehf-ehf_last;

		std::cout << std::setprecision(10) << "SCF STEP " << iter << "\tENERGY: " << ehf << "\tEDIFF: " << ediff <<std::endl;
	}while((std::abs(ediff)>converge));

	res.eelec=ehf;
	res.C=C;
	res.D=D;
	res.E=E;
	res.G=G;
	res.F=F;

	res.norbs=C.r;

	return res;
}

dmath::mat HF::computeGMat(dmath::mat & D, Integrals & ic){
	int n = ic.ERIs.size();
	
	dmath::mat G = dmath::mat::Zero(n,n);
	for(int i = 0; i < n; i++){
	for(int j = 0; j < n; j++){
	for(int s = 0; s < n; s++){
	for(int r = 0; r < n; r++){
		double J = ic.getERI(i,j,r,s);
		double K = ic.getERI(i,s,r,j);

		G(i,j)+=D(s,r)*(J-.5*K);
	}}}}

	return G;

}
