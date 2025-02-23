#include "myints.h"

//Integrals::Integrals():primes{2,3,5,7,11,13,17,19,23,29,31,37}{
//}
Integrals::Integrals(){

}

void Integrals::computeIntegrals(std::vector<atOrb> & orbitals){
	orbs=orbitals;

	computeOverlap();
	computeKinetic();
	computePotential();
	computeERIs();
}

void Integrals::computeOverlap(){
	S=dmath::mat(orbs.size(),orbs.size());
	for(int i = 0; i < S.r;i++){
		for(int j = 0; j < S.c;j++){
			precalc.clear();
			S(i,j)=0.0;
			for(int k=0; k < orbs[i].contraction.size();k++){
				for(int l=0; l < orbs[j].contraction.size();l++){
					S(i,j)+=orbs[i].contraction[k].co*orbs[j].contraction[l].co*
						computeOverlap_el(orbs[i].contraction[k],orbs[j].contraction[l]);
				}
			}
		}
	}
}

void Integrals::computeKinetic(){

}

void Integrals::computePotential(){

}

void Integrals::computeERIs(){

}

double Integrals::computeOverlap_el(primitive a, primitive b){
	if(a.l()<0||b.l()<0){
		return 0;
	}
	double sm=a.exp+b.exp;
	double red=a.exp*b.exp/sm;
	dmath::vec A=a.at.center,B=b.at.center;
	dmath::vec P=(a.exp*A+b.exp*B)*(1.0/sm);
	if(a.l()==0 && b.l()==0){
		dmath::vec diff=A-B;
		double dt=diff*diff;
		return std::pow(3.141592/sm,3.0/2.0)*std::exp(-sm*dt);
	}
	//recursive
	unsigned int hv=hash(std::vector<primitive>{a,b});
	if(precalc.find(hv)!=precalc.end()){
		return precalc[hv];
	}

	bool targeta=a.l()!=0;
	int dim=0;
	for(int i = 0; i < 3;i++){
		if((targeta?a:b).ang[i]>0){
			dim=i;
			break;
		}
	}
	double thexp=(targeta)?a.exp:b.exp;

	primitive an=a,bn=b;
	if(targeta) an.ang[dim]-=1;
	else bn.ang[dim]-=1;

	primitive ad=an,bd=bn;
	ad.ang[dim]-=1;
	bd.ang[dim]-=1;

	double val=(P(dim)-a.at.center(dim))*computeOverlap_el(an,bn)+
		1.0/(2.0*thexp)*ad.ang[dim]*computeOverlap_el(ad,bn)+
		1.0/(2.0*thexp)*bd.ang[dim]*computeOverlap_el(an,bd);

	precalc[hv]=val;

	return val;
}

unsigned int Integrals::hash(std::vector<primitive> orbs){
	int hashv=0;
	for(int i = 0; i < orbs.size(); i++){
		for(int j=0; j < 3; j++){
			hashv+=std::pow(primes[i*3+j],orbs[i].ang[j]);
		}
	}
	return hashv;
}

