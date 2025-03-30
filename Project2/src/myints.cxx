#include "myints.h"

Integrals::Integrals():primes{2,3,5,7,11,13,17,19,23,29,31,37}{
}

void Integrals::computeIntegrals(molecule ml){
	orbs=ml.orbitals;
	mol=ml;

	computeOverlap();
	computeKinetic();
	computePotential();
	computeERIs();
}

void Integrals::computeOverlap(){
	S=dmath::mat(orbs.size(),orbs.size());
	for(int i = 0; i < S.r;i++){
		for(int j = 0; j < S.c;j++){
			S(i,j)=0.0;
			for(int k=0; k < orbs[i].contraction.size();k++){
			for(int l=0; l < orbs[j].contraction.size();l++){
					precalcOverlap.clear();
					double n1=orbs[i].contraction[k].norm,
					       n2=orbs[j].contraction[l].norm;
					S(i,j)+=n1*n2*orbs[i].contraction[k].co*orbs[j].contraction[l].co*
						computeOverlap_el(orbs[i].contraction[k],orbs[j].contraction[l]);
			}}
		}
	}
}

void Integrals::computeKinetic(){
	T=dmath::mat(orbs.size(),orbs.size());
	for(int i = 0; i < S.r;i++){
		for(int j = 0; j < S.c;j++){
			T(i,j)=0.0;
			for(int k=0; k < orbs[i].contraction.size();k++){
			for(int l=0; l < orbs[j].contraction.size();l++){
				precalcKinetic.clear();
				precalcOverlap.clear();
				double n1=orbs[i].contraction[k].norm,
				       n2=orbs[j].contraction[l].norm;
				T(i,j)+=n1*n2*orbs[i].contraction[k].co*orbs[j].contraction[l].co*
					computeKinetic_el(orbs[i].contraction[k],orbs[j].contraction[l]);
			}}
		}
	}
}

void Integrals::computePotential(){
	std::vector<dmath::mat> nucattract;
	nucattract.resize(mol.atoms.size());
	for(auto & mt:nucattract) mt.cSize(orbs.size(),orbs.size());

	V=dmath::mat::Zero(orbs.size(),orbs.size());

	for(int a = 0; a < mol.atoms.size(); a++){
		for(int i = 0; i < S.r;i++){
			for(int j = 0; j < S.c;j++){
				nucattract[a](i,j)=0.0;
				for(int k=0; k < orbs[i].contraction.size();k++){
				for(int l=0; l < orbs[j].contraction.size();l++){
					precalcPotential.clear();
					precalcOverlap.clear();
					double n1=orbs[i].contraction[k].norm,
				       		n2=orbs[j].contraction[l].norm;
					nucattract[a](i,j)+=n1*n2*orbs[i].contraction[k].co*orbs[j].contraction[l].co*
						computePotential_el(orbs[i].contraction[k],orbs[j].contraction[l],0,mol.atoms[a]);
				}}
			}
		}
		V=V-(double)mol.atoms[a].an*nucattract[a];
		
	}

}

void Integrals::computeERIs(){
	std::cout << "COMPUTING ERIS\n";
	ERIs.resize(orbs.size());
	for(int i = 0; i < ERIs.size();i++){
		ERIs[i].resize(i+1);
		for(int j = 0; j < ERIs[i].size(); j++){
			ERIs[i][j].resize(i+1);
			for(int k = 0; k < ERIs[i][j].size();k++){
				int lsize=(i==k)?j:k;
				ERIs[i][j][k].resize(lsize+1);
				for(int l = 0; l < ERIs[i][j][k].size();l++){
					ERIs[i][j][k][l]=0.0;
				}
			}
		}
	}

	std::cout << "ATTEMPTING\n";
	for(int a = 0; a < ERIs.size(); a++){
	for(int b = 0; b < ERIs[a].size();b++){
	for(int c = 0; c < ERIs[a][b].size();c++){
	for(int d = 0; d < ERIs[a][b][c].size();d++){
		//std::cout << "ATTEMPTING: (" << a<<b<< "|" << c << d << ")"<<std::endl;
		for(int i = 0; i < orbs[a].contraction.size(); i++){
		for(int j = 0; j < orbs[b].contraction.size(); j++){
		for(int k = 0; k < orbs[c].contraction.size(); k++){
		for(int l = 0; l < orbs[d].contraction.size(); l++){
			precalcERI.clear();
			primitive ap=orbs[a].contraction[i],
				bp=orbs[b].contraction[j],
				cp=orbs[c].contraction[k],
				dp=orbs[d].contraction[l];
			double normfull=ap.norm*bp.norm*cp.norm*dp.norm;
			ERIs[a][b][c][d]+=normfull*ap.co*bp.co*cp.co*dp.co*computeERI_el(ap,bp,cp,dp,0);
		}}}}
	}}}}
}

double Integrals::getERI(int a, int b, int c, int d){
	std::array<int,4> crd=getcord(a,b,c,d);
	return ERIs[crd[0]][crd[1]][crd[2]][crd[3]];
}

double Integrals::computeOverlap_el(primitive a, primitive b){
	if(a.l()<0||b.l()<0){
		return 0;
	}
	double sm=a.exp+b.exp;
	double rd=a.exp*b.exp/sm;
	dmath::vec A=a.at.center,B=b.at.center;
	
	dmath::vec P=(a.exp*A+b.exp*B)*(1.0/sm);
	if(a.l()==0 && b.l()==0){
		dmath::vec diff=A-B;
		double dt=diff*diff;
		return std::pow(mypi/sm,3.0/2.0)*std::exp(-rd*dt);
	}
	//recursive
	unsigned long hv=hash(std::vector<primitive>{a,b});
	if(precalcOverlap.find(hv)!=precalcOverlap.end()){
		return precalcOverlap[hv];
	}

	bool targeta=a.l()!=0;
	int dim=0;
	for(int i = 0; i < 3;i++){
		if((targeta?a:b).ang[i]>0){
			dim=i;
			break;
		}
	}
	//double thexp=(targeta)?a.exp:b.exp;
	//double thexp=a.exp+b.exp;

	primitive an=a,bn=b;
	if(targeta) an.ang[dim]-=1;
	else bn.ang[dim]-=1;

	primitive ad=an,bd=bn;
	ad.ang[dim]-=1;
	bd.ang[dim]-=1;

	double val=(P(dim)-(targeta?A:B)(dim))*computeOverlap_el(an,bn)+
		1.0/(2.0*sm)*an.ang[dim]*computeOverlap_el(ad,bn)+
		1.0/(2.0*sm)*bn.ang[dim]*computeOverlap_el(an,bd);

	precalcOverlap[hv]=val;

	return val;
}

double Integrals::computeKinetic_el(primitive a, primitive b){
	if(a.l()<0||b.l()<0){
		return 0;
	}
	double sm=a.exp+b.exp;
	double rd=a.exp*b.exp/sm;
	dmath::vec A=a.at.center,B=b.at.center;
	dmath::vec P=(a.exp*A+b.exp*B)*(1.0/sm);
	
	if(a.l()==0 && b.l()==0){
		dmath::vec diff=A-B;
		double dt=diff*diff;
		double over=computeOverlap_el(a,b);
		return rd*(3-2*rd*dt)*over;
	}

	unsigned long hv=hash(std::vector<primitive>{a,b});
	if(precalcKinetic.find(hv)!=precalcKinetic.end()) return precalcKinetic[hv];

	bool targeta=a.l()!=0;
	int dim=0;
	for(int i = 0; i < 3; i++){
		if((targeta?a:b).ang[i]>0){
			dim=i;
			break;
		}
	}
	//double thexp=(targeta)?a.exp:b.exp;
	
	primitive an=a,bn=b;
	if(targeta) an.ang[dim]-=1;
	else bn.ang[dim]-=1;

	primitive ad=an,bd=bn;
	ad.ang[dim]-=1;
	bd.ang[dim]-=1;
	
	precalcOverlap.clear();
	double O1=computeOverlap_el((targeta)?a:an,(targeta)?bn:b);
	//precalcOverlap.clear();
	double O2=computeOverlap_el((targeta)?ad:an,(targeta)?bn:bd);

	double val=(P(dim)-(targeta?A:B)(dim))*computeKinetic_el(an,bn)+
			1.0/(2.0*sm)*an.ang[dim]*computeKinetic_el(ad,bn)+
			1.0/(2.0*sm)*bn.ang[dim]*computeKinetic_el(an,bd)+
			2.0*rd*(O1-1.0/(2.0*(targeta?a.exp:b.exp))*(targeta?an:bn).ang[dim]*O2);
	precalcKinetic[hv]=val;
	return val;
}

double Integrals::computePotential_el(primitive a, primitive b, int aux, Atom at){
	if(a.l()<0||b.l()<0){
		return 0;
	}
	double sm=a.exp+b.exp;
	double rd=a.exp*b.exp/sm;
	dmath::vec A=a.at.center,B=b.at.center;
	dmath::vec P=(a.exp*A+b.exp*B)*(1.0/sm);
	dmath::vec C=at.center;
	if(a.l()==0 && b.l()==0){
		dmath::vec diff=P-C;
		double dt=diff*diff;
		double U=sm*dt;
		return 2.0*std::sqrt(sm/mypi)*computeOverlap_el(a,b)*Fa(aux,U,0.000001);
	}
	unsigned long hv=hash(std::vector<primitive>{a,b},aux);
	if(precalcPotential.find(hv)!=precalcPotential.end())return precalcPotential[hv];

	bool targeta=a.l()!=0;
	int dim=0;
	for(int i = 0; i < 3; i++){
		if((targeta?a:b).ang[i]>0){
			dim=i;
			break;
		}
	}
	//double thexp=(targeta)?a.exp:b.exp;

	primitive an=a,bn=b;
	if(targeta) an.ang[dim]-=1;
	else bn.ang[dim]-=1;

	primitive ad=an,bd=bn;
	ad.ang[dim]-=1;
	bd.ang[dim]-=1;
	
	double val=(P(dim)-(targeta?A:B)(dim))*computePotential_el(an,bn,aux,at)-(P(dim)-C(dim))*computePotential_el(an,bn,aux+1,at)+
		1.0/(2.0*sm)*an.ang[dim]*(computePotential_el(ad,bn,aux,at)-computePotential_el(ad,bn,aux+1,at))+
		1.0/(2.0*sm)*bn.ang[dim]*(computePotential_el(an,bd,aux,at)-computePotential_el(an,bd,aux+1,at));
	precalcPotential[hv]=val;
	return val;
}

double Integrals::computeERI_el(primitive a, primitive b, primitive c, primitive d, int aux){
	if(a.l()<0 || b.l()<0 || c.l()<0||d.l()<0) return 0.0;

	std::vector<primitive> eriprims{a,b,c,d};

	std::vector<dmath::vec> vecvec{a.at.center,b.at.center,c.at.center,d.at.center};

	double sm=a.exp+b.exp;
	double rd=a.exp*b.exp/sm;
	double nu=c.exp+d.exp;
	double rho=sm*nu/(sm+nu);

	dmath::vec P=(1.0/sm)*(a.exp*vecvec[0]+b.exp*vecvec[1]);
	dmath::vec Q=(1.0/nu)*(c.exp*vecvec[2]+d.exp*vecvec[3]);
	dmath::vec W=(1.0/(sm+nu))*(sm*P+nu*Q);

	dmath::vec pqdiff=P-Q;
	double pqdot=pqdiff*pqdiff;

	double fullnorm=1.0;
	//for(int i = 0; i < 4; i++) fullnorm*=eriprims[i].norm;

	if(a.l()==0 && b.l()==0 && c.l()==0 && d.l()==0){
		return fullnorm*std::pow(sm+nu,-1.0/2.0)*K(a.exp,b.exp,vecvec[0],vecvec[1])*K(c.exp,d.exp,vecvec[2],vecvec[3])*Fa(aux,rho*pqdot,0.0000001);
	}

	unsigned long hv=hash(eriprims,aux);
	if(precalcERI.find(hv)!=precalcERI.end())return precalcERI[hv];

	int target=-1;
	for(int i=0; i < eriprims.size();i++){
		if(eriprims[i].l()!=0) {
			target=i;
			break;
		}
	}
	int dim =-1;
	for(int i = 0; i < 3; i++){
		if(eriprims[target].ang[i]!=0){
			dim=i;
			break;
		}
	}

	std::vector<primitive> pn=eriprims;
	pn[target].ang[dim]-=1;

	std::vector<primitive> pd=pn;
	for(int i = 0; i < pd.size();i++){
		pd[i].ang[dim]-=1;
	}

	dmath::vec v1=(target<2)?P:Q, v2=(target<2)?Q:P;
	double s1=(target<2)?sm:nu, s2=(target<2)?nu:sm;

	int offtarget=(target+1)%2+target/2;
	int opptarget=target%2+(target/2+1)%2;
	int oppofftarget=(opptarget+1)%2+opptarget/2;

	//Not finished yet
	double val=(v1(dim)-vecvec[target](dim))*computeERI_el(pn[0],pn[1],pn[2],pn[3],aux)+
		(W(dim)-v1(dim))*computeERI_el(pn[0],pn[1],pn[2],pn[3],aux+1)+
		1.0/(2.0*s1)*pn[(target<2)?0:2].ang[dim]*(
			((target<2)?computeERI_el(pd[0],pn[1],pn[2],pn[3],aux):computeERI_el(pn[0],pn[1],pd[2],pn[3],aux))-
			rho/s1*((target<2)?computeERI_el(pd[0],pn[1],pn[2],pn[3],aux+1):computeERI_el(pn[0],pn[1],pd[2],pn[3],aux+1))
			)+
		1.0/(2.0*s1)*pn[(target<2)?1:3].ang[dim]*(
			((target<2)?computeERI_el(pn[0],pd[1],pn[2],pn[3],aux):computeERI_el(pn[0],pn[1],pn[2],pd[3],aux))-
			rho/s1*((target<2)?computeERI_el(pn[0],pd[1],pn[2],pn[3],aux+1):computeERI_el(pn[0],pn[1],pn[2],pd[3],aux+1))
			)+
		1.0/(2*(sm+nu))*pn[(target<2)?2:0].ang[dim]*(
				(target<2)?computeERI_el(pn[0],pn[1],pd[2],pn[3],aux+1):computeERI_el(pd[0],pn[1],pn[2],pn[3],aux+1))+
		1.0/(2*(sm+nu))*pn[(target<2)?3:1].ang[dim]*(
				(target<2)?computeERI_el(pn[0],pn[1],pn[2],pd[3],aux+1):computeERI_el(pn[0],pd[1],pn[2],pn[3],aux+1));

	precalcERI[hv]=val;
	return val;

}

unsigned long Integrals::hash(std::vector<primitive> orbs,int aux){
	unsigned long hashv=0;
	for(int i = 0; i < orbs.size(); i++){
		for(int j=0; j < 3; j++){
			hashv+=std::pow(primes[i*3+j+1],orbs[i].ang[j]);
		}
	}
	hashv+=std::pow(primes[0],aux);
	return hashv;
}

double Integrals::Fa(int m, double a, double T){
	if(a<0.0000000001)return 1.0;
	double base=std::exp(-a);
	double diff=2.0*T;
	double sm=0.0;
	for(int k = 0; diff>T&&k<100;k++){
		double ism=base*sm;
		sm+=dmath::dfactorial2(2*m-1)*std::pow(2.0*a,(double)k)/dmath::dfactorial2(2*m+2*k+1);
		diff=std::abs(ism-sm*base);
	}
	return base*sm;
}

double Integrals::K(double sa, double sb, dmath::vec v1, dmath::vec v2){
	dmath::vec diff=v2-v1;
	double dot = diff*diff;

	return std::sqrt(2)*std::pow(mypi,5.0/4.0)/(sa+sb)*std::exp(-sa*sb/(sa+sb)*dot);
}

std::array<int,4> Integrals::getcord(int a, int b, int c, int d){
	std::array<int,4> cd{a,b,c,d};

	if(cd[1]>cd[0]) std::swap(cd[0],cd[1]);
	if(cd[3]>cd[2]) std::swap(cd[2],cd[3]);
	if(cd[2]>cd[0]){
		std::swap(cd[0],cd[2]);
		std::swap(cd[1],cd[3]);
	}
	if(cd[0]==cd[2]&&cd[3]>cd[1]){
		std::swap(cd[0],cd[2]);
		std::swap(cd[1],cd[3]);
	}
	else if(cd[3]>cd[2]) std::swap(cd[2],cd[3]);
	
	return cd;
}

