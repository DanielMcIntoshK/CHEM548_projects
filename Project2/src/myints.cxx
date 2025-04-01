#include "myints.h"


Integrals::Integrals():primes{2,3,5,7,11,13,17,19,23,29,31,37}{
	fcalc.load("tabulated.txt");
}

void Integrals::computeIntegrals(molecule ml){
	orbs=ml.orbitals;
	mol=ml;
	
	std::cout << "CALCULATING OVERLAP\n";
	computeOverlap();
	std::cout << "CALCULATING KINETIC\n";
	computeKinetic();
	std::cout << "CALCULATING POTENTIAL\n";
	computePotential();
	std::cout << "CALCULATING ERIs\n";
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
					double co1=orbs[i].contraction[k].co,
					       co2=orbs[j].contraction[l].co;
					double v = computeOverlap_el(orbs[i].contraction[k],orbs[j].contraction[l]);
					S(i,j)+=n1*n2*co1*co2*v;
					verbose=false;
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
				//std::cout << "\nATTEMPTING: " << i <<" "<< j << std::endl;
				for(int k=0; k < orbs[i].contraction.size();k++){
				for(int l=0; l < orbs[j].contraction.size();l++){
					precalcPotential.clear();
					precalcOverlap.clear();
					double n1=orbs[i].contraction[k].norm,
				       		n2=orbs[j].contraction[l].norm;
					double co1=orbs[i].contraction[k].co,
					       co2=orbs[j].contraction[l].co;
					double v=computePotential_el(orbs[i].contraction[k],
								    orbs[j].contraction[l],0,mol.atoms[a]);
					nucattract[a](i,j)+=n1*n2*co1*co2*v;
				}}
			}
		}
		V=V-(double)mol.atoms[a].an*nucattract[a];
		
	}

}

void Integrals::computeERIs(){
	precalcK();
	std::cout << "Ks PRECALCULATED\n";
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

	for(int a = 0; a < ERIs.size(); a++){
		std::cout << a << " " << ERIs.size()<<std::endl;
	for(int b = 0; b < ERIs[a].size();b++){
	for(int c = 0; c < ERIs[a][b].size();c++){
	for(int d = 0; d < ERIs[a][b][c].size();d++){
		for(int i = 0; i < orbs[a].contraction.size(); i++){
		for(int j = 0; j < orbs[b].contraction.size(); j++){
		for(int k = 0; k < orbs[c].contraction.size(); k++){
		for(int l = 0; l < orbs[d].contraction.size(); l++){
			std::vector<int> cr{a,b,c,d,i,j,k,l};
			precalcERI.clear();
			primitive ap=orbs[a].contraction[i],
				bp=orbs[b].contraction[j],
				cp=orbs[c].contraction[k],
				dp=orbs[d].contraction[l];
			double normfull=ap.norm*bp.norm*cp.norm*dp.norm;
			if(a==0&&b==0&&c==4&&d==4){
				verbose=true;
			}
			ERIs[a][b][c][d]+=normfull*ap.co*bp.co*cp.co*dp.co*computeERI_el(ap,bp,cp,dp,0,cr);
			verbose=false;
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

	dmath::vec A=a.at.center,B=b.at.center;

	double xi=a.exp+b.exp;
	double zeta=(1.0/xi)*a.exp*b.exp;

	if(a.l()==0 && b.l()==0){
		return std::pow(mypi/xi,3.0/2.0)*std::exp(-zeta*(A-B)*(A-B));
	}


	dmath::vec P=1.0/xi*(a.exp*A+b.exp*B);

	std::vector<primitive> np{a,b};
	if(a.l()==0){
		std::swap(np[0],np[1]);
		std::swap(A,B);
	}

	int anx=0;
	while(np[0].ang[anx]==0){anx+=1;}

	np[0]=np[0].r(anx);

	return (P(anx)-A(anx))*computeOverlap_el(np[0],np[1])
		+1.0/(2.0*xi)*(double)np[0].ang[anx]*computeOverlap_el(np[0].r(anx),np[1])
		+1.0/(2.0*xi)*(double)np[1].ang[anx]*computeOverlap_el(np[0],np[1].r(anx));
}

double Integrals::computeKinetic_el(primitive a, primitive b){
	if(a.l()<0||b.l()<0){
		return 0;
	}

	dmath::vec A=a.at.center,B=b.at.center;
	double xi=a.exp+b.exp;
	double zeta=(1.0/xi)*a.exp*b.exp;

	if(a.l()==0&&b.l()==0){
		return zeta*(3.0-2.0*zeta*(A-B)*(A-B))*computeOverlap_el(a,b);
	}

	dmath::vec P=(1.0/xi)*(a.exp*A+b.exp*B);

	std::vector<primitive> np{a,b};
	if(a.l()==0){
		std::swap(np[0],np[1]);
		std::swap(A,B);
	}

	int anx=0;
	while(np[0].ang[anx]==0){anx+=1;}

	np[0]=np[0].r(anx);

	return (P(anx)-A(anx))*computeKinetic_el(np[0],np[1])
		+1.0/(2.0*xi)*np[0].ang[anx]*(computeKinetic_el(np[0].r(anx),np[1]))
		+1.0/(2.0*xi)*np[1].ang[anx]*(computeKinetic_el(np[0],np[1].r(anx)))
		+2.0*zeta*(computeOverlap_el(np[0].u(anx),np[1])
				-1.0/(2.0*np[0].exp)*np[0].ang[anx]*computeOverlap_el(np[0].r(anx),np[1]));
}

double Integrals::computePotential_el(primitive a, primitive b, int aux, Atom at){
	if(a.l()<0||b.l()<0){
		return 0;
	}
	
	dmath::vec A=a.at.center, B=b.at.center;

	double xi = a.exp+b.exp;
	dmath::vec P=(1.0/xi)*(a.exp*A+b.exp*B),
		C=at.center;
	double U=xi*(P-C)*(P-C);

	//std::cout << "POTENTIAL: " << xi << " "<< aux << " " << fcalc.Fn(aux,U)<< " " << 
	//	a.exp << " " << b.exp << " : "<<
	//	A << "/////"<<B<<std::endl;
	//std::cout << a.exp << " " << b.exp << " (" << A << ") : (" << B << ") : (" << P <<")"<< std::endl;
	if(a.l()==0&&b.l()==0){
		return 2.0*std::sqrt(xi/mypi)*computeOverlap_el(a,b)*fcalc.Fn(aux,U);
	}

	std::vector<primitive> np{a,b};
	if(a.l()==0){
		std::swap(np[0],np[1]);
		std::swap(A,B);
	}

	int anx=0;
	while(np[0].ang[anx]==0){anx+=1;}

	np[0]=np[0].r(anx);

	return (P(anx)-A(anx))*computePotential_el(np[0],np[1],aux,at)
		-(P(anx)-C(anx))*computePotential_el(np[0],np[1],aux+1,at)
		+1.0/(2.0*xi)*np[0].ang[anx]*(computePotential_el(np[0].r(anx),np[1],aux,at)
				-computePotential_el(np[0].r(anx),np[1],aux+1,at))
		+1.0/(2.0*xi)*np[1].ang[anx]*(computePotential_el(np[0],np[1].r(anx),aux,at)
				-computePotential_el(np[0],np[1].r(anx),aux+1,at));
}

double Integrals::computeERI_el(primitive a, primitive b, primitive c, primitive d, int aux,std::vector<int> &cr){
	if(a.l()<0 || b.l()<0 || c.l()<0||d.l()<0) return 0.0;

	std::vector<primitive> prims{a,b,c,d};
		
	double xi=a.exp+b.exp;
	double zeta=a.exp*b.exp/xi;
	double eta=c.exp+d.exp;

	double rho=xi*eta/(xi+eta);

	dmath::vec P=(1.0/xi)*(a.exp*a.at.center+b.exp*b.at.center);
	dmath::vec Q=(1.0/eta)*(c.exp*c.at.center+d.exp*d.at.center);
	dmath::vec W=(1.0/(xi+eta))*(xi*P+eta*Q);

	double T=rho*(P-Q)*(P-Q);

	if(verbose){
		std::cout << a.l() << " " << b.l() << " " << c.l() << " " << d.l() <<std::endl;
	}

	if(a.l()==0 && b.l()==0 && c.l()==0 && d.l()==0){
		double k1=Ks[cr[0]][cr[1]][cr[4]][cr[5]],
		       k2=Ks[cr[2]][cr[3]][cr[6]][cr[7]];
		return 1.0/std::sqrt(xi+eta)*k1*k2*fcalc.Fn(aux,T);
	}

	int itx=0;
	while(prims[itx].l()==0){itx+=1;}
	int anx=0;
	while(prims[itx].ang[anx]==0){anx+=1;}

	//double nxi=(itx>=2)?eta:xi;
	//double neta=(itx>=2)?xi:eta;
	dmath::vec V1=P, V2=Q;
	std::vector<primitive> np=prims;
	if(itx>=2){
		std::swap(xi,eta);
		std::swap(np[0],np[2]);
		std::swap(np[1],np[3]);
		std::swap(V1,V2);
	}
	if(itx%2==1){
		std::swap(np[0],np[1]);
		//std::swap(np[2],np[3]);
	}
	dmath::vec Acen=np[0].at.center;

	np[0]=np[0].r(anx);

	return (V1(anx)-Acen(anx))*computeERI_el(np[0],np[1],np[2],np[3],aux,cr)
		+(W(anx)-V1(anx))*computeERI_el(np[0],np[1],np[2],np[3],aux+1,cr)
		+1.0/(2.0*xi)*np[0].ang[anx]*(computeERI_el(np[0].r(anx),np[1],np[2],np[3],aux,cr)
						-rho/xi*computeERI_el(np[0].r(anx),np[1],np[2],np[3],aux+1,cr))
		+1.0/(2.0*xi)*np[1].ang[anx]*(computeERI_el(np[0],np[1].r(anx),np[2],np[3],aux,cr)
						-rho/xi*computeERI_el(np[0],np[1].r(anx),np[2],np[3],aux+1,cr))
		+1.0/(2.0*(xi+eta))*np[2].ang[anx]*computeERI_el(np[0],np[1],np[2].r(anx),np[3],aux+1,cr)
		+1.0/(2.0*(xi+eta))*np[3].ang[anx]*computeERI_el(np[0],np[1],np[2],np[3].r(anx),aux+1,cr);
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

double Integrals::K(double sa, double sb, dmath::vec v1, dmath::vec v2){
	dmath::vec diff=v1-v2;
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

void Integrals::precalcK(){
	Ks.resize(orbs.size());
	for(int a = 0; a < orbs.size(); a++){
	Ks[a].resize(orbs.size());
	for(int b = 0; b < orbs.size(); b++){
		Ks[a][b].resize(orbs[a].contraction.size());
		for(int i = 0; i < orbs[a].contraction.size(); i++){
		Ks[a][b][i].resize(orbs[b].contraction.size());
		for(int j = 0; j < orbs[b].contraction.size(); j++){
			Ks[a][b][i][j]=K(orbs[a].contraction[i].exp,
					orbs[b].contraction[j].exp,
					orbs[a].at.center,orbs[b].at.center);
		}}
	}}	
}

