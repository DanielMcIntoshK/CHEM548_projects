#include "mymath.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>


void dmath::mat::cSize(int i, int j){
	r=i;
	c=j;
	vals.resize(r);
	for(int ri = 0; ri < i; ri++){
		vals[ri].resize(c);
	}
}

dmath::mat dmath::mat::Zero(int i, int j){
	dmath::mat zmat(i,j);
	for(int r = 0; r < i; r++){
		for(int c = 0; c < j; c++){
			zmat(r,c)=0.0;
		}
	}
	return zmat;
}

dmath::mat dmath::mat::Identity(int i, int j){
	dmath::mat imat(i,j);
	for(int r = 0; r < i; r++){
		for(int c = 0; c < j; c++){
			imat(r,c)=(r==c)?1.0:0.0;
		}
	}
	return imat;
}

std::ostream & dmath::operator<<(std::ostream & os,dmath::vec & v){
	for(int i = 0; i < v.vals.size(); i++) os<<std::fixed << std::setprecision(8)<<v(i)<<" ";
	return os;
}

std::ostream & dmath::operator<<(std::ostream & os,dmath::mat & m){
	for(int i = 0; i < m.r; i++){
		for(int j = 0; j < m.c; j++){
			if(m(i,j)>=0.0){
				os<< " ";
			}
			os << std::fixed << std::setprecision(8)<<m(i,j) << " ";
		}
		os<<std::endl;
	}
	return os;
}

int dmath::factorial(int n){
	int f=1.0;
	for(int i = n; i >1;i--){
		f*=i;
	}
	return f;
}

int dmath::factorial2(int n){
	int f =1.0;
	for(int i = n; i >1;i-=2){
		f*=i;
	}
	return f;
}

double dmath::dfactorial(int n){
	double v=1.0;
	for(int i =n; i >1;i--){
		v*=(double)i;
	}
	return v;
}

double dmath::dfactorial2(int n){
	double v=1.0;
	for(int i = n; i >1; i-=2){
		v*=(double)i;
	}
	return v;
}

dmath::vec dmath::operator*(dmath::vec v, double s){
	dmath::vec vn=v;
	for(int i = 0; i < vn.vals.size(); i++) vn.vals[i]*=s;
	return vn;
}

dmath::vec dmath::operator*(double s, dmath::vec v){
	dmath::vec vn=v;
	for(int i = 0; i < vn.vals.size(); i++) vn.vals[i]*=s;
	return vn;
}

dmath::mat dmath::operator*(dmath::mat m, double s){
	dmath::mat mn=m;
	for(int i = 0; i < mn.r; i++) {
	for(int j = 0; j < mn.c; j++){	
		mn.vals[i][j]*=s;
	}}
	return mn;
}

dmath::mat dmath::operator*(double s, dmath::mat m){
	return m*s;
}

dmath::mat dmath::mat::operator*(dmath::mat m){
	if(c!=m.r) return dmath::mat(0,0);

	dmath::mat mn(r,m.c);
	for(int i = 0; i < r; i++){
		dmath::vec v1(vals[i]);
		for(int j = 0; j < m.c; j++){
			dmath::vec v2(m.c);
			for(int k = 0; k < m.r;k++){
				v2(k)=m(k,j);
			}
			mn(i,j)=v1*v2;
		}
	}
	return mn;
}

dmath::mat dmath::JacobiDiagonalizer::diagonalize(dmath::mat m, double threshold, dmath::mat & U,bool sort){
	U=dmath::mat::Identity(m.r,m.c);
	dmath::mat d = dmath::mat::Zero(m.r,m.c);

	int citer=0;
	int maxiter=m.r*m.r*m.r;

	findMaxUpperTriangle(m);
	while(std::abs(m(maxr,maxc))>threshold && citer++<maxiter){
		jacobiRotate(m,U);
		findMaxUpperTriangle(m);

	}

	if(citer>=maxiter) std::cout << "TIMEOUT\n";
	//U=U.transpose();

	if(sort){
	for(int i = 0; i < m.r; i++){
		int li=i;
		for(int j = i+1; j < m.r; j++){
			if(m(j,j)<m(li,li)){li=j;}
		}
		//std::cout << "Lowest Eigen: " << li << " "<< m(li,li) << std::endl;
		if(li==i) continue;
		std::swap(m(i,i),m(li,li));
		for(int k = 0; k<U.r; k++){
			//std::cout << U(i,k) << " " << U(li,k) << std::endl;
			std::swap(U(k,i),U(k,li));
		}
	}
	}

	return m;
}

void dmath::JacobiDiagonalizer::findMaxUpperTriangle(dmath::mat & m){
	maxr=0; maxc=1;
	double maxval=-1;
	for(int i = 0; i < m.r-1; i++){
		for(int j=i+1;j<m.c;j++){
			if(std::abs(m(i,j))>maxval){
				maxr=i;
				maxc=j;
				maxval=std::abs(m(i,j));
			}
		}
	}
}

void dmath::JacobiDiagonalizer::jacobiRotate(dmath::mat & m, dmath::mat & U){
	double angle = mypi/4.0;
	if(std::abs(m(maxr,maxr)-m(maxc,maxc))>0.0000001){
		angle=.5*std::atan(2.0*m(maxr,maxc)/(m(maxc,maxc)-m(maxr,maxr)));
	}

	double c=std::cos(angle),
		s=std::sin(angle);

	dmath::mat G{dmath::mat::Identity(m.r,m.c)};
	G(maxr,maxr)=c;
	G(maxc,maxc)=c;
	G(maxr,maxc)=s;
	G(maxc,maxr)=-s;
	dmath::mat test=G.transpose()*G;
	//std::cout << "GMAT:\n" << G << std::endl << test<<std::endl;

	m=G.transpose()*m*G;
	U=U*G;
}

/*
//I'm confident this can be written 1000X better
void dmath::JacobiDiagonalizer::jacobiRotate(dmath::mat & m, dmath::mat & U){
	double c,s;

	if(m(maxc,maxr)!=0.0){
		double tau, t;
		tau=(m(maxc,maxc)-m(maxr,maxr))/(2.0*m(maxc,maxr));
		if(tau>0.00000001){
			t=1.0/(tau+std::sqrt(1.0+tau*tau));
		}
		else{
			t=-1.0/(-tau+std::sqrt(1.0+tau*tau));
		}

		c=1.0/std::sqrt(1+t*t);
		s=c*t;
	}
	else{
		c=1.0;
		s=0.0;
	}

	dmath::vec evIRow(m.r), evJRow(m.r);
	for(int j = 0; j < m.r; j++){
		evIRow(j)=U(maxr,j);
		evJRow(j)=U(maxc,j);

		U(maxr,j)=evIRow(j)*c-evJRow(j)*s;
		U(maxc,j)=evIRow(j)*s+evJRow(j)*c;
	}

	double e_ii=m(maxr,maxr),
	       e_jj=m(maxc,maxc),
	       e_il,e_jk;

	m(maxr,maxr)=c*c*e_ii-2.0*s*c*m(maxr,maxc)+s*s*e_jj;
	m(maxc,maxc)=s*s*e_ii+2.0*s*c*m(maxr,maxc)+c*c*e_jj;
	m(maxr,maxc)=0.0;
	m(maxc,maxr)=0.0;

	for(int l = 0; l < m.r; l++){
		if(l!=maxr&&l!=maxc){
			e_il=m(maxr,l);
			e_jk=m(maxc,l);

			m(maxr,l)=c*e_il-s*e_jk;
			m(l,maxr)=m(maxr,l);

			m(maxc,l)=s*e_il+c*e_jk;
			m(l,maxc)=m(maxc,l);
		}
	}

}
*/

void dmath::FnCalc::load(std::string tabfile){
	std::ifstream infile(tabfile);
	if(!infile){
		std::cout << "TABULATED FILE NOT FOUND!\n";
		return;
	}
	tabulated.resize(25);
	std::string linein;
	while(!infile.eof()){
		std::getline(infile,linein);
		if(linein=="")continue;

		std::stringstream ss(linein);
		double v=0.0;
		ss >> v;
		for(int i = 0; i < tabulated.size(); i++){
			ss>>v;
			tabulated[i].push_back(v);
		}
	}
	infile.close();

	/*
	for(int i = 0; i < tabulated[0].size(); i++){
		std::cout << 0.05*i << " ";
		for(int j = 0; j < tabulated.size(); j++){
			std::cout << tabulated[j][i] << " ";
		}
		std::cout << std::endl;
	}
	*/
}

double dmath::FnCalc::Fn(int n, double T){
	if(T<12.0){
		if(n>16) return 0.0;
		int Tstari=(int)((T+0.025)/0.05);

		//std::cout << Tstari << std::endl; 

		double Tstar=(double)Tstari*0.05;
		double sum = 0.0;
		for(int k = 0; k <= 6; k++){
			sum+=tabulated[k+n][Tstari]*std::pow(Tstar-T,k)/dmath::dfactorial(k);
		}
		return sum;
	}
	else{
		return dmath::dfactorial2(2*n-1)/(2.0*std::pow(2.0*T,(double)n))*std::sqrt(mypi/T);
	}
	return 0.0;
}


