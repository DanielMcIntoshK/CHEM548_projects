#include "mymath.h"


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
	for(int i = 0; i < v.vals.size(); i++) os<<v(i)<<" ";
	return os;
}

std::ostream & dmath::operator<<(std::ostream & os,dmath::mat & m){
	for(int i = 0; i < m.r; i++){
		for(int j = 0; j < m.c; j++){
			os << m(i,j) << " ";
		}
		os<<std::endl;
	}
	return os;
}

int factorial(int n){
	int f=1.0;
	for(int i = n; i >1;i--){
		f*=i;
	}
	return f;
}

int factorial2(int n){
	int f =1.0;
	for(int i = n; i >1;i-=2){
		f*=i;
	}
	return f;
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

