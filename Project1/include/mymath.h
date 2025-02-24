#ifndef MYMATH__H_
#define MYMATH__H_
#include <vector>
#include <iostream>
#define mypi 3.1415926

namespace dmath{

	struct vec{
		vec(){}
		vec(std::vector<double> init):vals{init}{n=vals.size();}
		vec(int s){vals.resize(s);n=s;}

		double & operator()(int i){return vals[i];}
		vec operator+(vec v2){
			vec vn=v2;
			for(int i =0; i < n;i++){vn.vals[i]+=vals[i];}
			return vn;
		}
		vec operator-(vec v2){
			vec vn(n);
			for(int i = 0; i < n; i++){vn(i)=vals[i]-v2.vals[i];}
			return vn;
		}
		double operator*(vec v2){
			if(v2.vals.size()!=vals.size())return 0.0;
			double dot=0.0;
			for(int i = 0; i < n;i++){dot+=vals[i]*v2.vals[i];}
			return dot;
		}
		vec operator-(){
			vec vn(n);
			for(int i =0; i < n;i++){vn(n)=-vals[n];}
			return vn;
		}
		int n;
		std::vector<double> vals;
	};

	class mat{
	public:
		mat(){}
		mat(int i, int j){cSize(i,j);}

		void cSize(int i, int j);

		double & operator()(int i, int j){return vals[i][j];}
		mat operator+(mat m){
			mat mn(r,c);
			for(int i = 0; i < r; i++){
			for(int j = 0; j < c; j++){
				mn(i,j)=vals[i][j]+m(i,j);
			}}
			return mn;
		}
		mat operator-(mat m){
			mat mn(r,c);
			for(int i = 0; i < r; i++){
			for(int j = 0; j < c; j++){
				mn(i,j)=vals[i][j]-m(i,j);
			}}
			return mn;
		}
		mat operator-(){
			mat mn(r,c);
			for(int i = 0; i < r; i++){
			for(int j = 0; j < c; j++){
				mn(i,j)=-vals[i][j];
			}}
			return mn;
		}
		mat operator*(mat m);

		mat transpose(){
			mat mn=*this;
			for(int i = 0; i < r-1; i++){
			for(int j = i+1;j<c;j++){
				std::swap(mn(i,j),mn(j,i));
			}}
			return mn;
		}

		mat diagonalize(mat & U, double threshold);
		
		int r,c;
		std::vector<std::vector<double>> vals;
	public:
		static mat Zero(int i, int j);
		static mat Identity(int i , int j);
	};

	int factorial(int n);
	int factorial2(int n);

	double dfactorial(int n);
	double dfactorial2(int n);

	class JacobiDiagonalizer{
	public:
		JacobiDiagonalizer(){}
		
		dmath::mat diagonalize(dmath::mat m, double threshold,dmath::mat & U);

		void findMaxUpperTriangle(dmath::mat & m);
		void jacobiRotate(dmath::mat & m, dmath::mat & U);

		void sort(dmath::mat & vals, dmath::mat & vecs);
	private:
		int maxr, maxc;
	};


	std::ostream & operator<<(std::ostream & os,dmath::vec & v);
	std::ostream & operator<<(std::ostream & os,dmath::mat & m);

	dmath::vec operator*(dmath::vec v, double s);
	dmath::vec operator*(double s, dmath::vec v);

	dmath::mat operator*(dmath::mat m, double s);
	dmath::mat operator*(double s, dmath::mat m);
}


#endif

