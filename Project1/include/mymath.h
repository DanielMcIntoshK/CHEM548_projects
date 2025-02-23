#ifndef MYMATH__H_
#define MYMATH__H_
#include <vector>
#include <iostream>

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
			double dot=0.0;
			for(int i = 0; i < n;i++){dot+=vals[i]*vals[i];}
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
		
		int r,c;
		std::vector<std::vector<double>> vals;
	public:
		static mat Zero(int i, int j);
		static mat Identity(int i , int j);
	};

	int factorial(int n);
	int factorial2(int n);

	std::ostream & operator<<(std::ostream & os,dmath::vec & v);
	std::ostream & operator<<(std::ostream & os,dmath::mat & m);

	dmath::vec operator*(dmath::vec v, double s);
	dmath::vec operator*(double s, dmath::vec v);
}


#endif

