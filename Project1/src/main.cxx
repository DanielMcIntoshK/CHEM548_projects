#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "mymath.h"
#include "potentials.h"
#include "fd.h"

int main(int argc, char ** argv){
	FD finitediff;

	int steps=100;
	dmath::vec pibV(steps), finV(steps),pbrecV(steps),harmV(steps),morse1(steps),morse2(steps);

	double x0=-5.0,xn=5.0;
	double dx=(xn-x0)/steps;
	for(int i = 0; i < steps; i++) {
		double x=x0+(double)i*dx;
		pibV(i)=finitewell(x,-4,4,100000.0);
		//pibV(i)=particleinabox(x,-3,3);
	}
	for(int i = 0; i < steps; i++) {
		double x=x0+(double)i*dx;
		finV(i)=finitewell(x,-4,4,2.5);
	}
	for(int i = 0; i < steps; i++) {
		double x=x0+(double)i*dx;
		double width=.2;
		double innerheight=5.0;
		pbrecV(i)=pbrectangular(x,-4,4,-width/2,width/2,innerheight);
	}
	for(int i = 0; i < steps; i++){
		double x=x0+(double)i*dx;
		harmV(i)=harmonic(x,1.0,0.0);
	}
	x0=-1.0,xn=5.0;
	steps=100;
	dx=(xn-x0)/steps;
	for(int i = 0; i < steps; i++){
		double x=x0+(double)i*dx;
		morse1(i)=morse(x,2,5,0,0);
		morse2(i)=morse(x,.5,3,7,1.5);
	}

	finitediff.compute(1.0,10.0/100.0,pibV);
	std::ofstream outfile1("Output/Pib");
	for(int i = 0; i < steps; i++){
		outfile1<<pibV(i)<<" ";
	}
	outfile1 << std::endl;
	for(int i = 0; i < steps; i++){
		outfile1<<finitediff.eigenvals(i,i)<<" ";
	}
	outfile1<<std::endl;
	outfile1<<finitediff.eigenvecs<<std::endl;
	
	finitediff.compute(1.0,10.0/100.0,finV);
	std::ofstream outfile2("Output/Fin");
	for(int i = 0; i < steps; i++){
		outfile2<<finV(i)<<" ";
	}
	outfile2 << std::endl;
	for(int i = 0; i < steps; i++){
		outfile2<<finitediff.eigenvals(i,i)<<" ";
	}
	outfile2<<std::endl;
	outfile2<<finitediff.eigenvecs<<std::endl;
	
	finitediff.compute(1.0,10.0/100.0,pbrecV);
	
	std::ofstream outfile3("Output/Pbrec");
	for(int i = 0; i < steps; i++){
		outfile3<<pbrecV(i)<<" ";
	}
	outfile3 << std::endl;
	for(int i = 0; i < steps; i++){
		outfile3<<finitediff.eigenvals(i,i)<<" ";
	}
	outfile3<<std::endl;
	outfile3<<finitediff.eigenvecs<<std::endl;
	
	finitediff.compute(1.0,10.0/100.0,harmV);
	
	std::ofstream outfile4("Output/Harm");
	for(int i = 0; i < steps; i++){
		outfile4<<harmV(i)<<" ";
	}
	outfile4 << std::endl;
	for(int i = 0; i < steps; i++){
		outfile4<<finitediff.eigenvals(i,i)<<" ";
	}
	outfile4<<std::endl;
	outfile4<<finitediff.eigenvecs<<std::endl;
	
	finitediff.compute(1.0,10.0/100.0,morse1);
	std::ofstream outfile5("Output/Morse1");
	for(int i = 0; i < steps; i++){
		outfile5<<morse1 <<" ";
	}
	outfile5<<std::endl;
	for(int i = 0; i < steps; i++){
		outfile5<<finitediff.eigenvals(i,i)<<" ";
	}
	outfile5<<std::endl;
	outfile5<<finitediff.eigenvecs<<std::endl;

	finitediff.compute(1.0,10.0/100.0,morse2);
	std::ofstream outfile6("Output/Morse2");
	for(int i = 0; i < steps; i++){
		outfile6<<morse2 <<" ";
	}
	outfile6<<std::endl;
	for(int i = 0; i < steps; i++){
		outfile6<<finitediff.eigenvals(i,i)<<" ";
	}
	outfile6<<std::endl;
	outfile6<<finitediff.eigenvecs<<std::endl;

	return 0;
}
