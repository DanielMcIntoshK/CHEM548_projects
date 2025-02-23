#include <iostream>
#include "mymath.h"

int main(int argc, char ** argv){

	dmath::vec v3(3);
	std::cout << v3<<std::endl;

	v3(0)=4.0;
	v3(1)=2.0;
	v3(2)=98.0;

	std::cout << v3<<std::endl;

	dmath::mat m(4,4);

	std::cout << m << std::endl;

	m=dmath::mat::Identity(4,4);
	m(2,1)=7.0;
	std::cout << m << std::endl;

	return 0;
}
