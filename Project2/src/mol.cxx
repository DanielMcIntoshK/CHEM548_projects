#include "mol.h"
#include <fstream>
#include <sstream>
#include <cmath>

void molecule::read(std::string filename){
	std::ifstream in(filename);
	if(!in){
		std::cout << "INPUT FILE DOES NOT EXIST\n";
		return;
	}

	std::string inputline;
	std::string basisfile;
	std::string basisdir;

	while(!in.eof()){
		std::getline(in,inputline);
		std::cout << "LINE: " << inputline << std::endl;
		if(in.eof())break;
		if(inputline.find('!')!=std::string::npos)continue;

		if(inputline=="GEOMETRY"){
			std::string a;
			double x,y,z;

			std::string mcountline;
			std::getline(in,mcountline);
			std::stringstream countunit(mcountline);
			int mcount =0;
			std::string unit;
			countunit>>mcount>>unit;

			double unitconvert=1.0;
			if(unit=="ANGSTROM") unitconvert=1.8897259886;

			for(int i = 0; i < mcount; i++){
				std::string gline;
				std::getline(in,gline);
				std::stringstream gs(gline);

				gs >> a >> x >> y >> z;
				
				Atom at;
				at.an=Atom::elements[a];
				at.center=dmath::vec(std::vector<double>{x,y,z});
				at.center=at.center*unitconvert;

				atoms.push_back(at);
			}
			std::cout << "GEO:\n";
			for(int i = 0; i < atoms.size(); i++){
				std::cout << atoms[i].an << " " << atoms[i].center << std::endl;
			}
		}
		if(inputline=="BASISSET"){
			std::string bsline;
			std::getline(in,basisfile);
		}
		if(inputline=="BASISDIR"){
			std::string bsline;
			std::getline(in,basisdir);
		}
		if(inputline=="CHARGE"){
			std::string chline;
			std::getline(in,chline);
			charge=std::stoi(chline);
		}
	}

	nelec=0;
	for(int i = 0; i < atoms.size(); i++){
		nelec+=atoms[i].an;
	}
	nelec-=charge;

	std::string fullbsname=basisdir+basisfile;
	buildbasis(fullbsname);

	/*
	std::cout << "BUILD BASIS\n";
	for(int i = 0; i < orbitals.size(); i++){
		std::cout << orbitals[i].at.an << " " << orbitals[i].at.center<< std::endl;
		for(int j = 0; j < orbitals[i].contraction.size(); j++){
			std::cout << orbitals[i].contraction[j].exp << " " << orbitals[i].contraction[j].co <<std::endl;
		}
	}
	*/
}

void molecule::buildbasis(std::string filename){
	std::ifstream in(filename);
	std::cout << "BUILDING BASIS\n";

	while(!in.eof()){
		std::string line;
		std::getline(in,line);
		
		if(line =="") continue;

		std::stringstream sline(line);

		std::string atname;
		int primscounttotal, orbscount;

		sline>>atname>>primscounttotal>>orbscount;

		bool reqel=false;
		int atnum=Atom::elements[atname];
		for(int i = 0; i < atoms.size();i++){
			if(atoms[i].an==atnum){
				reqel=true;
				break;
			}
		}

		if(reqel){
			std::string orbinfo;
			for(int i = 0; i < orbscount;i++){
				std::vector<atOrb> orbs;
				std::vector<primitive> prims;
				
				std::getline(in,orbinfo);
				std::stringstream sorbinfo(orbinfo);

				std::string otype;
				int primcount;

				sorbinfo>>otype>>primcount;

				double factor=1.0;
				std::string priminfo;
				for(int i = 0; i < primcount; i++){
					std::getline(in,priminfo);
					std::stringstream spinfo(priminfo);
					double c,a;
					spinfo>>a>>c;
										
					primitive prim;
					prim.exp=a;
					prim.co=c;
					prim.ang[0]=0;
					prim.ang[1]=0;
					prim.ang[2]=0;
					
					prims.push_back(prim);
				}
				switch(otype[0]){
				case 'S':{
					atOrb ob;
					ob.contraction=prims;
					orbs.push_back(ob);
				}break;
				case 'P':{
					for(int i = 0; i < 3; i++){
						atOrb ob;
						ob.contraction=prims;
						for(auto & a: ob.contraction){
							a.ang[i]+=1;
						}
						orbs.push_back(ob);
					}
				}break;
				case 'D':{
					for(int i = 0; i < 3; i++){
					for(int j = i; j < 3; j++){
						atOrb ob;
						ob.contraction=prims;
						for(auto & a: ob.contraction){
							a.ang[i]+=1;
							a.ang[j]+=1;
						}
						orbs.push_back(ob);
					}}
				};break;
				case 'F':{
					for(int i = 0; i < 3; i++){
					for(int j = i; j < 3; j++){
					for(int k = j; k < 3; k++){
						atOrb ob;
						ob.contraction=prims;
						for(auto & a: ob.contraction){
							a.ang[i]+=1;
							a.ang[j]+=1;
							a.ang[k]+=1;
						}
						orbs.push_back(ob);
					}}}
				};break;
				}
				for(auto &a: orbs){
					for(auto &o: a.contraction){
						o.norm=std::pow(2.0*o.exp/mypi,3.0/4.0)*std::pow(4.0*o.exp,(double)o.l()/2.0)*
						       std::pow(dmath::dfactorial2(2*o.ang[0]-1)*
								dmath::dfactorial2(2*o.ang[1]-1)*
								dmath::dfactorial2(2*o.ang[2]-1),-1.0/2.0);
						if(o.l()==2&&(o.ang[0]==0 ||o.ang[1]==0||o.ang[2]==0)){
							o.norm/=std::sqrt(3);
						}
					}
				}

				for(int i = 0; i < atoms.size(); i++){
					if(atnum==atoms[i].an){
						for(auto & a: orbs) {
							a.at=atoms[i];
							for(auto & p: a.contraction){
								p.at=atoms[i];
							}
						}
						orbitals.insert(orbitals.end(),orbs.begin(),orbs.end());
					}
				}
			}
		}
		else{
			for(int i = 0; i < primscounttotal+orbscount;i++){
				std::getline(in,line);
			}
		}
	}
}

void molecule::swapbs(std::vector<int> swapvec){
	if(swapvec.size()!=orbitals.size()) return;

	std::vector<atOrb> neworbs;
	for(int i = 0; i < swapvec.size(); i++){
		neworbs.push_back(orbitals[swapvec[i]]);
	}
	orbitals=neworbs;
}

std::map<std::string,int> Atom::elements=std::map<std::string,int>{{"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},
	{"N",7},{"O",8},{"F",9},{"Ne",10}};
