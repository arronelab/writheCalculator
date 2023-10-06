#include "writhe.h"
#include <string.h>
#include <string>
#include <iostream>



int main( int argc, const char* argv[] )
{
       if(argc >=2){
	bool check=false;
	std::ifstream myfile;
 	myfile.open(argv[1]);
	std::ofstream outfile;
 	outfile.open(argv[2]);
	std::string output;
	std::vector<point> points;
	if (myfile.is_open()) { 
	  while(std::getline(myfile,output)){
		points.push_back(point(output));
  		}
	}else{
	std::cout<<"Curve data file failed to open";
	}
	if(points.size()>3){
	    check =true;
	}
	myfile.close();
	if(check){
	  writhe w;
	  std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > > fps = w.getWritheFingerprints(points);
	  for(int i=0;i<fps.first.size();i++){
	    for(int j=0;j<fps.first[i].size();j++){
	      outfile<<fps.first[i][j].getX()+1<<" "<<fps.first[i][j].getY()+1<<" "<<fps.first[i][j].getZ()<<" "<<fps.second[i][j].getZ()<<"\n";
	    }
	  }
	  outfile.close(); 
	}else{
	  std::cout<<"need more than 3 points in the curve\n";
	}
       }else{
	  std::cout<<"must supply a file containing the curve data";
       }
       return 0;
}
