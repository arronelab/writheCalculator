#include "writhe.h"
#include <string.h>
#include <string>
#include <iostream>
#include <filesystem>


int main( int argc, const char* argv[] )
{
       if(argc >=2){
	bool check=false;
	//Open and read the main curve file
	std::ifstream myfile;
 	myfile.open(argv[1]);
 	std::string output;
	std::vector<point> points1;
	if (myfile.is_open()) { 
	  while(std::getline(myfile,output)){
	    points1.push_back(point(output));
	  }
	}else{
	std::cout<<"Curve data file 1 failed to open";
	}
	writhe w;
	int len1 = points1.size();
	 myfile.close();
	//calculate its finger print
	
	double cutOff = std::atof(argv[3]);
	std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > >  fp1=w.getWritheFingerprints(points1);

	// open and read the folder with all the other curves
	
	std::string directory(argv[2]);

	// open the output file to write to

	std::string outFileLoc = "comparisons/"+std::string(argv[4])+"_"+std::string(argv[2])+"_"+std::string(argv[3])+".dat";
	std::ofstream outFile;
	outFile.open(outFileLoc);
	int index=0;
	for(const auto & entry : std::filesystem::directory_iterator(directory)){
	  std::ifstream myfile2;
	  std::string fl = entry.path();
	  index++;
	  using std::filesystem::directory_iterator;
	  using fp = bool (*)( const std::filesystem::path&);
	  int size_t = std::count_if(directory_iterator(directory), directory_iterator{}, (fp)std::filesystem::is_regular_file);
	  ///int dir_size = std::count_if(std::filesystem::directory_iterator(directory), {}, std::filesystem::is_regular_file);
	  std::cout << "\r" <<" Compared "<< index <<" of "<<size_t;
	  myfile2.open(fl);
	  std::vector<point> points2;
	  if (myfile2.is_open()) { 
	    while(std::getline(myfile2,output)){
	      points2.push_back(point(output));
	    }
	  }else{
	    std::cout<<"Curve data file 2 failed to open";
	  }
	  if((points1.size()>3) && (points2.size()>3)){
	    check =true;
	  }
	  myfile2.close();
	  if(check){
	    int len2 = points2.size();
	    std::pair<std::vector<std::vector<point> >,std::vector<std::vector<point> > >  fp2=w.getWritheFingerprints(points2);
	    std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > fSetBest;
	    int len1Best=0;
	    for(int i=0;i<10;i++){
	      double cutofftest = cutOff - i*(cutOff-0.01)/9;
	      //std::cout<<cutofftest<<"\n";
	      std::vector<std::pair<std::pair<int,int>,std::pair<int,int> > > fSet =w.compareFingerPrints(fp1.first,len1,fp2.first,len2,cutofftest);
	      int lenTot1=0;int lenTot2=0;
	      for(int i=0;i<fSet.size();i++){
		if(lenTot1+fSet[i].first.second-fSet[i].first.first>0){
		  lenTot1 = lenTot1+fSet[i].first.second-fSet[i].first.first+1;
		}
		if(lenTot1+fSet[i].second.second-fSet[i].second.first>0){
		  lenTot2 = lenTot2+fSet[i].second.second-fSet[i].second.first+1;
		}
	      }
	      //std::cout<<lenTot1<<"\n";
	      if(lenTot1>len1Best || i==0){
		fSetBest = fSet;
		len1Best = lenTot1;
	      }
	    }
	    int lenTot1=0;int lenTot2=0;
	    for(int i=0;i<fSetBest.size();i++){
	      outFile<<fSetBest[i].first.first<<" "<<fSetBest[i].first.second<<" "<<fSetBest[i].second.first<<" "<<fSetBest[i].second.second<<" ";
	      if(lenTot1+fSetBest[i].first.second-fSetBest[i].first.first>0){
		lenTot1 = lenTot1+fSetBest[i].first.second-fSetBest[i].first.first+1;
	      }
	    if(lenTot1+fSetBest[i].second.second-fSetBest[i].second.first>0){
	      lenTot2 = lenTot2+fSetBest[i].second.second-fSetBest[i].second.first+1;
	    }
	    }
	    outFile<<double(lenTot1)/double(len1)<<" ";
	    outFile<<double(lenTot2)/double(len2)<<" ";
	    outFile<<fl<<"\n";
	  
	  }else{
	    std::cout<<"need more than 3 points in the curve\n";
	  }
	  
	}
	outFile.close();
       }else{
	 std::cout<<"must supply a file containing the curve data";
       }
       return 0;
}
