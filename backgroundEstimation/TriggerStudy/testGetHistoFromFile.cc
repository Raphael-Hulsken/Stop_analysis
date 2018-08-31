//to get histogram(s) from file and canvas
#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

#include <iostream>

using namespace std;

int main(){

        vector<string> processes;
	theDoctor::SonicScrewdriver sonic(false);
        sonic.LoadXMLConfig("config.xml");
        //sonic.GetProcessClassTagList(&processes);
        processes.push_back("rare");
        vector<TH1D*> h;
        h = sonic.Get1DHistoCloneFromFile("plotsTest", "1DStack","MET",processes,"SR1l2jMET250to350","lepChannel");
        cout<<"size of histos: " <<h.size()<<endl;
       	std::cout<<h.at(0)->GetEntries()<<std::endl;
}
