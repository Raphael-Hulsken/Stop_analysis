#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

#include <iostream>

using namespace std;

int main(){
	//One need to follow the sequence of methods in order to properly load the histograms

        vector<string> processes;
	cout<<" - Create an instance of SonicScrewdriver with a parameter = false in order to not create default yield plots"<<endl;
	theDoctor::SonicScrewdriver sonic(false);
	cout<<" - Load XML file "<<endl;
	sonic.LoadXMLConfig("config.xml");
	cout<<" - Create the Histos (here 1D)"<<endl;
    	sonic.Create1DHistos();
	cout<<" - Schedule plots"<<endl;
	sonic.SchedulePlots("1DSuperimposed");
	cout<<" - Import histos from the root-file"<<endl;
	sonic.ImportHistosFromFile("plotsTest/1DSuperimposed.root");		
	cout<<" - Here we could access the plots and change them ... "<<endl;
	//
	//
	// to be filled in needed
	//
	//
	cout<<" - Make plot " <<endl;
	cout<<" - Write plots" <<endl;
	sonic.MakePlots();
	sonic.WritePlots("plotsTest/new");
}
