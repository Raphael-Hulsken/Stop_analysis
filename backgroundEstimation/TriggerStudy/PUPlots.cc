//to get histogram(s) from file and canvas
#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

#include <iostream>

using namespace std;

int main(){

        TH1::SetDefaultSumw2();

        vector< vector<string> > regions;
        //
        //TO READ ALL SIGNAL REGIONS
        //
/*	string line;
        ifstream regfile ("datacards/signalReg.txt");
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                vector<string> names;
                names.push_back(line);
                names.push_back(line+"LowPU");
                names.push_back(line+"HighPU");
                regions.push_back(names);
            }
            regfile.close();
        }
*/
        vector<string> processes;

        vector<string> regionsBl = {"Baseline", "BaselineLowPU","BaselineHighPU",};
        regions.push_back(regionsBl);
        vector<string> variable;
	theDoctor::SonicScrewdriver sonic(false);
        sonic.LoadXMLConfig("config.xml");
        sonic.GetVariablesTagList(&variable);

        uint16_t num = 1;
        uint16_t den = 2;

        //signals
        processes.push_back("500.000000_400.000000");
        processes.push_back("1000.000000_1.000000");
        processes.push_back("1000.000000_650.000000");

        vector< vector<TH1D*> > yield;
        for(uint32_t v = 0; v<variable.size(); v++)
        {
        
        for(uint32_t m = 0; m<regions.size(); m++)
        {
        for(uint32_t r = 0; r < regions.at(m).size(); r++  )
        {
            yield.push_back(sonic.Get1DHistoCloneFromFile("plotsTest", "1DSuperimposedNoNorm",variable.at(v),processes, regions.at(m).at(r),"lepChannel"));
        }
        if(v==0)
            TFile fi2("PUplots.root","RECREATE");
        TFile fi2("PUplots.root","UPDATE");
        for(uint32_t p = 0; p < processes.size(); p++  )
        {
            //@MJ@ TODO could create new histogram copying the old one and than divide by other method
            yield.at(num).at(p)->Divide(yield.at(den).at(p));
            TString name = variable.at(v) + ":point:"+processes.at(p)+":"+regions.at(m).at(num)+"/"+regions.at(m).at(den);
            yield.at(num).at(p)->SetName(name);
            yield.at(num).at(p)->Write();
        }
        fi2.Close();
        yield.clear();
        }
        }

}
