#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"



using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 3)
            throw std::runtime_error("Bad number of arguments!");

        string inputTab = argv[1];
        TString inputFile = argv[2];

        vector<string> regions;
        vector<string> datasets = {"totalSM", "signal", "data"}; //@MJ@ TODO change the datasets to meaningful ones


	string line;
        ifstream regfile(inputFile);
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                regions.push_back(line);
            }
            regfile.close();
        }

       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

       vector<TH1F*> histo;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), regions.size(), 0, regions.size()));
       }
       for(uint32_t r=0; r<regions.size();r++)
       {
           for(uint32_t d=0; d<datasets.size();d++)
           {
               theDoctor::Figure res = tab.Get(datasets.at(d), regions.at(r));
               histo.at(d)->SetBinContent(r+1,res.value());
               histo.at(d)->SetBinError(r+1,res.error());
 
           }
       }

       TCanvas *can = new TCanvas("can","can");
       can->cd();
       for(uint32_t c=0; c<histo.size(); c++)
       {
           if(c==0)
              histo.at(c)->Draw("hist e");
           else
           {               
              histo.at(c)->SetLineColor(c+2);
              histo.at(c)->Draw("same hist e");
           }
       }
       can->SaveAs("plotTest.eps"); //@MJ@ TODO svae to root file
}
