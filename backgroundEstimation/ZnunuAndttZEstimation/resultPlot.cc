#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"


//./resultPlot yieldZnunuMorTTZNLOSR.tab realregions.txt 
using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 3)
            throw std::runtime_error("Bad number of arguments!");

        string inputTab = argv[1];
        TString inputFile = argv[2];

        vector<string> regions;
        vector<string> datasets = {"ttZ", "ttZNLO"}; //@MJ@ TODO change the datasets to meaningful ones


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

        vector<string> regLabels = regions;
        const std::string ss = "less";
        const std::string t = "<";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ss, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ss.size(), t );
                n += t.size();
            }
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
               theDoctor::Figure res = tab.Get(regions.at(r),datasets.at(d));
               histo.at(d)->SetBinContent(r+1,res.value());
               histo.at(d)->SetBinError(r+1,res.error());
               histo.at(d)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
 
           }
       }

       TCanvas *can = new TCanvas("can","can");
       can->cd();
       for(uint32_t c=0; c<histo.size(); c++)
       {
           if(c==0)
           {
              histo.at(c)->SetLineColor(kBlack);
              histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
              histo.at(c)->GetXaxis()->SetTitle("signal regions");
              histo.at(c)->GetYaxis()->SetTitle("yield");
           }
           else
           {               
              //histo.at(c)->SetLineColor(c+2);
              histo.at(c)->SetLineColor(kRed);
              histo.at(c)->Draw("same hist e");
           }
       }

       cout << "chi2 test 1: " << histo.at(0)->Chi2Test(histo.at(1), "WW P") << endl;

       TH1D* h1 = new TH1D("ttZ1", "ttZ1", regions.size()-6, 0, regions.size()-6);
       TH1D* h2 = new TH1D("ttZ2", "ttZ2", regions.size()-6, 0, regions.size()-6);

       uint32_t bin = 1;
       for(uint32_t b = 0; b<regions.size(); b++)
       {
           if(b+1==6 || b+1==7 || b+1==12 || b+1==15 || b+1==16 || b+1==27)
               continue;
           double value1 = histo.at(0)->GetBinContent(b+1); 
           double error1 = histo.at(0)->GetBinError(b+1); 
           h1->SetBinContent(bin,value1);
           h1->SetBinError(bin,error1);
           double value2 = histo.at(1)->GetBinContent(b+1); 
           double error2 = histo.at(1)->GetBinError(b+1); 
           h2->SetBinContent(bin,value2);
           h2->SetBinError(bin,error2);
           bin++;
       }

       cout << "chi2 test 2: " << h1->Chi2Test(h2, "WW P") << endl;
       can->SaveAs("plotTest.root"); //@MJ@ TODO svae to root file
}
