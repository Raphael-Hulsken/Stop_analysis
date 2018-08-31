#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//./PUYields signalRegPU.tab datacards/signalReg.txt

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 3)
            throw std::runtime_error("Bad number of arguments!");
         
        TFile *f = new TFile("/opt/sbg/scratch1/cms/mjansova/store/tmp/1309/T2tt_400to1200.root");
        f->ls();
        TH1D * h1 = (TH1D*)f->Get("hnvertex");
        double total = h1->Integral();
        double iLow = h1->Integral(1,15); //@MJ@ TODO integral and error
        //cout << "bin up to: " << h1->GetBinLowEdge(20) << endl;
        double iHigh = h1->Integral(32,50);
        double totalHigh = total/iHigh;
        double totalLow = total/iLow;
        cout << iLow << "totalLow " << totalLow << " : " << iHigh <<  " totalHigh " << totalHigh << " total: " << total << endl; 
         

        string inputTab = argv[1];
        TString inputFile = argv[2];

        vector<string> regions;
        vector<string> datasets = {"500.000000_400.000000", "1000.000000_1.000000", "1000.000000_650.000000"}; //@MJ@ TODO change the datasets to meaningful ones

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

       TH1::SetDefaultSumw2();
       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

       vector<TH1F*> histo;
       vector<TH1F*> histoL;
       vector<TH1F*> histoH;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), regions.size(), 0, regions.size()));
           histoL.push_back(new TH1F((datasets.at(h)+"Low").c_str(), (datasets.at(h)+"Low").c_str(), regions.size(), 0, regions.size()));
           histoH.push_back(new TH1F((datasets.at(h)+"High").c_str(), (datasets.at(h)+"Low").c_str(), regions.size(), 0, regions.size()));
       }
       for(uint32_t r=0; r<regions.size();r++)
       {
           for(uint32_t d=0; d<datasets.size();d++)
           {
               theDoctor::Figure resLPU = tab.Get(regions.at(r)+"LowPU", datasets.at(d));
               theDoctor::Figure resHPU = tab.Get(regions.at(r)+"HighPU", datasets.at(d));
               theDoctor::Figure resall = tab.Get(regions.at(r), datasets.at(d));
               Figure res = (resLPU*totalLow)/(resHPU*totalHigh);
               Figure resH = (resHPU*totalHigh)/resall;
               Figure resL = (resLPU*totalLow)/resall;
               histo.at(d)->SetBinContent(r+1,res.value());
               histo.at(d)->SetBinError(r+1,res.error());
               histoH.at(d)->SetBinContent(r+1,resH.value());
               histoH.at(d)->SetBinError(r+1,resH.error());
               histoL.at(d)->SetBinContent(r+1,resL.value());
               histoL.at(d)->SetBinError(r+1,resL.error());
 
           }
       }

       //TCanvas *can = new TCanvas("can","can");
       //can->cd();
       TFile fi2("PUYields.root","RECREATE");
       for(uint32_t c=0; c<histo.size(); c++)
       {
           histo.at(c)->SetTitle("LowPU/HighPU");
           histoH.at(c)->SetTitle("HighPU/all");
           histoL.at(c)->SetTitle("LowPU/all");
           for(uint32_t b=0; b<regions.size(); b++)
           {
               histo.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
               histoH.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
               histoL.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
           }
           //if(c==0)
           //{
           //   histo.at(c)->Draw("hist e");
           //}
           //else
           //{               
           //histo.at(c)->SetLineColor(c+4);
           //histoL.at(c)->SetLineColor(c+4);
           //histoH.at(c)->SetLineColor(c+4);
              //histo.at(c)->Draw("same hist e");
           //}
           histo.at(c)->Write();
           histoH.at(c)->Write();
           histoL.at(c)->Write();
       }
       fi2.Close();
       //can->SaveAs("plotPUyield.eps"); //@MJ@ TODO svae to root file
}
