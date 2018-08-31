#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//./PUYieldsMorePoints signalRegOneBinPU.tab

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 2)
            throw std::runtime_error("Bad number of arguments!");

        string region = "allRegions";
         
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

        vector<string> regions;
        vector<string> datasets = {
"500.000000_1.000000",
"500.000000_50.000000",
"500.000000_100.000000",
"500.000000_150.000000",
"500.000000_200.000000",
"500.000000_250.000000",
"500.000000_275.000000",
"500.000000_300.000000",
"500.000000_325.000000",
"500.000000_350.000000",
"500.000000_375.000000",
"500.000000_400.000000",
"1000.000000_1.000000",
"1000.000000_50.000000",
"1000.000000_100.000000",
"1000.000000_150.000000",
"1000.000000_200.000000",
"1000.000000_250.000000",
"1000.000000_300.000000",
"1000.000000_350.000000",
"1000.000000_400.000000",
"1000.000000_450.000000",
"1000.000000_500.000000",
"1000.000000_550.000000",
"1000.000000_600.000000",
"1000.000000_650.000000"};



       TH1::SetDefaultSumw2();
       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

           TH1F* histo = new TH1F("Low/High", "Low/High", datasets.size(), 0, datasets.size());
           TH1F* histoL = new TH1F("Low/all", "Low/all", datasets.size(), 0, datasets.size());
           TH1F* histoH = new TH1F("High/all", "High/all", datasets.size(), 0, datasets.size());
           for(uint32_t d=0; d<datasets.size();d++)
           {
               theDoctor::Figure resLPU = tab.Get(region+"LowPU", datasets.at(d));
               theDoctor::Figure resHPU = tab.Get(region+"HighPU", datasets.at(d));
               theDoctor::Figure resall = tab.Get(region, datasets.at(d));
               Figure res = (resLPU*totalLow)/(resHPU*totalHigh);
               Figure resH = (resHPU*totalHigh)/resall;
               Figure resL = (resLPU*totalLow)/resall;
               histo->SetBinContent(d+1,res.value());
               histo->SetBinError(d+1,res.error());
               histoH->SetBinContent(d+1,resH.value());
               histoH->SetBinError(d+1,resH.error());
               histoL->SetBinContent(d+1,resL.value());
               histoL->SetBinError(d+1,resL.error());
 
           }

       //TCanvas *can = new TCanvas("can","can");
       //can->cd();
       TFile fi2("PUYieldsMorePoints2.root","RECREATE");
           histo->SetTitle("LowPU/HighPU");
           histoH->SetTitle("HighPU/all");
           histoL->SetTitle("LowPU/all");
           for(uint32_t b=0; b<datasets.size(); b++)
           {
               histo->GetXaxis()->SetBinLabel(b+1,datasets.at(b).c_str());
               histoH->GetXaxis()->SetBinLabel(b+1,datasets.at(b).c_str());
               histoL->GetXaxis()->SetBinLabel(b+1,datasets.at(b).c_str());
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
           histo->Write();
           histoH->Write();
           histoL->Write();
       fi2.Close();
       //can->SaveAs("plotPUyield.eps"); //@MJ@ TODO svae to root file
}
