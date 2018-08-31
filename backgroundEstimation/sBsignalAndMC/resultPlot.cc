#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TStyle.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"


//./resultPlot card.tab signalRegMor.txt 
using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        gStyle->SetOptStat(0);
        gStyle->SetLegendBorderSize(0);
        gROOT->ForceStyle();

        if(argc != 3)
            throw std::runtime_error("Bad number of arguments!");

        string inputTab = argv[1];
        TString inputFile = argv[2];

        vector<string> regions;
        vector<string> datasets = {"totalSM", "(900,300)", "(400,300)", "(1000,100)"}; //@MJ@ TODO change the datasets to meaningful ones
        vector<string> datasetsLeg = {"SM background", "#tilde{b}#rightarrowt#tilde{#chi}^{#pm}_{1}(900,300)", "#tilde{b}#rightarrowt#tilde{#chi}^{#pm}_{1}(400,300)", "#tilde{b}#rightarrowt#tilde{#chi}^{#pm}_{1}(1000,100)"}; //@MJ@ TODO change the datasets to meaningful ones


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
        const std::string ss = "lessMETless";
        const std::string t = "-";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ss, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ss.size(), t );
                n += t.size();
            }
        }

        const std::string uu = "SR1l_";
        const std::string v ="";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( uu, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, uu.size(), v );
                n += v.size();
            }
        }

        const std::string ww = "_";
        const std::string x = ";";

        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ww, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ww.size(), x );
                n += x.size();
            }
        }

        /*for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            regLabels.at(reg) = regLabels.at(reg)+")";
        }*/


       TLegend *leg = new TLegend(0.6,0.65,0.85,0.85);

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

      vector<uint16_t> colors = {810, 603, 828};
      if(colors.size() < (histo.size()-1))
          throw std::runtime_error("Not enough colors for all histograms");

       TCanvas *can = new TCanvas("can","can");
       TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0);
       TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2);
       pad1->Draw();
       pad2->Draw();
       
        //can->SetLogy();

       pad1->cd();
       //pad1->GetFrame()->SetFillColor(0);
       pad1->SetFillStyle(0);
       //pad1->GetFrame()->SetFillStyle(0);
       pad1->SetLogy();
       for(uint32_t c=0; c<histo.size(); c++)
       {
           if(c==0)
           {
              histo.at(c)->SetTitle("");
              //histo.at(c)->SetLineColor(kGray+2);
              histo.at(c)->SetLineColor(kCyan+2);
              histo.at(c)->SetLineWidth(3);
              histo.at(c)->SetFillColor(kCyan+2);
              histo.at(c)->SetFillStyle(1001);
              float max = histo.at(c)->GetMaximum();
              histo.at(c)->SetMaximum(1.5*max);
              //histo.at(c)->GetXaxis()->LabelsOption("v");
              //histo.at(c)->GetXaxis()->SetTitleOffset(1.2);
              histo.at(c)->Draw("hist e"); //SET titles (per bin?!)
              histo.at(c)->GetXaxis()->SetTitle("");
              histo.at(c)->GetYaxis()->SetTitle("Events");
           }
           else
           {               
              //histo.at(c)->SetLineColor(c+2);
              histo.at(c)->SetLineColor(colors.at(c-1));
              histo.at(c)->SetLineWidth(3);
              histo.at(c)->SetLineStyle(2);
              histo.at(c)->Draw("same hist e");
           }
              leg->AddEntry(histo.at(c), (TString)datasetsLeg.at(c));
       }
       leg->Draw("l");

       pad2->cd();
              TH1F* ratio = (TH1F*) histo.at(0)->Clone();
              for(uint32_t b=0; b< ratio->GetNbinsX(); b++)
              {
                  ratio->SetBinContent(b+1, 1);
              }
              ratio->SetFillColor(0);
              ratio->SetLineColor(kBlack);
              ratio->SetMinimum(0);
              ratio->SetMaximum(2);
              ratio->GetXaxis()->SetTitle("MET");
              ratio->GetYaxis()->SetTitle("data/MC");
              ratio->Draw("hist");



       /*cout << "chi2 test 1: " << histo.at(0)->Chi2Test(histo.at(1), "WW P") << endl;

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

       cout << "chi2 test 2: " << h1->Chi2Test(h2, "WW P") << endl;*/
       can->Update();
       can->SaveAs("plotTest.root"); //@MJ@ TODO svae to root file
}
