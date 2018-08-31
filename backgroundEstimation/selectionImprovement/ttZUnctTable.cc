#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//./ttZUnctTable yield.tab signalReg.txt statNames.txt nUnc+1 yieldJECDown.tab yieldJECUp.tab(=nparameters per one signal region, nUnc+1 = 16 currently)
//./ttZUnctTable yieldMor.tab signalRegMor.txt statNames.txt nUnc+1 yieldJECDownMor.tab yieldJECUpMor.tab(=nparameters per one signal region, nUnc+1 = 15 currently)

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 7)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTab = argv[1];
        string inputFile = argv[2];
        string uncNames = argv[3];
        TString nUncs = argv[4];
        string JECDownTab = argv[5];
        string JECUpTab = argv[6];
        int nUnc = nUncs.Atoi();
        if(nUnc%2 !=1)
            cout << "wrong number of systematics!! tables will not be filled correctly" << endl;

        vector<string> regions;
        //specify dataset to be read
        vector<string> datasets = {"Znunu"}; //@MJ@ TODO change the datasets to meaningful ones

        //read signal regions 
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
      
       if(((regions.size())%nUnc) != 0)
           throw std::runtime_error("wrong number of regions");

       TH1::SetDefaultSumw2();
       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

       Table tJECdown(JECDownTab);
       Table tJECup(JECUpTab);
       

        //read names of sysematics 
        vector<string> colI;
	string line2;
        ifstream systfile(uncNames);
        ofstream myfile2;
        myfile2.open ("uncertainties.txt");
        if (systfile.is_open())
        {
            colI.push_back("yield");
            while ( getline (systfile,line2) )
            { 
                if(colI.size() == nUnc)
                    break;
                colI.push_back(line2);
                myfile2 << line2 << endl;
            }
            colI.push_back("JECdown");
            myfile2 << "JECdown" << endl;
            colI.push_back("JECup");
            myfile2 << "JECup" << endl;
            systfile.close();
        }
        //get nemaes of signal regions without regions for systematics
        uint32_t rowIdI = 0;
        vector<Double_t> error;
        Double_t result;
        vector<string> realReg;
        ofstream myfile;
        myfile.open ("realregions.txt");
        for(uint32_t e=0; e< regions.size(); e++)
        {
            if(e==0 || (e%nUnc)==0) //nUnct for 8 syst should be 9!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            {
                realReg.push_back(regions.at(e));
                myfile<< regions.at(e) << endl;
            }
        }
        myfile.close();
       //define output histogram
       vector<TH1F*> histo;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), realReg.size(), 0, realReg.size()));
       }

       //output syst names
       vector<string> systOutNames;
       for(uint32_t s=0; s<colI.size(); s++)
       {
           if(s%2 == 1)
               systOutNames.push_back(colI.at(s));
       }
       //tables for higher uncertainty and summary table
       Table thigher(systOutNames, realReg);
       vector<string> uncRange = {"syst low (%)", "syst high (%)"};
       Table tsummary( uncRange, systOutNames);

       
       //tables for all realtive and absolute uncertainties
       Table tI(colI, realReg);
       Table trel(colI, realReg);
       uint32_t uncLine = 0;
       float uncTot = 0;
       float uncBef = 0;      
       float uncHigher = 0;      
       float binValue = 0;
       uint32_t l = 0;
       theDoctor::Figure yieldVal;

       //loop over all regions
       for(uint32_t r=0; r<regions.size();r++)
       {
           theDoctor::Figure resall = tab.Get(regions.at(r), datasets.at(0));
           
           //first there is real region value
           if(r == 0 || r%nUnc == 0)
           {
               l=0;
               uncLine++;
               yieldVal = resall;
               histo.at(0)->SetBinContent(uncLine,resall.value());
               if(r != 0)
                   histo.at(0)->SetBinError(uncLine-1,uncTot);
               uncTot = 0;
               uncBef = 0;
 
           }
           //then all systematics
           else
           {
               if(uncBef == 0)
               {
                   uncBef = resall.value();
               }
               else
               { 
                   binValue = histo.at(0)->GetBinContent(uncLine);
                   uncHigher = abs(uncBef-binValue) > abs(resall.value() - binValue) ? abs(uncBef-binValue): abs(resall.value() - binValue);   
                   uncTot += uncHigher*uncHigher ;//@MJ@TODO now average unc, change to maximal?!
                   uncBef = 0;

                   //fill higher
                   thigher.Set(systOutNames.at(l), realReg.at(uncLine-1), Figure((uncHigher/binValue)*100, 0));
                   l++;


               }
           }
           //fill the tables
           //cout << "colI " << colI.at(r-((uncLine-1)*(nUnc))) << endl;
	   tI.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), resall);
           Figure g = (resall-yieldVal)/yieldVal;
           if(r == 0 || r%nUnc == 0)
           {
               trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(abs(yieldVal.error()/yieldVal.value())*100,abs(0)));
           }
           else
           {
               trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(abs(g.value())*100,abs(0)));
           }
           //cout << realReg.at(uncLine-1) << " value " <<  colI.at(r-((uncLine-1)*(nUnc)))<< endl;
       }
       //fill JEC
       for(uint32_t j = 0; j<realReg.size(); j++)
       {     
           theDoctor::Figure JECd = tJECdown.Get(realReg.at(j),"totalSM" );
           theDoctor::Figure JECu = tJECup.Get(realReg.at(j), "totalSM");
           cout << "JEC down value " << JECd.value() << endl;
 
           tI.Set("JECdown", realReg.at(j), JECd);
           tI.Set("JECup", realReg.at(j), JECu);
           Figure a = (JECd - tab.Get(realReg.at(j), datasets.at(0)))/(tab.Get(realReg.at(j), datasets.at(0)));
           trel.Set("JECdown", realReg.at(j), Figure(abs(a.value())*100,abs(0)));
           Figure c = (JECu - tab.Get(realReg.at(j), datasets.at(0)))/(tab.Get(realReg.at(j), datasets.at(0)));
           trel.Set("JECup", realReg.at(j), Figure(abs(c.value())*100,abs(0)));
           float higherJEC = a.value() > c.value()? a.value(): c.value();
           uint32_t lastElement = systOutNames.size() -1;
           thigher.Set(systOutNames.at(lastElement), realReg.at(j), Figure(higherJEC*100, 0));
       }

        //now find the extreme values
        for(uint32_t s=0; s<systOutNames.size(); s++)
        {
            vector<float> oneSysts;
            for(uint32_t t=0; t<realReg.size(); t++)
            {
                Figure oneSyst = thigher.Get(systOutNames.at(s), realReg.at(t));
                oneSysts.push_back(oneSyst.value());
            }
            sort(oneSysts.begin(), oneSysts.end());
            tsummary.Set(uncRange.at(0), systOutNames.at(s), oneSysts.at(0));
            tsummary.Set(uncRange.at(1), systOutNames.at(s), oneSysts.at(realReg.size()-1));
        }

        //write output tables
        tI.Print(static_cast<string>("tableUncZnunu.tab"),2);
        tI.PrintLatex(static_cast<string>("tableUncZnunu.tex"),2);
        trel.Print(static_cast<string>("tableUncZnunuRel.tab"),2);
        trel.PrintLatex(static_cast<string>("tableUncZunuRel.tex"),2);
        thigher.Print(static_cast<string>("tableUncZnunuRelHiger.tab"),2);
        thigher.PrintLatex(static_cast<string>("tableUncZunuRelHigher.tex"),2);
        tsummary.Print(static_cast<string>("tableUncZnunuRelSummaryTab.tab"),2);
        tsummary.PrintLatex(static_cast<string>("tableUncZunuRelSummaryTab.tex"),2);


       //TCanvas *can = new TCanvas("can","can");
       //can->cd();
       TFile fi2("ttZuncertainties.root","RECREATE");
       histo.at(0)->Write();
       /*for(uint32_t c=0; c<histo.size(); c++)
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
       }*/
       fi2.Close();
       //can->SaveAs("plotPUyield.eps"); //@MJ@ TODO svae to root file
}
