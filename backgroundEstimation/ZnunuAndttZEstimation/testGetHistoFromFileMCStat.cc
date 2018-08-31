//to get histogram(s) from file and canvas
#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

#include <iostream>
#include "TList.h"
#include "THStack.h"

using namespace std;

int main(){

        vector<string> processes;
	theDoctor::SonicScrewdriver sonic(false);
        sonic.LoadXMLConfig("config.xml");
        //sonic.GetProcessClassTagList(&processes);
        processes.push_back("ttZ");
        vector<string> systematicsDN = {""};
        //vector<string> systematicsUP = {"PDFup"};
        //vector<string> systematicsDN = {"PDFdown"};
        //vector<string> systematicsUP = {"alphaSup"};
        //vector<string> systematicsDN = {"alphaSdown"};
        //vector<string> systematicsUP = {"Q2up" };
        //vector<string> systematicsDN = {"Q2down"};
        //vector<string> regions = {"SR1l_AB_250lessMETlessInf", "SR1l_CD_250lessMETlessInf", "SR1l_EFGH_250lessMETlessInf", "SR1l_I_250lessMETlessInf"};
        //vector<string> regions = {"SR1l_AB_250lessMETlessInf", "SR1l_CDEFGH_250lessMETlessInf", "SR1l_I_250lessMETlessInf"};
        vector<string> regions = {"SR1l_AB_250lessMETlessInf", "SR1l_CDEFGH_250lessMETlessInf"};
        //vector<string> regions = {"SR1l_NJlowTM_250lessMETlessInf", "SR1l_NJmidTM_250lessMETlessInf", "SR1l_NJhighTM_250lessMETlessInf"};
        //vector<string> regions = {"SR1l_NJ_250lessMETlessInf"};
        //vector<string> variable = {"METAB", "METCD", "MET3EFGHI", "MET3EFGHI"};
        vector<string> variable = {"METN2", "METN2"};
        //vector<string> variable = {"Mlb", "Mlb"};
        //vector<string> variable = {"Njets", "Njets","Njets"};
        //vector<string> variable = {"Njets"};
        //vector<string> variable = {"MET"};
        //vector<string> variable = {"topnessMod", "topnessMod" };
        //string sample = "METplotsttZNLOplotsRatios";
        string sample = "METplotsttZplotsRatios";

        //per each region
        for(uint32_t r =0; r<regions.size(); r++)
        {

		vector<TH1D*> down;
	//	vector<TH1D*> up;
                TH1D* absUnc;
                //read systematics
                for(uint32_t u=0; u<systematicsDN.size(); u++)
                {
		    vector<TH1D*> h = sonic.Get1DHistoCloneFromFile(sample, "1DSuperimposedNoNorm", variable.at(r),processes,regions.at(r)+systematicsDN.at(u),"lepChannel");
		    down.push_back(h.at(0));
                    cout << "region dn " << regions.at(r)+systematicsDN.at(u) << endl;
                }
		
                absUnc = dynamic_cast<TH1D*>(down.at(0)->Clone());

                ///fil/ the histogram with rel Syst
                for(uint32_t b = 0; b<down.at(0)->GetNbinsX(); b++)
                {
                    Double_t Yield = 0;
                    Double_t dnYield = 0;
                    Double_t all = 0;
                    for(uint32_t hist = 0; hist<down.size(); hist++)
                    {
                         Yield = down.at(hist)->GetBinContent(b+1);
                         dnYield = down.at(hist)->GetBinContent(b+1) - down.at(hist)->GetBinError(b+1);
                         //if(upYield+dnYield != 0)
                         all+= (abs(Yield-dnYield))*(abs(Yield-dnYield)); //@MJ@ TODO protect for 0,, start here!!!!
                         //cout << "syst " <<  hist << " rel bin" << b << " " << (abs(upYield-dnYield))/(abs(upYield+dnYield)) << endl ; 
                    }
                    cout << "all " << sqrt(all) << endl;
                    absUnc->SetBinContent(b+1,sqrt(all));
                }


               //open file and read pads
	       TFile* f = NULL;
               f = new TFile((sample+"/1DDataMCComparison.root").c_str());
	       TCanvas * c = NULL;
	       c = dynamic_cast<TCanvas*>(f->Get(("lepChannel/"+regions.at(r)+"/"+variable.at(r)).c_str()));
	       TList *t = NULL;
	       t = c->GetListOfPrimitives();

               //main pad
	       TPad *main = NULL;
	       main = dynamic_cast<TPad*>(t->At(0));
	       TList *m = 0;
	       m = main->GetListOfPrimitives();
	       cout << "m " << m << " last " << m->Last() << endl;
	       TH1D* theHist = dynamic_cast<TH1D*>(m->At(2));
	       THStack* theHistStack = dynamic_cast<THStack*>(m->At(1));
               cout << "stack " << theHistStack << endl; 

               //ratio pad
	       TPad *ratio = NULL;
	       ratio = dynamic_cast<TPad*>(t->At(1));
	       TList* o = NULL;
	       o = ratio->GetListOfPrimitives();
	       TH1D* theRatio = dynamic_cast<TH1D*>(o->At(0));

               //set statistics and draw new histograms
               if(absUnc->GetNbinsX() != theHist->GetNbinsX())
                   cout << "somethin went terribly wrong" << endl;
               if(absUnc->GetNbinsX() != theRatio->GetNbinsX())
                   cout << "somethin went terribly wrong again" << endl;
               TH1D* theRatioDN = dynamic_cast<TH1D*>(theRatio->Clone());
               TH1D* theRatioUP = dynamic_cast<TH1D*>(theRatio->Clone());
               for(uint32_t b=0; b<absUnc->GetNbinsX(); b++)
               {
                   Double_t bValue = theHist->GetBinContent(b+1);
                   Double_t bAbs = absUnc->GetBinContent(b+1);
                   cout << "b value " << bValue << " b error " << bAbs << "b relative error " << bAbs/bValue << endl;
                   theHist->SetBinError(b+1, bAbs);
	           theRatio->SetBinContent(b+1,1);
                   if(bValue !=0)
                   {
	               theRatio->SetBinError(b+1,bAbs/bValue);
	               theRatioDN->SetBinContent(b+1,1-(bAbs/bValue));
	               theRatioUP->SetBinContent(b+1,1+(bAbs/bValue));
                   }
                   else
                   {
	               theRatio->SetBinError(b+1,0);
	               theRatioDN->SetBinContent(b+1,1);
	               theRatioUP->SetBinContent(b+1,1);

                   }

               }
               main->cd();
               theHist->SetFillColor(kMagenta-3);
               theHist->SetLineColor(kMagenta-3);
               theHist->SetFillStyle(1001);
               theHist->SetMarkerStyle(8);
               theHist->SetMarkerColor(kBlack);
               //TH1* theHistStackH= theHistStack->GetHistogram();
               //cout << "nstack hist" << theHistStack->GetNhists() << endl;
               //theHistStackH->SetFillColor(kBlack);
               //theHistStackH->SetFillStyle(1001);
               //theHistStackH->Draw("same");
               //theHistStack->Paint("noclear");
               theHist->Draw("same E2");

	       ratio->cd();
	       theRatio->SetMarkerStyle(1);
	       theRatio->SetFillColor(kMagenta-3);
	       theRatio->SetFillStyle(1001);
	       theRatio->Draw("E2");
	       theRatioDN->Draw("same");
	       theRatioUP->Draw("same");
	       cout << "c " << c << " p1 " << main << " p2 " << ratio << "TH1D " << theRatio << " TH1 main " << theHist << endl;
               

	       TCanvas *can = new TCanvas((regions.at(r)+"Can").c_str(),(regions.at(r)+"Can").c_str());
	       can->cd();
	       main->Draw();
	       ratio->Draw("same");
	       can->SaveAs((regions.at(r)+variable.at(r)+sample+"Statistics.root").c_str());
	       can->SaveAs((regions.at(r)+variable.at(r)+sample+"Statistics.eps").c_str());
               //delete main;
               //delete ratio;
               //delete can;
               //delete f;

       }
 
}
