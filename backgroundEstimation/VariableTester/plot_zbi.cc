#include "../../sonicScrewdriver/interface/Table.h"
#include "../../sonicScrewdriver/interface/Figure.h"

#include <iostream>

using namespace std;
using namespace theDoctor;

int main(){
   
    vector<string> signals = {"(175,0)","(200,25)","(225,50)","(250,75)"};
    vector<string> inputRegionTags = {"SR_2l_2j","SR_2l_3j","SR_2l_4j","SR_2l_ISR"};
    
    // vectors below should have the same size
    vector<string> files = {"isr_150/card_zbi.tab", "isr_200/card_zbi.tab", "isr_250/card_zbi.tab", "isr_300/card_zbi.tab", "isr_350/card_zbi.tab", "isr_400/card_zbi.tab"};
    vector<EColor> colors = {kBlack, kBlue, kGreen, kMagenta, kRed, kOrange};
    vector<string> labels = {"jet p_{T} > 150 GeV","jet p_{T} > 200 GeV","jet p_{T} > 250 GeV","jet p_{T} > 300 GeV","jet p_{T} > 350 GeV","jet p_{T} > 400 GeV",};


    TH1F***  plots = new TH1F**[files.size()];
    Table** tables = new Table*[files.size()];	
    for(int i=0 ;i<files.size();i++){
      plots[i] = new TH1F*[signals.size()];
      tables[i] = new Table(files[i]);
      for(int j=0 ;j<signals.size();j++){
      	 string name = files[i]+"_"+signals[j];
	 plots[i][j] = new TH1F(name.c_str(),"",4,1,4);
      }

    }
    
    //Analysis of the tables
    Figure fig;

    TLegend* leg = new TLegend(0.1,0.1,0.3,0.3);
    for(int i=0 ;i<files.size();i++){
	for(int j=0;j<signals.size();j++){
		plots[i][j]->SetLineWidth(2);
		plots[i][j]->SetLineColor(colors[i]);
		for(int k=0;k<inputRegionTags.size();k++){
			fig = tables[i]->Get(inputRegionTags[k],signals[j]);
			plots[i][j]->SetBinContent(k+1,fig.value());
		}
	}

    }

    for(int i=0;i<files.size();i++){
    	leg->AddEntry(plots[i][0],labels[i].c_str(),"l");
    }
    /*
    Table table("card_zbi.tab"); 
    Figure fig;
    fig = table.Get("SR_2l_2j","(200,0)");
    

    cout<<fig.Print()<<endl;
    */

    TFile fout("plots.root","RECREATE");
    fout.cd();

    gStyle->SetOptStat(0);

    TCanvas** canvas = new TCanvas*[signals.size()];
    for(int i=0;i<signals.size();i++){
    	string name = "c_"+signals[i];
    	canvas[i] = new TCanvas(name.c_str());
	canvas[i]->cd();
	// plots
	for(int j=0;j<files.size();j++){
		plots[j][i]->GetXaxis()->SetTitle("boxes");
		plots[j][i]->GetYaxis()->SetTitle("Zbi");
		plots[j][i]->GetYaxis()->SetRangeUser(0,6);
		
		if(j==0) plots[j][i]->Draw();
		else plots[j][i]->Draw("same");
	}
	leg->Draw("same");
        //fout.cd();	
	canvas[i]->Write();
	string filename = signals[i]+".eps";
    	canvas[i]->Print(filename.c_str());
    }

   fout.Write();

}
