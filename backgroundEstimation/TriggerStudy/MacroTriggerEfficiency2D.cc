#include <iostream>
#include <string>

using namespace std;

#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"

#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

bool f(){return true;};
	
void AddEfficiency(TFile* fout, theDoctor::SonicScrewdriver& sonic, string var, string processClass, string trigger, string refTrigger){
	string all = "all";
	string cname = var+string("_")+trigger+string("_")+refTrigger;
	TCanvas* c1 = new TCanvas(cname.c_str());
	TH1D* hNum = sonic.Get1DHistoClone(var.c_str(), processClass, trigger.c_str(), refTrigger.c_str());
	TH1D* hDenom = sonic.Get1DHistoClone(var.c_str(), processClass, all, refTrigger.c_str());
	// create a TEfficiency
	TEfficiency* hEff = new TEfficiency(*hNum,*hDenom);
	string title = string(";")+var+string(" [GeV] ; Efficiency");
	hEff->SetTitle(title.c_str());
	c1->cd();
	hEff->Draw("AP");
	c1->Update();
	//hEff->GetPaintedGraph()->GetXaxis()->SetTitle(var.c_str());
	hEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.);
	hEff->GetPaintedGraph()->Draw("AP");
	//hEff.GetXaxis()->SetTitle(hNum->GetXaxis()->GetTitle());
	//hEff->Draw("AP");
	c1->Update();
	fout->cd();
	c1->Write();
}

int main(){
	//style for error printing
        gStyle->SetPaintTextFormat("3.2f");

	//TFile* fin = TFile::Open("plotsTrigger/1DSuperimposed.root");
	//string ifilename = "plotsTrigger/1DSuperimposed.root";
	string ifilename = "plotsTrigger/2D.root";
	
	//string region = "CombinedMETEl";
	//string region = "METTrigger";
	//string region = "CombinedMETMHTMu";
	//string region = "CombinedMET";
	string region = "CombinedMETMu";
	//string region = "CombinedMETMHTEl";
	//string region = "CombinedMETEl";
	//string channel = "All";
	//string channel = "DoubleMuTrigger";
	//string channel = "DoubleMuTriggerAndMuSel";
	//string channel = "DoubleElTriggerAndElSel";
	//string channel = "DoubleMuTriggerAndDiMuSel";
	//string channel = "DoubleElTriggerAndDiElSel";
	//string channel = "DoubleMuTriggerAndDiMuSelBaseline";
	//string channel = "DoubleElTriggerAndDiElSelBaseline";
	string channel = "MuChannel";
	//string channel = "LepChannel";
	string variables = "LeptonPT[vs]MET";
	//string variables = "nJets[vs]MET";
	//string var_plotlabel = "vX:nJets|vY:MET";
	string var_plotlabel = "vX:LeptonPT|vY:MET";
	//string title = string("; Lepton PT [GeV]; MET [GeV]");
	string title = string("; nJets; MET [GeV]");
	string pname = "test";
	string cfilename = region+"_"+channel+".pdf";
	cout<<"toto"<<endl;
	TFile* fin = TFile::Open("plotsTrigger/2D.root");
	cout<<"toto"<<endl;
	TFile* fout = new TFile("eff2D.root","RECREATE");

	// global variables
	TDirectory* d;
	TCanvas* c;

        // --------------------------
        // Estimation of 
	// --------------------------

	// retrieve the numerator
	string dir = channel+string("/")+region+string("/")+variables;
	cout<<dir<<endl;
	d = (TDirectory*) fin->Get(dir.c_str());
	c = (TCanvas*) d->Get(pname.c_str());
	//string name2D = string("vX:LeptonPT|vY:MET|p:")+pname+string("|r:")+region+string("|c:")+channel+string("|t:2DEntries");
	string name2D = var_plotlabel+string("|p:")+pname+string("|r:")+region+string("|c:")+channel+string("|t:2DEntries");
	cout<<name2D<<endl;
	c->GetListOfPrimitives()->Print();
	TH2D* hNum = (TH2D*) c->GetPrimitive(name2D.c_str());

	cout<<"Numerator: "<<hNum<<endl;
	// retrive the denominator
	region = "all";
	//region = "ElChannel";
	dir = channel+string("/")+region+string("/")+variables;
	d = (TDirectory*) fin->Get(dir.c_str());
	cout<<dir<<" "<<d<<endl;
	c = (TCanvas*) d->Get(pname.c_str());
	cout<<"c = "<<c<<endl;
	//name2D = string("vX:LeptonPT|vY:MET|r:")+region+string("|c:")+channel+string("|t:2DEntries");
	//name2D = string("vX:LeptonPT|vY:MET|p:")+pname+string("|r:")+region+string("|c:")+channel+string("|t:2DEntries");
	name2D = var_plotlabel+string("|p:")+pname+string("|r:")+region+string("|c:")+channel+string("|t:2DEntries");
	TH2D* hDenom = (TH2D*) c->GetPrimitive(name2D.c_str());
	cout<<"Denomirator: "<<hNum<<endl;
	
	// create Efficiency
	TCanvas* c1  = new TCanvas("c1");
	//TEfficiency* hEff = new TEfficiency(*hNum,*hDenom);
	//hNum->GetZaxis()->SetRangeUser(0., 1.);
	TH2D* hEff = (TH2D*) hNum->Clone();
	hEff->Divide(hDenom);
	cout<<"compute efficiency !"<<endl;
	hEff->SetTitle(title.c_str());
	hEff->GetZaxis()->SetRangeUser(0., 1.);
	c1->cd();
	hEff->Draw("COLETEXT");
	c1->Update();
	//hEff->GetPaintedGraph()->GetXaxis()->SetTitle(var.c_str());
	//hEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.);
	//hEff->GetPaintedGraph()->Draw("AP");
	//hEff.GetXaxis()->SetTitle(hNum->GetXaxis()->GetTitle());
	//hEff->Draw("AP");
	c1->Update();
	fout->cd();
	c1->Write();
	//c1->Print("c1.pdf");
	c1->Print(cfilename.c_str());
	fout->Write();
	fout->Close();
}

