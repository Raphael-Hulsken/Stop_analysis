#include <iostream>
#include <string>

using namespace std;

#include "TFile.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

bool f(){return true;};
	
TEfficiency* AddEfficiency(TFile* fout, theDoctor::SonicScrewdriver& sonic, string var, string processClass, string trigger, string refTrigger){
	string all = "all";
	string cname = var+string("_")+trigger+string("_")+refTrigger;
	TCanvas* c1 = new TCanvas(cname.c_str());
	TH1D* hNum = sonic.Get1DHistoClone(var.c_str(), processClass, trigger.c_str(), refTrigger.c_str());
	TH1D* hDenom = sonic.Get1DHistoClone(var.c_str(), processClass, all, refTrigger.c_str());
	// create a TEfficiency
	hNum->GetYaxis()->SetRangeUser(0,1.);
	hDenom->GetYaxis()->SetRangeUser(0,1.);
	TEfficiency* hEff = new TEfficiency(*hNum,*hDenom);
	string title = string(";")+var+string(" [GeV] ; Efficiency");
	hEff->SetTitle(title.c_str());
	//hEff->GetPaintedHistogram()->GetYaxis()->SetRangeUser(0,1);
	c1->cd();
	//hEff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.);
	//hEff->GetPassedHistogram()->GetYaxis()->SetRangerUser(0,1);
	hEff->Draw("AP");
	//hEff->GetPaintedGraph()->GetXaxis()->SetTitle(var.c_str());
	//hEff->GetPaintedGraph()->Draw("AP");
	c1->Update();
	//hEff.GetXaxis()->SetTitle(hNum->GetXaxis()->GetTitle());
	//hEff->Draw("AP");
	c1->Update();
	string cfilename = trigger+string("_")+refTrigger+string("_")+var+string(".pdf");
	c1->Print(cfilename.c_str());
	fout->cd();
	//c1->Write();
	//return c1;
	return hEff;
}

//----------------------------------------------
// Superimpose 2 graphs and save that in a file
//----------------------------------------------
void SuperimposeEfficiency(TEfficiency* h1, TEfficiency* h2, string leg1, string leg2, string cname ){
	TCanvas c("c");
	c.cd();
	h1->SetLineColor(kBlue);
	h2->SetLineColor(kRed);
	h1->Draw("AP");
	h2->Draw("same");
	TLegend leg(0.15,0.7,0.3,0.85);
	leg.AddEntry(h1,leg1.c_str(),"l");
	leg.AddEntry(h2,leg2.c_str(),"l");
	leg.Draw("same");
	c.Print(cname.c_str());
}

int main(){

	//TFile* fin = TFile::Open("plotsTrigger/1DSuperimposed.root");
	string ifilename = "plotsTrigger/1DSuperimposed.root";
	//string processClass = "signal";
	string processClass = "test";
	//------------------------
	// Computate efficiency
	//------------------------
   	//numerator
	vector<string> AnaTrigger = {"METTrigger","METMHTTrigger","ElTrigger","MuTrigger","ORTrigger"}; //correspond to Regions in the file
 	//denominator
	string all = "all"; //no cut - correspond to a Region in the file
	
	//Several trigger used as reference - so several measurement (not independant)
	vector<string> RefTrigger = {"ORTriggerAndMuSel","ORTriggerAndElSel","ElTriggerAndElSel","MuTriggerAndMuSel"}; // corresponds to Channel in the file
	
	//create an output file
	TFile* fout = new TFile("eff.root","RECREATE");

	//import all histos
	theDoctor::SonicScrewdriver sonic;
        
	//sonic.AddProcessClass("test", "test");
        sonic.AddProcessClass("test", "test", "signal", kBlue);
	
        //sonic.AddVariable("MET", "MET");
        //sonic.AddVariable("LeptonPT", "LeptonPT");
        
        vector<float> METBins = {150,175,200,225,250,275,300,350,400,500,800};
        vector<float> LepPtBins = {20,25,30,35,40,50,80,100,125,150,200,300};
	float dummy;
	//sonic.AddVariable("MET", "MET",  "MET", 20,   20, 500,  &(dummy), "noUnderflowInFirstBin");
        sonic.AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", 20,   0, 200,  &(dummy), "noUnderflowInFirstBin");
        //sonic.AddVariable("MET", "MET",  "MET", (int) (METBins.size()-1), METBins.data(),  &(dummy), "");
        sonic.AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", (int) (LepPtBins.size()-1), LepPtBins.data(), &(dummy), "");
        sonic.AddVariable("MET", "MET",  "MET", 4 ,150,550,  &(dummy), "");
        sonic.AddVariable("nJets","nJets","nJets",5,1,5,&(dummy),"");
    
        sonic.AddChannel("ORTriggerAndMuSel","",&f,"");
        sonic.AddChannel("ORTriggerAndElSel","",&f,"");
        sonic.AddChannel("ElTriggerAndElSel","",&f,"");
        sonic.AddChannel("MuTriggerAndMuSel","",&f,"");
        sonic.AddChannel("METTriggerAndMuSel","",&f,"");
	sonic.AddChannel("METTriggerAndElSel","",&f,"");
	sonic.AddChannel("MuTriggerAndMETSel","",&f,"");
	sonic.AddChannel("ElTriggerAndMETSel","",&f,"");
	sonic.AddChannel("METTrigger","",&f,"");
	sonic.AddChannel("All","",&f,"");
        sonic.AddChannel("ElChannel","ElChannel", &f, "");
        sonic.AddChannel("MuChannel","MuChannel", &f, "");
        sonic.AddChannel("LepChannel","MuChannel", &f, "");

        sonic.AddChannel("DoubleMuTrigger","DoubleMuTrigger",&f,"");
        sonic.AddChannel("DoubleElTrigger","DoubleElTrigger",&f,"");
        sonic.AddChannel("MuElTrigger","MuElTrigger",&f,"");
        sonic.AddChannel("DoubleLeptonTrigger","DoubleLeptonTrigger",&f,"");
        sonic.AddChannel("DoubleMuTriggerAndMuSel","PassDoubleMuTriggerAndMuSel",&f);
        sonic.AddChannel("DoubleElTriggerAndElSel","PassDoubleElTriggerAndElSel",&f);
        sonic.AddChannel("DoubleMuTriggerAndDiMuSel","PassDoubleMuTriggerAndDiMuSel",&f);
        sonic.AddChannel("DoubleElTriggerAndDiElSel","PassDoubleElTriggerAndDiElSel",&f);
        sonic.AddChannel("DoubleMuTriggerAndDiMuSelBaseline","PassDoubleMuTriggerAndDiMuSelBaseline",&f);
        sonic.AddChannel("DoubleElTriggerAndDiElSelBaseline","PassDoubleElTriggerAndDiElSelBaseline",&f);


        sonic.AddRegion("all","all",&f,"");
        sonic.AddRegion("METTrigger","MET Trigger",&f,"");
        sonic.AddRegion("METMHTTrigger","MET Trigger",&f,"");
        sonic.AddRegion("ElTrigger","Electron Trigger",&f,"");
        sonic.AddRegion("MuTrigger","Muon Trigger",&f,"");
        sonic.AddRegion("ORTrigger","ORTrigger",&f,"");
        sonic.AddRegion("CombinedMET","CombinedMETTrigger", &f,"");
        sonic.AddRegion("CombinedMETMu","CombinedMETMuTrigger", &f,"");
        sonic.AddRegion("CombinedMETEl","CombinedMETElTrigger", &f, "");
        sonic.AddRegion("CombinedMETMHTMu","CombinedMETMHTMuTrigger", &f, "");
        sonic.AddRegion("CombinedMETMHTEl","CombinedMETMHTElTrigger", &f, "");
    
    
	sonic.Create1DHistos();

	cout<<"before import"<<endl;
	//sonic.ImportHistosEntries(ifilename);
	sonic.ImportHistosFromFile(ifilename);
	cout<<"after import"<<endl;
	
	vector<theDoctor::Histo1DEntries>* entries;
	entries = sonic.Get1DHistosEntries();
	cout<<"nof files:"<<entries->size()<<endl;
	for(unsigned int i=0;i<(*entries).size();i++){
		cout<<(*entries)[i].getEntriesHisto()->GetName()<<endl;
	}
	cout<<"there"<<endl;
	
	//Compute the efficiencies via several references
	/*
	for(unsigned int ref=0; ref<RefTrigger.size(); ref++){
		cout<<"RefTrigger = "<<RefTrigger[ref]<<endl;
		//Compute the efficiency for each trigger used
		for(unsigned int trig=0;trig<AnaTrigger.size();trig++){
			cout<<" > "<<AnaTrigger[trig]<<endl;
			// --------------------------------
      			//  Efficiency as function of MET
			// --------------------------------
			TH1D* hNum = sonic.Get1DHistoClone("MET", processClass, AnaTrigger[trig], RefTrigger[ref]);
			TH1D* hDenom = sonic.Get1DHistoClone("MET", processClass, all, RefTrigger[ref]);
			cout<<hNum<<" "<<hDenom<<" "<<hNum->GetEntries()<<" "<<hDenom->GetEntries()<<endl;
			cout<<hNum->GetEntries()<<" & "<<hDenom->GetEntries()<<endl;
			cout<<hNum->Integral()<<" && "<<hDenom->Integral()<<endl;
			//hNum->SetBinContent(0, 0);
			//hNum->SetBinContent(hNum->GetNbinsX()+1,0);
			//hDenom->SetBinContent(0, 0);
			//hDenom->SetBinContent(hNum->GetNbinsX()+1,0);
			for(int b=0;b<hNum->GetNbinsX()+2;b++)
				if(hNum->GetBinContent(b)>hDenom->GetBinContent(b))
					cout<<"bins: "<<b<<" > "<<hNum->GetBinContent(b)<<" "<<hDenom->GetBinContent(b)<<endl;
			        else
					cout<<"bins: "<<b<<" < "<<hNum->GetBinContent(b)<<" "<<hDenom->GetBinContent(b)<<endl;
			// create a TEfficiency
			TEfficiency hEff(*hNum,*hDenom);
			//hEff.GetXaxis()->SetTitle(hNum->GetXaxis()->GetTitle());
			//hEff.GetYaxis()->SetTitle("eff");
			cout<<" HERE "<<endl;
			//change name and title
			// save it on the file
			fout->cd(); hEff.Write();
			
			// --------------------------------
      			//  Efficiency as function of lepPT
			// --------------------------------
			TH1D* hNum2 = sonic.Get1DHistoClone("LeptonPT", processClass, AnaTrigger[trig], RefTrigger[ref]);
			TH1D* hDenom2 = sonic.Get1DHistoClone("LeptonPT", processClass, all, RefTrigger[ref]);
			cout<<hNum2<<" "<<hDenom2<<" "<<hNum2->GetEntries()<<" "<<hDenom->GetEntries()<<endl;
			// create a TEfficiency
			TEfficiency hEff2(*hNum2,*hDenom2);
			//hEff2.GetXaxis()->SetTitle(hNum2->GetXaxis()->GetTitle());
			//hEff2.GetYaxis()->SetTitle("eff");
			cout<<" HERE2 "<<endl;
			// save it on the file
			fout->cd(); hEff2.Write();
			
		
		}	
	}
	*/

	TCanvas* c1 = new TCanvas("c1");
	TH1D* hNum = sonic.Get1DHistoClone("MET", processClass, "ElTrigger", "METTriggerAndElSel");
	cout<<"numerator = "<<hNum->GetEntries()<<endl;
	TH1D* hDenom = sonic.Get1DHistoClone("MET", processClass, all,"METTriggerAndElSel"); 
	cout<<"denominator = "<<hDenom->GetEntries()<<endl;
	// create a TEfficiency
	TEfficiency* hEff = new TEfficiency(*hNum,*hDenom);
	//hEff.GetXaxis()->SetTitle(hNum->GetXaxis()->GetTitle());
	//hEff.GetYaxis()->SetTitle("eff");
	c1->cd();
	hEff->Draw("AP");
	//c1->Write();
	/*
        AddEfficiency(fout, sonic, "MET", processClass, "ElTrigger", "METTriggerAndElSel");
        AddEfficiency(fout, sonic, "MET", processClass, "MuTrigger", "METTriggerAndMuSel");
        AddEfficiency(fout, sonic, "LeptonPT", processClass, "ElTrigger", "METTrigger");
        AddEfficiency(fout, sonic, "LeptonPT", processClass, "MuTrigger", "METTrigger");
        AddEfficiency(fout, sonic, "LeptonPT", processClass, "ElTrigger", "METTriggerAndElSel");
        AddEfficiency(fout, sonic, "LeptonPT", processClass, "MuTrigger", "METTriggerAndMuSel");
        AddEfficiency(fout, sonic, "LeptonPT", processClass, "METTrigger", "ElTriggerAndElSel");
        AddEfficiency(fout, sonic, "LeptonPT", processClass, "METTrigger", "MuTriggerAndMuSel");
	AddEfficiency(fout, sonic, "MET", processClass, "METMHTTrigger", "ElTriggerAndElSel");
	AddEfficiency(fout, sonic, "MET", processClass, "METMHTTrigger", "MuTriggerAndMuSel");
	
	
        AddEfficiency(fout, sonic, "MET", processClass, "ElTrigger", "All");
        AddEfficiency(fout, sonic, "MET", processClass, "MuTrigger", "All");
	
	AddEfficiency(fout, sonic, "MET", processClass, "METTrigger", "MuTriggerAndMETSel" );
	AddEfficiency(fout, sonic, "MET", processClass, "METTrigger", "ElTriggerAndMETSel" );
	AddEfficiency(fout, sonic, "MET", processClass, "METMHTTrigger", "MuTriggerAndMETSel" );
	AddEfficiency(fout, sonic, "MET", processClass, "METMHTTrigger", "ElTriggerAndMETSel" );
	*/
	//AddEfficiency(fout, sonic, "MET",processClass,"ElTrigger","All");
	//AddEfficiency(fout, sonic, "MET",processClass,"ElTrigger","inclusive");
	//AddEfficiency(fout, sonic, "MET",processClass,"MuTrigger","inclusive");
	//AddEfficiency(fout, sonic, "LeptonPT",processClass,"METTrigger","inclusive");
	//*/

	/*
	AddEfficiency(fout,sonic,"MET", processClass, "METTrigger", "DoubleMuTrigger");
	AddEfficiency(fout,sonic,"MET", processClass, "METTrigger", "DoubleElTrigger");
	AddEfficiency(fout,sonic,"MET", processClass, "METTrigger", "MuElTrigger");
	AddEfficiency(fout,sonic,"MET", processClass, "METTrigger", "DoubleLeptonTrigger");
	*/
	//string region = "CombinedMETMHTMu";
	//string region = "CombinedMETMu";
	//string region = "CombinedMETMHTEl";
	string region = "CombinedMETEl";
	//string channel = "All";
	//string channel = "DoubleMuTrigger";
	//string channel = "DoubleMuTriggerAndMuSel";
	string channel = "DoubleElTriggerAndElSel";
	
	
	//Dilepton-trigger - single lepton selection
	TEfficiency* hMET_El_fMET = AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETEl", "DoubleElTriggerAndElSel");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTEl", "DoubleElTriggerAndElSel");
	TEfficiency* hMET_El_fnJets = AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETEl", "DoubleElTriggerAndElSel");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMHTEl", "DoubleElTriggerAndElSel");
	
	TEfficiency* hMET_Mu_fMET = AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMu","DoubleMuTriggerAndMuSel");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTMu","DoubleMuTriggerAndMuSel");
	TEfficiency* hMET_Mu_fnJets = AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMu","DoubleMuTriggerAndMuSel");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMHTMu","DoubleMuTriggerAndMuSel");
	

	//Dilepton-trigger - double lepton  selection
	TEfficiency* hMET_DiEl_fMET = AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETEl", "DoubleElTriggerAndDiElSel");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTEl", "DoubleElTriggerAndDiElSel");
	TEfficiency* hMET_DiEl_fnJets = AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETEl", "DoubleElTriggerAndDiElSel");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMHTEl", "DoubleElTriggerAndDiElSel");

	TEfficiency* hMET_DiMu_fMET = AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMu","DoubleMuTriggerAndDiMuSel");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTMu","DoubleMuTriggerAndDiMuSel");
	TEfficiency* hMET_DiMu_fnJets = AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMu","DoubleMuTriggerAndDiMuSel");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMHTMu","DoubleMuTriggerAndDiMuSel");
	
	//Dilepton-trigger - double lepton selection + MET>250+ Njet>=2
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETEl", "DoubleElTriggerAndDiElSelBaseline");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTEl", "DoubleElTriggerAndDiElSelBaseline");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETEl", "DoubleElTriggerAndDiElSelBaseline");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMHTEl", "DoubleElTriggerAndDiElSelBaseline");

	
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMu","DoubleMuTriggerAndDiMuSelBaseline");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTMu","DoubleMuTriggerAndDiMuSelBaseline");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMu","DoubleMuTriggerAndDiMuSelBaseline");
	AddEfficiency(fout,sonic,"nJets", processClass,  "CombinedMETMHTMu","DoubleMuTriggerAndDiMuSelBaseline");	

	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTEl","ElChannel");
	AddEfficiency(fout,sonic,"MET", processClass,  "CombinedMETMHTMu","MuChannel");



	//Superpose graph
	/*
	TCanvas* cEl_Diff = new TCanvas("c_ElDiff");
	cEl_Diff->cd();
	hMET_El_fMET->SetLineColor(kBlue);
	hMET_DiEl_fMET->SetLineColor(kRed);
	hMET_El_fMET->Draw("AP");
	hMET_DiEl_fMET->Draw("same");
	TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
	leg->AddEntry(hMET_El_fMET,"1-lep","l");
	leg->AddEntry(hMET_DiEl_fMET,"2-lep","l");
	leg->Draw("same");
	cEl_Diff->Print("Diff_El_CombMET.pdf");
	*/
        SuperimposeEfficiency(hMET_El_fMET,hMET_DiEl_fMET,"1-lep","2-lep","Diff_El_CombMET_MET.pdf");
        SuperimposeEfficiency(hMET_El_fnJets,hMET_DiEl_fnJets,"1-lep","2-lep","Diff_El_CombMET_nJets.pdf");
        SuperimposeEfficiency(hMET_El_fMET,hMET_DiMu_fMET,"1-lep","2-lep","Diff_Mu_CombMET_MET.pdf");
        SuperimposeEfficiency(hMET_El_fMET,hMET_DiMu_fnJets,"1-lep","2-lep","Diff_Mu_CombMET_nJets.pdf");


        fout->Write();
	fout->Close();
}




