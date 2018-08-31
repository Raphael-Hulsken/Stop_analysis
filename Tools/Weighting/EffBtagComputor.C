#include "TFile.h"
#include "TH2F.h"
#include "TChain.h"

#include <vector>
#include <string>
#include <iostream>

#include "../../sonicScrewdriver/test/common.h"

using namespace std;


#define BTAGCUT 0.80

//-----------
// The name of the plots are not yet the one used in WeighFactory
// ----------

int main(){
	
	cout<<"##################################"<<endl;
	cout<<"##################################"<<endl;
	cout<<" Btagging efficiency computor    "<<endl;
	cout<<"##################################"<<endl;

	//--------------------------------
	// Configuration
	//--------------------------------
	string path = "/opt/sbg/scratch1/cms/mjansova/store/tmp/0709/";
	string file_basename = "";
	vector<string> filenames = {"TTJetsSLtop.root","TTJetsSLatopv1.root","TTJetsSLatopext.root","TTJetsDLv0v4.root","TTJetsDLext.root"};
	cout<<"Files from the directory "<<path<<" will be used: "<<endl;
	for(unsigned int i=0;i<filenames.size();i++) cout<<filenames[i]<<" ";
	cout<<endl;
	

	//--------------------------------
	//--------------------------------
	// Create a chain
	//--------------------------------

	TChain* chain = new TChain("babyTuple");
	//Add the files
	for(unsigned int i=0;i<filenames.size();i++){
		string fname = path+filenames[i];
		chain->Add(fname.c_str());
	}

	//--------------------------------
	// Load only relevant variables
	//--------------------------------
	
	// desactivate everything by default
	chain->SetBranchStatus("*",0);
	// activate the desired branches
	// selection: njets, nlep, MET, MT
	cout<<"# set status for branches of global variables "<<endl;
	chain->SetBranchStatus("numberOfSelectedJets",1);
	chain->SetBranchStatus("numberOfSelectedLeptons",1);
	chain->SetBranchStatus("pfmet",1);
	chain->SetBranchStatus("mt_met_lep",1);
	int njets, nlep;
	float MET, MT;
	chain->SetBranchAddress("numberOfSelectedJets",&njets);
	chain->SetBranchAddress("numberOfSelectedLeptons",&nlep);
	chain->SetBranchAddress("pfmet",&MET);
	chain->SetBranchAddress("mt_met_lep",&MT);

	// access to the jet collection
	// jet pt/eta/CSV/flavor
	// remark: use hadronFlavour and not partonFlavour
	cout<<"# set status for branches of jets variables "<<endl;
	chain->SetBranchStatus("ak4pfjets_pt",1);
	chain->SetBranchStatus("ak4pfjets_eta",1);
	chain->SetBranchStatus("ak4pfjets_hadronFlavour",1);
	vector<float>* p_jet_pt = 0;
	vector<float>* p_jet_eta = 0;
	vector<float>* p_jet_CSV = 0;
	vector<float>* p_jet_hadronFlavour = 0;
	vector<float> jet_pt;
	vector<float> jet_eta;
	vector<float> jet_CSV;
	vector<float> jet_hadronFlavour;
	chain->SetBranchAddress("ak4pfjets_pt",&p_jet_pt);
	chain->SetBranchAddress("ak4pfjets_eta",&p_jet_eta);
	chain->SetBranchAddress("ak4pfjets_CSV",&p_jet_CSV);
	chain->SetBranchAddress("ak4pfjets_hadronFlavour",&p_jet_hadronFlavour);

	//--------------------------------
	// Define the 2D binning
	//--------------------------------
	cout<<"# Prepare the histograms"<<endl;
	vector<double> PtBins = {30,50,70,100,140,200,300,670,1000};
	vector<double> EtaBins = {0,0.4,0.8,1.2,1.6,2.0,2.4};
	

	//--------------------------------
	// Define the 2D plots
	//--------------------------------
	TH2F* h2_BTaggingEff_Denom_b = new TH2F("h2_BTaggingEff_Denom_b", ";p_{T} [GeV];#eta", PtBins.size()-1, PtBins.data(), EtaBins.size()-1, EtaBins.data());
	TH2F* h2_BTaggingEff_Denom_c = new TH2F("h2_BTaggingEff_Denom_c", ";p_{T} [GeV];#eta", PtBins.size()-1, PtBins.data(), EtaBins.size()-1, EtaBins.data());
	TH2F* h2_BTaggingEff_Denom_udsg = new TH2F("h2_BTaggingEff_Denom_udsg", ";p_{T} [GeV];#eta", PtBins.size()-1, PtBins.data(), EtaBins.size()-1, EtaBins.data());
	TH2F* h2_BTaggingEff_Num_b = new TH2F("h2_BTaggingEff_Num_b", ";p_{T} [GeV];#eta", PtBins.size()-1, PtBins.data(), EtaBins.size()-1, EtaBins.data());
	TH2F* h2_BTaggingEff_Num_c = new TH2F("h2_BTaggingEff_Num_c", ";p_{T} [GeV];#eta", PtBins.size()-1, PtBins.data(), EtaBins.size()-1, EtaBins.data());
	TH2F* h2_BTaggingEff_Num_udsg = new TH2F("h2_BTaggingEff_Num_udsg", ";p_{T} [GeV];#eta", PtBins.size()-1, PtBins.data(), EtaBins.size()-1, EtaBins.data());

	//--------------------------------
	// Loop over events
	//--------------------------------
	cout<<"# loop over the events "<<endl;
	long int nEntries = chain->GetEntries();
	for(unsigned i=0;i<nEntries;i++){
		chain->GetEntry(i);
		if(i%(nEntries/100)==0) 
			printProgressBar(i,nEntries);

		//special treatment for vector
		jet_pt = *p_jet_pt;
		jet_eta = *p_jet_eta;
		jet_CSV = *p_jet_CSV;
		jet_hadronFlavour = *p_jet_hadronFlavour;

		//  apply the selection
		if(njets>=2 && nlep>=1 && MET>250 && MT>150){
			// run over the jets collection
			for(unsigned int j=0;j<jet_pt.size();j++){
			   //jet should be selected
			   if(jet_pt[j]>30 && abs(jet_eta[j])<2.4){
				//check the flavor
				if( abs(jet_hadronFlavour[j]) == 4){
					h2_BTaggingEff_Denom_c->Fill(jet_pt[j],jet_eta[j]);
					if(jet_CSV[j]>BTAGCUT) h2_BTaggingEff_Num_c->Fill(jet_pt[j],jet_eta[j]);
						
				}
				else if ( abs(jet_hadronFlavour[j]) == 5){
					h2_BTaggingEff_Denom_b->Fill(jet_pt[j],jet_eta[j]);
					if(jet_CSV[j]>BTAGCUT) h2_BTaggingEff_Num_b->Fill(jet_pt[j],jet_eta[j]);
				}
				else{
					h2_BTaggingEff_Denom_udsg->Fill(jet_pt[j],jet_eta[j]);
					if(jet_CSV[j]>BTAGCUT) h2_BTaggingEff_Num_udsg->Fill(jet_pt[j],jet_eta[j]);
				}
			   }
			}
		}
	}
	cout<<endl;
	//--------------------------------
	// Compute the ratios and save
	//--------------------------------

	cout<<"Divide the histograms"<<endl;	
	h2_BTaggingEff_Num_b->Divide(h2_BTaggingEff_Denom_b);
	h2_BTaggingEff_Num_c->Divide(h2_BTaggingEff_Denom_c);
	h2_BTaggingEff_Num_udsg->Divide(h2_BTaggingEff_Denom_udsg);

	TFile* fout = new TFile("eff_btag.root","RECREATE");
	fout->cd();
	h2_BTaggingEff_Num_b->Write();
	h2_BTaggingEff_Num_c->Write();
	h2_BTaggingEff_Num_udsg->Write();
	fout->Close();

}

