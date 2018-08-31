#ifndef _WEIGHT_FACTORY_
#define _WEIGHT_FACTORY_

#include "../../sonicScrewdriver/test/common.h"
#include "BTagCalibrationStandalone.h"

#include "TFile.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <string>
#include <iostream>

using namespace std;

#define BTAG_DISCRI_CUT 0.80


//---------------------------
// To Be DONE
// -btagging weights: done as US group
//   does not treat SF for C ...
//   do not propage the uncertainties yet


class WeightFactory{


	public:
	  WeightFactory();
	  ~WeightFactory();

	//-- Initalialization

	void Initialize();
	void InitializeBtagSFTool(string CSVfileFullSim, string CSVfileFastSim, string BtagEffFullSimFilename, string BtagEffFastSimFilename);
	//many things are hard-coded in that function such as the WPoint etc ...
	void InitializeLeptonSFTool();
	//many files should be given as arguments ... to much. Find a solution ...
	//many things are hard-coded in that function such as the WPoint etc ...

	//-- Computation
	float GetBtagEff(const float& pt, const float& eta, const int& mcFlavour); // called by BtagWeighComputor
	void BtagWeighComputor (const vector<float>& jets_pt, const vector<float>& jets_eta, const vector<int>& jets_hadronFlavour, const vector<float>& jets_CSV);// should be called once per event
	//float GetLepEff();
	void LeptonWeightComputor(float lep1_pt, float lep1_eta, float lep1_pdgid, float lep2_pt, float lep2_eta, float lep2_pdgid, int nVetoLeptons, int nGoodLeptons, int NgenLeptons, vector<float> genLostLeptons_pt, vector<float> genLostLeptons_eta, vector<int> genLostLeptons_pdgid);
	double TopPTWeightComputor(float top_pt);
	
	//-- Accessors
	double GetGlobalW() {return Wglobal;}
	double GetMCW(){return Wmc;}
	double GetBtagW(){return Wbtag;}
	double GetLepW(){return Wlep;}
	double GetTopPtW(){return Wtop_pt;}

	//-- Mutators
	//-- should be called for each dataset
	void SetIsFastSim(bool isFS) {isFastSim = isFS;}
	void SetIsData(bool isD) {isData = isD;}

	private:
	  
	  //---- Global evnt info ----//
	  bool isFastSim;
	  bool isData;

	  //---- Global Weights ----//
	  double Wglobal;
	  double Wmc;
	  double Wbtag;
	  double Wlep;
	  double Wtop_pt; //top pt reweighing

	  //---- break down of the weights ----//


      	  double lepSF;
     	  double lepSF_Up;
     	  double lepSF_Dn;
      
      	  double vetoLepSF;
   	  double vetoLepSF_Up;
      	  double vetoLepSF_Dn;

    	  double lepSF_FS;
     	  double lepSF_FS_Up;
     	  double lepSF_FS_Dn;	

	  //----------------------------------------// 
	  //---- Btag weighting tools    	----//
	  //----------------------------------------// 

		//-- for fullsim --//
	  BTagCalibration* calib;
	  BTagCalibrationReader* reader_heavy;
	  BTagCalibrationReader* reader_heavy_UP;
	  BTagCalibrationReader* reader_heavy_DN;
	  BTagCalibrationReader* reader_c;
	  BTagCalibrationReader* reader_c_UP;
	  BTagCalibrationReader* reader_c_DN;
	  BTagCalibrationReader* reader_light;
	  BTagCalibrationReader* reader_light_UP;
	  BTagCalibrationReader* reader_light_DN;
		
		//-- fullsim plots
	  TH2D* h_btag_eff_b_fullSim;
	  TH2D* h_btag_eff_c_fullSim;
	  TH2D* h_btag_eff_udsg_fullSim;
	  
		//-- for fastsim --//
	  
	  BTagCalibration* calib_fastsim;
	  BTagCalibrationReader* reader_light_fastsim;
	  BTagCalibrationReader* reader_light_fastsim_UP;
	  BTagCalibrationReader* reader_light_fastsim_DN;
	  BTagCalibrationReader* reader_heavy_fastsim;
	  BTagCalibrationReader* reader_heavy_fastsim_UP;
	  BTagCalibrationReader* reader_heavy_fastsim_DN;

		//-- fastsim plots
	  TH2D* h_btag_eff_b_fastsim;
	  TH2D* h_btag_eff_c_fastsim;
	  TH2D* h_btag_eff_udsg_fastsim;;

	  //----------------------------------------// 
	  //---- Lepton  weighting tools 	----//
	  //----------------------------------------// 
	  // Fullsim Electron file
	  TFile *f_el_SF;
	  TFile *f_el_SF_tracking;

	  // Fullsim Muon files
	  TFile *f_mu_SF_id;
	  TFile *f_mu_SF_iso;
	  TFile *f_mu_SF_ip;
	  TFile *f_mu_SF_tracking;
	  TFile *f_mu_SF_veto_id;
	  TFile *f_mu_SF_veto_iso;
	  TFile *f_mu_SF_veto_ip;

	  // Fullsim/Fastsim Electron files
	  TFile *f_el_FS_ID;
	  TFile *f_el_FS_Iso;
	  TFile *f_el_veto_FS_ID;
	  TFile *f_el_veto_FS_Iso;

	  // Fullsim/Fastsim Muon files
	  TFile *f_mu_FS_ID;
	  TFile *f_mu_FS_Iso;
	  TFile *f_mu_FS_Ip;
	  TFile *f_mu_veto_FS_ID;
	  TFile *f_mu_veto_FS_Iso;
	  TFile *f_mu_veto_FS_Ip;

	  // Lepton MC reco efficiency for veto lep IDs
	  TFile *f_vetoLep_eff;
	 
	  // Final scale factor histograms for selected leptons
	  TH2D *h_el_SF = NULL;
	  TH2D *h_mu_SF = NULL;
	  
	  // Final scale factor histograms for veto leptons
	  TH2D *h_el_SF_veto = NULL;
	  TH2D *h_mu_SF_veto = NULL;

	  // Final scale factor histograms for tracking
	  TH2D *h_el_SF_tracking = NULL;
	  TH1D *h_mu_SF_tracking = NULL;

	  // Final scale factor histograms for fastim/fullsim for selected leptons
	  TH2D *h_el_FS = NULL;
	  TH2D *h_mu_FS = NULL;

	  // Final scale factor histograms for fastim/fullsim for veto leptons
	  TH2D *h_el_veto_FS = NULL;
	  TH2D *h_mu_veto_FS = NULL;

	  // Final scale factor histograms for lost leptons
	  TH2D *h_el_vetoLepEff = NULL;
	  TH2D *h_mu_vetoLepEff = NULL; 

};

#endif
