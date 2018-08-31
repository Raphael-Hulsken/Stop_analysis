#include "WeightFactory.h"

#ifndef PATH
  #error PATH need to be defined
#endif

WeightFactory::WeightFactory ()
{
	  //---- Global evnt info ----//
	  isFastSim = false;
	  isData = false;

	  //---- Global Weights ----//
	  Wglobal = 1.0;
	  Wmc = 1.0;
	  Wbtag = 1.0;
	  Wlep = 1.0;
	  Wtop_pt = 1.0;

          lepSF    = 1.0;
          lepSF_Up = 1.0;
          lepSF_Dn = 1.0;
      
          vetoLepSF    = 1.0;
          vetoLepSF_Up = 1.0;
          vetoLepSF_Dn = 1.0;


          lepSF_FS    = 1.0;
          lepSF_FS_Up = 1.0;
          lepSF_FS_Dn = 1.0;	
	  //----------------------------------------// 
	  //---- Btag weighting tools    	----//
	  //----------------------------------------// 

		//-- for fullsim --//
	  calib = NULL;
	  reader_heavy = NULL;
	  reader_heavy_UP = NULL;
	  reader_heavy_DN = NULL;
	  reader_c = NULL;
	  reader_c_UP = NULL;
	  reader_c_DN = NULL;
	  reader_light = NULL;
	  reader_light_UP = NULL;
	  reader_light_DN = NULL;
		
		//-- fullsim plots
	  h_btag_eff_b_fullSim = NULL;
	  h_btag_eff_c_fullSim = NULL;
	  h_btag_eff_udsg_fullSim = NULL;
	  
		//-- for fastsim --//
	  
	  calib_fastsim = NULL;
	  reader_light_fastsim = NULL;
	  reader_light_fastsim_UP = NULL;
	  reader_light_fastsim_DN = NULL;
	  reader_heavy_fastsim = NULL;
	  reader_heavy_fastsim_UP = NULL;
	  reader_heavy_fastsim_DN = NULL;

		//-- fastsim plots
	  h_btag_eff_b_fastsim = NULL;
	  h_btag_eff_c_fastsim = NULL;
	  h_btag_eff_udsg_fastsim = NULL;

	  //----------------------------------------// 
	  //---- Lepton  weighting tools 	----//
	  //----------------------------------------// 
	  // Fullsim Electron file
	  f_el_SF = NULL;
	  f_el_SF_tracking = NULL;

	  // Fullsim Muon files
	  f_mu_SF_id = NULL;
	  f_mu_SF_iso = NULL;
	  f_mu_SF_ip = NULL;
	  f_mu_SF_tracking = NULL;
	  f_mu_SF_veto_id = NULL;
	  f_mu_SF_veto_iso = NULL;
	  f_mu_SF_veto_ip = NULL;

	  // Fullsim/Fastsim Electron files
	  f_el_FS_ID = NULL;
	  f_el_FS_Iso = NULL;
	  f_el_veto_FS_ID = NULL;
	  f_el_veto_FS_Iso = NULL;

	  // Fullsim/Fastsim Muon files
	  f_mu_FS_ID = NULL;
	  f_mu_FS_Iso = NULL;
	  f_mu_FS_Ip = NULL;
	  f_mu_veto_FS_ID = NULL;
	  f_mu_veto_FS_Iso = NULL;
	  f_mu_veto_FS_Ip = NULL;

	  // Lepton MC reco efficiency for veto lep IDs
	  f_vetoLep_eff = NULL;
	 
	  // Final scale factor histograms for selected leptons
	  h_el_SF = NULL;
	  h_mu_SF = NULL;
	  
	  // Final scale factor histograms for veto leptons
	  h_el_SF_veto = NULL;
	  h_mu_SF_veto = NULL;

	  // Final scale factor histograms for tracking
	  h_el_SF_tracking = NULL;
	  h_mu_SF_tracking = NULL;

	  // Final scale factor histograms for fastim/fullsim for selected leptons
	  h_el_FS = NULL;
	  h_mu_FS = NULL;

	  // Final scale factor histograms for fastim/fullsim for veto leptons
	  h_el_veto_FS = NULL;
	  h_mu_veto_FS = NULL;

	  // Final scale factor histograms for lost leptons
	  h_el_vetoLepEff = NULL;
	  h_mu_vetoLepEff = NULL; 

}

WeightFactory::~WeightFactory ()
{
}

/*
void WeightFactory::Initialize ()
{
}
*/

//do we also need a file for Btag eff FastStim sample ?
void WeightFactory::InitializeBtagSFTool (string CSVfileFullSim, string CSVfileFastSim, string BtagEffFullSimFilename, string BtagEffFastSimFilename)
{
  //CSVfileFullSim = "btagsf/CSVv2_ichep_slimmed.csv";
  //CSVfileFastSim = "btagsf/CSV_13TEV_Combined_14_7_2016.csv"
  printBoxedMessage ("WeightFactory: initialization of BtagSFTool");
  printBoxedMessage (" > FullSim SF file: " + CSVfileFullSim);
  printBoxedMessage (" > FastSim SF file: " + CSVfileFastSim);
  printBoxedMessage (" > FullSim BtagEff file: " + BtagEffFullSimFilename);
  printBoxedMessage (" > FastSim BtagEff file: " + BtagEffFastSimFilename);

  //-- Loading SF --//
  //
  // fullSim
  calib = new BTagCalibration ("csvv2", CSVfileFullSim.c_str ());
  //- heavy (b) flavor
  //reader_heavy = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "central");	// central
  reader_heavy = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  reader_heavy->load(*calib,BTagEntry::FLAV_B,"comb");
  //reader_heavy_UP = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "up");	// sys up
  reader_heavy_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_heavy_UP->load(*calib,BTagEntry::FLAV_B,"comb");
  //reader_heavy_DN = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "down");	// sys down
  reader_heavy_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_heavy_DN->load(*calib,BTagEntry::FLAV_B,"comb");
  
  //- c flavor
  reader_c = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  reader_c->load(*calib,BTagEntry::FLAV_C,"comb");
  //reader_c_UP = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "up");	// sys up
  reader_c_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_c_UP->load(*calib,BTagEntry::FLAV_C,"comb");
  //reader_c_DN = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "down");	// sys down
  reader_c_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_c_DN->load(*calib,BTagEntry::FLAV_C,"comb");
  
  //- light flavor
  //reader_light = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "incl", "central");	// central
  reader_light = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  reader_light->load(*calib,BTagEntry::FLAV_UDSG,"incl");
  //reader_light_UP = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "incl", "up");	// sys up
  reader_light_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_light_UP->load(*calib,BTagEntry::FLAV_UDSG,"incl");
  //reader_light_DN = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "incl", "down");	// sys down
  reader_light_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_light_DN->load(*calib,BTagEntry::FLAV_UDSG,"incl");
  
  // fastSim
  calib_fastsim = new BTagCalibration("CSV", CSVfileFastSim.c_str()); 
  //reader_light_fastsim = new BTagCalibrationReader (calib_light_fastsim, BTagEntry::OP_MEDIUM, "light_fastsim", "central");	// central
  reader_light_fastsim = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  reader_light_fastsim->load(*calib_fastsim,BTagEntry::FLAV_UDSG,"fastsim");
  //reader_light_fastsim_UP = new BTagCalibrationReader (calib_light_fastsim, BTagEntry::OP_MEDIUM, "light_fastsim", "up");	// sys up
  reader_light_fastsim_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_light_fastsim_UP->load(*calib_fastsim,BTagEntry::FLAV_UDSG,"fastsim");
  //reader_light_fastsim_DN = new BTagCalibrationReader (calib_light_fastsim, BTagEntry::OP_MEDIUM, "light_fastsim", "down");	// sys down
  reader_light_fastsim_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_light_fastsim_DN->load(*calib_fastsim,BTagEntry::FLAV_UDSG,"fastsim");
  
  //reader_heavy_fastsim = new BTagCalibrationReader (calib_heavy_fastsim, BTagEntry::OP_MEDIUM, "heavy_fastsim", "central");	// central
  reader_heavy_fastsim = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  reader_heavy_fastsim->load(*calib_fastsim,BTagEntry::FLAV_B,"fastsim");
  //reader_heavy_fastsim_UP = new BTagCalibrationReader (calib_heavy_fastsim, BTagEntry::OP_MEDIUM, "heavy_fastsim", "up");	// sys up
  reader_heavy_fastsim_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_heavy_fastsim_UP->load(*calib_fastsim,BTagEntry::FLAV_B,"fastsim");
  //reader_heavy_fastsim_DN = new BTagCalibrationReader (calib_heavy_fastsim, BTagEntry::OP_MEDIUM, "heavy_fastsim", "down");	// sys down
  reader_heavy_fastsim_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_heavy_fastsim_DN->load(*calib_fastsim,BTagEntry::FLAV_B,"fastsim");

  //-- Loading efficiency --//
  //-- fullSim
  TFile *file_eff = TFile::Open (BtagEffFullSimFilename.c_str ());
  h_btag_eff_b_fullSim = (TH2D *) file_eff->Get ("h2_BTaggingEff_csv_med_Eff_b")->Clone ();
  h_btag_eff_c_fullSim = (TH2D *) file_eff->Get ("h2_BTaggingEff_csv_med_Eff_c")->Clone ();
  h_btag_eff_udsg_fullSim = (TH2D *) file_eff->Get ("h2_BTaggingEff_csv_med_Eff_udsg")->Clone ();
  //-- fastSim
  TFile *file_eff_FastSim = TFile::Open (BtagEffFastSimFilename.c_str ());
  h_btag_eff_b_fastsim = (TH2D *) file_eff_FastSim->Get ("h2_BTaggingEff_csv_med_Eff_b")->Clone ();
  h_btag_eff_c_fastsim = (TH2D *) file_eff_FastSim->Get ("h2_BTaggingEff_csv_med_Eff_c")->Clone ();
  h_btag_eff_udsg_fastsim = (TH2D *) file_eff_FastSim->Get ("h2_BTaggingEff_csv_med_Eff_udsg")->Clone ();

  if (!h_btag_eff_b_fullSim || !h_btag_eff_c_fullSim || !h_btag_eff_udsg_fullSim)
    std::cerr << " WeightFactory::InitializeBtagSFTool:\t fullSim efficiency plots were not retrieved .." << endl;
  if (!h_btag_eff_b_fastsim || !h_btag_eff_c_fastsim || !h_btag_eff_udsg_fastsim)
    std::cerr << " WeightFactory::InitializeBtagSFTool:\t fastSim efficiency plots were not retrieved .." << endl;

}

float WeightFactory::GetBtagEff (const float& pt, const float& eta, const int& mcFlavour)
{
  // only use pt bins up to 400 GeV for charm and udsg
  float pt_cutoff = std::max (20., std::min (399., double (pt)));
  TH2D *h (0);
  if (abs (mcFlavour) == 5) {
    if (isFastSim)
      h = h_btag_eff_b_fastsim;
    else
      h = h_btag_eff_b_fullSim;
    // use pt bins up to 600 GeV for b
    pt_cutoff = std::max (20., std::min (599., double (pt)));
  }
  else if (abs (mcFlavour) == 4) {
    if (isFastSim)
      h = h_btag_eff_c_fastsim;
    else
      h = h_btag_eff_c_fullSim;
  }
  else {
    if (isFastSim)
      h = h_btag_eff_udsg_fastsim;
    else
      h = h_btag_eff_udsg_fullSim;
  }

  int binx = h->GetXaxis ()->FindBin (pt_cutoff);
  int biny = h->GetYaxis ()->FindBin (fabs (eta));
  return h->GetBinContent (binx, biny);
}
void WeightFactory::BtagWeighComputor (const vector<float>& jets_pt, const vector<float>& jets_eta, const vector<int>& jets_hadronFlavour, const vector<float>& jets_CSV)
{
  double btagprob_data (1.), btagprob_mc (1.), btagprob_heavy_UP (1.), btagprob_heavy_DN (1.), btagprob_light_UP (1.), btagprob_light_DN (1.), btagprob_FS_UP (1.), btagprob_FS_DN (1.);

  for (size_t i = 0; i < jets_pt.size (); ++i) {
    // Check if the jet pass the medium btag cut
    if (jets_CSV[i] > BTAG_DISCRI_CUT) {
      if (!isData) {
	float eff = GetBtagEff (jets_pt[i], jets_eta[i], jets_hadronFlavour[i]);
	//if(eff==0) cout<<"b-tagg eff is null !!"<<endl;
	BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;

	if (abs (jets_hadronFlavour[i]) == 5)
	  flavor = BTagEntry::FLAV_B;
	else if (abs (jets_hadronFlavour[i]) == 4)
	  //flavor = BTagEntry::FLAV_C;
	  flavor = BTagEntry::FLAV_B;

	float pt_cutoff = std::max (30., std::min (669., double (jets_pt[i])));
	float eta_cutoff = std::min (2.39, fabs (double (jets_eta[i])));
	float weight_cent (1.), weight_UP (1.), weight_DN (1.), weight_FS_UP (1.), weight_FS_DN (1.);

	if (flavor == BTagEntry::FLAV_UDSG) {
	  weight_cent = reader_light->eval (flavor, eta_cutoff, pt_cutoff);
	  //if(weight_cent==0) cout<<"weight_cent = "<<weight_cent<<" UDSG "<<flavor<<" eta "<<jets_eta[i]<<" % "<< eta_cutoff<<" pt "<<jets_pt[i]<<" % "<< pt_cutoff<<endl;
	  weight_UP = reader_light_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_light_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	else if(flavor == BTagEntry::FLAV_B){
	  //if(weight_cent==0) cout<<"weight_cent = "<<weight_cent<<" BorC "<<flavor<<" eta "<<jets_eta[i]<<" % "<< eta_cutoff<<" pt "<<jets_pt[i]<<" % "<< pt_cutoff<<endl;
	  weight_cent = reader_heavy->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP = reader_heavy_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_heavy_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	else{
	  weight_cent = reader_c->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP = reader_c_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_c_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	if (isFastSim) {
	  if (flavor == BTagEntry::FLAV_UDSG) {
	    weight_FS_UP = weight_cent * reader_light_fastsim_UP->eval (flavor, eta_cutoff, pt_cutoff);
	    weight_FS_DN = weight_cent * reader_light_fastsim_DN->eval (flavor, eta_cutoff, pt_cutoff);
	    weight_cent *= reader_light_fastsim->eval (flavor, eta_cutoff, pt_cutoff);
	    //if(weight_cent==0) cout<<"FASTSIM effect "<<flavor<<endl;
	    weight_UP *= reader_light_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	    weight_DN *= reader_light_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	  }
	  else{// need to differentiate b and c
	    weight_FS_UP = weight_cent * reader_heavy_fastsim_UP->eval (flavor, eta_cutoff, pt_cutoff);
	    weight_FS_DN = weight_cent * reader_heavy_fastsim_DN->eval (flavor, eta_cutoff, pt_cutoff);
	    weight_cent *= reader_heavy_fastsim->eval (flavor, eta_cutoff, pt_cutoff);
	    //if(weight_cent==0) cout<<"FASTSIM effect "<<flavor<<endl;
	    weight_UP *= reader_heavy_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	    weight_DN *= reader_heavy_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	  }
	}
	//if(weight_cent==0 || eff==0)
	//  cout<<"weight_cent = "<<weight_cent<<" || eff = "<<eff<<endl;
	btagprob_data *= weight_cent * eff;
	btagprob_mc *= eff;
	if (flavor == BTagEntry::FLAV_UDSG) {
	  btagprob_light_UP *= weight_UP * eff;
	  btagprob_light_DN *= weight_DN * eff;
	  btagprob_heavy_UP *= weight_cent * eff;
	  btagprob_heavy_DN *= weight_cent * eff;
	}
	else {
	  btagprob_light_UP *= weight_cent * eff;
	  btagprob_light_DN *= weight_cent * eff;
	  btagprob_heavy_UP *= weight_UP * eff;
	  btagprob_heavy_DN *= weight_DN * eff;
	}
	if (isFastSim) {
	  btagprob_FS_UP *= weight_FS_UP * eff;
	  btagprob_FS_DN *= weight_FS_DN * eff;
	}
      }
    }
    else {			// event if it fails med btag -- it's  needed for SF event weights
      if (!isData) {
	float eff = GetBtagEff (jets_eta[i], jets_eta[i], jets_hadronFlavour[i]);
	BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;
	if (abs (jets_hadronFlavour[i]) == 5)
	  flavor = BTagEntry::FLAV_B;
	else if (abs (jets_hadronFlavour[i]) == 4)
	  //flavor = BTagEntry::FLAV_C;
	  // wrong .. but just for test ...
	  flavor = BTagEntry::FLAV_B;
	float pt_cutoff = std::max (30., std::min (669., double (jets_pt[i])));
	float eta_cutoff = std::min (2.39, fabs (double (jets_eta[i])));
	float weight_cent (1.), weight_UP (1.), weight_DN (1.), weight_FS_UP (1.), weight_FS_DN (1.);
	if (flavor == BTagEntry::FLAV_UDSG) {
	  weight_cent = reader_light->eval (flavor, eta_cutoff, pt_cutoff);
	  if(weight_cent==0) cout<<"Weight_cent = "<<weight_cent<<" UDSG "<<flavor<<" eta "<<jets_eta[i]<<" % "<< eta_cutoff<<" pt "<<jets_pt[i]<<" % "<< pt_cutoff<<endl;
	  weight_UP = reader_light_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_light_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	else if(flavor == BTagEntry::FLAV_B) {
	  weight_cent = reader_heavy->eval (flavor, eta_cutoff, pt_cutoff);
	  if(weight_cent==0) cout<<"else weight_cent = "<<weight_cent<<" UDSG "<<flavor<<" eta "<<jets_eta[i]<<" % "<< eta_cutoff<<" pt "<<jets_pt[i]<<" % "<< pt_cutoff<<endl;
	  weight_UP = reader_heavy_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_heavy_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	else{
	  weight_cent = reader_c->eval (flavor, eta_cutoff, pt_cutoff);
	  if(weight_cent==0) cout<<"else weight_cent = "<<weight_cent<<" UDSG "<<flavor<<" eta "<<jets_eta[i]<<" % "<< eta_cutoff<<" pt "<<jets_pt[i]<<" % "<< pt_cutoff<<endl;
	  weight_UP = reader_c_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_c_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	if (isFastSim) {
	  if (flavor == BTagEntry::FLAV_UDSG) {
	    weight_FS_UP = weight_cent * reader_light_fastsim_UP->eval (flavor, eta_cutoff, pt_cutoff);	//this is pure light_fastsimSF
	    weight_FS_DN = weight_cent * reader_light_fastsim_DN->eval (flavor, eta_cutoff, pt_cutoff);	//this is pure light_fastsimSF
	    weight_cent *= reader_light_fastsim->eval (flavor, eta_cutoff, pt_cutoff);
	    if(weight_cent==0) cout<<"FASTSIM effect"<<endl;
	    weight_UP *= reader_light_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	    weight_DN *= reader_light_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	  }
	  else{// need to differentiate btw b and c
	    weight_FS_UP = weight_cent * reader_heavy_fastsim_UP->eval (flavor, eta_cutoff, pt_cutoff);	//this is pure heavy_fastsimSF
	    weight_FS_DN = weight_cent * reader_heavy_fastsim_DN->eval (flavor, eta_cutoff, pt_cutoff);	//this is pure heavy_fastsimSF
	    weight_cent *= reader_heavy_fastsim->eval (flavor, eta_cutoff, pt_cutoff);
	    if(weight_cent==0) cout<<"FASTSIM effect"<<endl;
	    weight_UP *= reader_heavy_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	    weight_DN *= reader_heavy_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	  }
	}

	btagprob_data *= (1. - weight_cent * eff);
	btagprob_mc *= (1. - eff);
	if (flavor == BTagEntry::FLAV_UDSG) {
	  btagprob_light_UP *= (1. - weight_UP * eff);
	  btagprob_light_DN *= (1. - weight_DN * eff);
	  btagprob_heavy_UP *= (1. - weight_cent * eff);
	  btagprob_heavy_DN *= (1. - weight_cent * eff);
	}
	else {
	  btagprob_light_UP *= (1. - weight_cent * eff);
	  btagprob_light_DN *= (1. - weight_cent * eff);
	  btagprob_heavy_UP *= (1. - weight_UP * eff);
	  btagprob_heavy_DN *= (1. - weight_DN * eff);
	}
	if (isFastSim) {
	  btagprob_FS_UP *= (1. - weight_FS_UP * eff);
	  btagprob_FS_DN *= (1. - weight_FS_DN * eff);
	}
      }
    }
  }
  //Compute a weight as the ratio of data over mc
  //if(btagprob_data==0) cout<<"btagprob_data = "<<btagprob_data<<"/"<<"btagprob_mc = "<<btagprob_mc<<endl;
  Wbtag = btagprob_data / btagprob_mc;

}

void WeightFactory::InitializeLeptonSFTool(){
    // Electron file
    string path = PATH;
    cout<<path<<endl;
    f_el_SF          = new TFile((string(PATH)+string("/files/lepsf/scaleFactors.root")).c_str(), "read");
    f_el_SF_tracking = new TFile((string(PATH)+string("/files/lepsf/egammaEffi.txt_SF2D.root")).c_str(), "read");
    
    // Muon files
    f_mu_SF_id       = new TFile((string(PATH)+string("/files/lepsf/TnP_MuonID_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root")).c_str(), "read"); // double unc
    f_mu_SF_iso      = new TFile((string(PATH)+string("/files/lepsf/TnP_MuonID_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root")).c_str(), "read"); // double unc
    f_mu_SF_ip       = new TFile((string(PATH)+string("/files/lepsf/TnP_MuonID_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root")).c_str(), "read"); // double unc
    f_mu_SF_tracking = new TFile((string(PATH)+string("/files/lepsf/muons_tracking_sf.root")).c_str(), "read"); 
    
    f_mu_SF_veto_id  = new TFile((string(PATH)+string("/files/lepsf/TnP_MuonID_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta.root")).c_str(), "read");
    f_mu_SF_veto_iso = new TFile((string(PATH)+string("/files/lepsf/TnP_MuonID_NUM_MiniIsoTight_DENOM_LooseID_VAR_map_pt_eta.root")).c_str(), "read");
    f_mu_SF_veto_ip  = new TFile((string(PATH)+string("/files/lepsf/TnP_MuonID_NUM_MediumIP2D_DENOM_LooseID_VAR_map_pt_eta.root")).c_str(), "read"); // double unc for this
    
    // Fastsim/Fullsim el files
    f_el_FS_ID  = new TFile((string(PATH)+string("/files/lepsf/sf_el_mediumCB.root")).c_str(), "read");
    f_el_FS_Iso = new TFile((string(PATH)+string("/files/lepsf/sf_el_mini01.root")).c_str(), "read");
    
    f_el_veto_FS_ID  = new TFile((string(PATH)+string("/files/lepsf/sf_el_vetoCB.root")).c_str(), "read");
    f_el_veto_FS_Iso = new TFile((string(PATH)+string("/files/lepsf/sf_el_mini02.root")).c_str(), "read"); 
    
    // Fastsim/Fullsim mu files
    f_mu_FS_ID  = new TFile((string(PATH)+string("/files/lepsf/sf_mu_medium.root")).c_str(), "read"); // double unc for this
    f_mu_FS_Iso = new TFile((string(PATH)+string("/files/lepsf/sf_mu_mediumID_mini02.root")).c_str(), "read"); // double unc for this
    f_mu_FS_Ip  = new TFile((string(PATH)+string("/files/lepsf/sf_mu_tightIP2D.root")).c_str(), "read"); // double unc for this
    
    f_mu_veto_FS_ID  = new TFile((string(PATH)+string("/files/lepsf/sf_mu_loose.root")).c_str(), "read");
    f_mu_veto_FS_Iso = new TFile((string(PATH)+string("/files/lepsf/sf_mu_looseID_mini02.root")).c_str(), "read");
    f_mu_veto_FS_Ip  = new TFile((string(PATH)+string("/files/lepsf/sf_mu_looseIP2D.root")).c_str(), "read"); // double unc for this

    // Veto lepton reco efficiency files
    f_vetoLep_eff = new TFile((string(PATH)+string("/files/lepsf/lepeff__ttbar_diLept_madgraph_pythia8_ext1_25ns__80x.root")).c_str(), "read");


    // Grab selected el histos
    TH2D *h_el_SF_id_temp       = (TH2D*)f_el_SF->Get("GsfElectronToMedium");
    TH2D *h_el_SF_iso_temp      = (TH2D*)f_el_SF->Get("MVAVLooseElectronToMini");
    TH2D *h_el_SF_tracking_temp = (TH2D*)f_el_SF_tracking->Get("EGamma_SF2D");

    // Grab veto el histos
    TH2D *h_el_SF_veto_id_temp  = (TH2D*)f_el_SF->Get("GsfElectronToVeto");
    TH2D *h_el_SF_veto_iso_temp = (TH2D*)f_el_SF->Get("MVAVLooseElectronToMini2");


    // Grab selected mu histos
    TH2D *h_mu_SF_id_temp  = (TH2D*)f_mu_SF_id->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
    TH2D *h_mu_SF_iso_temp = (TH2D*)f_mu_SF_iso->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass");
    TH2D *h_mu_SF_ip_temp  = (TH2D*)f_mu_SF_ip->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass");
    TGraphAsymmErrors *h_mu_SF_tracking_temp  = (TGraphAsymmErrors*)f_mu_SF_tracking->Get("ratio_eta");

    // Grab veto mu histos
    TH2D *h_mu_SF_veto_id_temp  = (TH2D*)f_mu_SF_veto_id->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
    TH2D *h_mu_SF_veto_iso_temp = (TH2D*)f_mu_SF_veto_iso->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_PF_pass");
    TH2D *h_mu_SF_veto_ip_temp  = (TH2D*)f_mu_SF_veto_ip->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_PF_pass");
    

    // Grab fastsim/fullsim selected el histos
    TH2D *h_el_FS_ID_temp  = (TH2D*)f_el_FS_ID->Get("histo2D");
    TH2D *h_el_FS_Iso_temp = (TH2D*)f_el_FS_Iso->Get("histo2D");

    // Grab fastsim/fullsim veto el histos
    TH2D *h_el_veto_FS_ID_temp  = (TH2D*)f_el_veto_FS_ID->Get("histo2D");
    TH2D *h_el_veto_FS_Iso_temp = (TH2D*)f_el_veto_FS_Iso->Get("histo2D");
    
    // Grab fastsim/fullsim selected mu histos
    TH2D *h_mu_FS_ID_temp  = (TH2D*)f_mu_FS_ID->Get("histo2D");
    TH2D *h_mu_FS_Iso_temp = (TH2D*)f_mu_FS_Iso->Get("histo2D");
    TH2D *h_mu_FS_Ip_temp  = (TH2D*)f_mu_FS_Ip->Get("histo2D");
    
    // Grab fastsim/fullsim veto mu histos
    TH2D *h_mu_veto_FS_ID_temp  = (TH2D*)f_mu_veto_FS_ID->Get("histo2D");
    TH2D *h_mu_veto_FS_Iso_temp = (TH2D*)f_mu_veto_FS_Iso->Get("histo2D");
    TH2D *h_mu_veto_FS_Ip_temp  = (TH2D*)f_mu_veto_FS_Ip->Get("histo2D");

    
    // Grab mc eff for veto lepton (for lost lepto SFs) histos
    TH2D *h_el_vetoLepEff_temp = (TH2D*)f_vetoLep_eff->Get("h2_lepEff_vetoSel_Eff_el");
    TH2D *h_mu_vetoLepEff_temp = (TH2D*)f_vetoLep_eff->Get("h2_lepEff_vetoSel_Eff_mu");
   

    // Get final fullsim, selected el, sfs
    TH2D *h_el_SF_id  = (TH2D*)h_el_SF_id_temp->Clone("h_el_SF_id");
    TH2D *h_el_SF_iso = (TH2D*)h_el_SF_iso_temp->Clone("h_el_SF_iso");
    h_el_SF = (TH2D*)h_el_SF_id->Clone("h_el_SF");
    h_el_SF->Multiply(h_el_SF_iso);
    
    // Get final fullsim, selected el, tracking sfs
    h_el_SF_tracking = (TH2D*)h_el_SF_tracking_temp->Clone("h_el_SF_iso");
    
    // Get final fullsim, veto el, sfs
    TH2D *h_el_SF_veto_id  = (TH2D*)h_el_SF_veto_id_temp->Clone("h_el_SF_veto_id");
    TH2D *h_el_SF_veto_iso = (TH2D*)h_el_SF_veto_iso_temp->Clone("h_el_SF_veto_iso");
    h_el_SF_veto = (TH2D*)h_el_SF_veto_id->Clone("h_el_SF_veto");
    h_el_SF_veto->Multiply(h_el_SF_veto_iso);
    
    
    // Get final fullsim, selected mu, sfs
    TH2D *h_mu_SF_id  = (TH2D*)h_mu_SF_id_temp->Clone("h_mu_SF_id");
    TH2D *h_mu_SF_iso = (TH2D*)h_mu_SF_iso_temp->Clone("h_mu_SF_iso");
    TH2D *h_mu_SF_ip  = (TH2D*)h_mu_SF_ip_temp->Clone("h_mu_SF_ip");
    // Double unc. on selected muon id sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_SF_id->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_SF_id->GetNbinsY(); y++){
	h_mu_SF_id->SetBinError(x,y,h_mu_SF_id->GetBinError(x,y)*2.0);
      }
    }
    // Double unc. on selected muon iso sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_SF_iso->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_SF_iso->GetNbinsY(); y++){
	h_mu_SF_iso->SetBinError(x,y,h_mu_SF_iso->GetBinError(x,y)*2.0);
      }
    }
    // Double unc. on selected muon ip sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_SF_ip->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_SF_ip->GetNbinsY(); y++){
	h_mu_SF_ip->SetBinError(x,y,h_mu_SF_ip->GetBinError(x,y)*2.0);
      }
    }
    h_mu_SF = (TH2D*)h_mu_SF_id->Clone("h_mu_SF");
    h_mu_SF->Multiply(h_mu_SF_iso);
    h_mu_SF->Multiply(h_mu_SF_ip);
    

    // Get final fullsim, selected muon, tracking sfs, convert TGraphErrors
    int nX = h_mu_SF_tracking_temp->GetN();
    Double_t *x_val = h_mu_SF_tracking_temp->GetX();
    Double_t *y_val = h_mu_SF_tracking_temp->GetY();
    Double_t *y_err_up = h_mu_SF_tracking_temp->GetEYhigh();
    Double_t *y_err_low = h_mu_SF_tracking_temp->GetEYhigh();
    h_mu_SF_tracking = new TH1D("h_mu_SF_tracking","h_mu_SF_tracking",nX-1,x_val);
    for(int i=0; i<nX; i++){
      h_mu_SF_tracking->SetBinContent(i+1,y_val[i]);
      h_mu_SF_tracking->SetBinError(i+1,std::max(y_err_up[i],y_err_low[i]));
    }
    

    // Get final fullsim, veto mu, sfs
    TH2D *h_mu_SF_veto_id  = (TH2D*)h_mu_SF_veto_id_temp->Clone("h_mu_SF_veto_id");
    TH2D *h_mu_SF_veto_iso = (TH2D*)h_mu_SF_veto_iso_temp->Clone("h_mu_SF_veto_iso");
    TH2D *h_mu_SF_veto_ip  = (TH2D*)h_mu_SF_veto_ip_temp->Clone("h_mu_SF_veto_ip");
    // Double unc. on veto muon ip sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_SF_veto_ip->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_SF_veto_ip->GetNbinsY(); y++){
	h_mu_SF_ip->SetBinError(x,y,h_mu_SF_veto_ip->GetBinError(x,y)*2.0);
      }
    }
    h_mu_SF_veto = (TH2D*)h_mu_SF_veto_id->Clone("h_mu_SF_veto");
    h_mu_SF_veto->Multiply(h_mu_SF_veto_iso);
    h_mu_SF_veto->Multiply(h_mu_SF_veto_ip);
    

    // Get final fullsim/fastsim, selected el, sfs
    TH2D* h_el_FS_ID  = (TH2D*)h_el_FS_ID_temp->Clone("h_el_FS_ID");
    TH2D* h_el_FS_Iso = (TH2D*)h_el_FS_Iso_temp->Clone("h_el_FS_Iso");
    h_el_FS = (TH2D*)h_el_FS_ID->Clone("h_el_FS");
    h_el_FS->Multiply(h_el_FS_Iso);
    
    // Get final fullsim/fastsim, veto el, sfs
    TH2D* h_el_veto_FS_ID  = (TH2D*)h_el_veto_FS_ID_temp->Clone("h_el_veto_FS_ID");
    TH2D* h_el_veto_FS_Iso = (TH2D*)h_el_veto_FS_Iso_temp->Clone("h_el_veto_FS_Iso");
    h_el_veto_FS = (TH2D*)h_el_veto_FS_ID->Clone("h_el_FS");
    h_el_veto_FS->Multiply(h_el_veto_FS_Iso);
    
    // Get final fullsim/fastsim, selected mu, sfs
    TH2D* h_mu_FS_ID  = (TH2D*)h_mu_FS_ID_temp->Clone("h_mu_FS_ID");
    TH2D* h_mu_FS_Iso = (TH2D*)h_mu_FS_Iso_temp->Clone("h_mu_FS_Iso");
    TH2D* h_mu_FS_Ip  = (TH2D*)h_mu_FS_Ip_temp->Clone("h_mu_FS_Ip");
    // Double unc. on selected muon FS id sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_FS_ID->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_FS_ID->GetNbinsY(); y++){
	h_mu_FS_ID->SetBinError(x,y,h_mu_FS_ID->GetBinError(x,y)*2.0);
      }
    }
    // Double unc. on selected muon FS iso sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_FS_Iso->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_FS_Iso->GetNbinsY(); y++){
	h_mu_FS_Iso->SetBinError(x,y,h_mu_FS_Iso->GetBinError(x,y)*2.0);
      }
    }
    // Double unc. on selected muon FS ip sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_FS_Ip->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_FS_Ip->GetNbinsY(); y++){
	h_mu_FS_Ip->SetBinError(x,y,h_mu_FS_Ip->GetBinError(x,y)*2.0);
      }
    }
    h_mu_FS = (TH2D*)h_mu_FS_ID->Clone("h_mu_FS");
    h_mu_FS->Multiply(h_mu_FS_Iso);
    h_mu_FS->Multiply(h_mu_FS_Ip);
    
    // Get final fullsim/fastsim, veto mu, sfs
    TH2D* h_mu_veto_FS_ID  = (TH2D*)h_mu_veto_FS_ID_temp->Clone("h_mu_veto_FS_ID");
    TH2D* h_mu_veto_FS_Iso = (TH2D*)h_mu_veto_FS_Iso_temp->Clone("h_mu_veto_FS_Iso");
    TH2D* h_mu_veto_FS_Ip  = (TH2D*)h_mu_veto_FS_Ip_temp->Clone("h_mu_veto_FS_Ip");
    // Double unc. on selected muon FS ip sfs, since not for our exact wp
    for(int x=1; x<=(int)h_mu_veto_FS_Ip->GetNbinsX(); x++){
      for(int y=1; y<=(int)h_mu_veto_FS_Ip->GetNbinsY(); y++){
	h_mu_veto_FS_Ip->SetBinError(x,y,h_mu_veto_FS_Ip->GetBinError(x,y)*2.0);
      }
    }
    h_mu_veto_FS = (TH2D*)h_mu_veto_FS_ID->Clone("h_mu_veto_FS");
    h_mu_veto_FS->Multiply(h_mu_veto_FS_Iso);
    h_mu_veto_FS->Multiply(h_mu_veto_FS_Ip);
    


    // Lepton efficiencies for Lost Leptons
    h_el_vetoLepEff = (TH2D*)h_el_vetoLepEff_temp->Clone("h_el_vetoLepEff");
    h_mu_vetoLepEff = (TH2D*)h_mu_vetoLepEff_temp->Clone("h_mu_vetoLepEff");

}


void WeightFactory::LeptonWeightComputor(float lep1_pt, float lep1_eta, float lep1_pdgid, float lep2_pt, float lep2_eta, float lep2_pdgid, int nVetoLeptons, int nGoodLeptons, int NgenLeptons, vector<float> genLostLeptons_pt, vector<float> genLostLeptons_eta, vector<int> genLostLeptons_pdgid){
      // Lepton SFs
      float lepSF_pt_cutoff = 99.999;
      float lepSF_pt_min    = 10.001;
      float lepSF_FS_pt_cutoff = 199.999;
      float lepSF_FS_pt_min    = 10.001;
      
      
      // Lep1 SF
      if( !isData && nVetoLeptons>0){
	
	if(abs(lep1_pdgid) == 11){
	  int binX = h_el_SF->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep1_pt), lepSF_pt_min ) );
	  int binY = h_el_SF->GetYaxis()->FindBin( fabs(lep1_eta) );
	  lepSF    = h_el_SF->GetBinContent( binX, binY );
	  lepSF_Up = lepSF + h_el_SF->GetBinError( binX, binY );
	  lepSF_Dn = lepSF - h_el_SF->GetBinError( binX, binY );

	  binX = h_el_SF_tracking->GetXaxis()->FindBin( std::max( std::min(2.39,(double)lep1_eta),-2.39) );
	  binY = h_el_SF_tracking->GetYaxis()->FindBin( std::max( std::min(199.0,(double)lep1_pt), 21.0 ) );
	  lepSF *= h_el_SF_tracking->GetBinContent( binX, binY );
	  lepSF_Up *= ( h_el_SF_tracking->GetBinContent(binX,binY) + h_el_SF_tracking->GetBinError(binX,binY) );
	  lepSF_Dn *= ( h_el_SF_tracking->GetBinContent(binX,binY) - h_el_SF_tracking->GetBinError(binX,binY) );

	  if(isFastSim){
	    int bin_FS  = h_el_FS->FindBin( std::max( std::min(lepSF_FS_pt_cutoff, (float)lep1_pt), lepSF_FS_pt_min ), fabs(lep1_eta) );
	    lepSF_FS    = h_el_FS->GetBinContent(bin_FS);
	    lepSF_FS_Up = lepSF_FS + h_el_FS->GetBinError(bin_FS);
	    lepSF_FS_Dn = lepSF_FS + h_el_FS->GetBinError(bin_FS);
	  }
	  
	}
	
	if(abs(lep1_pdgid) == 13){
	  int binX = h_mu_SF->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep1_pt), lepSF_pt_min ) );
	  int binY = h_mu_SF->GetYaxis()->FindBin( fabs(lep1_eta) );
	  lepSF    = h_mu_SF->GetBinContent( binX, binY );
	  lepSF_Up = lepSF + h_mu_SF->GetBinError( binX, binY );
	  lepSF_Dn = lepSF - h_mu_SF->GetBinError( binX, binY );

	  binX = h_mu_SF_tracking->GetXaxis()->FindBin( std::max(-2.2,(double)lep1_eta) );
	  lepSF *= h_mu_SF_tracking->GetBinContent( binX );
	  lepSF_Up *= ( h_mu_SF_tracking->GetBinContent(binX) + h_mu_SF_tracking->GetBinError(binX) );
	  lepSF_Dn *= ( h_mu_SF_tracking->GetBinContent(binX) - h_mu_SF_tracking->GetBinError(binX) );
							
	  if(isFastSim){
	    int bin_FS  = h_mu_FS->FindBin( std::max( std::min(lepSF_FS_pt_cutoff, (float)lep1_pt), lepSF_FS_pt_min ), fabs(lep1_eta) );
	    lepSF_FS    = h_mu_FS->GetBinContent(bin_FS);
	    lepSF_FS_Up = lepSF_FS + h_mu_FS->GetBinError(bin_FS);
	    lepSF_FS_Dn = lepSF_FS + h_mu_FS->GetBinError(bin_FS);
	  }

	}
		    
      }
      
      // Lep2 SF
      if(!isData && nVetoLeptons>1 ){

	if(abs(lep2_pdgid) == 11){

	  if(nGoodLeptons>1){
	    int binX = h_el_SF->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep2_pt), lepSF_pt_min ) );
	    int binY = h_el_SF->GetYaxis()->FindBin( fabs(lep2_eta) );
	    lepSF    *= h_el_SF->GetBinContent( binX, binY );
	    lepSF_Up *= ( lepSF + h_el_SF->GetBinError( binX, binY ) );
	    lepSF_Dn *= ( lepSF - h_el_SF->GetBinError( binX, binY ) );
	    
	    if(isFastSim){
	      int bin_FS  = h_el_FS->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep2_pt), lepSF_pt_min ), fabs(lep2_eta) );
	      lepSF_FS    *= h_el_FS->GetBinContent(bin_FS);
	      lepSF_FS_Up *= (lepSF_FS + h_el_FS->GetBinError(bin_FS));
	      lepSF_FS_Dn *= (lepSF_FS + h_el_FS->GetBinError(bin_FS));
	    }
	  } // end if 2 good electrons
	  else{
	    int binX = h_el_SF_veto->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep2_pt), lepSF_pt_min ) );
	    int binY = h_el_SF_veto->GetYaxis()->FindBin( fabs(lep2_eta) );
	    lepSF    *= h_el_SF_veto->GetBinContent( binX, binY );
	    lepSF_Up *= ( lepSF + h_el_SF_veto->GetBinError( binX, binY ) );
	    lepSF_Dn *= ( lepSF - h_el_SF_veto->GetBinError( binX, binY ) );

	    if(isFastSim){
	      int bin_FS  = h_el_veto_FS->FindBin( std::max( std::min(lepSF_FS_pt_cutoff, (float)lep2_pt), lepSF_FS_pt_min ), fabs(lep2_eta) );
	      lepSF_FS    *= h_el_veto_FS->GetBinContent(bin_FS);
	      lepSF_FS_Up *= (lepSF_FS + h_el_veto_FS->GetBinError(bin_FS));
	      lepSF_FS_Dn *= (lepSF_FS + h_el_veto_FS->GetBinError(bin_FS));
	    }
	  }

	  int binX = h_el_SF_tracking->GetXaxis()->FindBin( std::max( std::min(2.39, (double)lep2_eta), -2.39) );
	  int binY = h_el_SF_tracking->GetYaxis()->FindBin( std::max( std::min(199.0, (double)lep2_pt), 21.0) );
	  lepSF *= h_el_SF_tracking->GetBinContent( binX, binY );
	  lepSF_Up *= ( h_el_SF_tracking->GetBinContent(binX,binY) + h_el_SF_tracking->GetBinError(binX,binY) );
	  lepSF_Dn *= ( h_el_SF_tracking->GetBinContent(binX,binY) - h_el_SF_tracking->GetBinError(binX,binY) );
	  
	} // end if 2nd lep if el
	

	if(abs(lep2_pdgid) == 13){

	  if(nGoodLeptons>1){
	    int binX = h_mu_SF->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep2_pt), lepSF_pt_min ) );
	    int binY = h_mu_SF->GetYaxis()->FindBin( fabs(lep2_eta) );
	    lepSF    *= h_mu_SF->GetBinContent( binX, binY );
	    lepSF_Up *= ( lepSF + h_mu_SF->GetBinError( binX, binY ) );
	    lepSF_Dn *= ( lepSF - h_mu_SF->GetBinError( binX, binY ) );
	  	
	    if(isFastSim){
	      int bin_FS  = h_mu_FS->FindBin( std::max( std::min(lepSF_FS_pt_cutoff, (float)lep2_pt), lepSF_FS_pt_min ), fabs(lep2_eta) );
	      lepSF_FS    *= h_mu_FS->GetBinContent(bin_FS);
	      lepSF_FS_Up *= lepSF_FS + h_mu_FS->GetBinError(bin_FS);
	      lepSF_FS_Dn *= lepSF_FS + h_mu_FS->GetBinError(bin_FS);
	    }
	  } // end if 2 good leptons
	  
	  else{
	    int binX = h_mu_SF_veto->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep2_pt), lepSF_pt_min ) );
	    int binY = h_mu_SF_veto->GetYaxis()->FindBin( fabs(lep2_eta) );
	    lepSF    *= h_mu_SF_veto->GetBinContent( binX, binY );
	    lepSF_Up *= ( lepSF + h_mu_SF_veto->GetBinError( binX, binY ) );
	    lepSF_Dn *= ( lepSF - h_mu_SF_veto->GetBinError( binX, binY ) );
	    
	    if(isFastSim){
	      int bin_FS  = h_mu_veto_FS->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)lep2_pt), lepSF_pt_min ), fabs(lep2_eta) );
	      lepSF_FS    *= h_mu_veto_FS->GetBinContent(bin_FS);
	      lepSF_FS_Up *= lepSF_FS + h_mu_veto_FS->GetBinError(bin_FS);
	      lepSF_FS_Dn *= lepSF_FS + h_mu_veto_FS->GetBinError(bin_FS);
	    }
	  }
	} // end if 2nd lep is mu

	int binX = h_mu_SF_tracking->GetXaxis()->FindBin( std::max(-2.23,(double)lep2_eta) );
	lepSF *= h_mu_SF_tracking->GetBinContent( binX );
	lepSF_Up *= ( h_mu_SF_tracking->GetBinContent(binX) + h_mu_SF_tracking->GetBinError(binX) );
	lepSF_Dn *= ( h_mu_SF_tracking->GetBinContent(binX) - h_mu_SF_tracking->GetBinError(binX) );

      } // end if 2nd lepton reco lepton exists

      
      // If only 1 reco lepton, and is2lep event, then find lost gen lepton
      // Not yet implemented ... need to be ported in pyROOF code ...
      ///*
      if( !isData && nVetoLeptons==1 && NgenLeptons ){

 	for(int i=0;i<NgenLeptons;i++){
	//this part commented was ported in PyROOT
	//for(int iGen=0; iGen<(int)gen_leps.p4.size(); iGen++){
	  /*
	  if( abs(gen_leps.id.at(iGen))!=11 && abs(gen_leps.id.at(iGen))!=13 ) continue;
	  if( !gen_leps.fromHardProcessFinalState.at(iGen) ) continue;
	  if( !gen_leps.isLastCopy.at(iGen) ) continue;
	  if( ROOT::Math::VectorUtil::DeltaR(gen_leps.p4.at(iGen), lep1.p4)<matched_dr ) continue;
	  if( gen_leps.p4.at(iGen).Pt()<5 || fabs(gen_leps.p4.at(iGen).Eta())>2.4 ) continue;
	  */

	  TH2D *h_vetoLep_eff = NULL;
	  if( abs(genLostLeptons_pdgid[i])==11 ) h_vetoLep_eff = h_el_vetoLepEff;
	  if( abs(genLostLeptons_pdgid[i])==13 ) h_vetoLep_eff = h_mu_vetoLepEff;
	  
	  int binX_eff = h_vetoLep_eff->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)genLostLeptons_pt[i]), lepSF_pt_min ) );
	  int binY_eff = h_vetoLep_eff->GetYaxis()->FindBin( fabs(genLostLeptons_eta[i]) );
	  double vetoEff = h_vetoLep_eff->GetBinContent( binX_eff, binY_eff );
	  
	  TH2D *h_lep_sf = NULL;
	  if( abs(genLostLeptons_pdgid[i])==11 ) h_lep_sf = h_el_SF_veto;
	  if( abs(genLostLeptons_pdgid[i])==13 ) h_lep_sf = h_mu_SF_veto;

	  int binX_sf = h_lep_sf->GetXaxis()->FindBin( std::max( std::min(lepSF_pt_cutoff, (float)genLostLeptons_pt[i]), lepSF_pt_min ) );
	  int binY_sf = h_lep_sf->GetYaxis()->FindBin( fabs(genLostLeptons_eta[i]) );
	  
	  double vetoLepSF_temp    = h_lep_sf->GetBinContent( binX_sf, binY_sf );
	  double vetoLepSF_temp_Up = vetoLepSF_temp + h_lep_sf->GetBinError( binX_sf, binY_sf );
	  double vetoLepSF_temp_Dn = vetoLepSF_temp - h_lep_sf->GetBinError( binX_sf, binY_sf );
	  
	  
	  if( vetoEff==1.0 ){
	    vetoLepSF    = 1.0;
	    vetoLepSF_Up = 1.0;
	    vetoLepSF_Dn = 1.0;
	  }
	  else{
	    vetoLepSF    = ( 1-(vetoEff*vetoLepSF_temp) )/( 1-vetoEff );
	    vetoLepSF_Up = ( 1-(vetoEff*vetoLepSF_temp_Up) )/( 1-vetoEff );
	    vetoLepSF_Dn = ( 1-(vetoEff*vetoLepSF_temp_Dn) )/( 1-vetoEff );	    
	  }
	  
	  break; // break after finding 2nd hard gen lepton

	} // end loop over gen leptons

      } // end if finding gen lost lepton for vetoEff SF
     //	*/

      Wlep = lepSF;
}

double WeightFactory::TopPTWeightComputor(float top_pt){
	double a = 0.156;
	double b = -0.00137;
	Wtop_pt =  exp(a+b*top_pt);
	return Wtop_pt;
}

