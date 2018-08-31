#include "Reader_CommonFormat.h"

babyEvent myEvent;


//------------------------------------
// Struct use to categorize the event
//------------------------------------


struct EventCategory{
    public:
	// One lepton events (CR or CR1l) - includes the vetoes
	bool OneLep;
	bool OneLep_e;
	bool OneLep_mu;

	// Two lepton events (CR-2l)
	bool TwoLep;
	bool TwoLep_ee;
	bool TwoLep_emu;
	bool TwoLep_mumu;

	// At least one b-tagged jets
	bool OneBtag;

        //Trigger result
        bool trigged;

	//Baseline: >=2 jets, MET, MT, DPhi
	bool PassBaselineCuts;
	
	// 2 jets bins (include mod_top>6.4
	bool TwoJets_MET_250;
	bool TwoJets_MET_350;
	bool TwoJets_MET_450;

	// 3 jets bins
	bool ThreeJets_hightMT2W_MET_250;
	bool ThreeJets_hightMT2W_MET_350;
	bool ThreeJets_hightMT2W_MET_450;
	bool ThreeJets_hightMT2W_MET_550;

	// >=4 jets bins
	
	// low MT2W
	bool FourJets_lowMT2W_MET_250;
	bool FourJets_lowMT2W_MET_350;
	bool FourJets_lowMT2W_MET_450;

	//- high MT2W
	bool FourJets_highMT2W_MET_250;
	bool FourJets_highMT2W_MET_350;
	bool FourJets_highMT2W_MET_450;
	bool FourJets_highMT2W_MET_550;
	bool FourJets_highMT2W_MET_650;

        //not MET divided
        bool TwoJets;
        bool ThreeJetsHighMT2W;
        bool FourJetsHighMT2W;
        bool FourJetsLowMT2W;

	void Reset();
	void Update(bool data, string dataset);

    private:
        bool passMETMHTTrigger;
        bool passElTrigger;
        bool passMuTrigger;

};

EventCategory eventCat;

void EventCategory::Reset(){
	// One lepton events (CR or CR1l) - includes the vetoes
	 OneLep = false;
	 OneLep_e = false;
	 OneLep_mu = false;
	// Two lepton events (CR-2l)
	 TwoLep = false;
	 TwoLep_ee = false;
	 TwoLep_emu = false;
	 TwoLep_mumu = false;
	// At least one b-tagged jets
	 OneBtag = false;
        //trigged
        trigged = false;
        passMETMHTTrigger = false;
        passElTrigger = false;
        passMuTrigger = false;
	//Baseline: >=2 jets, MET, MT, DPhi
	 PassBaselineCuts = false;
	// 2 jets bins (include mod_top>6.4
	 TwoJets_MET_250 = false;
	 TwoJets_MET_350 = false;
	 TwoJets_MET_450 = false;
	// 3 jets bins
	 ThreeJets_hightMT2W_MET_250 = false;
	 ThreeJets_hightMT2W_MET_350 = false;
	 ThreeJets_hightMT2W_MET_450 = false;
	 ThreeJets_hightMT2W_MET_550 = false;
	// >=4 jets bins
	// low MT2W
	 FourJets_lowMT2W_MET_250 = false;
	 FourJets_lowMT2W_MET_350 = false;
	 FourJets_lowMT2W_MET_450 = false;
	//- high MT2W
	 FourJets_highMT2W_MET_250 = false;
	 FourJets_highMT2W_MET_350 = false;
	 FourJets_highMT2W_MET_450 = false;
	 FourJets_highMT2W_MET_550 = false;
	 FourJets_highMT2W_MET_650 = false;
         //no MET binning
         TwoJets = false;
         ThreeJetsHighMT2W = false;
         FourJetsHighMT2W = false;
         FourJetsLowMT2W = false;
}

void EventCategory::Update(bool data, string dataset){
        if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
	// One lepton events (CR or CR1l) - includes the vetoes
	if(myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto){
		OneLep = true;
		if(abs(myEvent.lep1_pdgid) == 11) OneLep_e = true;
		if(abs(myEvent.lep1_pdgid) == 13) OneLep_mu = true;
	}

	if( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto){ //@MJ@ TODO tau veto should be here?!
                //@MJ@ TODO lep2 does not have to be veto lepton!!!
		TwoLep = true;
		if(abs(myEvent.lep1_pdgid) == 11 && abs(myEvent.lep2_pdgid) == 11) TwoLep_ee = true;
		if(abs(myEvent.lep1_pdgid) == 11 && abs(myEvent.lep2_pdgid) == 13) TwoLep_emu = true;
		if(abs(myEvent.lep1_pdgid) == 13 && abs(myEvent.lep2_pdgid) == 11) TwoLep_emu = true;
		if(abs(myEvent.lep1_pdgid) == 13 && abs(myEvent.lep2_pdgid) == 13) TwoLep_mumu = true;
	}

	// At least one b-tagged jets
	if( myEvent.ngoodbtags>=1 ) OneBtag = true;

        if(data)
        {
		// METMH, electron or muon trigger 
		//for(unsigned int i=0; i<myEvent.trigger_name.size();i++)
		//{
			// -- MET MHT trigger
			if(myEvent.HLT_PFMET170 == true && dataset.find("MET")!=std::string::npos)
                                passMETMHTTrigger = true;
			if(myEvent.HLT_PFMET100_PFMHT100_IDTight  == true && dataset.find("MET")!=std::string::npos)
				passMETMHTTrigger = true;
			// -- Electron trigger
			if(myEvent.HLT_Ele25_eta2p1_WPLoose  == true && dataset.find("SE")!=std::string::npos)
			        if(abs(myEvent.lep1_pdgid) == 11) passElTrigger = true;
			if(myEvent.HLT_Ele27_eta2p1_WPLoose  == true && dataset.find("SE")!=std::string::npos )		
				if(abs(myEvent.lep1_pdgid) == 11) passElTrigger = true;
			
			// -- Muon trigger
			if(myEvent.HLT_IsoMu20 == true && dataset.find("SM")!=std::string::npos)
				if(abs(myEvent.lep1_pdgid) == 13) passMuTrigger = true;
			if(myEvent.HLT_IsoMu22  == true && dataset.find("SM")!=std::string::npos)
				if(abs(myEvent.lep1_pdgid) == 13) passMuTrigger = true; //string in recompute that I have the right dataset name
		 //}

		 if(passMETMHTTrigger == true || passElTrigger == true || passMuTrigger == true) trigged = true;
         }
         else
         {
             trigged = true;
         }


	//Baseline: >=2 jets, MET, MT, DPhi
    	if (myEvent.pfmet >= 250 && myEvent.mt_met_lep > 150 && myEvent.ngoodjets >=2 && myEvent.dphi_ak4pfjets_met >= 0.8 && trigged == true )   PassBaselineCuts = true;
        //cout << "met " << myEvent.pfmet << " mtmet " << myEvent.mt_met_lep << " goodjets " << myEvent.ngoodjets << " dphi " << myEvent.dphi_ak4pfjets_met << " trigged: " << trigged << "blcutpasssed: " << PassBaselineCuts << endl;
	
	// 2 jets bins (include mod_top>6.4
	if(myEvent.ngoodjets == 2 && myEvent.topness>6.4){
                TwoJets = true;
		if(myEvent.pfmet >= 250 && myEvent.pfmet <350) 
                {
                 TwoJets_MET_250 = true;
                }
		if(myEvent.pfmet >= 350 && myEvent.pfmet <450)
                {
                 TwoJets_MET_350 = true;
                }
		if(myEvent.pfmet >= 450)
                {
                 TwoJets_MET_450 = true;
                }
	}

	// 3 jets bins
	if(myEvent.ngoodjets == 3 && myEvent.MT2W>200){
                ThreeJetsHighMT2W = true;
		if(myEvent.pfmet >= 250 && myEvent.pfmet <350) ThreeJets_hightMT2W_MET_250 = true;
		if(myEvent.pfmet >= 350 && myEvent.pfmet <450) ThreeJets_hightMT2W_MET_350 = true;
		if(myEvent.pfmet >= 450 && myEvent.pfmet <550) ThreeJets_hightMT2W_MET_450 = true;
		if(myEvent.pfmet >= 550) ThreeJets_hightMT2W_MET_550 = true;
	}

	// >=4 jets bins
	if(myEvent.ngoodjets >= 4){
		if(myEvent.MT2W<200){
                        FourJetsLowMT2W = true;
			if(myEvent.pfmet >= 250 && myEvent.pfmet <350) FourJets_lowMT2W_MET_250 = true;
			if(myEvent.pfmet >= 350 && myEvent.pfmet <450) FourJets_lowMT2W_MET_350 = true;
			if(myEvent.pfmet >= 450) FourJets_lowMT2W_MET_450 = true;
		}
		else{
                        FourJetsHighMT2W = true;
			if(myEvent.pfmet >= 250 && myEvent.pfmet <350) FourJets_highMT2W_MET_250 = true;
			if(myEvent.pfmet >= 350 && myEvent.pfmet <450) FourJets_highMT2W_MET_350 = true;
			if(myEvent.pfmet >= 450 && myEvent.pfmet <550) FourJets_highMT2W_MET_450 = true;
			if(myEvent.pfmet >= 550 && myEvent.pfmet <650) FourJets_highMT2W_MET_550 = true;
			if(myEvent.pfmet >= 650) FourJets_highMT2W_MET_650 = true;
		}
	}
	
}



// -- Function used for CR2l events where 2nd lepton is added to the MET
void RecomputeMETAndDerivedQuantities(){
}

void recompute(bool dataIn = false, string dataset = "")
{
    eventCat.Reset();
    eventCat.Update(dataIn, dataset);
}

//--------------
//   Channels 
//--------------
//
bool Baseline()
{return eventCat.PassBaselineCuts;}
//-- Signal regions
bool SR1l() {
//cout << "SR1l; one lep: " << eventCat.OneLep << " pass bl: " << eventCat.PassBaselineCuts << " oneB: " << eventCat.OneBtag << endl;
return eventCat.OneLep && eventCat.PassBaselineCuts && eventCat.OneBtag;}
bool SR1l_e() {return eventCat.OneLep_e && eventCat.PassBaselineCuts && eventCat.OneBtag;}
bool SR1l_mu() {return eventCat.OneLep_mu && eventCat.PassBaselineCuts && eventCat.OneBtag;}

//-- CR1l regions
bool CR1l() {
//cout << "CR1l: one lep: " << eventCat.OneLep << " pass bl: " << eventCat.PassBaselineCuts << " oneB: " << eventCat.OneBtag << endl;
return eventCat.OneLep && eventCat.PassBaselineCuts && !eventCat.OneBtag;}
bool CR1l_e() {return eventCat.OneLep_e && eventCat.PassBaselineCuts && !eventCat.OneBtag;}
bool CR1l_mu() {return eventCat.OneLep_mu && eventCat.PassBaselineCuts && !eventCat.OneBtag;}


// Should be add Z-veto and b-tagging requirement ??
bool CR2l() {
//cout << "CR2l: one lep: " << eventCat.TwoLep << " pass bl: " << eventCat.PassBaselineCuts << endl;
return eventCat.TwoLep && eventCat.PassBaselineCuts;};
bool CR2l_ee() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};
bool CR2l_emu() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};
bool CR2l_mumu() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};





//-----------------
// Selection  Bins
//-----------------


// 2 jets bins (include mod_top>6.4)
bool TwoJets_MET_250(){ return eventCat.TwoJets_MET_350;}
bool TwoJets_MET_350(){ return eventCat.TwoJets_MET_350;}
bool TwoJets_MET_450(){ return eventCat.TwoJets_MET_450;}

// 3 jets bins
bool ThreeJets_hightMT2W_MET_250(){ return eventCat.ThreeJets_hightMT2W_MET_250;}
bool ThreeJets_hightMT2W_MET_350(){ return eventCat.ThreeJets_hightMT2W_MET_350;}
bool ThreeJets_hightMT2W_MET_450(){ return eventCat.ThreeJets_hightMT2W_MET_450;}
bool ThreeJets_hightMT2W_MET_550(){ return eventCat.ThreeJets_hightMT2W_MET_550;}

// >=4 jets bins

// low MT2W
bool FourJets_lowMT2W_MET_250(){ return eventCat.FourJets_lowMT2W_MET_250;}
bool FourJets_lowMT2W_MET_350(){ return eventCat.FourJets_lowMT2W_MET_350;}
bool FourJets_lowMT2W_MET_450(){ return eventCat.FourJets_lowMT2W_MET_450;}

//- high MT2W
bool FourJets_highMT2W_MET_250(){ return eventCat.FourJets_highMT2W_MET_250;}
bool FourJets_highMT2W_MET_350(){ return eventCat.FourJets_highMT2W_MET_350;}
bool FourJets_highMT2W_MET_450(){ return eventCat.FourJets_highMT2W_MET_450;}
bool FourJets_highMT2W_MET_550(){ return eventCat.FourJets_highMT2W_MET_550;}
bool FourJets_highMT2W_MET_650(){ return eventCat.FourJets_highMT2W_MET_650;}

//no MET binning
bool TwoJets(){return eventCat.TwoJets;}
bool ThreeJetsHighMT2W(){return eventCat.ThreeJetsHighMT2W;}
bool FourJetsHighMT2W(){return eventCat.FourJetsHighMT2W;}
bool FourJetsLowMT2W(){return eventCat.FourJetsLowMT2W;}

//////////////
/// SR and CR
/////////////

//SR1l
// 2 jets bins (include mod_top>6.4)
bool SR1l2jMET250to350(){ return eventCat.TwoJets_MET_250 && SR1l();}
bool SR1l2jMET350to450(){ return eventCat.TwoJets_MET_350 && SR1l();}
bool SR1l2jMET450toInf(){ return eventCat.TwoJets_MET_450 && SR1l();}

// 3 jets bins
bool SR1l3jMET250to350(){ return eventCat.ThreeJets_hightMT2W_MET_250 && SR1l();}
bool SR1l3jMET350to450(){ return eventCat.ThreeJets_hightMT2W_MET_350 && SR1l();}
bool SR1l3jMET450to550(){ return eventCat.ThreeJets_hightMT2W_MET_450 && SR1l();}
bool SR1l3jMET550toInf(){ return eventCat.ThreeJets_hightMT2W_MET_550 && SR1l();}

// >=4 jets bins

// low MT2W
bool SR1l4jMET250to350lowMT2W(){ return eventCat.FourJets_lowMT2W_MET_250 && SR1l();}
bool SR1l4jMET350to450lowMT2W(){ return eventCat.FourJets_lowMT2W_MET_350 && SR1l();}
bool SR1l4jMET450toInflowMT2W(){ return eventCat.FourJets_lowMT2W_MET_450 && SR1l();}

//- high MT2W
bool SR1l4jMET250to350highMT2W(){ return eventCat.FourJets_highMT2W_MET_250 && SR1l();}
bool SR1l4jMET350to450highMT2W(){ return eventCat.FourJets_highMT2W_MET_350 && SR1l();}
bool SR1l4jMET450to550highMT2W(){ return eventCat.FourJets_highMT2W_MET_450 && SR1l();}
bool SR1l4jMET550to650highMT2W(){ return eventCat.FourJets_highMT2W_MET_550 && SR1l();}
bool SR1l4jMET650toInfhighMT2W(){ return eventCat.FourJets_highMT2W_MET_650 && SR1l();}

//no MET binning

bool SR1lTwoJets() {return eventCat.TwoJets && SR1l();}
bool SR1lThreeJetsHighMT2W() {return eventCat.ThreeJetsHighMT2W && SR1l();}
bool SR1lFourJetsHighMT2W() {return eventCat.FourJetsHighMT2W && SR1l();}
bool SR1lFourJetsLowMT2W() {return eventCat.FourJetsLowMT2W && SR1l();}

/////CR2l
// 2 jets bins (include mod_top>6.4)
bool CR2l2jMET250to350(){ return eventCat.TwoJets_MET_250 && CR2l();}
bool CR2l2jMET350to450(){ return eventCat.TwoJets_MET_350 && CR2l();}
bool CR2l2jMET450toInf(){ return eventCat.TwoJets_MET_450 && CR2l();}

// 3 jets bins
bool CR2l3jMET250to350(){ return eventCat.ThreeJets_hightMT2W_MET_250 && CR2l();}
bool CR2l3jMET350to450(){ return eventCat.ThreeJets_hightMT2W_MET_350 && CR2l();}
bool CR2l3jMET450to550(){ return eventCat.ThreeJets_hightMT2W_MET_450 && CR2l();}
bool CR2l3jMET550toInf(){ return eventCat.ThreeJets_hightMT2W_MET_550 && CR2l();}

// >=4 jets bins

// low MT2W
bool CR2l4jMET250to350lowMT2W(){ return eventCat.FourJets_lowMT2W_MET_250 && CR2l();}
bool CR2l4jMET350to450lowMT2W(){ return eventCat.FourJets_lowMT2W_MET_350 && CR2l();}
bool CR2l4jMET450toInflowMT2W(){ return eventCat.FourJets_lowMT2W_MET_450&& CR2l();}

//- high MT2W
bool CR2l4jMET250to350highMT2W(){ return eventCat.FourJets_highMT2W_MET_250 && CR2l();}
bool CR2l4jMET350to450highMT2W(){ return eventCat.FourJets_highMT2W_MET_350 && CR2l();}
bool CR2l4jMET450to550highMT2W(){ return eventCat.FourJets_highMT2W_MET_450 && CR2l();}
bool CR2l4jMET550to650highMT2W(){ return eventCat.FourJets_highMT2W_MET_550 && CR2l();}
bool CR2l4jMET650toInfhighMT2W(){ return eventCat.FourJets_highMT2W_MET_650 && CR2l();}

/////CR1l
// 2 jets bins (include mod_top>6.4)
bool CR1l2jMET250to350(){ return eventCat.TwoJets_MET_250 && CR1l();}
bool CR1l2jMET350to450(){ return eventCat.TwoJets_MET_350 && CR1l();}
bool CR1l2jMET450toInf(){ return eventCat.TwoJets_MET_450 && CR1l();}

// 3 jets bins
bool CR1l3jMET250to350(){ return eventCat.ThreeJets_hightMT2W_MET_250 && CR1l();}
bool CR1l3jMET350to450(){ return eventCat.ThreeJets_hightMT2W_MET_350 && CR1l();}
bool CR1l3jMET450to550(){ return eventCat.ThreeJets_hightMT2W_MET_450 && CR1l();}
bool CR1l3jMET550toInf(){ return eventCat.ThreeJets_hightMT2W_MET_550 && CR1l();}

// >=4 jets bins

// low MT2W
bool CR1l4jMET250to350lowMT2W(){ return eventCat.FourJets_lowMT2W_MET_250 && CR1l();}
bool CR1l4jMET350to450lowMT2W(){ return eventCat.FourJets_lowMT2W_MET_350 && CR1l();}
bool CR1l4jMET450toInflowMT2W(){ return eventCat.FourJets_lowMT2W_MET_450 && CR1l();}

//- high MT2W
bool CR1l4jMET250to350highMT2W(){ return eventCat.FourJets_highMT2W_MET_250 && CR1l();}
bool CR1l4jMET350to450highMT2W(){ return eventCat.FourJets_highMT2W_MET_350 && CR1l();}
bool CR1l4jMET450to550highMT2W(){ return eventCat.FourJets_highMT2W_MET_450 && CR1l();}
bool CR1l4jMET550to650highMT2W(){ return eventCat.FourJets_highMT2W_MET_550 && CR1l();}
bool CR1l4jMET650toInfhighMT2W(){ return eventCat.FourJets_highMT2W_MET_650 && CR1l();}

//all regions
bool allRegions()
{
return SR1l2jMET250to350() || SR1l2jMET350to450() || SR1l2jMET450toInf() || SR1l3jMET250to350() || SR1l3jMET350to450() || SR1l3jMET450to550() || SR1l3jMET550toInf() ||
SR1l4jMET250to350lowMT2W() || SR1l4jMET350to450lowMT2W() || SR1l4jMET450toInflowMT2W() || 
SR1l4jMET250to350highMT2W() || SR1l4jMET350to450highMT2W() || SR1l4jMET450to550highMT2W() || SR1l4jMET550to650highMT2W() || SR1l4jMET650toInfhighMT2W();
}
