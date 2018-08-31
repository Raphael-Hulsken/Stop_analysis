#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//#include "../../common/common.h"

#define _METCUT_ 50
#define _LEPTONPTCUT_ 40


using namespace std;


// --------------------------
// Define the input format 
// --------------------------

typedef struct 
{
    Int_t      totalNumberOfInitialEvent = -13;
    Float_t	   leadingLeptonPt;		// PT of the leading selected lepton
    Int_t	   leadingLeptonId;		// Id of the leading selected lepton
    Float_t	   secondLeptonPt;		// PT of the second leading selected lepton
    Int_t	   secondLeptonId;		// Id of the second leading selected lepton
    
    Int_t                nJets;               // Number of selected jets
    Float_t MET;                                // Type-1 + phi-corrected PF MET

    Int_t	ngoodleps;
    Int_t	nvetoleps;
    Int_t	ngoodjets;
    bool	PassTrackVeto;
    bool	PassTauVeto;

    vector<string>	trigger_name;	 	// Name of the trigger paths
    vector<bool>	trigger_pass;		// Boolean corresponding to trigger results
    
    // Use pointers for reading 
    vector<string>*	    pointerToTrigger_name;
    vector<bool>*	    pointerToTrigger_pass;

} babyEvent;

babyEvent myEvent;




// --------------------------
// Functions used to read the tree
// --------------------------
void InitializeBranchesForReading(TTree* theTree, babyEvent* myEvent)
{
    theTree->SetBranchAddress("totalNumberOfInitialEvent", &(myEvent->totalNumberOfInitialEvent));
    theTree->SetBranchAddress("leadingLeptonPt",                              &(myEvent->leadingLeptonPt));
    theTree->SetBranchAddress("leadingLeptonId",                              &(myEvent->leadingLeptonId));
    theTree->SetBranchAddress("secondLeptonPt",                               &(myEvent->secondLeptonPt));
    theTree->SetBranchAddress("secondLeptonId",                               &(myEvent->secondLeptonId));
    theTree->SetBranchAddress("numberOfSelectedJets",                         &(myEvent->nJets));
    theTree->SetBranchAddress("pfmet",                                          &(myEvent->MET));
    theTree->SetBranchAddress("trigger_name",				      &(myEvent->pointerToTrigger_name));
    theTree->SetBranchAddress("trigger_pass",				      &(myEvent->pointerToTrigger_pass));
    theTree->SetBranchAddress("ngoodleps", 	&myEvent->ngoodleps);
    theTree->SetBranchAddress("nvetoleps", 	&myEvent->nvetoleps);
    theTree->SetBranchAddress("ngoodjets", 	&myEvent->ngoodjets);
    theTree->SetBranchAddress("PassTrackVeto", 	&myEvent->PassTrackVeto);
    theTree->SetBranchAddress("PassTauVeto", 	&myEvent->PassTauVeto);
}


// --------------------------
// Function called @ each event
// --------------------------
void ReadEvent(TTree* theTree, long int i, babyEvent* myEvent)
{
      theTree->GetEntry(i);
      myEvent->trigger_name		    = *(myEvent->pointerToTrigger_name);
      myEvent->trigger_pass		    = *(myEvent->pointerToTrigger_pass);
}

// --------------------------
// Intermetidate struct to
// store trigger results 
// --------------------------
typedef struct{
	// triggers used in the analysis;
	bool passElTrigger;
	bool passMuTrigger;
	bool passMETTrigger;
	bool passMETMHTTrigger;
	bool passORTrigger;
  	bool passHTTrigger;

	bool passDoubleMuTrigger;
	bool passDoubleElTrigger;
	bool passMuElTrigger;

	// orthogonal triggers used for measurement
	bool passOrthogElTrigger;
	bool passOrthogMuTrigger;
	bool passOrthogORTrigger;

	//Reset
	void Reset(){
		passElTrigger = false; passMuTrigger = false; passMETTrigger = false; passMETMHTTrigger = false; passORTrigger = false;
		passDoubleMuTrigger = false; passDoubleElTrigger = false; passMuElTrigger = false;
		passOrthogElTrigger = false; passOrthogMuTrigger = false; passOrthogORTrigger = false;
		passHTTrigger = false;
	}
	void Print(){
			cout<<"Analysis triggers :\t electron:"<<passElTrigger<<"\t muon:"<<passMuTrigger<<"\t MET:"<<passMETTrigger<<"\tMETMHT:"<<passMETMHTTrigger<<"\tOR"<<passORTrigger<<endl;
		cout<<"Orthognal triggers:\t electron:"<<passOrthogElTrigger<<"\t muon:"<<passOrthogMuTrigger<<"\t OR:"<<passOrthogORTrigger<<endl;
	}
}trigResults;



trigResults myTrigResults;


// Produce Trigger efficiency as function of MET, Lepton pT, and a 2D version
// Trigger to be considered:
// - MET triggers
// - Lepton triggers
// - Their combination

// --------------------------
// Function to check 
// the trigger results 
// --------------------------

///*
bool CheckTriggerResults(){
	myTrigResults.Reset();
	for(unsigned int i=0; i<myEvent.trigger_name.size();i++){
		//cout<<myEvent.trigger_name[i]<<endl;
		//if(myEvent.trigger_pass[i]) cout<<myEvent.trigger_name[i]<<endl;
		//----------------------
		// ANALYSIS TRIGGERS
		//----------------------
		// -- MET trigger
		//HLT_PFMET170_NoiseCleaned_v4 HLT_PFMET170_HBHECleaned_v3 HLT_PFMET170_JetIdCleaned_v3 HLT_PFMET170_NotCleaned_v2 HLT_PFMET170_BeamHaloCleaned_v1
		if(myEvent.trigger_name[i].find("HLT_PFMET170_")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMETTrigger = true;
		// -- MET MHT trigger
		if(myEvent.trigger_name[i].find("HLT_PFMET100_PFMHT100_IDTight")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMETMHTTrigger = true;
		// -- HT Trigger
		if(myEvent.trigger_name[i].find("HLT_HT")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passHTTrigger = true;
		// -- Electron trigger
		//if(myEvent.trigger_name[i].find("HLT_Ele23_WPLoose_Gsf")!=std::string::npos)		
		if(myEvent.trigger_name[i].find("HLT_Ele25_eta2p1_WPLoose")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Ele27_eta2p1_WPLoose")!=std::string::npos )		
			if(myEvent.trigger_pass[i]) myTrigResults.passElTrigger = true;
		
		// -- Muon trigger
		if(myEvent.trigger_name[i].find("HLT_IsoMu20_")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_IsoTkMu20_")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuTrigger = true;
	
 		// -- Double lepton trigger
		if(myEvent.trigger_name[i].find("HLT_DoubleMu")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passDoubleMuTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_DoubleEl")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passDoubleElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v5")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v5")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v3")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v5")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passMuElTrigger = true;

		//----------------------
		// ORTHOGONAL TRIGGERS
		//----------------------
		if(myEvent.trigger_name[i].find("HLT_IsoMu8_")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passOrthogMuTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_IsoMu20_")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passOrthogMuTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_Ele22_eta2p1_WPLoose_Gsf_*")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passOrthogElTrigger = true;
		if(myEvent.trigger_name[i].find("HLT_PF_")!=std::string::npos)
			if(myEvent.trigger_pass[i]) myTrigResults.passOrthogORTrigger = true;
		
	}
	return false;
}

//*/
bool PassMETMHTTrigger(){ return myTrigResults.passMETMHTTrigger;}
bool PassMETTrigger(){ return myTrigResults.passMETTrigger;}
bool PassElTrigger(){return myTrigResults.passElTrigger;}
bool PassMuTrigger(){return myTrigResults.passMuTrigger;}
bool PassORTrigger(){ return myTrigResults.passMETTrigger || myTrigResults.passElTrigger || myTrigResults.passMuTrigger; }

bool PassDoubleMuTrigger() {return myTrigResults.passDoubleMuTrigger;}
bool PassDoubleElTrigger() {return myTrigResults.passDoubleElTrigger;}
bool PassMuElTrigger() {return myTrigResults.passMuElTrigger;}
bool PassDoubleLeptonTrigger() {return (myTrigResults.passDoubleMuTrigger || myTrigResults.passDoubleElTrigger || myTrigResults.passMuElTrigger);}

bool PassDoubleMuTriggerAndMuSel() {return myTrigResults.passDoubleMuTrigger && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassDoubleElTriggerAndElSel() {return myTrigResults.passDoubleElTrigger && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}

bool PassDoubleMuTriggerAndDiMuSel() {return myTrigResults.passDoubleMuTrigger && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && 
	abs(myEvent.secondLeptonId) == 13 && (myEvent.secondLeptonPt > _LEPTONPTCUT_);}
bool PassDoubleElTriggerAndDiElSel() {return myTrigResults.passDoubleElTrigger && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && 
	abs(myEvent.secondLeptonId) == 11 && (myEvent.secondLeptonPt > _LEPTONPTCUT_);}

bool PassDoubleMuTriggerAndDiMuSelBaseline() {return myTrigResults.passDoubleMuTrigger && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && 
	abs(myEvent.secondLeptonId) == 13 && (myEvent.secondLeptonPt > _LEPTONPTCUT_) && myEvent.MET>250 && myEvent.nJets>=2; }
bool PassDoubleElTriggerAndDiElSelBaseline() {return myTrigResults.passDoubleElTrigger && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && 
	abs(myEvent.secondLeptonId) == 11 && (myEvent.secondLeptonPt > _LEPTONPTCUT_) && myEvent.MET>250 && myEvent.nJets>=2 ;}

// Selection for orthogonal trigger:
// - the OR of all PF*obbject triggers
// - the presecaled lepton trigger with lower threshold

bool PassOrthogORTrigger(){ return myTrigResults.passOrthogORTrigger;}
bool PassOrthogElTrigger(){ return myTrigResults.passOrthogElTrigger;}
bool PassOrthogMuTrigger(){ return myTrigResults.passOrthogMuTrigger;}

//PFMET100 PFMHT100 IDTight
//PFMET170

// Define an offline selection
// - exactly one selected leptons (separated muon & electron)
// - MET > 250

bool OfflineElSelection(){
	if (abs(myEvent.leadingLeptonId) != 11) return false;
	if (myEvent.MET<_METCUT_) return false;
	if (myEvent.leadingLeptonPt < _LEPTONPTCUT_) return false;
	//if (myEvent.secondLeptonPt > 20) return false;
	return true;
}

bool OfflineMuSelection(){
	if (abs(myEvent.leadingLeptonId) != 13) return false;
	if (myEvent.MET<_METCUT_) return false;
	if (myEvent.leadingLeptonPt < _LEPTONPTCUT_) return false;
	//if (myEvent.secondLeptonPt > 20) return false;
	return true;
}

//bool ElChannel(){	 (abs(myEvent.leadingLepton.PDGId) == 13) ? return true: return false;	 }
//bool MuChannel(){	 (abs(myEvent.leadingLepton.PDGId) == 13) ? return true: return false;	 }


bool PassOrthogORTriggerAndMuSelection(){ return myTrigResults.passOrthogORTrigger && OfflineMuSelection();}
bool PassOrthogORTriggerAndElSelection(){ return myTrigResults.passOrthogORTrigger && OfflineElSelection();}
bool PassOrthogMuTriggerAndMuSelection(){ return myTrigResults.passOrthogMuTrigger && OfflineMuSelection();}
bool PassOrthogElTriggerAndElSelection(){ return myTrigResults.passOrthogElTrigger && OfflineElSelection();}

//bool PassMETTriggerAndMuSelection() { return myTrigResults.passMETTrigger && OfflineMuSelection();}
//bool PassMETTriggerAndElSelection() { return myTrigResults.passMETTrigger && OfflineElSelection();}

bool PassElTriggerAndElSelection(){ return myTrigResults.passElTrigger && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassMuTriggerAndMuSelection(){ return myTrigResults.passMuTrigger && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}

bool PassMETTriggerAndElSelection() { return myTrigResults.passMETTrigger && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassMETTriggerAndMuSelection() { return myTrigResults.passMETTrigger && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}

bool PassMuTriggerAndMETSelection(){ return myTrigResults.passMuTrigger && (myEvent.MET>_METCUT_);}
bool PassElTriggerAndMETSelection(){ return myTrigResults.passElTrigger && (myEvent.MET>_METCUT_);}

bool PassCombinedMETTrigger(){ return (myTrigResults.passMETTrigger || myTrigResults.passElTrigger || myTrigResults.passMuTrigger) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassCombinedMETElTrigger(){ return (myTrigResults.passMETTrigger || myTrigResults.passElTrigger) && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassCombinedMETMuTrigger(){ return (myTrigResults.passMETTrigger || myTrigResults.passMuTrigger) && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
//bool PassCombinedMETMuTrigger(){ return (myTrigResults.passMETTrigger && myTrigResults.passMuTrigger) && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassCombinedMETMHTElTrigger(){ return (myTrigResults.passMETMHTTrigger || myTrigResults.passElTrigger) && (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}
bool PassCombinedMETMHTMuTrigger(){ return (myTrigResults.passMETMHTTrigger || myTrigResults.passMuTrigger) && (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_);}

//bool ElChannel(){ return (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && myEvent.nJets>=3;}
//bool MuChannel(){ return (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && myEvent.nJets>=3;}

bool BaselineSel(){ return myEvent.ngoodleps==1 && myEvent.nvetoleps==0 && myEvent.PassTauVeto && myEvent.PassTrackVeto && myEvent.ngoodjets>1;}
bool ElChannel(){ return (abs(myEvent.leadingLeptonId) == 11) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && BaselineSel() && myTrigResults.passHTTrigger;}
bool MuChannel(){ return (abs(myEvent.leadingLeptonId) == 13) && (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && BaselineSel() && myTrigResults.passHTTrigger;}
bool LepChannel(){ return (myEvent.leadingLeptonPt > _LEPTONPTCUT_) && BaselineSel() && myTrigResults.passHTTrigger;}

bool all(){return true;}

// ----------------------------------------------
// Should be called only here because many
// struct and fuctions have to be declare first
// ----------------------------------------------
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    //babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/2902/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/Synchro/CMSSW_8_0_5/src/PyROOF/test/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/Synchro/CMSSW_8_0_5/src/PyROOF/";
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/Synchro/CMSSW_8_0_5/src/store/babyTuples/TriggerStudy/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/CMSSW_7_4_12_patch4/src/StopAF/backgroundEstimation/TriggerStudy/";
    totalNumberOfWorkers = 10;


    // ------------------
    // Histograms 
    // ------------------
    //vector<float> METBins = {150,175,200,225,250,275,300,350,400,500,800};
    vector<float> METBins = {150,200,250,350,450,600,800};
    vector<float> LepPtBins = {20,25,30,35,40,50,80,100,125,150,200,300};
    //vector<float> METBins = {250,350,450,800};
    //vector<float> LepPtBins = {30,40,50,100,200,300};

    // ------------------
    // Variables
    // ------------------
    // fixed binning
    //AddVariable("MET", "MET",  "MET", 20,   20, 500,  &(myEvent.MET), "noUnderflowInFirstBin,noOverflowInLastBin");
    //AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", 20,   0, 200,  &(myEvent.leadingLeptonPt), "noUnderflowInFirstBin,noOverflowInLastBin");
    //AddVariable("MET", "MET",  "MET", 20,   20, 500,  &(myEvent.MET), "");
    //AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", 20,   0, 200,  &(myEvent.leadingLeptonPt), "");
    // defined binning
    //AddVariable("MET", "MET",  "MET", (int) (METBins.size()-1), METBins.data(),  &(myEvent.MET), "");
    AddVariable("MET", "MET",  "MET", 4 ,150,550,  &(myEvent.MET), "");
    AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", (int) (LepPtBins.size()-1), LepPtBins.data(), &(myEvent.leadingLeptonPt), "");
    AddVariable("nJets","nJets","nJets",5,1,5,&(myEvent.nJets),"");

    // ------------------
    // Datasets
    // ------------------
    //AddProcessClass("T2tt_600-950_1to450", "T2tt_600-950_1to450", "signal", kBlue);
    // 	AddDataset("600-950_1to450", "T2tt_600-950_1to450", 0, 0 );
    
    AddProcessClass("test", "test", "signal", kBlue);
    	//AddDataset("test","test",0,0);
	AddDataset("SingleElectron", "test", 0, 0 );
	AddDataset("SingleMuon", "test", 0, 0 );
	//AddDataset("DoubleMuon", "test", 0, 0 );
	//AddDataset("DoubleElectron", "test", 0, 0 );
	AddDataset("JetHT", "test", 0, 0 );
    	//AddDataset("SElectron_1", "test", 0, 0 );
    
    
    // ------------------
    // Regions
    // ------------------
    // Referecne
    AddRegion("all","all",&all);
    // triggers used in the analysis
    AddRegion("METTrigger","MET Trigger", &PassMETTrigger);
    AddRegion("METMHTTrigger","MET Trigger", &PassMETMHTTrigger);
    AddRegion("ElTrigger","Electron Trigger", &PassElTrigger);
    AddRegion("MuTrigger","Muon Trigger", &PassMuTrigger);
    AddRegion("ORTrigger","ORTrigger", &PassORTrigger);
    AddRegion("CombinedMET","CombinedMETTrigger", & PassCombinedMETTrigger);
    AddRegion("CombinedMETMu","CombinedMETMuTrigger", & PassCombinedMETMuTrigger);
    AddRegion("CombinedMETEl","CombinedMETElTrigger", & PassCombinedMETElTrigger);
    AddRegion("CombinedMETMHTMu","CombinedMETMHTMuTrigger", & PassCombinedMETMHTMuTrigger);
    AddRegion("CombinedMETMHTEl","CombinedMETMHTElTrigger", & PassCombinedMETMHTElTrigger);


    // ------------------
    // Channels
    // ------------------
    
    AddChannel("ORTriggerAndMuSel","ORTrigger - Muon Selection",&PassOrthogORTriggerAndMuSelection);
    AddChannel("ORTriggerAndElSel","ORTrigger - Electron Selection",&PassOrthogORTriggerAndElSelection);
    //AddChannel("ElTriggerAndElSel","ElTrigger - Electron Selection",&PassOrthogElTriggerAndElSelection);
    //AddChannel("MuTriggerAndMuSel","MuTrigger - Muon Selection",&PassOrthogMuTriggerAndMuSelection);
    AddChannel("ElTriggerAndElSel","ElTrigger - Electron Selection",&PassElTriggerAndElSelection);
    AddChannel("MuTriggerAndMuSel","MuTrigger - Muon Selection",&PassMuTriggerAndMuSelection);
    AddChannel("METTriggerAndMuSel","METTrigger - Muon Selection", &PassMETTriggerAndMuSelection);
    AddChannel("METTriggerAndElSel","METTrigger - Electron Selection", &PassMETTriggerAndElSelection);
    AddChannel("MuTriggerAndMETSel","MuTrigger - MET Selection",&PassMuTriggerAndMETSelection);
    AddChannel("ElTriggerAndMETSel","ElTrigger - MET Selection",&PassElTriggerAndMETSelection);
    AddChannel("METTrigger","METTrigger",&PassMETTrigger);
    AddChannel("ElTrigger","ElectronTrigger",&PassElTrigger);
    AddChannel("MuTrigger","MuonTrigger",&PassMuTrigger);
    //AddChannel("inclusive","inclusive",&all);
    AddChannel("All","All",&all);
    
    AddChannel("ElChannel","ElChannel", &ElChannel);
    AddChannel("MuChannel","MuChannel", &MuChannel);
    AddChannel("LepChannel","MuChannel", &MuChannel);
    // ...
    AddChannel("DoubleMuTrigger","DoubleMuTrigger",&PassDoubleMuTrigger);
    AddChannel("DoubleElTrigger","DoubleElTrigger",&PassDoubleElTrigger);
    AddChannel("MuElTrigger","MuElTrigger",&PassMuElTrigger);
    AddChannel("DoubleLeptonTrigger","DoubleLeptonTrigger",&PassDoubleLeptonTrigger);

    AddChannel("DoubleMuTriggerAndMuSel","PassDoubleMuTriggerAndMuSel",&PassDoubleMuTriggerAndMuSel);
    AddChannel("DoubleElTriggerAndElSel","PassDoubleElTriggerAndElSel",&PassDoubleElTriggerAndElSel);
    AddChannel("DoubleMuTriggerAndDiMuSel","PassDoubleMuTriggerAndDiMuSel",&PassDoubleMuTriggerAndDiMuSel);
    AddChannel("DoubleElTriggerAndDiElSel","PassDoubleElTriggerAndDiElSel",&PassDoubleElTriggerAndDiElSel);
    AddChannel("DoubleMuTriggerAndDiMuSelBaseline","PassDoubleMuTriggerAndDiMuSelBaseline",&PassDoubleMuTriggerAndDiMuSelBaseline);
    AddChannel("DoubleElTriggerAndDiElSelBaseline","PassDoubleElTriggerAndDiElSelBaseline",&PassDoubleElTriggerAndDiElSelBaseline);

    SetLumi(2440.);

    Create1DHistos();
    Add2DHisto("LeptonPT","MET");
    Add2DHisto("nJets","MET");

    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    //cout<<"MET = "<<myEvent.MET<<endl;
    //Check the output of the trigger results
    CheckTriggerResults();
    //myTrigResults.Print();
    //exit(1);
    float weight = 1.0;
    AutoFillProcessClass(currentProcessClass, weight);
    /* 
    cout<<"Regions:  "<<PassMETTrigger()<<" "<<PassMETMHTTrigger()<<" "<<PassElTrigger()<<" "<<PassMuTrigger()<<" "<<PassORTrigger()<<endl;
    cout<<"Channels: "<<PassOrthogORTriggerAndMuSelection()<<" "<<PassOrthogORTriggerAndElSelection()<<" "<<PassOrthogElTriggerAndElSelection()<<" "<<PassOrthogMuTriggerAndMuSelection();
    cout<<" "<<PassMETTriggerAndMuSelection()<<" "<<PassMETTriggerAndElSelection()<<endl;
    cout<<"Selection: "<<OfflineElSelection()<<" "<<OfflineMuSelection()<<endl;
    */
}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    // ######################
    //  Plot configuration and production
    // ######################

    // Schedule plots

    SchedulePlots("1DSuperimposed");
    //SchedulePlots("1DStack");
    SchedulePlots("2D");

    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", true);
    SetGlobalBoolOption("Plot", "exportEps", false);
    SetGlobalBoolOption("Plot", "exportPng", false);

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plotsTrigger/");
 }
