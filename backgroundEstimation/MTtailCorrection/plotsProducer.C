#include <vector>
#include <iostream>
using namespace std;

#include "TTree.h"
//babyEvent myEvent;




//-----------------------------------
//-- Load only necessary branches
//-----------------------------------

#define USE_GEN_INFO
//#define USE_GEN_INFO_EXT
//#define USE_SKIMMING_VAR
#define USE_JETS
//#define USE_JETS_EXT
//#define USE_GEN_INFO
//#define USE_GEN_INFO_EXT
#define USE_LEP1
//#define USE_LEP1_EXT
#define USE_LEP2
//#define USE_LEP2_EXT
#define USE_JETS
//#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_VAR_BASELINE
#define USE_GLOBAL_VAR
//#define USE_OLD_VAR


//-----------------------------------
// -- load Babytuple Format --
// -- Need first to activate the branches
//-----------------------------------
//#include "../common/BabyTupleFormat.h"
#include "../../common/Reader_CommonFormat.h"

#include "../../common/selectionDefinitions13TeV.h"


//---------------------------------------
// -- Add variables computed on the flyt
// -- Need first to load the format
//--------------------------------------
//#include "OnTheFlyVariables.h"

// Should be called after selectionDefinition 'cause it need myEvent varible ...
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"

//------------------------------
//  All selection needed
//------------------------------


float topPtWeight(){
        float weight = 1;
	float TopPt = 0.0;
        for (unsigned int i = 0 ; i < myEvent.gen_status.size() ; i++)
            if (myEvent.gen_id[i] == 6 && myEvent.gen_status[i] == 22) { 
		TopPt = myEvent.gen_pt[i]; 
        	weight = exp(0.156 - 0.00137 * TopPt);
		break; 
	}
	cout<<weight<<endl;
    	return weight;
}

bool genInfo(){
        for (unsigned int i = 0 ; i < myEvent.gen_status.size() ; i++)
            if ( (abs(myEvent.gen_id[i]) == 4 || abs(myEvent.gen_id[i])  == 5) && myEvent.gen_status[i] == 23) { 
		//cout<<myEvent.gen_id[i]<<" "<<myEvent.gen_status[i]<<" "<<myEvent.gen_pt[i]<<endl;
		return true;
	}
	return false;
}

float DeltaRlj(){
	TLorentzVector lep;
	lep.SetPtEtaPhiM(myEvent.lep1_pt,myEvent.lep1_eta,myEvent.lep1_phi,myEvent.lep1_mass);
	float minDR = 9999;
	for(unsigned int i=0;i<myEvent.jet_pt.size();i++){
		TLorentzVector p4;
		p4.SetPtEtaPhiM(myEvent.jet_pt[i],myEvent.jet_eta[i],myEvent.jet_phi[i],myEvent.jet_mass[i]);
		float dr = lep.DeltaR(p4);
		if(dr<minDR) minDR = dr; 
	}
	return minDR;
}

void dumpSelection(){
	cout<<"MET = " <<myEvent.pfmet<<"\t";
	cout<<"MT = " <<myEvent.mt_met_lep<<"\t";
	cout<<"Njets = " <<myEvent.ngoodjets<<"\t";
	cout<<"Nbjets = " <<myEvent.ngoodbtags<<"\t";
	cout<<"Nleps = " <<myEvent.ngoodleps<<"\t";
	cout<<"lep-id = " <<myEvent.lep1_pdgid<<"\t";
	cout<<"TauVeto = " <<myEvent.PassTauVeto<<"\t";
	cout<<"TauVeto = " <<myEvent.PassTrackVeto<<"\t";
	cout<<endl;
}


// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");
    cout << "non correct version" << endl;
    //babyTuplePath = "./";
    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_all_3j_lowpt_dphi_v2/";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_missingInfo/";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_all/";

    totalNumberOfWorkers = 8;
    
    
    
    AddVariable("Mlb_leadb", "M_lb [GeV]",  "", 12,   0, 600,  &(myEvent.Mlb_leadb)); // binning 12 or 24
    float Mlbbins[] = {0,50,100,150,200,250,350};
    AddVariable("Mlb_leadb_bin1", "M_lb [GeV]",  "", 6,  Mlbbins,  &(myEvent.Mlb_leadb)); // binning 12 or 24
    float Mlbbins2[] = {0,50,150,250,350};
    AddVariable("Mlb_leadb_bin2", "M_lb [GeV]",  "", 4,   Mlbbins2,  &(myEvent.Mlb_leadb)); // binning 12 or 24
    
    AddVariable("MT", "M_T [GeV]",  "", 30,   0, 300,  &(myEvent.mt_met_lep)); // use less bins
    float MTbins[] = {0,10,20,30,40,50,60,70,80,90,100,120,150,200,250,300};
    AddVariable("MT_bin", "M_T [GeV]",  "", 15, MTbins,  &(myEvent.mt_met_lep)); // use less bins
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 300,  &(myEvent.pfmet));
    AddVariable("MET_bin", "MET [GeV]",  "", 6,   0, 300,  &(myEvent.pfmet));
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W));
    AddVariable("MT2W_bin", "MT2W",  "", 6,   0, 300,  &(myEvent.MT2W));
    
    AddVariable("DeltaRlj", "Min(#DeltaR(j,l))",  "", 60,   0, 6,  &(myEvent.DeltaRlj));
    /*
    AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets));
    AddVariable("Nbjets", "N_{b-jets}",  "", 6,   0, 5,  &(myEvent.ngoodbtags));
    */
    //AddVariable("SRBins", "Signal Region box",  "", 11,   -5, 5,  &(onTheFlyVariables.SRbins));
    
    
    //AddVariable("x", "whatever2", "", 20, -40,  40,  &(onTheFlyVariables.someCoolVariable));
    // ...

    ///*
    AddProcessClass("OTHER", "other processes", "background", kGreen);
    	AddDataset("ST_s", "OTHER" ,-1, 10.32*0.3272);
    	AddDataset("ST_t-atop", "OTHER", -1, 80.95*0.3272 );
    	AddDataset("ST_t-top", "OTHER", -1, 136.02*0.3272 );
    	AddDataset("ST_tW-atop", "OTHER", -1, 35.85 );
    	AddDataset("ST_tW-top", "OTHER", -1, 35.85 );
    	AddDataset("TTW_ln", "OTHER", -1, 0.70*0.3);
    	AddDataset("TTW_qq", "OTHER", -1, 0.70*0.675);
    	AddDataset("TTZ_ll", "OTHER", -1, 0.62*0.2);
    	AddDataset("TTZ_qq", "OTHER", -1, 0.62*0.7);
    	AddDataset("WW_aMC", "OTHER", -1, 48.4);
    	AddDataset("ZZ", "OTHER", -1, 5.4);
    	AddDataset("WZ", "OTHER", -1, 48.4);
    	AddDataset("QCD_HT200to300", "OTHER", -1,  1717000);
    	AddDataset("QCD_HT300to500", "OTHER", -1, 351300);
    	AddDataset("QCD_HT500to700", "OTHER", -1, 31630 );
    	AddDataset("QCD_HT700to1000", "OTHER", -1, 6802);
    	AddDataset("QCD_HT1000to1500", "OTHER", -1, 1206);
    	AddDataset("QCD_HT1500to2000", "OTHER", -1,  120.4 );
    	AddDataset("QCD_HT2000toInf", "OTHER", -1,  25.24 );
    
    AddProcessClass("TT_1l", "TT+jets - 1l", "background", kRed);
    	AddDataset("TTjets_M5", "TT_1l", -1, 831.76);
    AddProcessClass("TT_2l", "TT+jets - 2l", "background", kCyan);
//*/
   
    //AddProcessClass("DY", "Drell-Yan", "background", kYellow);
    //	AddDataset("DYJetsToNuNu", "DY", -1, 6025.2);
    
    AddProcessClass("WJets", "W+jets", "background", kBlue);
    AddProcessClass("WJets_HF", "W+jets", "background", kBlue);
    //AddDataset("Wjets_aMC", "WJets", -1, 61466);
    ///*
	AddDataset("WJetsToLNu_HT-100To200", "WJets", -1, 1345*1.21);
    	AddDataset("WJetsToLNu_HT-200To400", "WJets", -1, 359.7*1.21);
    	AddDataset("WJetsToLNu_HT-400To600", "WJets", -1, 48.91*1.21);
    	AddDataset("WJetsToLNu_HT-600To800", "WJets", -1, 12.05*1.21);
    	AddDataset("WJetsToLNu_HT-800To1200", "WJets", -1, 5.501*1.21);
    	AddDataset("WJetsToLNu_HT-1200To2500", "WJets", -1, 1.329*1.21);
    	AddDataset("WJetsToLNu_HT-2500ToInf", "WJets", -1, 0.03216*1.21);
    	//*/
    
   /*
    AddProcessClass("QCD", "QCD", "background", kYellow);
    	AddDataset("QCD_HT200to300", "QCD", -1,  1717000);
    	AddDataset("QCD_HT300to500", "QCD", -1, 351300);
    	AddDataset("QCD_HT500to700", "QCD", -1, 31630 );
    	AddDataset("QCD_HT700to1000", "QCD", -1, 6802);
    	AddDataset("QCD_HT1000to1500", "QCD", -1, 1206);
    	AddDataset("QCD_HT1500to2000", "QCD", -1,  120.4 );
    	AddDataset("QCD_HT2000toInf", "QCD", -1,  25.24 );
    */
    
    AddProcessClass("DATA", "DATA", "data", kRed);
    	AddDataset("SingleMuon_RunD", "DATA");
    	//AddDataset("SingleElectron_RunD", "DATA");
    	AddDataset("SingleMuon_RunD_prompt", "DATA");
    	//AddDataset("SingleElectron_RunD_prompt", "DATA");
    
    // ----------------------------
    // Region definition
    // ----------------------------


    //AddRegion("preselection", "Preselection", &goesInPreselection);
    //AddRegion("BaselineSelection", "BaselineSelection", &goesInBaselineSearchSR);
    // ...
    AddRegion("CR_MET250_lowMT","",&CR_MET250_lowMT);
    AddRegion("CR0b_MET250_lowMT","",&CR0b_MET250_lowMT);
    AddRegion("CR_MET250_lowMT_3j","",&CR_MET250_lowMT_3j);
    AddRegion("CR0b_MET250_lowMT_3j","",&CR0b_MET250_lowMT_3j);
    AddRegion("CR0b_presel", "CR - == 0b - >=4j", &CR0b_presel);
    AddRegion("CR0b_presel_MET50", "CR - == 0b - >=4j - MET >50", &CR0b_presel_MET50);
    AddRegion("CR0b_presel_MET100", "CR - == 0b - >=4j - MET >100", &CR0b_presel_MET100);
    AddRegion("CR0b_presel_MET150", "CR - == 0b - >=4j - MET >150", &CR0b_presel_MET150);
    AddRegion("CR0b_presel_MET200", "CR - == 0b - >=4j - MET >200", &CR0b_presel_MET200);
    AddRegion("CR0b_presel_MET250", "CR - == 0b - >=4j - MET >250", &CR0b_presel_MET250);
    AddRegion("CR0b_presel_MET300", "CR - == 0b - >=4j - MET >250", &CR0b_presel_MET300);
    AddRegion("CR0b_presel_MT2Wtail", "CR - == 0b - >=4j - MT2W>200", &CR0b_presel_MT2Wtail);
    AddRegion("CR0b_presel_MT2W200", "CR - == 0b - >=4j - 200<MT2W<250 ", &CR0b_presel_MT2W200);
    AddRegion("CR0b_presel_MT2W250", "CR - == 0b - >=4j - MT2W>250", &CR0b_presel_MT2W250);
    AddRegion("CR0b_presel_MTpeak", "CR - == 0b - >=4j - MTpeak", &CR0b_presel_MTpeak);
    AddRegion("CR0b_presel_MTtail", "CR - == 0b - >=4j - MTtail", &CR0b_presel_MTtail);
    AddRegion("CR0b_presel_MTtail_80", "CR - == 0b - >=4j - MT>80", &CR0b_presel_MTtail_80);
    AddRegion("CR0b_presel_MTtail_90", "CR - == 0b - >=4j - MT>90", &CR0b_presel_MTtail_90);
    AddRegion("CR0b_presel_MTtail_100", "CR - == 0b - >=4j - MT>100", &CR0b_presel_MTtail_100);
    AddRegion("CR0b_presel_MTtail_110", "CR - == 0b - >=4j - MT>110", &CR0b_presel_MTtail_110);
    AddRegion("CR0b_presel_MTtail_120", "CR - == 0b - >=4j - MT>120", &CR0b_presel_MTtail_120);
    AddRegion("CR0b_presel_MTtail_130", "CR - == 0b - >=4j - MT>130", &CR0b_presel_MTtail_130);
    AddRegion("CR0b_presel_MTtail_140", "CR - == 0b - >=4j - MT>140", &CR0b_presel_MTtail_140);
    AddRegion("CR0b_presel_MTtail_150", "CR - == 0b - >=4j - MT>150", &CR0b_presel_MTtail_150);
    AddRegion("CR0b_presel_MTtail_80_ex", "CR - == 0b - >=4j - MT>80", &CR0b_presel_MTtail_80_ex);
    AddRegion("CR0b_presel_MTtail_90_ex", "CR - == 0b - >=4j - MT>90", &CR0b_presel_MTtail_90_ex);
    AddRegion("CR0b_presel_MTtail_100_ex", "CR - == 0b - >=4j - MT>100", &CR0b_presel_MTtail_100_ex);
    AddRegion("CR0b_presel_MTtail_110_ex", "CR - == 0b - >=4j - MT>110", &CR0b_presel_MTtail_110_ex);
    AddRegion("CR0b_presel_MTtail_120_ex", "CR - == 0b - >=4j - MT>120", &CR0b_presel_MTtail_120_ex);
    AddRegion("CR0b_presel_MTtail_130_ex", "CR - == 0b - >=4j - MT>130", &CR0b_presel_MTtail_130_ex);
    AddRegion("CR0b_presel_2j_MTpeak", "CR - == 0b - >=4j - MTpeak", &CR0b_presel_2j_MTpeak);
    AddRegion("CR0b_presel_2j_MTtail", "CR - == 0b - >=4j - MTtail", &CR0b_presel_2j_MTtail);
    AddRegion("CR0b_presel_3j_MTpeak", "CR - == 0b - >=4j - MTpeak", &CR0b_presel_3j_MTpeak);
    AddRegion("CR0b_presel_3j_MTtail", "CR - == 0b - >=4j - MTtail", &CR0b_presel_3j_MTtail);
    AddRegion("CR0b_presel_4j_MTpeak", "CR - == 0b - >=4j - MTpeak", &CR0b_presel_4j_MTpeak);
    AddRegion("CR0b_presel_4j_MTtail", "CR - == 0b - >=4j - MTtail", &CR0b_presel_4j_MTtail);

    AddChannel("muon", "#mu channel", &goesInSingleMuonChannel);
    /*
    AddChannel("electron", "electron channel", &goesInSingleElecChannel);
    AddChannel("electronEC", "electron (EC) channel", &goesInSingleElecECChannel);
    AddChannel("singleLepton", "single lepton channel", &goesInSingleLeptonChannel);
    */
    // ...

    SetLumi(1226);
    //SetLumi(578.3);

    Create1DHistos();
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{

   /*
   cout << "testing run" << endl;
   cout << "nr. of selected leptons: " << myEvent.numberOfSelectedLeptons << endl;
   cout << "nr. of init events: " << myEvent.totalNumberOfInitialEvent << endl;


    cout<<myEvent.mt_met_lep<<endl;
    cout<<"sel: "<<goesInPreselection()<<endl;
    */
    //myEvent.crossSection  = CrossSection(currentDataset);
    myEvent.crossSection = GetDatasetCrossSection(currentDataset);
    //myEvent.totalNumberOfInitialEvent = 
    /*
    cout<<"xsextion : "<<myEvent.crossSection<<endl;
    cout<<"tot weight: "<<myEvent.totalNumberOfInitialEvent<<endl;
    cout<<"weight: "<<myEvent.mc_weight<<endl;
    cout<<"met: "<<myEvent.pfmet<<endl;
    */
    // Compute on the fly variables if needed

    //ComputeOnTheFlyVariables();

    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    // Compute weight for current event

    double weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;
    //dumpSelection();
    //cout<<"in presel: "<<goesInPreselection()<<endl;

    double weight     = weightLumi;
    //cout<<myEvent.crossSection<<" "<<GetLumi()<<" " <<myEvent.mc_weight<<" "<<myEvent.totalNumberOfInitialEvent<<endl;
    //if(goesInPreselection() && goesInSingleMuonChannel()) cout<<weight<<endl;
    if (currentProcessType == "data") weight = 1.0;
    //if (currentProcessClass == "WJets") weight *= 2.5;

    // Fill this event in the histo collections
    //cout<<"weight = "<<weight<<" "<<myEvent.crossSection<<" "<<GetLumi()<<" "<< myEvent.mc_weight<<" "<<myEvent.totalNumberOfInitialEvent<<endl;
    // Split 1-lepton ttbar and 2-lepton ttbar
    //if (ttbarDatasetToBeSplitted && (myEvent.numberOfGeneratedLeptons == 2))
    if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons >= 2)){
	currentProcessClass = "TT_2l";
	/*
	if(myEvent.hasGenInfo && currentProcessType != "data"){
		weight*=topPtWeight();
		cout<<"weight modified"<<endl;
		}
		*/
    }
    myEvent.DeltaRlj = DeltaRlj();
  		


    //cout<<"evts"<<endl;
    if (currentProcessClass == "WJets" && genInfo()){
   	 currentProcessClass = "WJets_HF";
    }

   AutoFillProcessClass(currentProcessClass, weight);
    //AutoFillProcessClass(currentProcessClass, 0.001);
}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    // ######################
    //  Plot configuration and production
    // ######################

    // Schedule plots

    SchedulePlots("1DSuperimposed");
    SchedulePlots("1DStack");
    SchedulePlots("1DDataMCComparison");
    
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
    WritePlots("./plots/");

    // ######################
    //  Tables and other stuff
    // ######################

    TableDataMC(this,{"preselection"},"muon", "").Print("table.tab",4);

}



