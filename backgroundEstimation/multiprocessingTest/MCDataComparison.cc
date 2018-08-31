#include <vector>
#include <iostream>
using namespace std;

#include "TTree.h"
//babyEvent myEvent;




//-----------------------------------
//-- Load only necessary branches
//-----------------------------------

//#define USE_GEN_INFO
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
#include "OnTheFlyVariables.h"

// Should be called after selectionDefinition 'cause it need myEvent varible ...
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"



bool muonChannelSelector() { return true; }

// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");
    cout << "non correct version" << endl;
    //babyTuplePath = "./";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_missingInfo/";
    // babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_all_3j_lowpt_dphi_v2/"; 
    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_all_2j/";

    totalNumberOfWorkers = 8;

    AddVariable("MT", "M_{T} [GeV]",  "", 23,   0, 460,  &(myEvent.mt_met_lep), "logY");
    AddVariable("MET", "MET [GeV]",  "", 30,   0, 300,  &(myEvent.pfmet), "logY");
    AddVariable("MT2W", "MT2W",  "", 25,   80, 380,  &(myEvent.MT2W), "logY");
    ///*
    AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets));
    AddVariable("Nbjets", "N_{b-jets}",  "", 6,   0, 5,  &(myEvent.ngoodbtags));
    AddVariable("Mlb_leadb", "M_lb [GeV]",  "", 40,   0, 600,  &(myEvent.Mlb_leadb));
    AddVariable("QGtag_leadJet","qg-tag (lead jet)", "", 20, 0, 1, &(onTheFlyVariables.QGtag_leadJet));
    AddVariable("met_sig", "MET significance", "", 100, 0, 50, &(myEvent.met_sig));
    //*/
    //AddVariable("SRBins", "Signal Region box",  "", 11,   -5, 5,  &(onTheFlyVariables.SRbins));
    
    
    //AddVariable("x", "whatever2", "", 20, -40,  40,  &(onTheFlyVariables.someCoolVariable));
    // ...
    
    AddProcessClass("TT_1l", "TT+jets - 1l", "background", kRed);
    	AddDataset("TTjets_M5", "TT_1l", -1, 831.76);
    AddProcessClass("TT_2l", "TT+jets - 2l", "background", kCyan);

    AddProcessClass("DY", "Drell Yan", "background", kGreen+1);
    	AddDataset("DYJetsToNuNu","DY", -1, 11433 );

    AddProcessClass("WJets", "W+jets", "background", kBlue);
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
    

    AddProcessClass("SingleTop", "Single top", "background", kGreen);
    	AddDataset("ST_s", "SingleTop" ,-1, 10.32*0.3272);
    	AddDataset("ST_t-atop", "SingleTop", -1, 80.95*0.3272 );
    	AddDataset("ST_t-top", "SingleTop", -1, 136.02*0.3272 );
    	AddDataset("ST_tW-atop", "SingleTop", -1, 35.85 );
    	AddDataset("ST_tW-top", "SingleTop", -1, 35.85 );
    /*
    AddProcessClass("ttV", "tt+boson", "background", kMagenta);
    	AddDataset("TTW_ln", "ttV", -1, 0.70*0.3);
    	AddDataset("TTW_qq", "ttV", -1, 0.70*0.675);
    	AddDataset("TTZ_ll", "ttV", -1, 0.62*0.2);
    	AddDataset("TTZ_qq", "ttV", -1, 0.62*0.7);

    AddProcessClass("VV", "di-bosons", "background", kMagenta+1);
    	AddDataset("WW_aMC", "VV", -1, 48.4);
    	AddDataset("ZZ", "VV", -1, 5.4);
    	AddDataset("WZ", "VV", -1, 48.4);
    */

    AddProcessClass("DATA", "DATA", "data", kRed);
    	AddDataset("SingleMuon_RunD", "DATA");
    	AddDataset("SingleElectron_RunD", "DATA");
    	AddDataset("SingleMuon_RunD_prompt", "DATA");
    	AddDataset("SingleElectron_RunD_prompt", "DATA");
    
    // ...

   // /*
    AddRegion("CR0b_presel_2j","0 bjet && == 2j", &CR0b_presel_2j);
    AddRegion("CR0b_presel_3j","0 bjet && >= 3j", &CR0b_presel_3j);
    AddRegion("CR0b_presel_4j","0 bjet && >= 4j", &CR0b_presel_4j);
    AddRegion("CR0b_presel_2j_MET100","0 bjet && == 2j && MET>100", &CR0b_presel_2j_MET100);
    AddRegion("CR0b_presel_3j_MET100","0 bjet && >= 3j && MET>100", &CR0b_presel_3j_MET100);
    AddRegion("CR0b_presel_4j_MET100","0 bjet && >= 4j && MET>100", &CR0b_presel_4j_MET100);
    //*/
    AddRegion("preselection", "Preselection", &goesInPreselection);
    //AddRegion("BaselineSelection", "BaselineSelection", &goesInBaselineSearchSR);
    // ...

    //AddChannel("def", "#mu channel", &muonChannelSelector);
    AddChannel("muon", "muon", &goesInSingleMuonChannel);
    //AddChannel("electron", "electron channel", &goesInSingleElecChannel);
    //AddChannel("electronEC", "endcap electron", &goesInSingleElecECChannel);
    //AddChannel("electronBarrel", "barrel electron", &goesInSingleElecBarrelChannel);
    //AddChannel("singleLepton", "single lepton channel", &goesInSingleLeptonChannel);
    // ...

    SetLumi(1226); // 1285.201
    //SetLumi(578.3); //595.182

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

    ComputeOnTheFlyVariables();

    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    // Compute weight for current event

    double weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;

    double weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    //if (currentProcessClass == "WJets") weight *= 2.5;
   
    if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons >= 2))
	currentProcessClass = "TT_2l";

    // Fill this event in the histo collections
    //cout<<"weight = "<<weight<<" "<<myEvent.crossSection<<" "<<GetLumi()<<" "<< myEvent.mc_weight<<" "<<myEvent.totalNumberOfInitialEvent<<endl;
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

/*
    TableDataMC(this,{"CR0b_presel_2j", "CR0b_presel_3j", "CR0b_presel_4j", "CR0b_presel_2j_MET100", "CR0b_presel_3j_MET100", "CR0b_presel_4j_MET100", "preselection"},"muon", "").Print("table.tab",4);
    TableDataMC(this,{"CR0b_presel_2j", "CR0b_presel_3j", "CR0b_presel_4j", "CR0b_presel_2j_MET100", "CR0b_presel_3j_MET100", "CR0b_presel_4j_MET100"},"muon", "").PrintLatex("table.tex",4);
*/
}



