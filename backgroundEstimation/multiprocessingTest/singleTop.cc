#include <vector>
#include <iostream>
using namespace std;

#include "TTree.h"
//babyEvent myEvent;




//-----------------------------------
//-- Load only necessary branches
//-----------------------------------

#define USE_GEN_INFO
#define USE_GEN_INFO_EXT
//#define USE_SKIMMING_VAR
#define USE_JETS
#define USE_JETS_EXT
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
#define USE_AK8_JETS


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



bool muonChannelSelector() 
{ 
    return (abs(myEvent.leadingLeptonId) == 13) ? true : false; 
}
bool electronChannelSelector() 
{   
    return (abs(myEvent.leadingLeptonId) == 11) ? true : false; 
}
bool tauChannelSelector() 
{   
    return (abs(myEvent.leadingLeptonId) == 15) ? true : false; 
}
bool combinedChannelSelector() 
{ 
    return true; 
}

uint64_t counting;
float MET2;
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_all_3j_lowpt_dphi_v2";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Oct21_postpatch8";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 5;
 


    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 500,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 500,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 50,   0, 500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("MET2", "MET2 [GeV]",  "", 10,   0, 500,  &(MET2), "noUnderflowInFirstBin");
    //AddVariable("genMET", "genMET",  "", 50,   0, 500,  &(myEvent.genMET), "noUnderflowInFirstBin");
    //AddVariable("genRecoDiff", "genRecoDiff",  "", 50,   0, 500,  &(onTheFlyVariables.m_genRecoDiff), "noUnderflowInFirstBin");
    AddVariable("lepMumId", "lepMumId",  "", 100,   -50, 50,  &(onTheFlyVariables.m_lepMumId), "noUnderflowInFirstBin");
    AddVariable("lepMumPt", "lepMumPt",  "", 100,   0, 500,  &(onTheFlyVariables.m_lepMumPt), "noUnderflowInFirstBin");
    AddVariable("lepMumMumId", "lepMumMumId",  "", 100,   -50, 50,  &(onTheFlyVariables.m_lepMumMumId), "noUnderflowInFirstBin");
    AddVariable("lepMumMumPt", "lepMumMumPt",  "", 50,   0, 500,  &(onTheFlyVariables.m_lepMumMumPt), "noUnderflowInFirstBin");
    AddVariable("lepPt", "lepPt",  "", 50,   0, 500,  &(onTheFlyVariables.m_lepPt), "noUnderflowInFirstBin");
    AddVariable("neutrinoPt", "neutrinoPt",  "", 50,   0, 500,  &(onTheFlyVariables.m_neutrinoPt), "noUnderflowInFirstBin");
    


    //AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets), "noUnderflowInFirstBin");
    //AddVariable("Nbjets", "N_{b-jets}",  "", 5,   0, 5,  &(myEvent.ngoodbtags), "noUnderflowInFirstBin");
// ...

    // background: cat 1
    AddProcessClass("sChannel", "sChannel", "background", kBlue-3);
    	AddDataset("ST_s", "sChannel", 0, 3.38 );

    AddProcessClass("tChannel", "tChannel", "background", kViolet+5);
    	AddDataset("ST_t-atop", "tChannel", 0, 26.49 );
    	AddDataset("ST_t-top", "tChannel", 0, 44.51 );

    AddProcessClass("tW", "tW", "background", kPink+7 );
    	AddDataset("ST_tW-atop", "tW", 0, 35.85 );
        AddDataset("ST_tW-top", "tW", 0 , 35.85 );

   //background cat 2
    AddProcessClass("st0Lep", "st0Lep", "background", kGray);
    	AddDataset("ST_s_2", "st0Lep", 0, 3.38 );
        AddDataset("ST_t-atop_2", "st0Lep", 0, 26.49 );
    	AddDataset("ST_t-top_2", "st0Lep", 0, 44.51 );
    	AddDataset("ST_tW-atop_2", "st0Lep", 0, 35.85 );
    	AddDataset("ST_tW-top_2", "st0Lep", 0 , 35.85 );

    AddProcessClass("st1Lep", "st1Lep", "background", kGray+2);
    AddProcessClass("st2Lep", "st2Lep", "background", kGray+3);
   
    //background cat 3
    AddProcessClass("tW0Lep", "tW0Lep", "background", kOrange+6);
    	AddDataset("ST_tW-atop_2_2", "tW0Lep", 0, 35.85 );
    	AddDataset("ST_tW-top_2_2", "tW0Lep", 0 , 35.85 );

    AddProcessClass("tW1Lep", "tW1Lep", "background", kOrange+2);
    AddProcessClass("tW2Lep", "tW2Lep", "background", kYellow-2);

   //ttJets
   AddProcessClass("TTJets", "TT+jets", "background", kSpring-8);
        AddDataset("TTjets_M5", "TTJets", 0, 831.76);
   AddProcessClass("TT_2l", "TT+jets - 2l", "background", kSpring-9);

    AddRegion("defaultSearch", "defaultSearch", &defaultSearch);
    AddRegion("defaultSearchAllCutsLowMT2W", "defaultSearchAllCutsLowMT2W", &defaultSearchAllCutsLowMT2W);
    AddRegion("defaultSearchMETCut", "defaultSearchMETCut", &defaultSearchMETCut);
    AddRegion("defaultSearchMTCut", "defaultSearchMTCut", &defaultSearchMTCut);
    AddRegion("defaultSearchMT2WCut", "defaultSearchMT2WCut", &defaultSearchMT2WCut);
    AddRegion("defaultSearchAllCuts", "defaultSearchAllCuts", &defaultSearchAllCuts);

    /*AddChannel("muon", "#mu channel", &muonChannelSelector);
    AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("tau", "#tau channel", &tauChannelSelector);*/
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    // ...

    SetLumi(2000.);

    Create1DHistos();

}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{

    counting++;
    
    //cout<<myEvent.mt_met_lep<<endl;
    //cout << " sel: "<< goesInPreselection() << endl;
    myEvent.crossSection  = GetDatasetCrossSection(currentDataset);
    //cout<<"xsextion : "<<myEvent.crossSection<<endl;
    //cout<<"tot weight: "<<myEvent.totalNumberOfInitialEvent<<endl;
    //cout<<"weight: "<<myEvent.mc_weight<<endl;
    // Compute on the fly variables if needed

     
    ComputeOnTheFlyVariables();
    
    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    //cout << "process class: " << currentProcessClass << " number of gen. leptons:  " << myEvent.numberOfGeneratedLeptons << endl;
 
    if (currentProcessClass == "TTJets" && (myEvent.numberOfGeneratedLeptons == 2))
    {
        currentProcessClass = "TT_2l";
    }
    
    if (currentProcessClass == "st0Lep" && (myEvent.numberOfGeneratedLeptons > 0))
    {
         currentProcessClass = myEvent.numberOfGeneratedLeptons == 1 ? "st1Lep" : "st2Lep";
        // cout << "now process is: " << currentProcessClass << endl;
    }
    
    if (currentProcessClass == "tW0Lep" && (myEvent.numberOfGeneratedLeptons > 0))
    {
         currentProcessClass = myEvent.numberOfGeneratedLeptons == 1 ? "tW1Lep" : "tW2Lep";
        // cout << "now process is: " << currentProcessClass << endl;
    }


    MET2 = myEvent.pfmet;
    //calculate variables
    findLeptonOrigin();
 
     // Compute weight for current event

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;
    float weight     = weightLumi;

    if (currentProcessType == "data") weight = 1.0;

    // Fill this event in the histo collections

    AutoFillProcessClass(currentProcessClass, weight);
    //cout << "process " << currentProcessClass << " weight " << weight << endl;
    //AutoFillProcessClass(currentProcessClass, 0.001);
    
    if(counting % 1000 == 0)
    {
        cout << counting << endl;

    }
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

   //vector<string> newBLregionsAk8 = {"Default2OneAk8Bin1", "Default2OneAk8Bin2", "Default2OneAk8Bin3"};
   //TableBackgroundSignal(this, newBLregionsAk8,"combinedChannel" ).Print("newBLak82fb.tab", 4);
   //TableBackgroundSignal(this, newBLregionsAk8,"combinedChannel" ).PrintLatex("newBLak82fb.tex", 4);
}

