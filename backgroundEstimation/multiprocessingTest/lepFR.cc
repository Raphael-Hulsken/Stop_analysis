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
    return (abs(myEvent.lep1_pdgid) == 13) ? true : false; 
}
bool electronChannelSelector() 
{   
    return (abs(myEvent.lep1_pdgid) == 11) ? true : false; 
}
bool tauChannelSelector() 
{   
    return (abs(myEvent.lep1_pdgid) == 15) ? true : false; 
}
bool combinedChannelSelector() 
{ 
    return true; 
}

uint64_t counting;
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/singleTopEstim";
    //babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/0801/";
    //babyTuplePath = "/home/mjansova/PyROOF/";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 6;

    float customBinning[3] = {0.0, 200.0, 600.0};
    float customBinningMET[5] = {50.0, 250.0, 350.0, 450.0, 600.0};

    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 500,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 50, 0, 600,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 50, 0, 600,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("FakesReals", "FakesReals",  "", 3, 0, 3,  &(onTheFlyVariables.m_fakeLep), "noUnderflowInFirstBin");
    AddVariable("lepDr", "lepDr",  "", 10, 0, 1,  &(onTheFlyVariables.m_lepDr), "noUnderflowInFirstBin");
    // @MJ@ TODO number of verices
    //AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets), "noUnderflowInFirstBin");
    //AddVariable("Nbjets", "N_{b-jets}",  "", 5,   0, 5,  &(myEvent.ngoodbtags), "noUnderflowInFirstBin");
// ...

    // background: cat 1
    //AddProcessClass("tW", "tW", "background", kPink+7 );
    //	AddDataset("ST_tW-atop", "tW", 0, 35.85 );
    //    AddDataset("ST_tW-top", "tW", 0 , 35.85 );
   //ttJets
   AddProcessClass("TTJets", "TT+jets - 1l", "background", kSpring-8);
        AddDataset("TTjets_M5", "TTJets", 0, 831.76);
   AddProcessClass("TT_2l", "TT+jets - 2l", "background", kSpring-9);
   AddProcessClass("TT_0l", "TT+jets - 0l", "background", kSpring-6);

//   AddProcessClass("TTJetsi", "TT+jets - 1li", "background", kSpring-8);
//        AddDataset("debuggingOutput.tmp", "TTJetsi", 0, 831.76);

    AddRegion("defaultSearchAllCuts", "defaultSearchAllCuts", &defaultSearchAllCuts);   
    AddRegion("defaultSearchAllCutsFakes", "defaultSearchAllCutsFakes", &defaultSearchAllCutsFakes);   
    AddRegion("noCuts", "noCuts", &noCuts);   
    AddRegion("defaultSearchAllCutsLooser", "defaultSearchAllCutsLooser", &defaultSearchAllCutsLooser);   

    /*AddChannel("muon", "#mu channel", &muonChannelSelector);
    AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("tau", "#tau channel", &tauChannelSelector);*/
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    AddChannel("electronChannel","e-channel",&electronChannelSelector);
    AddChannel("muonChannel","#mu-channel",&muonChannelSelector);
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

     
    //ComputeOnTheFlyVariables();
    
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
    
    if (currentProcessClass == "TTJets" && (myEvent.numberOfGeneratedLeptons == 0))
    {
        currentProcessClass = "TT_0l";
    }
    

    //calculate variables
    //findLeptonOrigin();    
    //countLepDeltas(); 
    //closeLepAndJet();
    //countAllPt();
    misidentifiedLepton(); 

 // Compute weight for current event
 //
    //if(stBackgroundMET501jPt20())
    // cout << "MT2W: " << myEvent.MT2W << endl;

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
    //SchedulePlots("2D");

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

