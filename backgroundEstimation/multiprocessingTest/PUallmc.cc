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
//#define USE_AK8_JETS


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
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/0801/";
    //babyTuplePath = "/home/mjansova/PyROOF/";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 5;

    float customBinning[3] = {0.0, 200.0, 600.0};
    float customBinningMET[5] = {50.0, 250.0, 350.0, 450.0, 600.0};

    AddVariable("MT", "M_T [GeV]",  "", 20,   0, 400,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 25, 0, 500,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 25, 0, 500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("nvertex", "nvertex",  "", 40, 0, 40,  &(myEvent.nvertex), "noUnderflowInFirstBin");
// ...

   //ttJets
   AddProcessClass("mc", "mc", "background", kPink);
        AddDataset("TTjets_M5", "mc", 0, 831.76);
        AddDataset("ST_s", "mc", 0, 3.38 );
        AddDataset("ST_t-atop", "mc", 0, 26.49 );
        AddDataset("ST_t-top", "mc", 0, 44.51 );
        AddDataset("ST_tW-atop", "mc", 0, 35.85 );
        AddDataset("ST_tW-top", "mc", 0 , 35.85 );
        AddDataset("WJetsToLNu_HT-100To200", "mc", 0, 1345*1.21);
        AddDataset("WJetsToLNu_HT-200To400", "mc", 0, 359.7*1.21 );
        AddDataset("WJetsToLNu_HT-400To600", "mc", 0, 48.91*1.21 );
        AddDataset("WJetsToLNu_HT-600To800", "mc", 0, 12.05*1.21);
        AddDataset("WJetsToLNu_HT-800To1200", "mc", 0, 5.501*1.21 );
        AddDataset("WJetsToLNu_HT-1200To2500", "mc", 0, 1.329*1.21);
        AddDataset("WJetsToLNu_HT-2500ToInf", "mc", 0, 0.03216*1.21);

   AddProcessClass("data", "data", "data", kPink);
        AddDataset("data", "data", -1, -1);


    AddRegion("defaultSearchAllCuts", "defaultSearchAllCuts", &defaultSearchAllCuts);   
    AddRegion("defaultSearchAllCutsFakes", "defaultSearchAllCutsFakes", &defaultSearchAllCutsFakes);  //preselection 
    AddRegion("noCuts", "noCuts", &noCuts);   
    AddRegion("preselectionPUBin1", "preselectionPUBin1", &preselectionPUBin1);   
    AddRegion("preselectionPUBin2", "preselectionPUBin2", &preselectionPUBin2);   
    AddRegion("preselectionPUBin3", "preselectionPUBin3", &preselectionPUBin3);   
    AddRegion("preselectionPUBin4", "preselectionPUBin4", &preselectionPUBin4);   
    AddRegion("preselectionPUBin5", "preselectionPUBin5", &preselectionPUBin5);   
    AddRegion("preselectionPUBin6", "preselectionPUBin6", &preselectionPUBin6);   
    AddRegion("preselectionPUBin7", "preselectionPUBin7", &preselectionPUBin7);   
    AddRegion("preselectionPUBin8", "preselectionPUBin8", &preselectionPUBin8);   

    /*AddChannel("muon", "#mu channel", &muonChannelSelector);
    AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("tau", "#tau channel", &tauChannelSelector);*/
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    // ...



    SetLumi(2440.);

    Create1DHistos();
    Add2DHisto("MET","nvertex"); //@MJ@ TODO
    Add2DHisto("MT","nvertex"); //@MJ@ TODO
    Add2DHisto("MT2W","nvertex"); //@MJ@ TODO

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
    //misidentifiedLepton(); 

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
    //

    //MC PU histo
    //
    //TH1F *h1 = new TH1F("mc_pu", "mc_pu", 53, 0, 4.4);

    // Schedule plots

    SchedulePlots("1DSuperimposed");
    SchedulePlots("1DStack");
    SchedulePlots("1DDataMCComparison");
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
    WritePlots("./plots6/");

   //vector<string> newBLregionsAk8 = {"Default2OneAk8Bin1", "Default2OneAk8Bin2", "Default2OneAk8Bin3"};
   //TableBackgroundSignal(this, newBLregionsAk8,"combinedChannel" ).Print("newBLak82fb.tab", 4);
   //TableBackgroundSignal(this, newBLregionsAk8,"combinedChannel" ).PrintLatex("newBLak82fb.tex", 4);
}

