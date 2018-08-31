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
Double_t Z;
Double_t MT_sig;
Double_t MT_sigToMT;
Double_t met_sigTomet;

TFile *fle = new TFile("ndivisionMTmcprofile.root");
//fle->ls();
TH1D * h1 = (TH1D*)fle->Get("pileup");

// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/1701";
    //babyTuplePath = "/home/mjansova/PyROOF/";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 5;

    float customBinning[3] = {0.0, 200.0, 600.0};
    float customBinningMET[5] = {50.0, 250.0, 350.0, 450.0, 600.0};

    AddVariable("MT", "M_T [GeV]",  "", 15,   0, 375,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 20, 0, 500,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("met_sig", "MET significance",  "", 20, 0, 160,  &(myEvent.met_sig), "noUnderflowInFirstBin");
    AddVariable("MT_sig", "MT significance",  "", 20, 0, 160,  &(MT_sig), "noUnderflowInFirstBin");
    AddVariable("Z", "Z",  "", 20, 0, 5,  &(Z), "noUnderflowInFirstBin");
    AddVariable("MT_sigToMT", "MT_sig / MT",  "", 20, 0, 1,  &(MT_sigToMT), "noUnderflowInFirstBin");
    AddVariable("met_sigTomet", "met_sig / met",  "", 20, 0, 1,  &(met_sigTomet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 20, 0, 500,  &(myEvent.MT2W), "noUnderflowInFirstBin");

   //ttJets
   AddProcessClass("TTJets", "TT+jets - 1l", "background", kBlue+2);
        AddDataset("TTjets_M5", "TTJets", 0, 831.76);
   AddProcessClass("TT_2l", "TT+jets - 2l", "background", kAzure+1);
   AddProcessClass("TT_0l", "TT+jets - 0l", "background", kViolet+6);

    AddProcessClass("SingleTop", "Single top", "background",kPink+7 );
        AddDataset("ST_s", "SingleTop", 0, 3.38 );
        AddDataset("ST_t-atop", "SingleTop", 0, 26.49 );
        AddDataset("ST_t-top", "SingleTop", 0, 44.51 );
        AddDataset("ST_tW-atop", "SingleTop", 0, 35.85 );
        AddDataset("ST_tW-top", "SingleTop", 0 , 35.85 );

   AddProcessClass("WJets", "W+jets", "background", kGray+3);
        AddDataset("WJetsToLNu_HT-100To200", "WJets", 0, 1345*1.21);
        AddDataset("WJetsToLNu_HT-200To400", "WJets", 0, 359.7*1.21 );
        AddDataset("WJetsToLNu_HT-400To600", "WJets", 0, 48.91*1.21 );
        AddDataset("WJetsToLNu_HT-600To800", "WJets", 0, 12.05*1.21);
        AddDataset("WJetsToLNu_HT-800To1200", "WJets", 0, 5.501*1.21 );
        AddDataset("WJetsToLNu_HT-1200To2500", "WJets", 0, 1.329*1.21);
        AddDataset("WJetsToLNu_HT-2500ToInf", "WJets", 0, 0.03216*1.21);

   AddProcessClass("TTZ_lep", "TTZ to leps", "background", kRed+3);
        AddDataset("TTZ_ll", "TTZ_lep", 0, 0.62*0.2);


   AddProcessClass("bkg", "bkg", "background", kBlack);
        AddDataset("TTjets_M52", "bkg", 0, 831.76);
        AddDataset("ST_s2", "bkg", 0, 3.38 );
        AddDataset("ST_t-atop2", "bkg", 0, 26.49 );
        AddDataset("ST_t-top2", "bkg", 0, 44.51 );
        AddDataset("ST_tW-atop2", "bkg", 0, 35.85 );
        AddDataset("ST_tW-top2", "bkg", 0 , 35.85 );
        AddDataset("WJetsToLNu_HT-100To2002", "bkg", 0, 1345*1.21);
        AddDataset("WJetsToLNu_HT-200To4002", "bkg", 0, 359.7*1.21 );
        AddDataset("WJetsToLNu_HT-400To6002", "bkg", 0, 48.91*1.21 );
        AddDataset("WJetsToLNu_HT-600To8002", "bkg", 0, 12.05*1.21);
        AddDataset("WJetsToLNu_HT-800To12002", "bkg", 0, 5.501*1.21 );
        AddDataset("WJetsToLNu_HT-1200To25002", "bkg", 0, 1.329*1.21);
        AddDataset("WJetsToLNu_HT-2500ToInf2", "bkg", 0, 0.03216*1.21);
        AddDataset("TTZ_ll2", "bkg", 0, 0.62*0.2);

    AddProcessClass("225_25to150", "225_25to150", "background",  kViolet+6); //@MJ@ TODO -extract stop mass and compute cross section and use for weight, aleso first th tot number of events for one stop mass must be known
        AddDataset("225_25to150","225_25to150", 0, 36.3818);

    AddProcessClass("500-525-550_1to425-300to450", "500-525-550_1to425-300to450", "background",  kBlue+6); //@MJ@ TODO -extract stop mass and compute cross section and use for weight, aleso first th tot number of events for one stop mass must be known
        AddDataset("500-525-550_1to425-300to450","500-525-550_1to425-300to450", 0, 0.390303);

   AddProcessClass("600-950_1to450", "600-950_1to450", "background",  kBlue); //@MJ@ TODO -extract stop mass and compute cross section and use for weight, aleso first th tot number of events for one stop mass must be known
        AddDataset("600-950_1to450","600-950_1to450", 0, 0.0431418);

/*    AddProcessClass("T2tt_850_325", "T2tt_850_325", "background",  kMagenta+3);
        AddDataset("T2tt_850_325","T2tt_850_325", 0, 0.01896);

    AddProcessClass("T2tt_650_325", "T2tt_650_325", "background",  kViolet+6);
        AddDataset("T2tt_650_325","T2tt_650_325", 0, 0.107045);

    AddProcessClass("T2tt_500_325", "T2tt_500_325", "background",  kOrange+7);
        AddDataset("T2tt_500_325","T2tt_500_325", 0, 0.51848);

    AddProcessClass("T2tt_425_325", "T2tt_425_325", "background",  kGray);
        AddDataset("T2tt_425_325","T2tt_425_325", 0, 1.31169);
   
   //AddProcessClass("data", "data", "data", kBlack);
   //     AddDataset("data", "data", -1, -1);
*/

    AddRegion("defaultSearchAllCuts", "defaultSearchAllCuts", &defaultSearchAllCuts);   
    AddRegion("defaultSearchAllCutsMS", "defaultSearchAllCutsMS", &defaultSearchAllCutsMS);   
    AddRegion("defaultSearchAllCutsMSToMET", "defaultSearchAllCutsMSToMET", &defaultSearchAllCutsMSToMET);   
    AddRegion("defaultSearchAllCutsLooseMTMS", "defaultSearchAllCutsLooseMTMS", &defaultSearchAllCutsLooseMTMS);   
    AddRegion("defaultSearchAllCutsMTS", "defaultSearchAllCutsMTS", &defaultSearchAllCutsMTS);   
    AddRegion("defaultSearchAllCutsMTSToMT", "defaultSearchAllCutsMTSToMT", &defaultSearchAllCutsMTSToMT);   
    AddRegion("defaultSearchAllCutsExcMET", "defaultSearchAllCutsExcMET", &defaultSearchAllCutsExcMET); 
    AddRegion("defaultSearchAllCutsExcMT", "defaultSearchAllCutsExcMT", &defaultSearchAllCutsExcMT); 
    AddRegion("defaultSearchAllCutsExcMETLooseMT", "defaultSearchAllCutsExcMETLooseMT", &defaultSearchAllCutsExcMETLooseMT); 
    AddRegion("defaultSearchMT2WCut", "defaultSearchMT2WCut", &defaultSearchMT2WCut);   
    AddRegion("defaultSearchAllCutsZ", "defaultSearchAllCutsZ", &defaultSearchAllCutsZ);   
    AddRegion("defaultSearch", "defaultSearch", &defaultSearch);   

    /*AddChannel("muon", "#mu channel", &muonChannelSelector);
    AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("tau", "#tau channel", &tauChannelSelector);*/
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    //AddChannel("muon","#mu-channel",&muonChannelSelector);
   //AddChannel("electron","e-channel",&electronChannelSelector);
    // ...



    SetLumi(2440.);

    Create1DHistos();
    //Add2DHisto("MET","nvertex"); //@MJ@ TODO
    //Add2DHisto("MT","nvertex"); //@MJ@ TODO
    //Add2DHisto("MT2W","nvertex"); //@MJ@ TODO

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
    
    Double_t deltaPhi = myEvent.lep1_phi - myEvent.pfmet_phi;
    met_sigTomet = myEvent.met_sig / myEvent.pfmet;
    MT_sig = sqrt(2 * myEvent.lep1_pt * myEvent.met_sig * (1 - cos(deltaPhi) ));
    Z = (myEvent.mt_met_lep - 80.385)/(myEvent.mt_met_lep*met_sigTomet);
    MT_sigToMT = MT_sig / myEvent.mt_met_lep;
    met_sigTomet = myEvent.met_sig / myEvent.pfmet;
   

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
    //cout << "weight " << weightLumi << endl;
 
    float weightRew;
    if(myEvent.puIntime <= h1->GetNbinsX())
    {
        //cout << "it weight " << myEvent.puIntime << endl;
        //cout << "reweighting by: " << h1->GetBinContent(myEvent.puIntime) << endl;
        weightRew = (h1->GetBinContent(myEvent.puIntime)) * myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;
        //cout << "rew weight " << weightRew << endl;
    }
    else
       weightRew = weightLumi;
    
    float weight     = weightLumi;
    //float weight     = weightRew;

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
    WritePlots("./plots11/");

   vector<string> regions = {"defaultSearchAllCuts", "defaultSearchAllCutsMS", "defaultSearchAllCutsMSToMET", "defaultSearchAllCutsLooseMTMS", "defaultSearchAllCutsMTS", "defaultSearchAllCutsMTSToMT", "defaultSearchAllCutsZ"};
   TableBackgroundSignal(this, regions,"combinedChannel" ).Print("metsigall2.tab", 4);
   TableBackgroundSignal(this, regions,"combinedChannel" ).PrintLatex("metsigall2.tex", 4);
}

