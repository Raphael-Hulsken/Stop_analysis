
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

#define USE_VAR_BASELINE

#define USE_LEP1
#define USE_LEP2
#define USE_JETS
#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_GLOBAL_VAR
#define USE_EVT_VAR

#include "../common/TFFactory.h"
// #include "goodrun.cc"
#include "../Selection/moriond.h"

using namespace std;


// ----------------------------------------------
// Should be called only here because many
// struct and fuctions have to be declare first
// ----------------------------------------------
#include "../sonicScrewdriver/interface/BabyScrewdriver.h"

uint32_t counter = 0;
string empty = "";
string storedDataset = "";
TH2D *h2 = NULL;
TH3D *h3 = NULL;
TAxis *xaxis = NULL;
TAxis *yaxis = NULL;
bool checkNegativeYields = false;
uint32_t nthentry = 0;
string outputName = "";
double scale1fbS2 =1;

//float getWeight(string currentProcessType, float lumi, float s1fb2=0);
double getWeight(string currentProcessType, float lumi, double s1fb2=1,bool mediumbtag=true);
map< pair<uint32_t,uint32_t>, string > scanMap;

TFile *fileX = new TFile("../common/xsec_stop_13TeV.root");
TH1D* stopXSEC = (TH1D*)fileX->Get("stop")->Clone();

bool lepChannel() 
{ 
    return true; 
}
    
//Add this as a global variable 

void BabyScrewdriver::Init()
{

    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/cms/ui3_data1/mjansova/Stop1lSharedBabies/v22/skim/";
//    totalNumberOfWorkers = 1;
    totalNumberOfWorkers = 5;

    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    TString fNameTmp =  babyTuplePath+"signal_T2tt_forBinning.root"; //stop
    ftmp = new TFile(fNameTmp);
    htmp = (TH2D*)ftmp->Get("histNEvts")->Clone();
    TH1D* pX = htmp->ProjectionX(); 
    TH1D* pY = htmp->ProjectionY();


    AddVariable("MET", "MET",  "MET", 100 ,0,1000,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "MT2W", 100 ,0,500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("MT", "MT",  "MT", 100 ,0,1000,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("nJets","nJets","nJets",10,1,10,&(myEvent.ngoodjets),"noUnderflowInFirstBin");
    AddVariable("nBJets","nBJets","nBJets",5,1,5,&(myEvent.ngoodbtags),"noUnderflowInFirstBin");
    AddVariable("topnessMod","topnessMod","topnessMod",20,-20,20,&(myEvent.topnessMod),"noUnderflowInFirstBin");
    AddVariable("dphi","dphi","dphi", 100,0,3.5,&(myEvent.dphi_ak4pfjets_met),"noUnderflowInFirstBin");
    AddVariable("Mlb","Mlb","Mlb", 100,0, 500,&(myEvent.Mlb),"noUnderflowInFirstBin");
    AddVariable("StopMass","stop mass", "GeV", pX->GetNbinsX(), pX->GetBinLowEdge(1) ,  (pX->GetBinLowEdge( pX->GetNbinsX())) + (pX->GetBinWidth(pX->GetNbinsX())) ,&myEvent.mass_stop);
    AddVariable("NeutralinoMass","lsp mass", "GeV", pY->GetNbinsX(), pY->GetBinLowEdge(1), (pY->GetBinLowEdge( pY->GetNbinsX())) + (pY->GetBinWidth(pY->GetNbinsX())) , &myEvent.mass_lsp);


    AddProcessClass
