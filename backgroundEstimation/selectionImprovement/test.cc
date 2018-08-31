#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

#define USE_VAR_BASELINE
//#define USE_VAR_BASELINE_UP
//#define USE_VAR_BASELINE_DOWN
#define USE_LEP1
#define USE_LEP2
//#define USE_SKIMMING_VAR
#define USE_JETS
#define USE_JETS_EXT
//#define USE_GEN_INFO
//#define USE_GEN_INFO_EXT
////#define USE_LEP1_EXT
////#define USE_LEP2_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_GLOBAL_VAR

//#include "../../common/common.h"
#include "../../common/TFFactory.h"
//#include "../../common/Reader_CommonFormat_CommonBabies.h" //@MJ@ TODO new selection do not forget to have all requirements (mainly leton), include the new selection
#include "../../Selection/test.h"
//#include "../../Tools/Weighting/WeightFactory.h" //@MJ@ TODO new function for all weights (b tagging, lepton sf)
#define _METCUT_ 50
#define _LEPTONPTCUT_ 40

using namespace std;

#include "signalCS.h"

// ----------------------------------------------
// Should be called only here because many
// struct and fuctions have to be declare first
// ----------------------------------------------
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"

uint32_t counter = 0;
string empty = "";
string storedDataset;
TH2D *h2 = NULL;
TAxis *xaxis = NULL;
TAxis *yaxis = NULL;

    string distingusihClassBkg(string currentProcessClass, vector<string> classLabels); // class, type, list of labels, ...
    void linkTheNameWithProcessClass(string PC);
    void distingusihcClassSig( string currentprocessclass, string currentprocesstype, vector<string> classlabels); // class, type, list of labels, ...
    float getWeight(string currentProcessType, float lumi);
map< pair<uint32_t,uint32_t>, string > scanMap;

bool lepChannel() 
{ 
    return true; 
}
    
ofstream statnames("statNames.txt");
//Add this as a global variable 

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11";
    totalNumberOfWorkers = 10;



    // ------------------
    // Histograms 
    // ------------------
    //vector<float> METBins = {150,200,250,350,450,600,800};
    //vector<float> LepPtBins = {20,25,30,35,40,50,80,100,125,150,200,300};

    // ------------------
    // Variables
    // ------------------
    // fixed binning
    //AddVariable("MET", "MET",  "MET", 100 ,200,1000,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    //AddVariable("MT2W", "MT2W",  "MT2W", 100 ,0,500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    //AddVariable("MT", "MT",  "MT", 100 ,100,1000,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    //AddVariable("nJets","nJets","nJets",5,1,5,&(myEvent.ngoodjets),"noUnderflowInFirstBin");
    //AddVariable("dphi","dphi","dphi", 100,0,3.5,&(myEvent.dphi_ak4pfjets_met),"noUnderflowInFirstBin");
    //AddVariable("nvertex","nvertex","nvertex",50,0,50,&(myEvent.nvertex),"noUnderflowInFirstBin");

    // ------------------
    // Datasets
    // ------------------
    //AddProcessClass("rare", "rare", "background", kBlue);//@MJ@ TODO K-factor?
    AddProcessClass("Znunu", "Znunu", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("ttZJets_13TeV_madgraphMLM","Znunu",0,0);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","Znunu",0,0);
    	//AddDataset("ZZTo2L2Nu_powheg_pythia8_25ns","ttZ",0,0);
    	//AddDataset("WWTo2l2Nu_powheg_25ns","ttZ",0,0);
    //	AddDataset("tZq","rare",0,0.0758);
  ///  	AddDataset("TTZToLLNuNu_M-10_amcnlo_pythia8_25ns","ttZ", 0, 0);
    //	AddDataset("WZ","rare",0,3.06);

    //AddProcessClass("ttZ", "ttZ", "background", kBlue);//@MJ@ TODO K-factor?
    //AddProcessClass("throw", "throw", "signal", kBlue);
     	//AddDataset("T2tt_400to1200", "throw", 0, 0 );
     	//AddDataset("T2tt_mStop_850_mLSP_100_25ns", "throw", 0, 0 );
    //
    //
    
    //signal examples
    //AddProcessClass( "850_100", "850_100", "signal", kBlue);
    //AddProcessClass( "1000_1", "1000_1", "signal", kBlue);

    //AddProcessClass( "grouped", "grouped", "signal", kBlue);
    //AddProcessClass("data", "data", "data", kViolet);
    	//AddDataset("SE_0", "data", 0, 0 );
    	//AddDataset("SE_1", "data", 0, 0 );
       // AddDataset("SM_0", "data", 0, 0 );
       // AddDataset("SM_1", "data", 0, 0 );
        //AddDataset("MET_0", "data", 0, 0 );
        //AddDataset("MET_1", "data", 0, 0 );
    
    //AddProcessClass("ST", "ST", "background", kRed);
    	//AddDataset("t_sch_4f_amcnlo_pythia8_25ns","ST",0,10.11*0.364176);
    	//AddDataset("t_tW_5f_powheg_pythia8_noHadDecays_25ns","ST",0,38.09*0.5135);
    	//AddDataset("t_tbarW_5f_powheg_pythia8_noHadDecays_25ns","ST",0,38.09*0.5135);
    	//AddDataset("t_tch_4f_powheg_pythia8_25ns","ST",0,80.95*0.324);
    	
	//AddDataset("TTJetsSLtop", "test", 0, 114.6*1.594 );
    	//AddDataset("TTJetsSLatopv1","test",0,114.6*1.594);
    	//AddDataset("TTJetsDLv0v4","test",0, 57.35*1.5225);
    	//AddDataset("WJetsToLNuTune","test",0,60781.5*1.01);
    //	AddDataset("W1JetsToLNuTune","test",0, 9493*1.238);
    //	AddDataset("W2JetsToLNuTune","test",0, 3120*1.231);
    //	AddDataset("W3JetsToLNuTune","test",0, 942.3*1.231);
    //	AddDataset("W4JetsToLNuTune","test",0, 524.2*1.114);
   // 	AddDataset("TTWtoQQ","test",0,0.4062);
    //	AddDataset("TTWtoLNu","test",0,0.2043);
    //	AddDataset("TTT","test",0,1.0);
    //	AddDataset("VV","test",0,12.05*0.9917);

    
    //AddProcessClass("lostLepton", "lostLepton", "background", kPink);
    //AddProcessClass("singleLepton", "singleLepton", "background", kGreen);
    //AddProcessClass("singleLeptonFromT", "singleLeptonFromT", "background", kGreen);
    
    // ------------------
    // Regions
    // ------------------
    
//    AddRegion("SR1l","SR1l",&SR1l);
//    AddRegion("CR1l","CR1l",&CR1l);
//    AddRegion("CR2l","CR2l",&CR2l);


//240 regions
AddRegion("SR1l2j_MET250to350","SR1l2j_MET250to350",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350PUdown","SR1l2j_MET250to350PUdown",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350PUup","SR1l2j_MET250to350PUup",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350LSFdown","SR1l2j_MET250to350LSFdown",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350LSFup","SR1l2j_MET250to350LSFup",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350BTlightDown","SR1l2j_MET250to350BTlightDown",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350BTlightUp","SR1l2j_MET250to350BTlightUp",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350BTheavyDown","SR1l2j_MET250to350BTheavyDown",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350BTheavyUp","SR1l2j_MET250to350BTheavyUp",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350PDFdown","SR1l2j_MET250to350PDFdown",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350PDFup","SR1l2j_MET250to350PDFup",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350alphaSdown","SR1l2j_MET250to350alphaSdown",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350alphaSup","SR1l2j_MET250to350alphaSup",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350Q2down","SR1l2j_MET250to350Q2down",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350Q2up","SR1l2j_MET250to350Q2up",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET250to350lumi","SR1l2j_MET250to350lumi",&SR1l2j_MET250to350);
AddRegion("SR1l2j_MET350to450","SR1l2j_MET350to450",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450PUdown","SR1l2j_MET350to450PUdown",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450PUup","SR1l2j_MET350to450PUup",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450LSFdown","SR1l2j_MET350to450LSFdown",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450LSFup","SR1l2j_MET350to450LSFup",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450BTlightDown","SR1l2j_MET350to450BTlightDown",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450BTlightUp","SR1l2j_MET350to450BTlightUp",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450BTheavyDown","SR1l2j_MET350to450BTheavyDown",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450BTheavyUp","SR1l2j_MET350to450BTheavyUp",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450PDFdown","SR1l2j_MET350to450PDFdown",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450PDFup","SR1l2j_MET350to450PDFup",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450alphaSdown","SR1l2j_MET350to450alphaSdown",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450alphaSup","SR1l2j_MET350to450alphaSup",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450Q2down","SR1l2j_MET350to450Q2down",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450Q2up","SR1l2j_MET350to450Q2up",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET350to450lumi","SR1l2j_MET350to450lumi",&SR1l2j_MET350to450);
AddRegion("SR1l2j_MET450toInf","SR1l2j_MET450toInf",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfPUdown","SR1l2j_MET450toInfPUdown",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfPUup","SR1l2j_MET450toInfPUup",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfLSFdown","SR1l2j_MET450toInfLSFdown",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfLSFup","SR1l2j_MET450toInfLSFup",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfBTlightDown","SR1l2j_MET450toInfBTlightDown",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfBTlightUp","SR1l2j_MET450toInfBTlightUp",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfBTheavyDown","SR1l2j_MET450toInfBTheavyDown",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfBTheavyUp","SR1l2j_MET450toInfBTheavyUp",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfPDFdown","SR1l2j_MET450toInfPDFdown",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfPDFup","SR1l2j_MET450toInfPDFup",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfalphaSdown","SR1l2j_MET450toInfalphaSdown",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfalphaSup","SR1l2j_MET450toInfalphaSup",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfQ2down","SR1l2j_MET450toInfQ2down",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInfQ2up","SR1l2j_MET450toInfQ2up",&SR1l2j_MET450toInf);
AddRegion("SR1l2j_MET450toInflumi","SR1l2j_MET450toInflumi",&SR1l2j_MET450toInf);
AddRegion("SR1l3j_MET250to350","SR1l3j_MET250to350",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350PUdown","SR1l3j_MET250to350PUdown",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350PUup","SR1l3j_MET250to350PUup",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350LSFdown","SR1l3j_MET250to350LSFdown",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350LSFup","SR1l3j_MET250to350LSFup",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350BTlightDown","SR1l3j_MET250to350BTlightDown",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350BTlightUp","SR1l3j_MET250to350BTlightUp",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350BTheavyDown","SR1l3j_MET250to350BTheavyDown",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350BTheavyUp","SR1l3j_MET250to350BTheavyUp",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350PDFdown","SR1l3j_MET250to350PDFdown",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350PDFup","SR1l3j_MET250to350PDFup",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350alphaSdown","SR1l3j_MET250to350alphaSdown",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350alphaSup","SR1l3j_MET250to350alphaSup",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350Q2down","SR1l3j_MET250to350Q2down",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350Q2up","SR1l3j_MET250to350Q2up",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET250to350lumi","SR1l3j_MET250to350lumi",&SR1l3j_MET250to350);
AddRegion("SR1l3j_MET350to450","SR1l3j_MET350to450",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450PUdown","SR1l3j_MET350to450PUdown",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450PUup","SR1l3j_MET350to450PUup",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450LSFdown","SR1l3j_MET350to450LSFdown",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450LSFup","SR1l3j_MET350to450LSFup",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450BTlightDown","SR1l3j_MET350to450BTlightDown",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450BTlightUp","SR1l3j_MET350to450BTlightUp",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450BTheavyDown","SR1l3j_MET350to450BTheavyDown",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450BTheavyUp","SR1l3j_MET350to450BTheavyUp",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450PDFdown","SR1l3j_MET350to450PDFdown",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450PDFup","SR1l3j_MET350to450PDFup",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450alphaSdown","SR1l3j_MET350to450alphaSdown",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450alphaSup","SR1l3j_MET350to450alphaSup",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450Q2down","SR1l3j_MET350to450Q2down",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450Q2up","SR1l3j_MET350to450Q2up",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET350to450lumi","SR1l3j_MET350to450lumi",&SR1l3j_MET350to450);
AddRegion("SR1l3j_MET450to550","SR1l3j_MET450to550",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550PUdown","SR1l3j_MET450to550PUdown",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550PUup","SR1l3j_MET450to550PUup",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550LSFdown","SR1l3j_MET450to550LSFdown",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550LSFup","SR1l3j_MET450to550LSFup",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550BTlightDown","SR1l3j_MET450to550BTlightDown",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550BTlightUp","SR1l3j_MET450to550BTlightUp",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550BTheavyDown","SR1l3j_MET450to550BTheavyDown",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550BTheavyUp","SR1l3j_MET450to550BTheavyUp",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550PDFdown","SR1l3j_MET450to550PDFdown",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550PDFup","SR1l3j_MET450to550PDFup",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550alphaSdown","SR1l3j_MET450to550alphaSdown",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550alphaSup","SR1l3j_MET450to550alphaSup",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550Q2down","SR1l3j_MET450to550Q2down",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550Q2up","SR1l3j_MET450to550Q2up",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET450to550lumi","SR1l3j_MET450to550lumi",&SR1l3j_MET450to550);
AddRegion("SR1l3j_MET550toInf","SR1l3j_MET550toInf",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfPUdown","SR1l3j_MET550toInfPUdown",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfPUup","SR1l3j_MET550toInfPUup",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfLSFdown","SR1l3j_MET550toInfLSFdown",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfLSFup","SR1l3j_MET550toInfLSFup",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfBTlightDown","SR1l3j_MET550toInfBTlightDown",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfBTlightUp","SR1l3j_MET550toInfBTlightUp",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfBTheavyDown","SR1l3j_MET550toInfBTheavyDown",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfBTheavyUp","SR1l3j_MET550toInfBTheavyUp",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfPDFdown","SR1l3j_MET550toInfPDFdown",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfPDFup","SR1l3j_MET550toInfPDFup",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfalphaSdown","SR1l3j_MET550toInfalphaSdown",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfalphaSup","SR1l3j_MET550toInfalphaSup",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfQ2down","SR1l3j_MET550toInfQ2down",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInfQ2up","SR1l3j_MET550toInfQ2up",&SR1l3j_MET550toInf);
AddRegion("SR1l3j_MET550toInflumi","SR1l3j_MET550toInflumi",&SR1l3j_MET550toInf);
AddRegion("SR1l4jLow_MET250to350","SR1l4jLow_MET250to350",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350PUdown","SR1l4jLow_MET250to350PUdown",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350PUup","SR1l4jLow_MET250to350PUup",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350LSFdown","SR1l4jLow_MET250to350LSFdown",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350LSFup","SR1l4jLow_MET250to350LSFup",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350BTlightDown","SR1l4jLow_MET250to350BTlightDown",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350BTlightUp","SR1l4jLow_MET250to350BTlightUp",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350BTheavyDown","SR1l4jLow_MET250to350BTheavyDown",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350BTheavyUp","SR1l4jLow_MET250to350BTheavyUp",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350PDFdown","SR1l4jLow_MET250to350PDFdown",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350PDFup","SR1l4jLow_MET250to350PDFup",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350alphaSdown","SR1l4jLow_MET250to350alphaSdown",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350alphaSup","SR1l4jLow_MET250to350alphaSup",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350Q2down","SR1l4jLow_MET250to350Q2down",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350Q2up","SR1l4jLow_MET250to350Q2up",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET250to350lumi","SR1l4jLow_MET250to350lumi",&SR1l4jLow_MET250to350);
AddRegion("SR1l4jLow_MET350to450","SR1l4jLow_MET350to450",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450PUdown","SR1l4jLow_MET350to450PUdown",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450PUup","SR1l4jLow_MET350to450PUup",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450LSFdown","SR1l4jLow_MET350to450LSFdown",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450LSFup","SR1l4jLow_MET350to450LSFup",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450BTlightDown","SR1l4jLow_MET350to450BTlightDown",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450BTlightUp","SR1l4jLow_MET350to450BTlightUp",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450BTheavyDown","SR1l4jLow_MET350to450BTheavyDown",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450BTheavyUp","SR1l4jLow_MET350to450BTheavyUp",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450PDFdown","SR1l4jLow_MET350to450PDFdown",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450PDFup","SR1l4jLow_MET350to450PDFup",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450alphaSdown","SR1l4jLow_MET350to450alphaSdown",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450alphaSup","SR1l4jLow_MET350to450alphaSup",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450Q2down","SR1l4jLow_MET350to450Q2down",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450Q2up","SR1l4jLow_MET350to450Q2up",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET350to450lumi","SR1l4jLow_MET350to450lumi",&SR1l4jLow_MET350to450);
AddRegion("SR1l4jLow_MET450toInf","SR1l4jLow_MET450toInf",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfPUdown","SR1l4jLow_MET450toInfPUdown",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfPUup","SR1l4jLow_MET450toInfPUup",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfLSFdown","SR1l4jLow_MET450toInfLSFdown",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfLSFup","SR1l4jLow_MET450toInfLSFup",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfBTlightDown","SR1l4jLow_MET450toInfBTlightDown",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfBTlightUp","SR1l4jLow_MET450toInfBTlightUp",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfBTheavyDown","SR1l4jLow_MET450toInfBTheavyDown",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfBTheavyUp","SR1l4jLow_MET450toInfBTheavyUp",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfPDFdown","SR1l4jLow_MET450toInfPDFdown",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfPDFup","SR1l4jLow_MET450toInfPDFup",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfalphaSdown","SR1l4jLow_MET450toInfalphaSdown",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfalphaSup","SR1l4jLow_MET450toInfalphaSup",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfQ2down","SR1l4jLow_MET450toInfQ2down",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInfQ2up","SR1l4jLow_MET450toInfQ2up",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jLow_MET450toInflumi","SR1l4jLow_MET450toInflumi",&SR1l4jLow_MET450toInf);
AddRegion("SR1l4jHighMT2W_MET250to350","SR1l4jHighMT2W_MET250to350",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350PUdown","SR1l4jHighMT2W_MET250to350PUdown",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350PUup","SR1l4jHighMT2W_MET250to350PUup",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350LSFdown","SR1l4jHighMT2W_MET250to350LSFdown",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350LSFup","SR1l4jHighMT2W_MET250to350LSFup",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350BTlightDown","SR1l4jHighMT2W_MET250to350BTlightDown",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350BTlightUp","SR1l4jHighMT2W_MET250to350BTlightUp",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350BTheavyDown","SR1l4jHighMT2W_MET250to350BTheavyDown",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350BTheavyUp","SR1l4jHighMT2W_MET250to350BTheavyUp",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350PDFdown","SR1l4jHighMT2W_MET250to350PDFdown",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350PDFup","SR1l4jHighMT2W_MET250to350PDFup",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350alphaSdown","SR1l4jHighMT2W_MET250to350alphaSdown",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350alphaSup","SR1l4jHighMT2W_MET250to350alphaSup",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350Q2down","SR1l4jHighMT2W_MET250to350Q2down",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350Q2up","SR1l4jHighMT2W_MET250to350Q2up",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET250to350lumi","SR1l4jHighMT2W_MET250to350lumi",&SR1l4jHighMT2W_MET250to350);
AddRegion("SR1l4jHighMT2W_MET350to450","SR1l4jHighMT2W_MET350to450",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450PUdown","SR1l4jHighMT2W_MET350to450PUdown",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450PUup","SR1l4jHighMT2W_MET350to450PUup",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450LSFdown","SR1l4jHighMT2W_MET350to450LSFdown",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450LSFup","SR1l4jHighMT2W_MET350to450LSFup",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450BTlightDown","SR1l4jHighMT2W_MET350to450BTlightDown",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450BTlightUp","SR1l4jHighMT2W_MET350to450BTlightUp",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450BTheavyDown","SR1l4jHighMT2W_MET350to450BTheavyDown",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450BTheavyUp","SR1l4jHighMT2W_MET350to450BTheavyUp",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450PDFdown","SR1l4jHighMT2W_MET350to450PDFdown",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450PDFup","SR1l4jHighMT2W_MET350to450PDFup",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450alphaSdown","SR1l4jHighMT2W_MET350to450alphaSdown",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450alphaSup","SR1l4jHighMT2W_MET350to450alphaSup",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450Q2down","SR1l4jHighMT2W_MET350to450Q2down",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450Q2up","SR1l4jHighMT2W_MET350to450Q2up",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET350to450lumi","SR1l4jHighMT2W_MET350to450lumi",&SR1l4jHighMT2W_MET350to450);
AddRegion("SR1l4jHighMT2W_MET450to550","SR1l4jHighMT2W_MET450to550",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550PUdown","SR1l4jHighMT2W_MET450to550PUdown",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550PUup","SR1l4jHighMT2W_MET450to550PUup",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550LSFdown","SR1l4jHighMT2W_MET450to550LSFdown",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550LSFup","SR1l4jHighMT2W_MET450to550LSFup",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550BTlightDown","SR1l4jHighMT2W_MET450to550BTlightDown",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550BTlightUp","SR1l4jHighMT2W_MET450to550BTlightUp",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550BTheavyDown","SR1l4jHighMT2W_MET450to550BTheavyDown",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550BTheavyUp","SR1l4jHighMT2W_MET450to550BTheavyUp",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550PDFdown","SR1l4jHighMT2W_MET450to550PDFdown",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550PDFup","SR1l4jHighMT2W_MET450to550PDFup",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550alphaSdown","SR1l4jHighMT2W_MET450to550alphaSdown",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550alphaSup","SR1l4jHighMT2W_MET450to550alphaSup",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550Q2down","SR1l4jHighMT2W_MET450to550Q2down",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550Q2up","SR1l4jHighMT2W_MET450to550Q2up",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET450to550lumi","SR1l4jHighMT2W_MET450to550lumi",&SR1l4jHighMT2W_MET450to550);
AddRegion("SR1l4jHighMT2W_MET550to650","SR1l4jHighMT2W_MET550to650",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650PUdown","SR1l4jHighMT2W_MET550to650PUdown",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650PUup","SR1l4jHighMT2W_MET550to650PUup",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650LSFdown","SR1l4jHighMT2W_MET550to650LSFdown",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650LSFup","SR1l4jHighMT2W_MET550to650LSFup",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650BTlightDown","SR1l4jHighMT2W_MET550to650BTlightDown",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650BTlightUp","SR1l4jHighMT2W_MET550to650BTlightUp",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650BTheavyDown","SR1l4jHighMT2W_MET550to650BTheavyDown",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650BTheavyUp","SR1l4jHighMT2W_MET550to650BTheavyUp",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650PDFdown","SR1l4jHighMT2W_MET550to650PDFdown",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650PDFup","SR1l4jHighMT2W_MET550to650PDFup",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650alphaSdown","SR1l4jHighMT2W_MET550to650alphaSdown",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650alphaSup","SR1l4jHighMT2W_MET550to650alphaSup",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650Q2down","SR1l4jHighMT2W_MET550to650Q2down",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650Q2up","SR1l4jHighMT2W_MET550to650Q2up",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET550to650lumi","SR1l4jHighMT2W_MET550to650lumi",&SR1l4jHighMT2W_MET550to650);
AddRegion("SR1l4jHighMT2W_MET650toInf","SR1l4jHighMT2W_MET650toInf",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfPUdown","SR1l4jHighMT2W_MET650toInfPUdown",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfPUup","SR1l4jHighMT2W_MET650toInfPUup",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfLSFdown","SR1l4jHighMT2W_MET650toInfLSFdown",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfLSFup","SR1l4jHighMT2W_MET650toInfLSFup",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfBTlightDown","SR1l4jHighMT2W_MET650toInfBTlightDown",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfBTlightUp","SR1l4jHighMT2W_MET650toInfBTlightUp",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfBTheavyDown","SR1l4jHighMT2W_MET650toInfBTheavyDown",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfBTheavyUp","SR1l4jHighMT2W_MET650toInfBTheavyUp",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfPDFdown","SR1l4jHighMT2W_MET650toInfPDFdown",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfPDFup","SR1l4jHighMT2W_MET650toInfPDFup",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfalphaSdown","SR1l4jHighMT2W_MET650toInfalphaSdown",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfalphaSup","SR1l4jHighMT2W_MET650toInfalphaSup",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfQ2down","SR1l4jHighMT2W_MET650toInfQ2down",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInfQ2up","SR1l4jHighMT2W_MET650toInfQ2up",&SR1l4jHighMT2W_MET650toInf);
AddRegion("SR1l4jHighMT2W_MET650toInflumi","SR1l4jHighMT2W_MET650toInflumi",&SR1l4jHighMT2W_MET650toInf);
                                                                                                                       

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(12.9);

    Create1DHistos();
    //Add2DHisto("nJets","MET");

    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);
    bool useTriggerInfo = currentProcessType == "data" ? true: false; //@MJ@ TODO fix the trigger


    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);

    /*if(currentProcessType == "blablabla" )//normally background
    {
    string PC = currentProcessClass;
    PC = distingusihClassBkg(currentProcessClass, classLabels);
    if(PC != "")
    {
        if (std::find(classLabels.begin(), classLabels.end(), PC) != classLabels.end())
        {  
            currentProcessClass = PC;
        }
        else
        {
            cout << "ProcessClass " << PC << " was not found in between the existing process classes, the existing classes are: " << endl;
            for(uint32_t c = 0; c<classLabels.size(); c++)
            {
                cout << classLabels.at(c);
            }            
            throw std::runtime_error("no class to atribute process to was not found");
        }
    }
    }

     //cout << " stop " << myEvent.mass_stop << " lsp " << myEvent.mass_lsp << endl; //@MJ@ TODO why these MASSES ARE 0

    for(uint32_t s= 0; s<myEvent.gensusy_id->size(); s++)
    {
        if(currentProcessType == "signal" && ( myEvent.gensusy_id->at(s) == 1000006 && sqrt(abs(myEvent.gensusy_p4->at(s).M2())) == 850))
        {
            for(uint32_t n= 0; n<myEvent.gensusy_id->size(); n++)
            {
                 if( ( myEvent.gensusy_id->at(n) == 1000022  && sqrt(abs(myEvent.gensusy_p4->at(n).M2())) == 100))
                 {
	             currentProcessClass = "850_100";
                     break;
                 }
            }
        }
    }
    */
    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
    if( currentProcessClass == "Znunu" && !(myEvent.isZtoNuNu) )
         currentProcessClass = "";

    float weightLumi = getWeight(currentProcessType, GetLumi()); //@MJ@ TODO cross section form file?!
    //@MJ@ TODO I hate myself for this, but no better solution foud
    //computation of up/down weights
    //NOTICE, important is to fill only weight histo and only have 1 process class

    vector<float> weightV;
    weightV.clear();
    float nEvents =  myEvent.wNormalization.at(22);
    //for number of SR
    for(uint32_t SR=0; SR<15; SR++) //@MJ@ TODO nr of sig regions changes
    {
        float w = 0;
        float btagmax = 0;
        //normal
        weightV.push_back(weightLumi);
        //PUdown
        if(counter == 1) statnames << "PUdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.weight_PUdown; //@MJ@ TODO PU without any normalization?!
        weightV.push_back(w);
        //PUup
        if(counter == 1) statnames << "PUup"<< endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.weight_PUup;
        weightV.push_back(w);
        //LSFdown
        if(counter == 1) statnames << "LSFdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF_down*( nEvents / myEvent.wNormalization.at(30) );
        weightV.push_back(w);
        //LSFup
        if(counter == 1) statnames << "LSFup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF_up*( nEvents / myEvent.wNormalization.at(29) );
        weightV.push_back(w);
        //BTlightdown
        if(counter == 1) statnames << "BTligtdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_light_DN*( nEvents / myEvent.wNormalization.at(18) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) );
        weightV.push_back(w);
        //BTlightup
        if(counter == 1) statnames << "BTligtup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_light_UP*( nEvents / myEvent.wNormalization.at(16) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) ;
        weightV.push_back(w);
        //BTheabydown
        if(counter == 1) statnames << "BTheavydown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_heavy_DN*( nEvents / myEvent.wNormalization.at(17) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) );
        weightV.push_back(w);
        //BTheavyup
        if(counter == 1) statnames << "BTheavyup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_heavy_UP*( nEvents / myEvent.wNormalization.at(15) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) ;
        weightV.push_back(w);
        //PDFdown
        if(counter == 1) statnames << "PDFdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.pdf_down_weight*( nEvents / myEvent.wNormalization.at(11) );
        weightV.push_back(w);
        //PDFup
        if(counter == 1) statnames << "PDFup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.pdf_up_weight*( nEvents / myEvent.wNormalization.at(10) );
        weightV.push_back(w);
        //alphaSdown
        if(counter == 1) statnames << "alphaSdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_alphas_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(13)) ; //TODO
        weightV.push_back(w);
        //alphaSup
        if(counter == 1) statnames << "alphaSup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_alphas_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(12))  ; //TODO
        weightV.push_back(w);
        //Q2down
        if(counter == 1) statnames << "Q2down" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_q2_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(9))  ; //TODO
        weightV.push_back(w);
        //Q2up
        if(counter == 1) statnames << "Q2up" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28)) * myEvent.weight_q2_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(5))  ; //TODO
        weightV.push_back(w);
        //lumi //for ICHEP now
        if(counter == 1) statnames << "lumi" << endl;
        w = GetLumi()*1.062 *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  ; //TODO
        weightV.push_back(w);
    }

    vector<string> theReg;
    GetRegionTagList(&theReg);
    if( weightV.size() != theReg.size())
        throw std::runtime_error("vector of weights does not have same size as regions, the weights will not be correctly assesed");

    
    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    //if (currentProcessType == "signal") weight = 1;
    AutoFillProcessClass(currentProcessClass, weight, weightV, true);

    //cout << "weight for process " << currentProcessType << " is " << weight << endl;

    if(counter % 10000 == 0)
    {
        cout << counter << endl;
    }

    //fle->Delete();
}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    // ######################
    //  Plot configuration and production
    // ######################

    // Schedule plots
    //

    //SchedulePlots("1DSuperimposed");
    //SchedulePlots("1DSuperimposedNoNorm");
    //SchedulePlots("1DStack");
    //SchedulePlots("2D");

    // Config plots
    statnames.close();

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", false);
    SetGlobalBoolOption("Plot", "exportEps", true);
    SetGlobalBoolOption("Plot", "exportPng", false);

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plotsTest/");

    // ######################
    //  Tables and other stuff
    // ######################


vector<string> totYield = { "SR1l2j_MET250to350" , "SR1l2j_MET250to350PUdown" , "SR1l2j_MET250to350PUup" , "SR1l2j_MET250to350LSFdown" , "SR1l2j_MET250to350LSFup" , "SR1l2j_MET250to350BTlightDown" , "SR1l2j_MET250to350BTlightUp" , "SR1l2j_MET250to350BTheavyDown" , "SR1l2j_MET250to350BTheavyUp" , "SR1l2j_MET250to350PDFdown" , "SR1l2j_MET250to350PDFup" , "SR1l2j_MET250to350alphaSdown" , "SR1l2j_MET250to350alphaSup" , "SR1l2j_MET250to350Q2down" , "SR1l2j_MET250to350Q2up" , "SR1l2j_MET250to350lumi" , "SR1l2j_MET350to450" , "SR1l2j_MET350to450PUdown" , "SR1l2j_MET350to450PUup" , "SR1l2j_MET350to450LSFdown" , "SR1l2j_MET350to450LSFup" , "SR1l2j_MET350to450BTlightDown" , "SR1l2j_MET350to450BTlightUp" , "SR1l2j_MET350to450BTheavyDown" , "SR1l2j_MET350to450BTheavyUp" , "SR1l2j_MET350to450PDFdown" , "SR1l2j_MET350to450PDFup" , "SR1l2j_MET350to450alphaSdown" , "SR1l2j_MET350to450alphaSup" , "SR1l2j_MET350to450Q2down" , "SR1l2j_MET350to450Q2up" , "SR1l2j_MET350to450lumi" , "SR1l2j_MET450toInf" , "SR1l2j_MET450toInfPUdown" , "SR1l2j_MET450toInfPUup" , "SR1l2j_MET450toInfLSFdown" , "SR1l2j_MET450toInfLSFup" , "SR1l2j_MET450toInfBTlightDown" , "SR1l2j_MET450toInfBTlightUp" , "SR1l2j_MET450toInfBTheavyDown" , "SR1l2j_MET450toInfBTheavyUp" , "SR1l2j_MET450toInfPDFdown" , "SR1l2j_MET450toInfPDFup" , "SR1l2j_MET450toInfalphaSdown" , "SR1l2j_MET450toInfalphaSup" , "SR1l2j_MET450toInfQ2down" , "SR1l2j_MET450toInfQ2up" , "SR1l2j_MET450toInflumi" , "SR1l3j_MET250to350" , "SR1l3j_MET250to350PUdown" , "SR1l3j_MET250to350PUup" , "SR1l3j_MET250to350LSFdown" , "SR1l3j_MET250to350LSFup" , "SR1l3j_MET250to350BTlightDown" , "SR1l3j_MET250to350BTlightUp" , "SR1l3j_MET250to350BTheavyDown" , "SR1l3j_MET250to350BTheavyUp" , "SR1l3j_MET250to350PDFdown" , "SR1l3j_MET250to350PDFup" , "SR1l3j_MET250to350alphaSdown" , "SR1l3j_MET250to350alphaSup" , "SR1l3j_MET250to350Q2down" , "SR1l3j_MET250to350Q2up" , "SR1l3j_MET250to350lumi" , "SR1l3j_MET350to450" , "SR1l3j_MET350to450PUdown" , "SR1l3j_MET350to450PUup" , "SR1l3j_MET350to450LSFdown" , "SR1l3j_MET350to450LSFup" , "SR1l3j_MET350to450BTlightDown" , "SR1l3j_MET350to450BTlightUp" , "SR1l3j_MET350to450BTheavyDown" , "SR1l3j_MET350to450BTheavyUp" , "SR1l3j_MET350to450PDFdown" , "SR1l3j_MET350to450PDFup" , "SR1l3j_MET350to450alphaSdown" , "SR1l3j_MET350to450alphaSup" , "SR1l3j_MET350to450Q2down" , "SR1l3j_MET350to450Q2up" , "SR1l3j_MET350to450lumi" , "SR1l3j_MET450to550" , "SR1l3j_MET450to550PUdown" , "SR1l3j_MET450to550PUup" , "SR1l3j_MET450to550LSFdown" , "SR1l3j_MET450to550LSFup" , "SR1l3j_MET450to550BTlightDown" , "SR1l3j_MET450to550BTlightUp" , "SR1l3j_MET450to550BTheavyDown" , "SR1l3j_MET450to550BTheavyUp" , "SR1l3j_MET450to550PDFdown" , "SR1l3j_MET450to550PDFup" , "SR1l3j_MET450to550alphaSdown" , "SR1l3j_MET450to550alphaSup" , "SR1l3j_MET450to550Q2down" , "SR1l3j_MET450to550Q2up" , "SR1l3j_MET450to550lumi" , "SR1l3j_MET550toInf" , "SR1l3j_MET550toInfPUdown" , "SR1l3j_MET550toInfPUup" , "SR1l3j_MET550toInfLSFdown" , "SR1l3j_MET550toInfLSFup" , "SR1l3j_MET550toInfBTlightDown" , "SR1l3j_MET550toInfBTlightUp" , "SR1l3j_MET550toInfBTheavyDown" , "SR1l3j_MET550toInfBTheavyUp" , "SR1l3j_MET550toInfPDFdown" , "SR1l3j_MET550toInfPDFup" , "SR1l3j_MET550toInfalphaSdown" , "SR1l3j_MET550toInfalphaSup" , "SR1l3j_MET550toInfQ2down" , "SR1l3j_MET550toInfQ2up" , "SR1l3j_MET550toInflumi" , "SR1l4jLow_MET250to350" , "SR1l4jLow_MET250to350PUdown" , "SR1l4jLow_MET250to350PUup" , "SR1l4jLow_MET250to350LSFdown" , "SR1l4jLow_MET250to350LSFup" , "SR1l4jLow_MET250to350BTlightDown" , "SR1l4jLow_MET250to350BTlightUp" , "SR1l4jLow_MET250to350BTheavyDown" , "SR1l4jLow_MET250to350BTheavyUp" , "SR1l4jLow_MET250to350PDFdown" , "SR1l4jLow_MET250to350PDFup" , "SR1l4jLow_MET250to350alphaSdown" , "SR1l4jLow_MET250to350alphaSup" , "SR1l4jLow_MET250to350Q2down" , "SR1l4jLow_MET250to350Q2up" , "SR1l4jLow_MET250to350lumi" , "SR1l4jLow_MET350to450" , "SR1l4jLow_MET350to450PUdown" , "SR1l4jLow_MET350to450PUup" , "SR1l4jLow_MET350to450LSFdown" , "SR1l4jLow_MET350to450LSFup" , "SR1l4jLow_MET350to450BTlightDown" , "SR1l4jLow_MET350to450BTlightUp" , "SR1l4jLow_MET350to450BTheavyDown" , "SR1l4jLow_MET350to450BTheavyUp" , "SR1l4jLow_MET350to450PDFdown" , "SR1l4jLow_MET350to450PDFup" , "SR1l4jLow_MET350to450alphaSdown" , "SR1l4jLow_MET350to450alphaSup" , "SR1l4jLow_MET350to450Q2down" , "SR1l4jLow_MET350to450Q2up" , "SR1l4jLow_MET350to450lumi" , "SR1l4jLow_MET450toInf" , "SR1l4jLow_MET450toInfPUdown" , "SR1l4jLow_MET450toInfPUup" , "SR1l4jLow_MET450toInfLSFdown" , "SR1l4jLow_MET450toInfLSFup" , "SR1l4jLow_MET450toInfBTlightDown" , "SR1l4jLow_MET450toInfBTlightUp" , "SR1l4jLow_MET450toInfBTheavyDown" , "SR1l4jLow_MET450toInfBTheavyUp" , "SR1l4jLow_MET450toInfPDFdown" , "SR1l4jLow_MET450toInfPDFup" , "SR1l4jLow_MET450toInfalphaSdown" , "SR1l4jLow_MET450toInfalphaSup" , "SR1l4jLow_MET450toInfQ2down" , "SR1l4jLow_MET450toInfQ2up" , "SR1l4jLow_MET450toInflumi" , "SR1l4jHighMT2W_MET250to350" , "SR1l4jHighMT2W_MET250to350PUdown" , "SR1l4jHighMT2W_MET250to350PUup" , "SR1l4jHighMT2W_MET250to350LSFdown" , "SR1l4jHighMT2W_MET250to350LSFup" , "SR1l4jHighMT2W_MET250to350BTlightDown" , "SR1l4jHighMT2W_MET250to350BTlightUp" , "SR1l4jHighMT2W_MET250to350BTheavyDown" , "SR1l4jHighMT2W_MET250to350BTheavyUp" , "SR1l4jHighMT2W_MET250to350PDFdown" , "SR1l4jHighMT2W_MET250to350PDFup" , "SR1l4jHighMT2W_MET250to350alphaSdown" , "SR1l4jHighMT2W_MET250to350alphaSup" , "SR1l4jHighMT2W_MET250to350Q2down" , "SR1l4jHighMT2W_MET250to350Q2up" , "SR1l4jHighMT2W_MET250to350lumi" , "SR1l4jHighMT2W_MET350to450" , "SR1l4jHighMT2W_MET350to450PUdown" , "SR1l4jHighMT2W_MET350to450PUup" , "SR1l4jHighMT2W_MET350to450LSFdown" , "SR1l4jHighMT2W_MET350to450LSFup" , "SR1l4jHighMT2W_MET350to450BTlightDown" , "SR1l4jHighMT2W_MET350to450BTlightUp" , "SR1l4jHighMT2W_MET350to450BTheavyDown" , "SR1l4jHighMT2W_MET350to450BTheavyUp" , "SR1l4jHighMT2W_MET350to450PDFdown" , "SR1l4jHighMT2W_MET350to450PDFup" , "SR1l4jHighMT2W_MET350to450alphaSdown" , "SR1l4jHighMT2W_MET350to450alphaSup" , "SR1l4jHighMT2W_MET350to450Q2down" , "SR1l4jHighMT2W_MET350to450Q2up" , "SR1l4jHighMT2W_MET350to450lumi" , "SR1l4jHighMT2W_MET450to550" , "SR1l4jHighMT2W_MET450to550PUdown" , "SR1l4jHighMT2W_MET450to550PUup" , "SR1l4jHighMT2W_MET450to550LSFdown" , "SR1l4jHighMT2W_MET450to550LSFup" , "SR1l4jHighMT2W_MET450to550BTlightDown" , "SR1l4jHighMT2W_MET450to550BTlightUp" , "SR1l4jHighMT2W_MET450to550BTheavyDown" , "SR1l4jHighMT2W_MET450to550BTheavyUp" , "SR1l4jHighMT2W_MET450to550PDFdown" , "SR1l4jHighMT2W_MET450to550PDFup" , "SR1l4jHighMT2W_MET450to550alphaSdown" , "SR1l4jHighMT2W_MET450to550alphaSup" , "SR1l4jHighMT2W_MET450to550Q2down" , "SR1l4jHighMT2W_MET450to550Q2up" , "SR1l4jHighMT2W_MET450to550lumi" , "SR1l4jHighMT2W_MET550to650" , "SR1l4jHighMT2W_MET550to650PUdown" , "SR1l4jHighMT2W_MET550to650PUup" , "SR1l4jHighMT2W_MET550to650LSFdown" , "SR1l4jHighMT2W_MET550to650LSFup" , "SR1l4jHighMT2W_MET550to650BTlightDown" , "SR1l4jHighMT2W_MET550to650BTlightUp" , "SR1l4jHighMT2W_MET550to650BTheavyDown" , "SR1l4jHighMT2W_MET550to650BTheavyUp" , "SR1l4jHighMT2W_MET550to650PDFdown" , "SR1l4jHighMT2W_MET550to650PDFup" , "SR1l4jHighMT2W_MET550to650alphaSdown" , "SR1l4jHighMT2W_MET550to650alphaSup" , "SR1l4jHighMT2W_MET550to650Q2down" , "SR1l4jHighMT2W_MET550to650Q2up" , "SR1l4jHighMT2W_MET550to650lumi" , "SR1l4jHighMT2W_MET650toInf" , "SR1l4jHighMT2W_MET650toInfPUdown" , "SR1l4jHighMT2W_MET650toInfPUup" , "SR1l4jHighMT2W_MET650toInfLSFdown" , "SR1l4jHighMT2W_MET650toInfLSFup" , "SR1l4jHighMT2W_MET650toInfBTlightDown" , "SR1l4jHighMT2W_MET650toInfBTlightUp" , "SR1l4jHighMT2W_MET650toInfBTheavyDown" , "SR1l4jHighMT2W_MET650toInfBTheavyUp" , "SR1l4jHighMT2W_MET650toInfPDFdown" , "SR1l4jHighMT2W_MET650toInfPDFup" , "SR1l4jHighMT2W_MET650toInfalphaSdown" , "SR1l4jHighMT2W_MET650toInfalphaSup" , "SR1l4jHighMT2W_MET650toInfQ2down" , "SR1l4jHighMT2W_MET650toInfQ2up" , "SR1l4jHighMT2W_MET650toInflumi"};

    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print("yield.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex("yield.tex", 4);

    ofstream sigfile("signalReg.txt");
    if (sigfile.is_open())
    {
        for(uint32_t r=0; r<totYield.size(); r++)
        {
            sigfile << totYield.at(r) << endl;
        }
            sigfile.close();
    }


    vector <string> tfreg{};


    TFProducer prod(tfreg, "lostLepton", "CR2l");
    prod.produceTFTable("yield.tab", "lostLeptonTF");
    TFProducer prod2(tfreg, "singleLepton", "CR1l");
    prod2.produceTFTable("yield.tab", "singleLeptonTF");
    TFProducer prod3(tfreg, "singleLeptonFromT", "CR1l");
    prod3.produceTFTable("yield.tab", "singleLeptonFromTTF");

    cout << "end of processing" << endl;
 }




    string distingusihClassBkg(string currentProcessClass, vector<string> classLabels) // class, type, list of labels, ...
    {
        string PC = "";
        if(currentProcessClass == "test" && myEvent.is2lep)
        {
	    PC = "lostLepton";
            //linkTheNameWithProcessClass(PC);
        }
        else if(currentProcessClass == "test" && (myEvent.is1lep || myEvent.is0lep ) && myEvent.is1lepFromTop)
        {
	    PC = "singleLeptonFromT";
            //linkTheNameWithProcessClass(PC);
        }
        else if(currentProcessClass == "test" && (myEvent.is1lep || myEvent.is0lep ))
        {
	    PC = "singleLepton";
            //linkTheNameWithProcessClass(PC);
        }
        else if(currentProcessClass == "test") 
        {
	    //cout << "nr of leptons " << myEvent.numberOfGeneratedLeptons <<endl;
	throw std::runtime_error("This should not happen");
        }
        else
        {}
        return PC;
    }

    void linkTheNameWithProcessClass(string PC)
    {/*
        if (std::find(classLabels.begin(), classLabels.end(), PC) != classLabels.end())
        {  
            currentProcessClass = PC;
        }
        else
        {
            count << "ProcessClass " << PC << " was not found in between the existing process classes, the existing classes are: " << endl;
            for(uint32_t c = 0; c<classLabels.size(); c++)
            {
                cout << classLabels.at(c);
            }            
            throw std::runtime_error("no class to atribute process to was not found");
        }
   */ }

    void distingusihcClassSig( string currentprocessclass, string currentprocesstype, vector<string> classlabels) // class, type, list of labels, ...
    {/*
        if(currentProcessType == "signal" && myEvent.mass_stop == 900 && myEvent.mass_lsp == 50)
        {
	    currentProcessClass = "900_50";
        }
        if(currentProcessType == "signal" && myEvent.mass_stop == 1000 && myEvent.mass_lsp == 1)
        {
	    currentProcessClass = "1000_1";
        }
       */
     }

    float getWeight(string currentProcessType, float lumi)
    {
        float nEvents =  myEvent.wNormalization.at(22);
        float all_weights = lumi*  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        //cout << "cs: " << myEvent.crossSection << "scale 1 fb " << myEvent.scale1fb << " weight total: " << all_weights << endl;
        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
        return all_weights;
    }


