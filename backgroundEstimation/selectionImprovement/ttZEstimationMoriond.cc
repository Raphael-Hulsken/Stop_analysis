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
#include "../../Selection/moriond.h"
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
bool checkNegativeYields = false;
uint32_t nthentry = 0;

    string distingusihClassBkg(string currentProcessClass, vector<string> classLabels); // class, type, list of labels, ...
    void linkTheNameWithProcessClass(string PC);
    void distingusihcClassSig( string currentprocessclass, string currentprocesstype, vector<string> classlabels); // class, type, list of labels, ...
    float getWeight(string currentProcessType, float lumi);
    float reweightTop(float topPt, float atopPt);
map< pair<uint32_t,uint32_t>, string > scanMap;

bool lepChannel() 
{ 
    return true; 
}
    
ofstream statnames("statNames.txt");
//Add this as a global variable 

TH1F * w_pdfnu = new TH1F("w_pdfnu", "w_pdfnu",30000 , -100000,100000);
TH1F * w_pdfnd = new TH1F("w_pdfnd", "w_pdfnd",30000 , -100000,100000);
TH1F * w_pdfu = new TH1F("w_pdfu", "w_pdfu", 30000 , -0.01,0.01);
TH1F * w_pdfd = new TH1F("w_pdfd", "w_pdfd", 3000 , -0.01,0.01);

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11";
    totalNumberOfWorkers = 1;



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
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","Znunu",0,0);
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


//405 regions
AddRegion("SR1l_A_250lessMETless350","SR1l_A_250lessMETless350",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350LSFdown","SR1l_A_250lessMETless350LSFdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350LSFup","SR1l_A_250lessMETless350LSFup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTlightDown","SR1l_A_250lessMETless350BTlightDown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTlightUp","SR1l_A_250lessMETless350BTlightUp",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTheavyDown","SR1l_A_250lessMETless350BTheavyDown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTheavyUp","SR1l_A_250lessMETless350BTheavyUp",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350PDFdown","SR1l_A_250lessMETless350PDFdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350PDFup","SR1l_A_250lessMETless350PDFup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350alphaSdown","SR1l_A_250lessMETless350alphaSdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350alphaSup","SR1l_A_250lessMETless350alphaSup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350Q2down","SR1l_A_250lessMETless350Q2down",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350Q2up","SR1l_A_250lessMETless350Q2up",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350topPtModelingdown","SR1l_A_250lessMETless350topPtModelingdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350topPtModelingup","SR1l_A_250lessMETless350topPtModelingup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_350lessMETless450","SR1l_A_350lessMETless450",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450LSFdown","SR1l_A_350lessMETless450LSFdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450LSFup","SR1l_A_350lessMETless450LSFup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTlightDown","SR1l_A_350lessMETless450BTlightDown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTlightUp","SR1l_A_350lessMETless450BTlightUp",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTheavyDown","SR1l_A_350lessMETless450BTheavyDown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTheavyUp","SR1l_A_350lessMETless450BTheavyUp",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450PDFdown","SR1l_A_350lessMETless450PDFdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450PDFup","SR1l_A_350lessMETless450PDFup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450alphaSdown","SR1l_A_350lessMETless450alphaSdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450alphaSup","SR1l_A_350lessMETless450alphaSup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450Q2down","SR1l_A_350lessMETless450Q2down",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450Q2up","SR1l_A_350lessMETless450Q2up",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450topPtModelingdown","SR1l_A_350lessMETless450topPtModelingdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450topPtModelingup","SR1l_A_350lessMETless450topPtModelingup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_450lessMETless600","SR1l_A_450lessMETless600",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600LSFdown","SR1l_A_450lessMETless600LSFdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600LSFup","SR1l_A_450lessMETless600LSFup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTlightDown","SR1l_A_450lessMETless600BTlightDown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTlightUp","SR1l_A_450lessMETless600BTlightUp",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTheavyDown","SR1l_A_450lessMETless600BTheavyDown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTheavyUp","SR1l_A_450lessMETless600BTheavyUp",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600PDFdown","SR1l_A_450lessMETless600PDFdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600PDFup","SR1l_A_450lessMETless600PDFup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600alphaSdown","SR1l_A_450lessMETless600alphaSdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600alphaSup","SR1l_A_450lessMETless600alphaSup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600Q2down","SR1l_A_450lessMETless600Q2down",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600Q2up","SR1l_A_450lessMETless600Q2up",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600topPtModelingdown","SR1l_A_450lessMETless600topPtModelingdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600topPtModelingup","SR1l_A_450lessMETless600topPtModelingup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_600lessMETlessInf","SR1l_A_600lessMETlessInf",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfLSFdown","SR1l_A_600lessMETlessInfLSFdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfLSFup","SR1l_A_600lessMETlessInfLSFup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTlightDown","SR1l_A_600lessMETlessInfBTlightDown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTlightUp","SR1l_A_600lessMETlessInfBTlightUp",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTheavyDown","SR1l_A_600lessMETlessInfBTheavyDown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTheavyUp","SR1l_A_600lessMETlessInfBTheavyUp",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfPDFdown","SR1l_A_600lessMETlessInfPDFdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfPDFup","SR1l_A_600lessMETlessInfPDFup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfalphaSdown","SR1l_A_600lessMETlessInfalphaSdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfalphaSup","SR1l_A_600lessMETlessInfalphaSup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfQ2down","SR1l_A_600lessMETlessInfQ2down",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfQ2up","SR1l_A_600lessMETlessInfQ2up",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInftopPtModelingdown","SR1l_A_600lessMETlessInftopPtModelingdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInftopPtModelingup","SR1l_A_600lessMETlessInftopPtModelingup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_B_250lessMETless450","SR1l_B_250lessMETless450",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450LSFdown","SR1l_B_250lessMETless450LSFdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450LSFup","SR1l_B_250lessMETless450LSFup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTlightDown","SR1l_B_250lessMETless450BTlightDown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTlightUp","SR1l_B_250lessMETless450BTlightUp",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTheavyDown","SR1l_B_250lessMETless450BTheavyDown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTheavyUp","SR1l_B_250lessMETless450BTheavyUp",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450PDFdown","SR1l_B_250lessMETless450PDFdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450PDFup","SR1l_B_250lessMETless450PDFup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450alphaSdown","SR1l_B_250lessMETless450alphaSdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450alphaSup","SR1l_B_250lessMETless450alphaSup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450Q2down","SR1l_B_250lessMETless450Q2down",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450Q2up","SR1l_B_250lessMETless450Q2up",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450topPtModelingdown","SR1l_B_250lessMETless450topPtModelingdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450topPtModelingup","SR1l_B_250lessMETless450topPtModelingup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_450lessMETless600","SR1l_B_450lessMETless600",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600LSFdown","SR1l_B_450lessMETless600LSFdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600LSFup","SR1l_B_450lessMETless600LSFup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTlightDown","SR1l_B_450lessMETless600BTlightDown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTlightUp","SR1l_B_450lessMETless600BTlightUp",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTheavyDown","SR1l_B_450lessMETless600BTheavyDown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTheavyUp","SR1l_B_450lessMETless600BTheavyUp",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600PDFdown","SR1l_B_450lessMETless600PDFdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600PDFup","SR1l_B_450lessMETless600PDFup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600alphaSdown","SR1l_B_450lessMETless600alphaSdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600alphaSup","SR1l_B_450lessMETless600alphaSup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600Q2down","SR1l_B_450lessMETless600Q2down",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600Q2up","SR1l_B_450lessMETless600Q2up",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600topPtModelingdown","SR1l_B_450lessMETless600topPtModelingdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600topPtModelingup","SR1l_B_450lessMETless600topPtModelingup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_600lessMETlessInf","SR1l_B_600lessMETlessInf",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfLSFdown","SR1l_B_600lessMETlessInfLSFdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfLSFup","SR1l_B_600lessMETlessInfLSFup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTlightDown","SR1l_B_600lessMETlessInfBTlightDown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTlightUp","SR1l_B_600lessMETlessInfBTlightUp",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTheavyDown","SR1l_B_600lessMETlessInfBTheavyDown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTheavyUp","SR1l_B_600lessMETlessInfBTheavyUp",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfPDFdown","SR1l_B_600lessMETlessInfPDFdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfPDFup","SR1l_B_600lessMETlessInfPDFup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfalphaSdown","SR1l_B_600lessMETlessInfalphaSdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfalphaSup","SR1l_B_600lessMETlessInfalphaSup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfQ2down","SR1l_B_600lessMETlessInfQ2down",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfQ2up","SR1l_B_600lessMETlessInfQ2up",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInftopPtModelingdown","SR1l_B_600lessMETlessInftopPtModelingdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInftopPtModelingup","SR1l_B_600lessMETlessInftopPtModelingup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_C_250lessMETless350","SR1l_C_250lessMETless350",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350LSFdown","SR1l_C_250lessMETless350LSFdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350LSFup","SR1l_C_250lessMETless350LSFup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTlightDown","SR1l_C_250lessMETless350BTlightDown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTlightUp","SR1l_C_250lessMETless350BTlightUp",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTheavyDown","SR1l_C_250lessMETless350BTheavyDown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTheavyUp","SR1l_C_250lessMETless350BTheavyUp",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350PDFdown","SR1l_C_250lessMETless350PDFdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350PDFup","SR1l_C_250lessMETless350PDFup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350alphaSdown","SR1l_C_250lessMETless350alphaSdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350alphaSup","SR1l_C_250lessMETless350alphaSup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350Q2down","SR1l_C_250lessMETless350Q2down",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350Q2up","SR1l_C_250lessMETless350Q2up",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350topPtModelingdown","SR1l_C_250lessMETless350topPtModelingdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350topPtModelingup","SR1l_C_250lessMETless350topPtModelingup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_350lessMETless450","SR1l_C_350lessMETless450",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450LSFdown","SR1l_C_350lessMETless450LSFdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450LSFup","SR1l_C_350lessMETless450LSFup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTlightDown","SR1l_C_350lessMETless450BTlightDown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTlightUp","SR1l_C_350lessMETless450BTlightUp",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTheavyDown","SR1l_C_350lessMETless450BTheavyDown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTheavyUp","SR1l_C_350lessMETless450BTheavyUp",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450PDFdown","SR1l_C_350lessMETless450PDFdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450PDFup","SR1l_C_350lessMETless450PDFup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450alphaSdown","SR1l_C_350lessMETless450alphaSdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450alphaSup","SR1l_C_350lessMETless450alphaSup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450Q2down","SR1l_C_350lessMETless450Q2down",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450Q2up","SR1l_C_350lessMETless450Q2up",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450topPtModelingdown","SR1l_C_350lessMETless450topPtModelingdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450topPtModelingup","SR1l_C_350lessMETless450topPtModelingup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_450lessMETless550","SR1l_C_450lessMETless550",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550LSFdown","SR1l_C_450lessMETless550LSFdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550LSFup","SR1l_C_450lessMETless550LSFup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTlightDown","SR1l_C_450lessMETless550BTlightDown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTlightUp","SR1l_C_450lessMETless550BTlightUp",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTheavyDown","SR1l_C_450lessMETless550BTheavyDown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTheavyUp","SR1l_C_450lessMETless550BTheavyUp",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550PDFdown","SR1l_C_450lessMETless550PDFdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550PDFup","SR1l_C_450lessMETless550PDFup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550alphaSdown","SR1l_C_450lessMETless550alphaSdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550alphaSup","SR1l_C_450lessMETless550alphaSup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550Q2down","SR1l_C_450lessMETless550Q2down",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550Q2up","SR1l_C_450lessMETless550Q2up",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550topPtModelingdown","SR1l_C_450lessMETless550topPtModelingdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550topPtModelingup","SR1l_C_450lessMETless550topPtModelingup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_550lessMETless650","SR1l_C_550lessMETless650",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650LSFdown","SR1l_C_550lessMETless650LSFdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650LSFup","SR1l_C_550lessMETless650LSFup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTlightDown","SR1l_C_550lessMETless650BTlightDown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTlightUp","SR1l_C_550lessMETless650BTlightUp",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTheavyDown","SR1l_C_550lessMETless650BTheavyDown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTheavyUp","SR1l_C_550lessMETless650BTheavyUp",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650PDFdown","SR1l_C_550lessMETless650PDFdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650PDFup","SR1l_C_550lessMETless650PDFup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650alphaSdown","SR1l_C_550lessMETless650alphaSdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650alphaSup","SR1l_C_550lessMETless650alphaSup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650Q2down","SR1l_C_550lessMETless650Q2down",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650Q2up","SR1l_C_550lessMETless650Q2up",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650topPtModelingdown","SR1l_C_550lessMETless650topPtModelingdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650topPtModelingup","SR1l_C_550lessMETless650topPtModelingup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_650lessMETlessInf","SR1l_C_650lessMETlessInf",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfLSFdown","SR1l_C_650lessMETlessInfLSFdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfLSFup","SR1l_C_650lessMETlessInfLSFup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTlightDown","SR1l_C_650lessMETlessInfBTlightDown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTlightUp","SR1l_C_650lessMETlessInfBTlightUp",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTheavyDown","SR1l_C_650lessMETlessInfBTheavyDown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTheavyUp","SR1l_C_650lessMETlessInfBTheavyUp",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfPDFdown","SR1l_C_650lessMETlessInfPDFdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfPDFup","SR1l_C_650lessMETlessInfPDFup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfalphaSdown","SR1l_C_650lessMETlessInfalphaSdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfalphaSup","SR1l_C_650lessMETlessInfalphaSup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfQ2down","SR1l_C_650lessMETlessInfQ2down",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfQ2up","SR1l_C_650lessMETlessInfQ2up",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInftopPtModelingdown","SR1l_C_650lessMETlessInftopPtModelingdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInftopPtModelingup","SR1l_C_650lessMETlessInftopPtModelingup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_D_250lessMETless350","SR1l_D_250lessMETless350",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350LSFdown","SR1l_D_250lessMETless350LSFdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350LSFup","SR1l_D_250lessMETless350LSFup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTlightDown","SR1l_D_250lessMETless350BTlightDown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTlightUp","SR1l_D_250lessMETless350BTlightUp",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTheavyDown","SR1l_D_250lessMETless350BTheavyDown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTheavyUp","SR1l_D_250lessMETless350BTheavyUp",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350PDFdown","SR1l_D_250lessMETless350PDFdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350PDFup","SR1l_D_250lessMETless350PDFup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350alphaSdown","SR1l_D_250lessMETless350alphaSdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350alphaSup","SR1l_D_250lessMETless350alphaSup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350Q2down","SR1l_D_250lessMETless350Q2down",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350Q2up","SR1l_D_250lessMETless350Q2up",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350topPtModelingdown","SR1l_D_250lessMETless350topPtModelingdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350topPtModelingup","SR1l_D_250lessMETless350topPtModelingup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_350lessMETless450","SR1l_D_350lessMETless450",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450LSFdown","SR1l_D_350lessMETless450LSFdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450LSFup","SR1l_D_350lessMETless450LSFup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTlightDown","SR1l_D_350lessMETless450BTlightDown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTlightUp","SR1l_D_350lessMETless450BTlightUp",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTheavyDown","SR1l_D_350lessMETless450BTheavyDown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTheavyUp","SR1l_D_350lessMETless450BTheavyUp",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450PDFdown","SR1l_D_350lessMETless450PDFdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450PDFup","SR1l_D_350lessMETless450PDFup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450alphaSdown","SR1l_D_350lessMETless450alphaSdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450alphaSup","SR1l_D_350lessMETless450alphaSup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450Q2down","SR1l_D_350lessMETless450Q2down",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450Q2up","SR1l_D_350lessMETless450Q2up",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450topPtModelingdown","SR1l_D_350lessMETless450topPtModelingdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450topPtModelingup","SR1l_D_350lessMETless450topPtModelingup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_450lessMETless550","SR1l_D_450lessMETless550",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550LSFdown","SR1l_D_450lessMETless550LSFdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550LSFup","SR1l_D_450lessMETless550LSFup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTlightDown","SR1l_D_450lessMETless550BTlightDown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTlightUp","SR1l_D_450lessMETless550BTlightUp",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTheavyDown","SR1l_D_450lessMETless550BTheavyDown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTheavyUp","SR1l_D_450lessMETless550BTheavyUp",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550PDFdown","SR1l_D_450lessMETless550PDFdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550PDFup","SR1l_D_450lessMETless550PDFup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550alphaSdown","SR1l_D_450lessMETless550alphaSdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550alphaSup","SR1l_D_450lessMETless550alphaSup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550Q2down","SR1l_D_450lessMETless550Q2down",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550Q2up","SR1l_D_450lessMETless550Q2up",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550topPtModelingdown","SR1l_D_450lessMETless550topPtModelingdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550topPtModelingup","SR1l_D_450lessMETless550topPtModelingup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_550lessMETlessInf","SR1l_D_550lessMETlessInf",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfLSFdown","SR1l_D_550lessMETlessInfLSFdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfLSFup","SR1l_D_550lessMETlessInfLSFup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTlightDown","SR1l_D_550lessMETlessInfBTlightDown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTlightUp","SR1l_D_550lessMETlessInfBTlightUp",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTheavyDown","SR1l_D_550lessMETlessInfBTheavyDown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTheavyUp","SR1l_D_550lessMETlessInfBTheavyUp",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfPDFdown","SR1l_D_550lessMETlessInfPDFdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfPDFup","SR1l_D_550lessMETlessInfPDFup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfalphaSdown","SR1l_D_550lessMETlessInfalphaSdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfalphaSup","SR1l_D_550lessMETlessInfalphaSup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfQ2down","SR1l_D_550lessMETlessInfQ2down",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfQ2up","SR1l_D_550lessMETlessInfQ2up",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInftopPtModelingdown","SR1l_D_550lessMETlessInftopPtModelingdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInftopPtModelingup","SR1l_D_550lessMETlessInftopPtModelingup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_E_250lessMETless350","SR1l_E_250lessMETless350",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350LSFdown","SR1l_E_250lessMETless350LSFdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350LSFup","SR1l_E_250lessMETless350LSFup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTlightDown","SR1l_E_250lessMETless350BTlightDown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTlightUp","SR1l_E_250lessMETless350BTlightUp",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTheavyDown","SR1l_E_250lessMETless350BTheavyDown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTheavyUp","SR1l_E_250lessMETless350BTheavyUp",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350PDFdown","SR1l_E_250lessMETless350PDFdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350PDFup","SR1l_E_250lessMETless350PDFup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350alphaSdown","SR1l_E_250lessMETless350alphaSdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350alphaSup","SR1l_E_250lessMETless350alphaSup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350Q2down","SR1l_E_250lessMETless350Q2down",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350Q2up","SR1l_E_250lessMETless350Q2up",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350topPtModelingdown","SR1l_E_250lessMETless350topPtModelingdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350topPtModelingup","SR1l_E_250lessMETless350topPtModelingup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_350lessMETless550","SR1l_E_350lessMETless550",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550LSFdown","SR1l_E_350lessMETless550LSFdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550LSFup","SR1l_E_350lessMETless550LSFup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTlightDown","SR1l_E_350lessMETless550BTlightDown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTlightUp","SR1l_E_350lessMETless550BTlightUp",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTheavyDown","SR1l_E_350lessMETless550BTheavyDown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTheavyUp","SR1l_E_350lessMETless550BTheavyUp",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550PDFdown","SR1l_E_350lessMETless550PDFdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550PDFup","SR1l_E_350lessMETless550PDFup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550alphaSdown","SR1l_E_350lessMETless550alphaSdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550alphaSup","SR1l_E_350lessMETless550alphaSup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550Q2down","SR1l_E_350lessMETless550Q2down",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550Q2up","SR1l_E_350lessMETless550Q2up",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550topPtModelingdown","SR1l_E_350lessMETless550topPtModelingdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550topPtModelingup","SR1l_E_350lessMETless550topPtModelingup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_550lessMETlessInf","SR1l_E_550lessMETlessInf",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfLSFdown","SR1l_E_550lessMETlessInfLSFdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfLSFup","SR1l_E_550lessMETlessInfLSFup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTlightDown","SR1l_E_550lessMETlessInfBTlightDown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTlightUp","SR1l_E_550lessMETlessInfBTlightUp",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTheavyDown","SR1l_E_550lessMETlessInfBTheavyDown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTheavyUp","SR1l_E_550lessMETlessInfBTheavyUp",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfPDFdown","SR1l_E_550lessMETlessInfPDFdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfPDFup","SR1l_E_550lessMETlessInfPDFup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfalphaSdown","SR1l_E_550lessMETlessInfalphaSdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfalphaSup","SR1l_E_550lessMETlessInfalphaSup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfQ2down","SR1l_E_550lessMETlessInfQ2down",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfQ2up","SR1l_E_550lessMETlessInfQ2up",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInftopPtModelingdown","SR1l_E_550lessMETlessInftopPtModelingdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInftopPtModelingup","SR1l_E_550lessMETlessInftopPtModelingup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_F_250lessMETless450","SR1l_F_250lessMETless450",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450LSFdown","SR1l_F_250lessMETless450LSFdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450LSFup","SR1l_F_250lessMETless450LSFup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTlightDown","SR1l_F_250lessMETless450BTlightDown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTlightUp","SR1l_F_250lessMETless450BTlightUp",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTheavyDown","SR1l_F_250lessMETless450BTheavyDown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTheavyUp","SR1l_F_250lessMETless450BTheavyUp",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450PDFdown","SR1l_F_250lessMETless450PDFdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450PDFup","SR1l_F_250lessMETless450PDFup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450alphaSdown","SR1l_F_250lessMETless450alphaSdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450alphaSup","SR1l_F_250lessMETless450alphaSup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450Q2down","SR1l_F_250lessMETless450Q2down",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450Q2up","SR1l_F_250lessMETless450Q2up",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450topPtModelingdown","SR1l_F_250lessMETless450topPtModelingdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450topPtModelingup","SR1l_F_250lessMETless450topPtModelingup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_450lessMETlessInf","SR1l_F_450lessMETlessInf",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfLSFdown","SR1l_F_450lessMETlessInfLSFdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfLSFup","SR1l_F_450lessMETlessInfLSFup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTlightDown","SR1l_F_450lessMETlessInfBTlightDown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTlightUp","SR1l_F_450lessMETlessInfBTlightUp",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTheavyDown","SR1l_F_450lessMETlessInfBTheavyDown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTheavyUp","SR1l_F_450lessMETlessInfBTheavyUp",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfPDFdown","SR1l_F_450lessMETlessInfPDFdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfPDFup","SR1l_F_450lessMETlessInfPDFup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfalphaSdown","SR1l_F_450lessMETlessInfalphaSdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfalphaSup","SR1l_F_450lessMETlessInfalphaSup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfQ2down","SR1l_F_450lessMETlessInfQ2down",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfQ2up","SR1l_F_450lessMETlessInfQ2up",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInftopPtModelingdown","SR1l_F_450lessMETlessInftopPtModelingdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInftopPtModelingup","SR1l_F_450lessMETlessInftopPtModelingup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_G_250lessMETless350","SR1l_G_250lessMETless350",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350LSFdown","SR1l_G_250lessMETless350LSFdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350LSFup","SR1l_G_250lessMETless350LSFup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTlightDown","SR1l_G_250lessMETless350BTlightDown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTlightUp","SR1l_G_250lessMETless350BTlightUp",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTheavyDown","SR1l_G_250lessMETless350BTheavyDown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTheavyUp","SR1l_G_250lessMETless350BTheavyUp",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350PDFdown","SR1l_G_250lessMETless350PDFdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350PDFup","SR1l_G_250lessMETless350PDFup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350alphaSdown","SR1l_G_250lessMETless350alphaSdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350alphaSup","SR1l_G_250lessMETless350alphaSup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350Q2down","SR1l_G_250lessMETless350Q2down",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350Q2up","SR1l_G_250lessMETless350Q2up",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350topPtModelingdown","SR1l_G_250lessMETless350topPtModelingdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350topPtModelingup","SR1l_G_250lessMETless350topPtModelingup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_350lessMETless450","SR1l_G_350lessMETless450",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450LSFdown","SR1l_G_350lessMETless450LSFdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450LSFup","SR1l_G_350lessMETless450LSFup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTlightDown","SR1l_G_350lessMETless450BTlightDown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTlightUp","SR1l_G_350lessMETless450BTlightUp",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTheavyDown","SR1l_G_350lessMETless450BTheavyDown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTheavyUp","SR1l_G_350lessMETless450BTheavyUp",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450PDFdown","SR1l_G_350lessMETless450PDFdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450PDFup","SR1l_G_350lessMETless450PDFup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450alphaSdown","SR1l_G_350lessMETless450alphaSdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450alphaSup","SR1l_G_350lessMETless450alphaSup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450Q2down","SR1l_G_350lessMETless450Q2down",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450Q2up","SR1l_G_350lessMETless450Q2up",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450topPtModelingdown","SR1l_G_350lessMETless450topPtModelingdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450topPtModelingup","SR1l_G_350lessMETless450topPtModelingup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_450lessMETless600","SR1l_G_450lessMETless600",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600LSFdown","SR1l_G_450lessMETless600LSFdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600LSFup","SR1l_G_450lessMETless600LSFup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTlightDown","SR1l_G_450lessMETless600BTlightDown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTlightUp","SR1l_G_450lessMETless600BTlightUp",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTheavyDown","SR1l_G_450lessMETless600BTheavyDown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTheavyUp","SR1l_G_450lessMETless600BTheavyUp",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600PDFdown","SR1l_G_450lessMETless600PDFdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600PDFup","SR1l_G_450lessMETless600PDFup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600alphaSdown","SR1l_G_450lessMETless600alphaSdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600alphaSup","SR1l_G_450lessMETless600alphaSup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600Q2down","SR1l_G_450lessMETless600Q2down",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600Q2up","SR1l_G_450lessMETless600Q2up",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600topPtModelingdown","SR1l_G_450lessMETless600topPtModelingdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600topPtModelingup","SR1l_G_450lessMETless600topPtModelingup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_600lessMETlessInf","SR1l_G_600lessMETlessInf",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfLSFdown","SR1l_G_600lessMETlessInfLSFdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfLSFup","SR1l_G_600lessMETlessInfLSFup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTlightDown","SR1l_G_600lessMETlessInfBTlightDown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTlightUp","SR1l_G_600lessMETlessInfBTlightUp",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTheavyDown","SR1l_G_600lessMETlessInfBTheavyDown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTheavyUp","SR1l_G_600lessMETlessInfBTheavyUp",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfPDFdown","SR1l_G_600lessMETlessInfPDFdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfPDFup","SR1l_G_600lessMETlessInfPDFup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfalphaSdown","SR1l_G_600lessMETlessInfalphaSdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfalphaSup","SR1l_G_600lessMETlessInfalphaSup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfQ2down","SR1l_G_600lessMETlessInfQ2down",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfQ2up","SR1l_G_600lessMETlessInfQ2up",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInftopPtModelingdown","SR1l_G_600lessMETlessInftopPtModelingdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInftopPtModelingup","SR1l_G_600lessMETlessInftopPtModelingup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_H_250lessMETless450","SR1l_H_250lessMETless450",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450LSFdown","SR1l_H_250lessMETless450LSFdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450LSFup","SR1l_H_250lessMETless450LSFup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTlightDown","SR1l_H_250lessMETless450BTlightDown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTlightUp","SR1l_H_250lessMETless450BTlightUp",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTheavyDown","SR1l_H_250lessMETless450BTheavyDown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTheavyUp","SR1l_H_250lessMETless450BTheavyUp",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450PDFdown","SR1l_H_250lessMETless450PDFdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450PDFup","SR1l_H_250lessMETless450PDFup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450alphaSdown","SR1l_H_250lessMETless450alphaSdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450alphaSup","SR1l_H_250lessMETless450alphaSup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450Q2down","SR1l_H_250lessMETless450Q2down",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450Q2up","SR1l_H_250lessMETless450Q2up",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450topPtModelingdown","SR1l_H_250lessMETless450topPtModelingdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450topPtModelingup","SR1l_H_250lessMETless450topPtModelingup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_450lessMETlessInf","SR1l_H_450lessMETlessInf",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfLSFdown","SR1l_H_450lessMETlessInfLSFdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfLSFup","SR1l_H_450lessMETlessInfLSFup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTlightDown","SR1l_H_450lessMETlessInfBTlightDown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTlightUp","SR1l_H_450lessMETlessInfBTlightUp",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTheavyDown","SR1l_H_450lessMETlessInfBTheavyDown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTheavyUp","SR1l_H_450lessMETlessInfBTheavyUp",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfPDFdown","SR1l_H_450lessMETlessInfPDFdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfPDFup","SR1l_H_450lessMETlessInfPDFup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfalphaSdown","SR1l_H_450lessMETlessInfalphaSdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfalphaSup","SR1l_H_450lessMETlessInfalphaSup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfQ2down","SR1l_H_450lessMETlessInfQ2down",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfQ2up","SR1l_H_450lessMETlessInfQ2up",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInftopPtModelingdown","SR1l_H_450lessMETlessInftopPtModelingdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInftopPtModelingup","SR1l_H_450lessMETlessInftopPtModelingup",&SR1l_H_450lessMETlessInf);

                                                                                             

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(36.46);

    Create1DHistos();
    //Add2DHisto("nJets","MET");

fillYieldsVector();
    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;
    nthentry++;

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);


    //cout << "nentries " << myEvent.nentries << " nthentry "<< nthentry << endl;
    checkNegativeYields = false;

    if(nthentry == myEvent.nentries)
    {
        nthentry=0;
        cout << "seting true to check negative " << endl;
        checkNegativeYields = true;
    }
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

    //cout << "pdf down weight " << myEvent.pdf_down_weight << " up weight " << myEvent.pdf_up_weight << endl;
    /*if(  myEvent.pdf_down_weight<0.0 ||  myEvent.pdf_up_weight<0.0 ){
          myEvent.pdf_down_weight = 1.0;
          myEvent.pdf_up_weight = 1.0;
    }*/  
    float weight_pt = reweightTop(myEvent.top_pt, myEvent.atop_pt);
 
    vector<float> weightV;
    weightV.clear();
    float nEvents =  myEvent.wNormalization.at(22);
    //for number of SR
    for(uint32_t SR=0; SR<27; SR++) //@MJ@ TODO nr of sig regions changes
    {

        float w = 0;
        float btagmax = 0;
        //normal
        weightV.push_back(weightLumi);
        //PUdown
        /*if(counter == 1) statnames << "PUdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.weight_PUdown; //@MJ@ TODO PU without any normalization?!
        weightV.push_back(w);
        //PUup
        if(counter == 1) statnames << "PUup"<< endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.weight_PUup;
        weightV.push_back(w);*/
        //LSFdown
        if(counter == 1) statnames << "LSFdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF_down*( nEvents / myEvent.wNormalization.at(30) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ));
        weightV.push_back(w);
        //LSFup
        if(counter == 1) statnames << "LSFup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF_up*( nEvents / myEvent.wNormalization.at(29) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ));
        weightV.push_back(w);
        //BTlightdown
        if(counter == 1) statnames << "BTligtdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_light_DN*( nEvents / myEvent.wNormalization.at(18) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTlightup
        if(counter == 1) statnames << "BTligtup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_light_UP*( nEvents / myEvent.wNormalization.at(16) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTheabydown
        if(counter == 1) statnames << "BTheavydown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_heavy_DN*( nEvents / myEvent.wNormalization.at(17) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTheavyup
        if(counter == 1) statnames << "BTheavyup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_heavy_UP*( nEvents / myEvent.wNormalization.at(15) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //PDFdown
        if(counter == 1) statnames << "PDFdown" << endl;
        w = GetLumi() *  myEvent.scale1fb  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs(myEvent.pdf_down_weight* ( nEvents / myEvent.wNormalization.at(11) ))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        //cout << "pdfd  " << pdfd << " w " << w << endl;
        weightV.push_back(w);
        //PDFup
        if(counter == 1) statnames << "PDFup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * abs(myEvent.pdf_up_weight*( nEvents / myEvent.wNormalization.at(10) ))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //alphaSdown
        if(counter == 1) statnames << "alphaSdown" << endl;
        w = GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_alphas_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(13))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) ; //TODO
        weightV.push_back(w);
        //alphaSup
        if(counter == 1) statnames << "alphaSup" << endl;
        w = GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_alphas_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(12))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) ; //TODO
        weightV.push_back(w);
        //Q2down
        if(counter == 1) statnames << "Q2down" << endl;
        w = GetLumi() *  myEvent.scale1fb  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_q2_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(9))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )  ; //TODO
        weightV.push_back(w);
        //Q2up
        if(counter == 1) statnames << "Q2up" << endl;
        w = GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28)) * myEvent.weight_q2_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(5))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ); //TODO
        weightV.push_back(w);
    /*    //ISRNjetsdown
        if(counter == 1) statnames << "ISRNjetsdown" << endl;
        w = GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_ISRnjets_DN*( nEvents / myEvent.wNormalization.at(27))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) ; //TODO
        weightV.push_back(w);
        //ISRnjetsup
        if(counter == 1) statnames << "ISRNjetsup" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_ISRnjets_UP*( nEvents / myEvent.wNormalization.at(26))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )  ; //TODO
        weightV.push_back(w);*/
       //topptmodelling
        if(counter == 1) statnames << "topPtmodeling" << endl;
        w = GetLumi() *   myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) * weight_pt;
        weightV.push_back(w);
        //topptmodelling obsolete
        if(counter == 1) statnames << "topPtmodeling2" << endl;
        w = GetLumi() *   myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ); //not ready now
        weightV.push_back(w);
    }

    vector<string> theReg;
    GetRegionTagList(&theReg);
    if( weightV.size() != theReg.size())
        throw std::runtime_error("vector of weights does not have same size as regions, the weights will not be correctly assesed");

   
        cout << "pdf weight nu" << nEvents / myEvent.wNormalization.at(10) << " nd " << nEvents / myEvent.wNormalization.at(11) << " down " << myEvent.pdf_down_weight << " up " << myEvent.pdf_up_weight << " weihgt " << weightLumi << endl;
 
    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    //if (currentProcessType == "signal") weight = 1;
    AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields, weightV, true);

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

  /*     TFile fi2("genw.root","RECREATE");
       w_pdfnu->Write();
       w_pdfnd->Write();
       w_pdfu->Write();
       w_pdfd->Write();
       fi2.Close();
*/

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

vector<string> totYield = { "SR1l_A_250lessMETless350" , "SR1l_A_250lessMETless350LSFdown" , "SR1l_A_250lessMETless350LSFup" , "SR1l_A_250lessMETless350BTlightDown" , "SR1l_A_250lessMETless350BTlightUp" , "SR1l_A_250lessMETless350BTheavyDown" , "SR1l_A_250lessMETless350BTheavyUp" , "SR1l_A_250lessMETless350PDFdown" , "SR1l_A_250lessMETless350PDFup" , "SR1l_A_250lessMETless350alphaSdown" , "SR1l_A_250lessMETless350alphaSup" , "SR1l_A_250lessMETless350Q2down" , "SR1l_A_250lessMETless350Q2up" , "SR1l_A_250lessMETless350topPtModelingdown" , "SR1l_A_250lessMETless350topPtModelingup" , "SR1l_A_350lessMETless450" , "SR1l_A_350lessMETless450LSFdown" , "SR1l_A_350lessMETless450LSFup" , "SR1l_A_350lessMETless450BTlightDown" , "SR1l_A_350lessMETless450BTlightUp" , "SR1l_A_350lessMETless450BTheavyDown" , "SR1l_A_350lessMETless450BTheavyUp" , "SR1l_A_350lessMETless450PDFdown" , "SR1l_A_350lessMETless450PDFup" , "SR1l_A_350lessMETless450alphaSdown" , "SR1l_A_350lessMETless450alphaSup" , "SR1l_A_350lessMETless450Q2down" , "SR1l_A_350lessMETless450Q2up" , "SR1l_A_350lessMETless450topPtModelingdown" , "SR1l_A_350lessMETless450topPtModelingup" , "SR1l_A_450lessMETless600" , "SR1l_A_450lessMETless600LSFdown" , "SR1l_A_450lessMETless600LSFup" , "SR1l_A_450lessMETless600BTlightDown" , "SR1l_A_450lessMETless600BTlightUp" , "SR1l_A_450lessMETless600BTheavyDown" , "SR1l_A_450lessMETless600BTheavyUp" , "SR1l_A_450lessMETless600PDFdown" , "SR1l_A_450lessMETless600PDFup" , "SR1l_A_450lessMETless600alphaSdown" , "SR1l_A_450lessMETless600alphaSup" , "SR1l_A_450lessMETless600Q2down" , "SR1l_A_450lessMETless600Q2up" , "SR1l_A_450lessMETless600topPtModelingdown" , "SR1l_A_450lessMETless600topPtModelingup" , "SR1l_A_600lessMETlessInf" , "SR1l_A_600lessMETlessInfLSFdown" , "SR1l_A_600lessMETlessInfLSFup" , "SR1l_A_600lessMETlessInfBTlightDown" , "SR1l_A_600lessMETlessInfBTlightUp" , "SR1l_A_600lessMETlessInfBTheavyDown" , "SR1l_A_600lessMETlessInfBTheavyUp" , "SR1l_A_600lessMETlessInfPDFdown" , "SR1l_A_600lessMETlessInfPDFup" , "SR1l_A_600lessMETlessInfalphaSdown" , "SR1l_A_600lessMETlessInfalphaSup" , "SR1l_A_600lessMETlessInfQ2down" , "SR1l_A_600lessMETlessInfQ2up" , "SR1l_A_600lessMETlessInftopPtModelingdown" , "SR1l_A_600lessMETlessInftopPtModelingup" , "SR1l_B_250lessMETless450" , "SR1l_B_250lessMETless450LSFdown" , "SR1l_B_250lessMETless450LSFup" , "SR1l_B_250lessMETless450BTlightDown" , "SR1l_B_250lessMETless450BTlightUp" , "SR1l_B_250lessMETless450BTheavyDown" , "SR1l_B_250lessMETless450BTheavyUp" , "SR1l_B_250lessMETless450PDFdown" , "SR1l_B_250lessMETless450PDFup" , "SR1l_B_250lessMETless450alphaSdown" , "SR1l_B_250lessMETless450alphaSup" , "SR1l_B_250lessMETless450Q2down" , "SR1l_B_250lessMETless450Q2up" , "SR1l_B_250lessMETless450topPtModelingdown" , "SR1l_B_250lessMETless450topPtModelingup" , "SR1l_B_450lessMETless600" , "SR1l_B_450lessMETless600LSFdown" , "SR1l_B_450lessMETless600LSFup" , "SR1l_B_450lessMETless600BTlightDown" , "SR1l_B_450lessMETless600BTlightUp" , "SR1l_B_450lessMETless600BTheavyDown" , "SR1l_B_450lessMETless600BTheavyUp" , "SR1l_B_450lessMETless600PDFdown" , "SR1l_B_450lessMETless600PDFup" , "SR1l_B_450lessMETless600alphaSdown" , "SR1l_B_450lessMETless600alphaSup" , "SR1l_B_450lessMETless600Q2down" , "SR1l_B_450lessMETless600Q2up" , "SR1l_B_450lessMETless600topPtModelingdown" , "SR1l_B_450lessMETless600topPtModelingup" , "SR1l_B_600lessMETlessInf" , "SR1l_B_600lessMETlessInfLSFdown" , "SR1l_B_600lessMETlessInfLSFup" , "SR1l_B_600lessMETlessInfBTlightDown" , "SR1l_B_600lessMETlessInfBTlightUp" , "SR1l_B_600lessMETlessInfBTheavyDown" , "SR1l_B_600lessMETlessInfBTheavyUp" , "SR1l_B_600lessMETlessInfPDFdown" , "SR1l_B_600lessMETlessInfPDFup" , "SR1l_B_600lessMETlessInfalphaSdown" , "SR1l_B_600lessMETlessInfalphaSup" , "SR1l_B_600lessMETlessInfQ2down" , "SR1l_B_600lessMETlessInfQ2up" , "SR1l_B_600lessMETlessInftopPtModelingdown" , "SR1l_B_600lessMETlessInftopPtModelingup" , "SR1l_C_250lessMETless350" , "SR1l_C_250lessMETless350LSFdown" , "SR1l_C_250lessMETless350LSFup" , "SR1l_C_250lessMETless350BTlightDown" , "SR1l_C_250lessMETless350BTlightUp" , "SR1l_C_250lessMETless350BTheavyDown" , "SR1l_C_250lessMETless350BTheavyUp" , "SR1l_C_250lessMETless350PDFdown" , "SR1l_C_250lessMETless350PDFup" , "SR1l_C_250lessMETless350alphaSdown" , "SR1l_C_250lessMETless350alphaSup" , "SR1l_C_250lessMETless350Q2down" , "SR1l_C_250lessMETless350Q2up" , "SR1l_C_250lessMETless350topPtModelingdown" , "SR1l_C_250lessMETless350topPtModelingup" , "SR1l_C_350lessMETless450" , "SR1l_C_350lessMETless450LSFdown" , "SR1l_C_350lessMETless450LSFup" , "SR1l_C_350lessMETless450BTlightDown" , "SR1l_C_350lessMETless450BTlightUp" , "SR1l_C_350lessMETless450BTheavyDown" , "SR1l_C_350lessMETless450BTheavyUp" , "SR1l_C_350lessMETless450PDFdown" , "SR1l_C_350lessMETless450PDFup" , "SR1l_C_350lessMETless450alphaSdown" , "SR1l_C_350lessMETless450alphaSup" , "SR1l_C_350lessMETless450Q2down" , "SR1l_C_350lessMETless450Q2up" , "SR1l_C_350lessMETless450topPtModelingdown" , "SR1l_C_350lessMETless450topPtModelingup" , "SR1l_C_450lessMETless550" , "SR1l_C_450lessMETless550LSFdown" , "SR1l_C_450lessMETless550LSFup" , "SR1l_C_450lessMETless550BTlightDown" , "SR1l_C_450lessMETless550BTlightUp" , "SR1l_C_450lessMETless550BTheavyDown" , "SR1l_C_450lessMETless550BTheavyUp" , "SR1l_C_450lessMETless550PDFdown" , "SR1l_C_450lessMETless550PDFup" , "SR1l_C_450lessMETless550alphaSdown" , "SR1l_C_450lessMETless550alphaSup" , "SR1l_C_450lessMETless550Q2down" , "SR1l_C_450lessMETless550Q2up" , "SR1l_C_450lessMETless550topPtModelingdown" , "SR1l_C_450lessMETless550topPtModelingup" , "SR1l_C_550lessMETless650" , "SR1l_C_550lessMETless650LSFdown" , "SR1l_C_550lessMETless650LSFup" , "SR1l_C_550lessMETless650BTlightDown" , "SR1l_C_550lessMETless650BTlightUp" , "SR1l_C_550lessMETless650BTheavyDown" , "SR1l_C_550lessMETless650BTheavyUp" , "SR1l_C_550lessMETless650PDFdown" , "SR1l_C_550lessMETless650PDFup" , "SR1l_C_550lessMETless650alphaSdown" , "SR1l_C_550lessMETless650alphaSup" , "SR1l_C_550lessMETless650Q2down" , "SR1l_C_550lessMETless650Q2up" , "SR1l_C_550lessMETless650topPtModelingdown" , "SR1l_C_550lessMETless650topPtModelingup" , "SR1l_C_650lessMETlessInf" , "SR1l_C_650lessMETlessInfLSFdown" , "SR1l_C_650lessMETlessInfLSFup" , "SR1l_C_650lessMETlessInfBTlightDown" , "SR1l_C_650lessMETlessInfBTlightUp" , "SR1l_C_650lessMETlessInfBTheavyDown" , "SR1l_C_650lessMETlessInfBTheavyUp" , "SR1l_C_650lessMETlessInfPDFdown" , "SR1l_C_650lessMETlessInfPDFup" , "SR1l_C_650lessMETlessInfalphaSdown" , "SR1l_C_650lessMETlessInfalphaSup" , "SR1l_C_650lessMETlessInfQ2down" , "SR1l_C_650lessMETlessInfQ2up" , "SR1l_C_650lessMETlessInftopPtModelingdown" , "SR1l_C_650lessMETlessInftopPtModelingup" , "SR1l_D_250lessMETless350" , "SR1l_D_250lessMETless350LSFdown" , "SR1l_D_250lessMETless350LSFup" , "SR1l_D_250lessMETless350BTlightDown" , "SR1l_D_250lessMETless350BTlightUp" , "SR1l_D_250lessMETless350BTheavyDown" , "SR1l_D_250lessMETless350BTheavyUp" , "SR1l_D_250lessMETless350PDFdown" , "SR1l_D_250lessMETless350PDFup" , "SR1l_D_250lessMETless350alphaSdown" , "SR1l_D_250lessMETless350alphaSup" , "SR1l_D_250lessMETless350Q2down" , "SR1l_D_250lessMETless350Q2up" , "SR1l_D_250lessMETless350topPtModelingdown" , "SR1l_D_250lessMETless350topPtModelingup" , "SR1l_D_350lessMETless450" , "SR1l_D_350lessMETless450LSFdown" , "SR1l_D_350lessMETless450LSFup" , "SR1l_D_350lessMETless450BTlightDown" , "SR1l_D_350lessMETless450BTlightUp" , "SR1l_D_350lessMETless450BTheavyDown" , "SR1l_D_350lessMETless450BTheavyUp" , "SR1l_D_350lessMETless450PDFdown" , "SR1l_D_350lessMETless450PDFup" , "SR1l_D_350lessMETless450alphaSdown" , "SR1l_D_350lessMETless450alphaSup" , "SR1l_D_350lessMETless450Q2down" , "SR1l_D_350lessMETless450Q2up" , "SR1l_D_350lessMETless450topPtModelingdown" , "SR1l_D_350lessMETless450topPtModelingup" , "SR1l_D_450lessMETless550" , "SR1l_D_450lessMETless550LSFdown" , "SR1l_D_450lessMETless550LSFup" , "SR1l_D_450lessMETless550BTlightDown" , "SR1l_D_450lessMETless550BTlightUp" , "SR1l_D_450lessMETless550BTheavyDown" , "SR1l_D_450lessMETless550BTheavyUp" , "SR1l_D_450lessMETless550PDFdown" , "SR1l_D_450lessMETless550PDFup" , "SR1l_D_450lessMETless550alphaSdown" , "SR1l_D_450lessMETless550alphaSup" , "SR1l_D_450lessMETless550Q2down" , "SR1l_D_450lessMETless550Q2up" , "SR1l_D_450lessMETless550topPtModelingdown" , "SR1l_D_450lessMETless550topPtModelingup" , "SR1l_D_550lessMETlessInf" , "SR1l_D_550lessMETlessInfLSFdown" , "SR1l_D_550lessMETlessInfLSFup" , "SR1l_D_550lessMETlessInfBTlightDown" , "SR1l_D_550lessMETlessInfBTlightUp" , "SR1l_D_550lessMETlessInfBTheavyDown" , "SR1l_D_550lessMETlessInfBTheavyUp" , "SR1l_D_550lessMETlessInfPDFdown" , "SR1l_D_550lessMETlessInfPDFup" , "SR1l_D_550lessMETlessInfalphaSdown" , "SR1l_D_550lessMETlessInfalphaSup" , "SR1l_D_550lessMETlessInfQ2down" , "SR1l_D_550lessMETlessInfQ2up" , "SR1l_D_550lessMETlessInftopPtModelingdown" , "SR1l_D_550lessMETlessInftopPtModelingup" , "SR1l_E_250lessMETless350" , "SR1l_E_250lessMETless350LSFdown" , "SR1l_E_250lessMETless350LSFup" , "SR1l_E_250lessMETless350BTlightDown" , "SR1l_E_250lessMETless350BTlightUp" , "SR1l_E_250lessMETless350BTheavyDown" , "SR1l_E_250lessMETless350BTheavyUp" , "SR1l_E_250lessMETless350PDFdown" , "SR1l_E_250lessMETless350PDFup" , "SR1l_E_250lessMETless350alphaSdown" , "SR1l_E_250lessMETless350alphaSup" , "SR1l_E_250lessMETless350Q2down" , "SR1l_E_250lessMETless350Q2up" , "SR1l_E_250lessMETless350topPtModelingdown" , "SR1l_E_250lessMETless350topPtModelingup" , "SR1l_E_350lessMETless550" , "SR1l_E_350lessMETless550LSFdown" , "SR1l_E_350lessMETless550LSFup" , "SR1l_E_350lessMETless550BTlightDown" , "SR1l_E_350lessMETless550BTlightUp" , "SR1l_E_350lessMETless550BTheavyDown" , "SR1l_E_350lessMETless550BTheavyUp" , "SR1l_E_350lessMETless550PDFdown" , "SR1l_E_350lessMETless550PDFup" , "SR1l_E_350lessMETless550alphaSdown" , "SR1l_E_350lessMETless550alphaSup" , "SR1l_E_350lessMETless550Q2down" , "SR1l_E_350lessMETless550Q2up" , "SR1l_E_350lessMETless550topPtModelingdown" , "SR1l_E_350lessMETless550topPtModelingup" , "SR1l_E_550lessMETlessInf" , "SR1l_E_550lessMETlessInfLSFdown" , "SR1l_E_550lessMETlessInfLSFup" , "SR1l_E_550lessMETlessInfBTlightDown" , "SR1l_E_550lessMETlessInfBTlightUp" , "SR1l_E_550lessMETlessInfBTheavyDown" , "SR1l_E_550lessMETlessInfBTheavyUp" , "SR1l_E_550lessMETlessInfPDFdown" , "SR1l_E_550lessMETlessInfPDFup" , "SR1l_E_550lessMETlessInfalphaSdown" , "SR1l_E_550lessMETlessInfalphaSup" , "SR1l_E_550lessMETlessInfQ2down" , "SR1l_E_550lessMETlessInfQ2up" , "SR1l_E_550lessMETlessInftopPtModelingdown" , "SR1l_E_550lessMETlessInftopPtModelingup" , "SR1l_F_250lessMETless450" , "SR1l_F_250lessMETless450LSFdown" , "SR1l_F_250lessMETless450LSFup" , "SR1l_F_250lessMETless450BTlightDown" , "SR1l_F_250lessMETless450BTlightUp" , "SR1l_F_250lessMETless450BTheavyDown" , "SR1l_F_250lessMETless450BTheavyUp" , "SR1l_F_250lessMETless450PDFdown" , "SR1l_F_250lessMETless450PDFup" , "SR1l_F_250lessMETless450alphaSdown" , "SR1l_F_250lessMETless450alphaSup" , "SR1l_F_250lessMETless450Q2down" , "SR1l_F_250lessMETless450Q2up" , "SR1l_F_250lessMETless450topPtModelingdown" , "SR1l_F_250lessMETless450topPtModelingup" , "SR1l_F_450lessMETlessInf" , "SR1l_F_450lessMETlessInfLSFdown" , "SR1l_F_450lessMETlessInfLSFup" , "SR1l_F_450lessMETlessInfBTlightDown" , "SR1l_F_450lessMETlessInfBTlightUp" , "SR1l_F_450lessMETlessInfBTheavyDown" , "SR1l_F_450lessMETlessInfBTheavyUp" , "SR1l_F_450lessMETlessInfPDFdown" , "SR1l_F_450lessMETlessInfPDFup" , "SR1l_F_450lessMETlessInfalphaSdown" , "SR1l_F_450lessMETlessInfalphaSup" , "SR1l_F_450lessMETlessInfQ2down" , "SR1l_F_450lessMETlessInfQ2up" , "SR1l_F_450lessMETlessInftopPtModelingdown" , "SR1l_F_450lessMETlessInftopPtModelingup" , "SR1l_G_250lessMETless350" , "SR1l_G_250lessMETless350LSFdown" , "SR1l_G_250lessMETless350LSFup" , "SR1l_G_250lessMETless350BTlightDown" , "SR1l_G_250lessMETless350BTlightUp" , "SR1l_G_250lessMETless350BTheavyDown" , "SR1l_G_250lessMETless350BTheavyUp" , "SR1l_G_250lessMETless350PDFdown" , "SR1l_G_250lessMETless350PDFup" , "SR1l_G_250lessMETless350alphaSdown" , "SR1l_G_250lessMETless350alphaSup" , "SR1l_G_250lessMETless350Q2down" , "SR1l_G_250lessMETless350Q2up" , "SR1l_G_250lessMETless350topPtModelingdown" , "SR1l_G_250lessMETless350topPtModelingup" , "SR1l_G_350lessMETless450" , "SR1l_G_350lessMETless450LSFdown" , "SR1l_G_350lessMETless450LSFup" , "SR1l_G_350lessMETless450BTlightDown" , "SR1l_G_350lessMETless450BTlightUp" , "SR1l_G_350lessMETless450BTheavyDown" , "SR1l_G_350lessMETless450BTheavyUp" , "SR1l_G_350lessMETless450PDFdown" , "SR1l_G_350lessMETless450PDFup" , "SR1l_G_350lessMETless450alphaSdown" , "SR1l_G_350lessMETless450alphaSup" , "SR1l_G_350lessMETless450Q2down" , "SR1l_G_350lessMETless450Q2up" , "SR1l_G_350lessMETless450topPtModelingdown" , "SR1l_G_350lessMETless450topPtModelingup" , "SR1l_G_450lessMETless600" , "SR1l_G_450lessMETless600LSFdown" , "SR1l_G_450lessMETless600LSFup" , "SR1l_G_450lessMETless600BTlightDown" , "SR1l_G_450lessMETless600BTlightUp" , "SR1l_G_450lessMETless600BTheavyDown" , "SR1l_G_450lessMETless600BTheavyUp" , "SR1l_G_450lessMETless600PDFdown" , "SR1l_G_450lessMETless600PDFup" , "SR1l_G_450lessMETless600alphaSdown" , "SR1l_G_450lessMETless600alphaSup" , "SR1l_G_450lessMETless600Q2down" , "SR1l_G_450lessMETless600Q2up" , "SR1l_G_450lessMETless600topPtModelingdown" , "SR1l_G_450lessMETless600topPtModelingup" , "SR1l_G_600lessMETlessInf" , "SR1l_G_600lessMETlessInfLSFdown" , "SR1l_G_600lessMETlessInfLSFup" , "SR1l_G_600lessMETlessInfBTlightDown" , "SR1l_G_600lessMETlessInfBTlightUp" , "SR1l_G_600lessMETlessInfBTheavyDown" , "SR1l_G_600lessMETlessInfBTheavyUp" , "SR1l_G_600lessMETlessInfPDFdown" , "SR1l_G_600lessMETlessInfPDFup" , "SR1l_G_600lessMETlessInfalphaSdown" , "SR1l_G_600lessMETlessInfalphaSup" , "SR1l_G_600lessMETlessInfQ2down" , "SR1l_G_600lessMETlessInfQ2up" , "SR1l_G_600lessMETlessInftopPtModelingdown" , "SR1l_G_600lessMETlessInftopPtModelingup" , "SR1l_H_250lessMETless450" , "SR1l_H_250lessMETless450LSFdown" , "SR1l_H_250lessMETless450LSFup" , "SR1l_H_250lessMETless450BTlightDown" , "SR1l_H_250lessMETless450BTlightUp" , "SR1l_H_250lessMETless450BTheavyDown" , "SR1l_H_250lessMETless450BTheavyUp" , "SR1l_H_250lessMETless450PDFdown" , "SR1l_H_250lessMETless450PDFup" , "SR1l_H_250lessMETless450alphaSdown" , "SR1l_H_250lessMETless450alphaSup" , "SR1l_H_250lessMETless450Q2down" , "SR1l_H_250lessMETless450Q2up" , "SR1l_H_250lessMETless450topPtModelingdown" , "SR1l_H_250lessMETless450topPtModelingup" , "SR1l_H_450lessMETlessInf" , "SR1l_H_450lessMETlessInfLSFdown" , "SR1l_H_450lessMETlessInfLSFup" , "SR1l_H_450lessMETlessInfBTlightDown" , "SR1l_H_450lessMETlessInfBTlightUp" , "SR1l_H_450lessMETlessInfBTheavyDown" , "SR1l_H_450lessMETlessInfBTheavyUp" , "SR1l_H_450lessMETlessInfPDFdown" , "SR1l_H_450lessMETlessInfPDFup" , "SR1l_H_450lessMETlessInfalphaSdown" , "SR1l_H_450lessMETlessInfalphaSup" , "SR1l_H_450lessMETlessInfQ2down" , "SR1l_H_450lessMETlessInfQ2up" , "SR1l_H_450lessMETlessInftopPtModelingdown" , "SR1l_H_450lessMETlessInftopPtModelingup"};


    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print("yieldMor.tab", 6);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex("yieldMor.tex", 6);

    ofstream sigfile("signalRegMor.txt");
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
    prod.produceTFTable("yieldMor.tab", "lostLeptonTF");
    TFProducer prod2(tfreg, "singleLepton", "CR1l");
    prod2.produceTFTable("yieldMor.tab", "singleLeptonTF");
    TFProducer prod3(tfreg, "singleLeptonFromT", "CR1l");
    prod3.produceTFTable("yieldMor.tab", "singleLeptonFromTTF");

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
        float all_weights = lumi*  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        //cout << "cs: " << myEvent.crossSection << "scale 1 fb " << myEvent.scale1fb << " weight total: " << all_weights << endl;
        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
        return all_weights;
    }

    
    float reweightTop(float topPt, float atopPt)
    {
        return sqrt(exp(0.0615-(0.0005*topPt))*exp(0.0615-(0.0005*atopPt))); 
    }
