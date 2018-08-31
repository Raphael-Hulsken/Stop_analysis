#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

//#define USE_TTZ
#define USE_WZ
//#define USE_ZZ

#define USE_VAR_BASELINE
#define USE_LEP1
#define USE_LEP2
#define USE_JETS
#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_GLOBAL_VAR

#include "../../common/TFFactory.h"
#include "../../Selection/moriond.h"

using namespace std;


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
string outputName = "";


float getWeight(string currentProcessType, float lumi);
float reweightTop(float topPt, float atopPt);
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

    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11/NotSkimmed";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11";
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/haweber_v22/";
    totalNumberOfWorkers = 1;

    #ifdef USE_TTZ
    AddProcessClass("ttZ", "ttZ", "background", kBlue);
    	AddDataset("ttZJets_13TeV_madgraphMLM","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_1","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_2","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_3","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_4","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_5","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_6","ttZ",0,0);
        outputName = "yieldMorttZ";
    #endif
    	
    #ifdef USE_ZZWZ
    AddProcessClass("ZZWZ", "ZZWZ", "background", kBlue);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","ZZWZ",0,0);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZWZ",0,0);
        outputName = "yieldMorZZWZ";
    #endif

    #ifdef USE_ZZ
    AddProcessClass("ZZ", "ZZ", "background", kBlue);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZ",0,0);
        outputName = "yieldMorZZ";
    #endif

    #ifdef USE_WZ
    AddProcessClass("WZ", "WZ", "background", kBlue);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","WZ",0,0);
        outputName = "yieldMorWZ";
    #endif

	//465 regions
AddRegion("SR1l_A_250lessMETless350","SR1l_A_250lessMETless350",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350LSFdown","SR1l_A_250lessMETless350LSFdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350LSFup","SR1l_A_250lessMETless350LSFup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTlightDown","SR1l_A_250lessMETless350BTlightDown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTlightUp","SR1l_A_250lessMETless350BTlightUp",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTheavyDown","SR1l_A_250lessMETless350BTheavyDown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350BTheavyUp","SR1l_A_250lessMETless350BTheavyUp",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350PUdown","SR1l_A_250lessMETless350PUdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350PUup","SR1l_A_250lessMETless350PUup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350PDFdown","SR1l_A_250lessMETless350PDFdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350PDFup","SR1l_A_250lessMETless350PDFup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350alphaSdown","SR1l_A_250lessMETless350alphaSdown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350alphaSup","SR1l_A_250lessMETless350alphaSup",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350Q2down","SR1l_A_250lessMETless350Q2down",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350Q2up","SR1l_A_250lessMETless350Q2up",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350ISRnjetsDown","SR1l_A_250lessMETless350ISRnjetsDown",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_250lessMETless350ISRnjetsUp","SR1l_A_250lessMETless350ISRnjetsUp",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_350lessMETless450","SR1l_A_350lessMETless450",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450LSFdown","SR1l_A_350lessMETless450LSFdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450LSFup","SR1l_A_350lessMETless450LSFup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTlightDown","SR1l_A_350lessMETless450BTlightDown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTlightUp","SR1l_A_350lessMETless450BTlightUp",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTheavyDown","SR1l_A_350lessMETless450BTheavyDown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450BTheavyUp","SR1l_A_350lessMETless450BTheavyUp",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450PUdown","SR1l_A_350lessMETless450PUdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450PUup","SR1l_A_350lessMETless450PUup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450PDFdown","SR1l_A_350lessMETless450PDFdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450PDFup","SR1l_A_350lessMETless450PDFup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450alphaSdown","SR1l_A_350lessMETless450alphaSdown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450alphaSup","SR1l_A_350lessMETless450alphaSup",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450Q2down","SR1l_A_350lessMETless450Q2down",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450Q2up","SR1l_A_350lessMETless450Q2up",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450ISRnjetsDown","SR1l_A_350lessMETless450ISRnjetsDown",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_350lessMETless450ISRnjetsUp","SR1l_A_350lessMETless450ISRnjetsUp",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_450lessMETless600","SR1l_A_450lessMETless600",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600LSFdown","SR1l_A_450lessMETless600LSFdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600LSFup","SR1l_A_450lessMETless600LSFup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTlightDown","SR1l_A_450lessMETless600BTlightDown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTlightUp","SR1l_A_450lessMETless600BTlightUp",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTheavyDown","SR1l_A_450lessMETless600BTheavyDown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600BTheavyUp","SR1l_A_450lessMETless600BTheavyUp",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600PUdown","SR1l_A_450lessMETless600PUdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600PUup","SR1l_A_450lessMETless600PUup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600PDFdown","SR1l_A_450lessMETless600PDFdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600PDFup","SR1l_A_450lessMETless600PDFup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600alphaSdown","SR1l_A_450lessMETless600alphaSdown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600alphaSup","SR1l_A_450lessMETless600alphaSup",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600Q2down","SR1l_A_450lessMETless600Q2down",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600Q2up","SR1l_A_450lessMETless600Q2up",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600ISRnjetsDown","SR1l_A_450lessMETless600ISRnjetsDown",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_450lessMETless600ISRnjetsUp","SR1l_A_450lessMETless600ISRnjetsUp",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_600lessMETlessInf","SR1l_A_600lessMETlessInf",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfLSFdown","SR1l_A_600lessMETlessInfLSFdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfLSFup","SR1l_A_600lessMETlessInfLSFup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTlightDown","SR1l_A_600lessMETlessInfBTlightDown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTlightUp","SR1l_A_600lessMETlessInfBTlightUp",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTheavyDown","SR1l_A_600lessMETlessInfBTheavyDown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfBTheavyUp","SR1l_A_600lessMETlessInfBTheavyUp",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfPUdown","SR1l_A_600lessMETlessInfPUdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfPUup","SR1l_A_600lessMETlessInfPUup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfPDFdown","SR1l_A_600lessMETlessInfPDFdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfPDFup","SR1l_A_600lessMETlessInfPDFup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfalphaSdown","SR1l_A_600lessMETlessInfalphaSdown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfalphaSup","SR1l_A_600lessMETlessInfalphaSup",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfQ2down","SR1l_A_600lessMETlessInfQ2down",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfQ2up","SR1l_A_600lessMETlessInfQ2up",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfISRnjetsDown","SR1l_A_600lessMETlessInfISRnjetsDown",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_A_600lessMETlessInfISRnjetsUp","SR1l_A_600lessMETlessInfISRnjetsUp",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_B_250lessMETless450","SR1l_B_250lessMETless450",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450LSFdown","SR1l_B_250lessMETless450LSFdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450LSFup","SR1l_B_250lessMETless450LSFup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTlightDown","SR1l_B_250lessMETless450BTlightDown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTlightUp","SR1l_B_250lessMETless450BTlightUp",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTheavyDown","SR1l_B_250lessMETless450BTheavyDown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450BTheavyUp","SR1l_B_250lessMETless450BTheavyUp",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450PUdown","SR1l_B_250lessMETless450PUdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450PUup","SR1l_B_250lessMETless450PUup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450PDFdown","SR1l_B_250lessMETless450PDFdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450PDFup","SR1l_B_250lessMETless450PDFup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450alphaSdown","SR1l_B_250lessMETless450alphaSdown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450alphaSup","SR1l_B_250lessMETless450alphaSup",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450Q2down","SR1l_B_250lessMETless450Q2down",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450Q2up","SR1l_B_250lessMETless450Q2up",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450ISRnjetsDown","SR1l_B_250lessMETless450ISRnjetsDown",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_250lessMETless450ISRnjetsUp","SR1l_B_250lessMETless450ISRnjetsUp",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_450lessMETless600","SR1l_B_450lessMETless600",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600LSFdown","SR1l_B_450lessMETless600LSFdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600LSFup","SR1l_B_450lessMETless600LSFup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTlightDown","SR1l_B_450lessMETless600BTlightDown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTlightUp","SR1l_B_450lessMETless600BTlightUp",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTheavyDown","SR1l_B_450lessMETless600BTheavyDown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600BTheavyUp","SR1l_B_450lessMETless600BTheavyUp",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600PUdown","SR1l_B_450lessMETless600PUdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600PUup","SR1l_B_450lessMETless600PUup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600PDFdown","SR1l_B_450lessMETless600PDFdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600PDFup","SR1l_B_450lessMETless600PDFup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600alphaSdown","SR1l_B_450lessMETless600alphaSdown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600alphaSup","SR1l_B_450lessMETless600alphaSup",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600Q2down","SR1l_B_450lessMETless600Q2down",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600Q2up","SR1l_B_450lessMETless600Q2up",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600ISRnjetsDown","SR1l_B_450lessMETless600ISRnjetsDown",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_450lessMETless600ISRnjetsUp","SR1l_B_450lessMETless600ISRnjetsUp",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_600lessMETlessInf","SR1l_B_600lessMETlessInf",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfLSFdown","SR1l_B_600lessMETlessInfLSFdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfLSFup","SR1l_B_600lessMETlessInfLSFup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTlightDown","SR1l_B_600lessMETlessInfBTlightDown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTlightUp","SR1l_B_600lessMETlessInfBTlightUp",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTheavyDown","SR1l_B_600lessMETlessInfBTheavyDown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfBTheavyUp","SR1l_B_600lessMETlessInfBTheavyUp",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfPUdown","SR1l_B_600lessMETlessInfPUdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfPUup","SR1l_B_600lessMETlessInfPUup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfPDFdown","SR1l_B_600lessMETlessInfPDFdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfPDFup","SR1l_B_600lessMETlessInfPDFup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfalphaSdown","SR1l_B_600lessMETlessInfalphaSdown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfalphaSup","SR1l_B_600lessMETlessInfalphaSup",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfQ2down","SR1l_B_600lessMETlessInfQ2down",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfQ2up","SR1l_B_600lessMETlessInfQ2up",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfISRnjetsDown","SR1l_B_600lessMETlessInfISRnjetsDown",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_B_600lessMETlessInfISRnjetsUp","SR1l_B_600lessMETlessInfISRnjetsUp",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_C_250lessMETless350","SR1l_C_250lessMETless350",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350LSFdown","SR1l_C_250lessMETless350LSFdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350LSFup","SR1l_C_250lessMETless350LSFup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTlightDown","SR1l_C_250lessMETless350BTlightDown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTlightUp","SR1l_C_250lessMETless350BTlightUp",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTheavyDown","SR1l_C_250lessMETless350BTheavyDown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350BTheavyUp","SR1l_C_250lessMETless350BTheavyUp",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350PUdown","SR1l_C_250lessMETless350PUdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350PUup","SR1l_C_250lessMETless350PUup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350PDFdown","SR1l_C_250lessMETless350PDFdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350PDFup","SR1l_C_250lessMETless350PDFup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350alphaSdown","SR1l_C_250lessMETless350alphaSdown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350alphaSup","SR1l_C_250lessMETless350alphaSup",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350Q2down","SR1l_C_250lessMETless350Q2down",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350Q2up","SR1l_C_250lessMETless350Q2up",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350ISRnjetsDown","SR1l_C_250lessMETless350ISRnjetsDown",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_250lessMETless350ISRnjetsUp","SR1l_C_250lessMETless350ISRnjetsUp",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_350lessMETless450","SR1l_C_350lessMETless450",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450LSFdown","SR1l_C_350lessMETless450LSFdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450LSFup","SR1l_C_350lessMETless450LSFup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTlightDown","SR1l_C_350lessMETless450BTlightDown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTlightUp","SR1l_C_350lessMETless450BTlightUp",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTheavyDown","SR1l_C_350lessMETless450BTheavyDown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450BTheavyUp","SR1l_C_350lessMETless450BTheavyUp",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450PUdown","SR1l_C_350lessMETless450PUdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450PUup","SR1l_C_350lessMETless450PUup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450PDFdown","SR1l_C_350lessMETless450PDFdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450PDFup","SR1l_C_350lessMETless450PDFup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450alphaSdown","SR1l_C_350lessMETless450alphaSdown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450alphaSup","SR1l_C_350lessMETless450alphaSup",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450Q2down","SR1l_C_350lessMETless450Q2down",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450Q2up","SR1l_C_350lessMETless450Q2up",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450ISRnjetsDown","SR1l_C_350lessMETless450ISRnjetsDown",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_350lessMETless450ISRnjetsUp","SR1l_C_350lessMETless450ISRnjetsUp",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_450lessMETless550","SR1l_C_450lessMETless550",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550LSFdown","SR1l_C_450lessMETless550LSFdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550LSFup","SR1l_C_450lessMETless550LSFup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTlightDown","SR1l_C_450lessMETless550BTlightDown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTlightUp","SR1l_C_450lessMETless550BTlightUp",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTheavyDown","SR1l_C_450lessMETless550BTheavyDown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550BTheavyUp","SR1l_C_450lessMETless550BTheavyUp",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550PUdown","SR1l_C_450lessMETless550PUdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550PUup","SR1l_C_450lessMETless550PUup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550PDFdown","SR1l_C_450lessMETless550PDFdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550PDFup","SR1l_C_450lessMETless550PDFup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550alphaSdown","SR1l_C_450lessMETless550alphaSdown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550alphaSup","SR1l_C_450lessMETless550alphaSup",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550Q2down","SR1l_C_450lessMETless550Q2down",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550Q2up","SR1l_C_450lessMETless550Q2up",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550ISRnjetsDown","SR1l_C_450lessMETless550ISRnjetsDown",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_450lessMETless550ISRnjetsUp","SR1l_C_450lessMETless550ISRnjetsUp",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_550lessMETless650","SR1l_C_550lessMETless650",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650LSFdown","SR1l_C_550lessMETless650LSFdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650LSFup","SR1l_C_550lessMETless650LSFup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTlightDown","SR1l_C_550lessMETless650BTlightDown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTlightUp","SR1l_C_550lessMETless650BTlightUp",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTheavyDown","SR1l_C_550lessMETless650BTheavyDown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650BTheavyUp","SR1l_C_550lessMETless650BTheavyUp",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650PUdown","SR1l_C_550lessMETless650PUdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650PUup","SR1l_C_550lessMETless650PUup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650PDFdown","SR1l_C_550lessMETless650PDFdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650PDFup","SR1l_C_550lessMETless650PDFup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650alphaSdown","SR1l_C_550lessMETless650alphaSdown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650alphaSup","SR1l_C_550lessMETless650alphaSup",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650Q2down","SR1l_C_550lessMETless650Q2down",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650Q2up","SR1l_C_550lessMETless650Q2up",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650ISRnjetsDown","SR1l_C_550lessMETless650ISRnjetsDown",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_550lessMETless650ISRnjetsUp","SR1l_C_550lessMETless650ISRnjetsUp",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_650lessMETlessInf","SR1l_C_650lessMETlessInf",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfLSFdown","SR1l_C_650lessMETlessInfLSFdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfLSFup","SR1l_C_650lessMETlessInfLSFup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTlightDown","SR1l_C_650lessMETlessInfBTlightDown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTlightUp","SR1l_C_650lessMETlessInfBTlightUp",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTheavyDown","SR1l_C_650lessMETlessInfBTheavyDown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfBTheavyUp","SR1l_C_650lessMETlessInfBTheavyUp",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfPUdown","SR1l_C_650lessMETlessInfPUdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfPUup","SR1l_C_650lessMETlessInfPUup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfPDFdown","SR1l_C_650lessMETlessInfPDFdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfPDFup","SR1l_C_650lessMETlessInfPDFup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfalphaSdown","SR1l_C_650lessMETlessInfalphaSdown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfalphaSup","SR1l_C_650lessMETlessInfalphaSup",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfQ2down","SR1l_C_650lessMETlessInfQ2down",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfQ2up","SR1l_C_650lessMETlessInfQ2up",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfISRnjetsDown","SR1l_C_650lessMETlessInfISRnjetsDown",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_C_650lessMETlessInfISRnjetsUp","SR1l_C_650lessMETlessInfISRnjetsUp",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_D_250lessMETless350","SR1l_D_250lessMETless350",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350LSFdown","SR1l_D_250lessMETless350LSFdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350LSFup","SR1l_D_250lessMETless350LSFup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTlightDown","SR1l_D_250lessMETless350BTlightDown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTlightUp","SR1l_D_250lessMETless350BTlightUp",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTheavyDown","SR1l_D_250lessMETless350BTheavyDown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350BTheavyUp","SR1l_D_250lessMETless350BTheavyUp",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350PUdown","SR1l_D_250lessMETless350PUdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350PUup","SR1l_D_250lessMETless350PUup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350PDFdown","SR1l_D_250lessMETless350PDFdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350PDFup","SR1l_D_250lessMETless350PDFup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350alphaSdown","SR1l_D_250lessMETless350alphaSdown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350alphaSup","SR1l_D_250lessMETless350alphaSup",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350Q2down","SR1l_D_250lessMETless350Q2down",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350Q2up","SR1l_D_250lessMETless350Q2up",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350ISRnjetsDown","SR1l_D_250lessMETless350ISRnjetsDown",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_250lessMETless350ISRnjetsUp","SR1l_D_250lessMETless350ISRnjetsUp",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_350lessMETless450","SR1l_D_350lessMETless450",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450LSFdown","SR1l_D_350lessMETless450LSFdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450LSFup","SR1l_D_350lessMETless450LSFup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTlightDown","SR1l_D_350lessMETless450BTlightDown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTlightUp","SR1l_D_350lessMETless450BTlightUp",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTheavyDown","SR1l_D_350lessMETless450BTheavyDown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450BTheavyUp","SR1l_D_350lessMETless450BTheavyUp",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450PUdown","SR1l_D_350lessMETless450PUdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450PUup","SR1l_D_350lessMETless450PUup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450PDFdown","SR1l_D_350lessMETless450PDFdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450PDFup","SR1l_D_350lessMETless450PDFup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450alphaSdown","SR1l_D_350lessMETless450alphaSdown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450alphaSup","SR1l_D_350lessMETless450alphaSup",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450Q2down","SR1l_D_350lessMETless450Q2down",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450Q2up","SR1l_D_350lessMETless450Q2up",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450ISRnjetsDown","SR1l_D_350lessMETless450ISRnjetsDown",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_350lessMETless450ISRnjetsUp","SR1l_D_350lessMETless450ISRnjetsUp",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_450lessMETless550","SR1l_D_450lessMETless550",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550LSFdown","SR1l_D_450lessMETless550LSFdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550LSFup","SR1l_D_450lessMETless550LSFup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTlightDown","SR1l_D_450lessMETless550BTlightDown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTlightUp","SR1l_D_450lessMETless550BTlightUp",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTheavyDown","SR1l_D_450lessMETless550BTheavyDown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550BTheavyUp","SR1l_D_450lessMETless550BTheavyUp",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550PUdown","SR1l_D_450lessMETless550PUdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550PUup","SR1l_D_450lessMETless550PUup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550PDFdown","SR1l_D_450lessMETless550PDFdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550PDFup","SR1l_D_450lessMETless550PDFup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550alphaSdown","SR1l_D_450lessMETless550alphaSdown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550alphaSup","SR1l_D_450lessMETless550alphaSup",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550Q2down","SR1l_D_450lessMETless550Q2down",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550Q2up","SR1l_D_450lessMETless550Q2up",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550ISRnjetsDown","SR1l_D_450lessMETless550ISRnjetsDown",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_450lessMETless550ISRnjetsUp","SR1l_D_450lessMETless550ISRnjetsUp",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_550lessMETlessInf","SR1l_D_550lessMETlessInf",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfLSFdown","SR1l_D_550lessMETlessInfLSFdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfLSFup","SR1l_D_550lessMETlessInfLSFup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTlightDown","SR1l_D_550lessMETlessInfBTlightDown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTlightUp","SR1l_D_550lessMETlessInfBTlightUp",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTheavyDown","SR1l_D_550lessMETlessInfBTheavyDown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfBTheavyUp","SR1l_D_550lessMETlessInfBTheavyUp",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfPUdown","SR1l_D_550lessMETlessInfPUdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfPUup","SR1l_D_550lessMETlessInfPUup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfPDFdown","SR1l_D_550lessMETlessInfPDFdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfPDFup","SR1l_D_550lessMETlessInfPDFup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfalphaSdown","SR1l_D_550lessMETlessInfalphaSdown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfalphaSup","SR1l_D_550lessMETlessInfalphaSup",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfQ2down","SR1l_D_550lessMETlessInfQ2down",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfQ2up","SR1l_D_550lessMETlessInfQ2up",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfISRnjetsDown","SR1l_D_550lessMETlessInfISRnjetsDown",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_D_550lessMETlessInfISRnjetsUp","SR1l_D_550lessMETlessInfISRnjetsUp",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_E_250lessMETless350","SR1l_E_250lessMETless350",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350LSFdown","SR1l_E_250lessMETless350LSFdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350LSFup","SR1l_E_250lessMETless350LSFup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTlightDown","SR1l_E_250lessMETless350BTlightDown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTlightUp","SR1l_E_250lessMETless350BTlightUp",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTheavyDown","SR1l_E_250lessMETless350BTheavyDown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350BTheavyUp","SR1l_E_250lessMETless350BTheavyUp",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350PUdown","SR1l_E_250lessMETless350PUdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350PUup","SR1l_E_250lessMETless350PUup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350PDFdown","SR1l_E_250lessMETless350PDFdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350PDFup","SR1l_E_250lessMETless350PDFup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350alphaSdown","SR1l_E_250lessMETless350alphaSdown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350alphaSup","SR1l_E_250lessMETless350alphaSup",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350Q2down","SR1l_E_250lessMETless350Q2down",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350Q2up","SR1l_E_250lessMETless350Q2up",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350ISRnjetsDown","SR1l_E_250lessMETless350ISRnjetsDown",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_250lessMETless350ISRnjetsUp","SR1l_E_250lessMETless350ISRnjetsUp",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_350lessMETless550","SR1l_E_350lessMETless550",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550LSFdown","SR1l_E_350lessMETless550LSFdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550LSFup","SR1l_E_350lessMETless550LSFup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTlightDown","SR1l_E_350lessMETless550BTlightDown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTlightUp","SR1l_E_350lessMETless550BTlightUp",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTheavyDown","SR1l_E_350lessMETless550BTheavyDown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550BTheavyUp","SR1l_E_350lessMETless550BTheavyUp",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550PUdown","SR1l_E_350lessMETless550PUdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550PUup","SR1l_E_350lessMETless550PUup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550PDFdown","SR1l_E_350lessMETless550PDFdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550PDFup","SR1l_E_350lessMETless550PDFup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550alphaSdown","SR1l_E_350lessMETless550alphaSdown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550alphaSup","SR1l_E_350lessMETless550alphaSup",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550Q2down","SR1l_E_350lessMETless550Q2down",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550Q2up","SR1l_E_350lessMETless550Q2up",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550ISRnjetsDown","SR1l_E_350lessMETless550ISRnjetsDown",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_350lessMETless550ISRnjetsUp","SR1l_E_350lessMETless550ISRnjetsUp",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_550lessMETlessInf","SR1l_E_550lessMETlessInf",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfLSFdown","SR1l_E_550lessMETlessInfLSFdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfLSFup","SR1l_E_550lessMETlessInfLSFup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTlightDown","SR1l_E_550lessMETlessInfBTlightDown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTlightUp","SR1l_E_550lessMETlessInfBTlightUp",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTheavyDown","SR1l_E_550lessMETlessInfBTheavyDown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfBTheavyUp","SR1l_E_550lessMETlessInfBTheavyUp",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfPUdown","SR1l_E_550lessMETlessInfPUdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfPUup","SR1l_E_550lessMETlessInfPUup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfPDFdown","SR1l_E_550lessMETlessInfPDFdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfPDFup","SR1l_E_550lessMETlessInfPDFup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfalphaSdown","SR1l_E_550lessMETlessInfalphaSdown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfalphaSup","SR1l_E_550lessMETlessInfalphaSup",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfQ2down","SR1l_E_550lessMETlessInfQ2down",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfQ2up","SR1l_E_550lessMETlessInfQ2up",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfISRnjetsDown","SR1l_E_550lessMETlessInfISRnjetsDown",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_E_550lessMETlessInfISRnjetsUp","SR1l_E_550lessMETlessInfISRnjetsUp",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_F_250lessMETless450","SR1l_F_250lessMETless450",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450LSFdown","SR1l_F_250lessMETless450LSFdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450LSFup","SR1l_F_250lessMETless450LSFup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTlightDown","SR1l_F_250lessMETless450BTlightDown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTlightUp","SR1l_F_250lessMETless450BTlightUp",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTheavyDown","SR1l_F_250lessMETless450BTheavyDown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450BTheavyUp","SR1l_F_250lessMETless450BTheavyUp",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450PUdown","SR1l_F_250lessMETless450PUdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450PUup","SR1l_F_250lessMETless450PUup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450PDFdown","SR1l_F_250lessMETless450PDFdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450PDFup","SR1l_F_250lessMETless450PDFup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450alphaSdown","SR1l_F_250lessMETless450alphaSdown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450alphaSup","SR1l_F_250lessMETless450alphaSup",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450Q2down","SR1l_F_250lessMETless450Q2down",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450Q2up","SR1l_F_250lessMETless450Q2up",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450ISRnjetsDown","SR1l_F_250lessMETless450ISRnjetsDown",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_250lessMETless450ISRnjetsUp","SR1l_F_250lessMETless450ISRnjetsUp",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_450lessMETlessInf","SR1l_F_450lessMETlessInf",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfLSFdown","SR1l_F_450lessMETlessInfLSFdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfLSFup","SR1l_F_450lessMETlessInfLSFup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTlightDown","SR1l_F_450lessMETlessInfBTlightDown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTlightUp","SR1l_F_450lessMETlessInfBTlightUp",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTheavyDown","SR1l_F_450lessMETlessInfBTheavyDown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfBTheavyUp","SR1l_F_450lessMETlessInfBTheavyUp",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfPUdown","SR1l_F_450lessMETlessInfPUdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfPUup","SR1l_F_450lessMETlessInfPUup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfPDFdown","SR1l_F_450lessMETlessInfPDFdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfPDFup","SR1l_F_450lessMETlessInfPDFup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfalphaSdown","SR1l_F_450lessMETlessInfalphaSdown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfalphaSup","SR1l_F_450lessMETlessInfalphaSup",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfQ2down","SR1l_F_450lessMETlessInfQ2down",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfQ2up","SR1l_F_450lessMETlessInfQ2up",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfISRnjetsDown","SR1l_F_450lessMETlessInfISRnjetsDown",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_F_450lessMETlessInfISRnjetsUp","SR1l_F_450lessMETlessInfISRnjetsUp",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_G_250lessMETless350","SR1l_G_250lessMETless350",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350LSFdown","SR1l_G_250lessMETless350LSFdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350LSFup","SR1l_G_250lessMETless350LSFup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTlightDown","SR1l_G_250lessMETless350BTlightDown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTlightUp","SR1l_G_250lessMETless350BTlightUp",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTheavyDown","SR1l_G_250lessMETless350BTheavyDown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350BTheavyUp","SR1l_G_250lessMETless350BTheavyUp",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350PUdown","SR1l_G_250lessMETless350PUdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350PUup","SR1l_G_250lessMETless350PUup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350PDFdown","SR1l_G_250lessMETless350PDFdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350PDFup","SR1l_G_250lessMETless350PDFup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350alphaSdown","SR1l_G_250lessMETless350alphaSdown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350alphaSup","SR1l_G_250lessMETless350alphaSup",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350Q2down","SR1l_G_250lessMETless350Q2down",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350Q2up","SR1l_G_250lessMETless350Q2up",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350ISRnjetsDown","SR1l_G_250lessMETless350ISRnjetsDown",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_250lessMETless350ISRnjetsUp","SR1l_G_250lessMETless350ISRnjetsUp",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_350lessMETless450","SR1l_G_350lessMETless450",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450LSFdown","SR1l_G_350lessMETless450LSFdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450LSFup","SR1l_G_350lessMETless450LSFup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTlightDown","SR1l_G_350lessMETless450BTlightDown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTlightUp","SR1l_G_350lessMETless450BTlightUp",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTheavyDown","SR1l_G_350lessMETless450BTheavyDown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450BTheavyUp","SR1l_G_350lessMETless450BTheavyUp",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450PUdown","SR1l_G_350lessMETless450PUdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450PUup","SR1l_G_350lessMETless450PUup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450PDFdown","SR1l_G_350lessMETless450PDFdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450PDFup","SR1l_G_350lessMETless450PDFup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450alphaSdown","SR1l_G_350lessMETless450alphaSdown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450alphaSup","SR1l_G_350lessMETless450alphaSup",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450Q2down","SR1l_G_350lessMETless450Q2down",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450Q2up","SR1l_G_350lessMETless450Q2up",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450ISRnjetsDown","SR1l_G_350lessMETless450ISRnjetsDown",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_350lessMETless450ISRnjetsUp","SR1l_G_350lessMETless450ISRnjetsUp",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_450lessMETless600","SR1l_G_450lessMETless600",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600LSFdown","SR1l_G_450lessMETless600LSFdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600LSFup","SR1l_G_450lessMETless600LSFup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTlightDown","SR1l_G_450lessMETless600BTlightDown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTlightUp","SR1l_G_450lessMETless600BTlightUp",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTheavyDown","SR1l_G_450lessMETless600BTheavyDown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600BTheavyUp","SR1l_G_450lessMETless600BTheavyUp",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600PUdown","SR1l_G_450lessMETless600PUdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600PUup","SR1l_G_450lessMETless600PUup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600PDFdown","SR1l_G_450lessMETless600PDFdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600PDFup","SR1l_G_450lessMETless600PDFup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600alphaSdown","SR1l_G_450lessMETless600alphaSdown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600alphaSup","SR1l_G_450lessMETless600alphaSup",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600Q2down","SR1l_G_450lessMETless600Q2down",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600Q2up","SR1l_G_450lessMETless600Q2up",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600ISRnjetsDown","SR1l_G_450lessMETless600ISRnjetsDown",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_450lessMETless600ISRnjetsUp","SR1l_G_450lessMETless600ISRnjetsUp",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_600lessMETlessInf","SR1l_G_600lessMETlessInf",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfLSFdown","SR1l_G_600lessMETlessInfLSFdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfLSFup","SR1l_G_600lessMETlessInfLSFup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTlightDown","SR1l_G_600lessMETlessInfBTlightDown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTlightUp","SR1l_G_600lessMETlessInfBTlightUp",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTheavyDown","SR1l_G_600lessMETlessInfBTheavyDown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfBTheavyUp","SR1l_G_600lessMETlessInfBTheavyUp",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfPUdown","SR1l_G_600lessMETlessInfPUdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfPUup","SR1l_G_600lessMETlessInfPUup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfPDFdown","SR1l_G_600lessMETlessInfPDFdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfPDFup","SR1l_G_600lessMETlessInfPDFup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfalphaSdown","SR1l_G_600lessMETlessInfalphaSdown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfalphaSup","SR1l_G_600lessMETlessInfalphaSup",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfQ2down","SR1l_G_600lessMETlessInfQ2down",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfQ2up","SR1l_G_600lessMETlessInfQ2up",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfISRnjetsDown","SR1l_G_600lessMETlessInfISRnjetsDown",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_G_600lessMETlessInfISRnjetsUp","SR1l_G_600lessMETlessInfISRnjetsUp",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_H_250lessMETless450","SR1l_H_250lessMETless450",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450LSFdown","SR1l_H_250lessMETless450LSFdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450LSFup","SR1l_H_250lessMETless450LSFup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTlightDown","SR1l_H_250lessMETless450BTlightDown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTlightUp","SR1l_H_250lessMETless450BTlightUp",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTheavyDown","SR1l_H_250lessMETless450BTheavyDown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450BTheavyUp","SR1l_H_250lessMETless450BTheavyUp",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450PUdown","SR1l_H_250lessMETless450PUdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450PUup","SR1l_H_250lessMETless450PUup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450PDFdown","SR1l_H_250lessMETless450PDFdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450PDFup","SR1l_H_250lessMETless450PDFup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450alphaSdown","SR1l_H_250lessMETless450alphaSdown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450alphaSup","SR1l_H_250lessMETless450alphaSup",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450Q2down","SR1l_H_250lessMETless450Q2down",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450Q2up","SR1l_H_250lessMETless450Q2up",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450ISRnjetsDown","SR1l_H_250lessMETless450ISRnjetsDown",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_250lessMETless450ISRnjetsUp","SR1l_H_250lessMETless450ISRnjetsUp",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_450lessMETlessInf","SR1l_H_450lessMETlessInf",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfLSFdown","SR1l_H_450lessMETlessInfLSFdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfLSFup","SR1l_H_450lessMETlessInfLSFup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTlightDown","SR1l_H_450lessMETlessInfBTlightDown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTlightUp","SR1l_H_450lessMETlessInfBTlightUp",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTheavyDown","SR1l_H_450lessMETlessInfBTheavyDown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfBTheavyUp","SR1l_H_450lessMETlessInfBTheavyUp",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfPUdown","SR1l_H_450lessMETlessInfPUdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfPUup","SR1l_H_450lessMETlessInfPUup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfPDFdown","SR1l_H_450lessMETlessInfPDFdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfPDFup","SR1l_H_450lessMETlessInfPDFup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfalphaSdown","SR1l_H_450lessMETlessInfalphaSdown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfalphaSup","SR1l_H_450lessMETlessInfalphaSup",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfQ2down","SR1l_H_450lessMETlessInfQ2down",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfQ2up","SR1l_H_450lessMETlessInfQ2up",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfISRnjetsDown","SR1l_H_450lessMETlessInfISRnjetsDown",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_H_450lessMETlessInfISRnjetsUp","SR1l_H_450lessMETlessInfISRnjetsUp",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_I_250lessMETless350","SR1l_I_250lessMETless350",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350LSFdown","SR1l_I_250lessMETless350LSFdown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350LSFup","SR1l_I_250lessMETless350LSFup",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350BTlightDown","SR1l_I_250lessMETless350BTlightDown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350BTlightUp","SR1l_I_250lessMETless350BTlightUp",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350BTheavyDown","SR1l_I_250lessMETless350BTheavyDown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350BTheavyUp","SR1l_I_250lessMETless350BTheavyUp",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350PUdown","SR1l_I_250lessMETless350PUdown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350PUup","SR1l_I_250lessMETless350PUup",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350PDFdown","SR1l_I_250lessMETless350PDFdown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350PDFup","SR1l_I_250lessMETless350PDFup",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350alphaSdown","SR1l_I_250lessMETless350alphaSdown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350alphaSup","SR1l_I_250lessMETless350alphaSup",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350Q2down","SR1l_I_250lessMETless350Q2down",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350Q2up","SR1l_I_250lessMETless350Q2up",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350ISRnjetsDown","SR1l_I_250lessMETless350ISRnjetsDown",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_250lessMETless350ISRnjetsUp","SR1l_I_250lessMETless350ISRnjetsUp",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_350lessMETless450","SR1l_I_350lessMETless450",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450LSFdown","SR1l_I_350lessMETless450LSFdown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450LSFup","SR1l_I_350lessMETless450LSFup",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450BTlightDown","SR1l_I_350lessMETless450BTlightDown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450BTlightUp","SR1l_I_350lessMETless450BTlightUp",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450BTheavyDown","SR1l_I_350lessMETless450BTheavyDown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450BTheavyUp","SR1l_I_350lessMETless450BTheavyUp",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450PUdown","SR1l_I_350lessMETless450PUdown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450PUup","SR1l_I_350lessMETless450PUup",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450PDFdown","SR1l_I_350lessMETless450PDFdown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450PDFup","SR1l_I_350lessMETless450PDFup",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450alphaSdown","SR1l_I_350lessMETless450alphaSdown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450alphaSup","SR1l_I_350lessMETless450alphaSup",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450Q2down","SR1l_I_350lessMETless450Q2down",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450Q2up","SR1l_I_350lessMETless450Q2up",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450ISRnjetsDown","SR1l_I_350lessMETless450ISRnjetsDown",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_350lessMETless450ISRnjetsUp","SR1l_I_350lessMETless450ISRnjetsUp",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_450lessMETless550","SR1l_I_450lessMETless550",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550LSFdown","SR1l_I_450lessMETless550LSFdown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550LSFup","SR1l_I_450lessMETless550LSFup",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550BTlightDown","SR1l_I_450lessMETless550BTlightDown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550BTlightUp","SR1l_I_450lessMETless550BTlightUp",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550BTheavyDown","SR1l_I_450lessMETless550BTheavyDown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550BTheavyUp","SR1l_I_450lessMETless550BTheavyUp",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550PUdown","SR1l_I_450lessMETless550PUdown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550PUup","SR1l_I_450lessMETless550PUup",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550PDFdown","SR1l_I_450lessMETless550PDFdown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550PDFup","SR1l_I_450lessMETless550PDFup",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550alphaSdown","SR1l_I_450lessMETless550alphaSdown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550alphaSup","SR1l_I_450lessMETless550alphaSup",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550Q2down","SR1l_I_450lessMETless550Q2down",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550Q2up","SR1l_I_450lessMETless550Q2up",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550ISRnjetsDown","SR1l_I_450lessMETless550ISRnjetsDown",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_450lessMETless550ISRnjetsUp","SR1l_I_450lessMETless550ISRnjetsUp",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_550lessMETlessInf","SR1l_I_550lessMETlessInf",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfLSFdown","SR1l_I_550lessMETlessInfLSFdown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfLSFup","SR1l_I_550lessMETlessInfLSFup",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfBTlightDown","SR1l_I_550lessMETlessInfBTlightDown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfBTlightUp","SR1l_I_550lessMETlessInfBTlightUp",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfBTheavyDown","SR1l_I_550lessMETlessInfBTheavyDown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfBTheavyUp","SR1l_I_550lessMETlessInfBTheavyUp",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfPUdown","SR1l_I_550lessMETlessInfPUdown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfPUup","SR1l_I_550lessMETlessInfPUup",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfPDFdown","SR1l_I_550lessMETlessInfPDFdown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfPDFup","SR1l_I_550lessMETlessInfPDFup",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfalphaSdown","SR1l_I_550lessMETlessInfalphaSdown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfalphaSup","SR1l_I_550lessMETlessInfalphaSup",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfQ2down","SR1l_I_550lessMETlessInfQ2down",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfQ2up","SR1l_I_550lessMETlessInfQ2up",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfISRnjetsDown","SR1l_I_550lessMETlessInfISRnjetsDown",&SR1l_I_550lessMETlessInf);
AddRegion("SR1l_I_550lessMETlessInfISRnjetsUp","SR1l_I_550lessMETlessInfISRnjetsUp",&SR1l_I_550lessMETlessInf);

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(35.867);

    Create1DHistos();

    fillYieldsVector();
    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;
    nthentry++;

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);


    checkNegativeYields = false;

    if(nthentry == myEvent.nentries)
    {
        nthentry=0;
        cout << "seting true to check negative " << endl;
        checkNegativeYields = true;
    }
    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);

    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
    if( ( currentProcessClass == "ttZ" || currentProcessClass == "ZZWZ" || currentProcessClass == "ZZ"|| currentProcessClass == "WZ" ) && !(myEvent.isZtoNuNu) )
         currentProcessClass = "";


    //@MJ@ TODO I hate myself for this, but no better solution foud
    //computation of up/down weights
    //NOTICE, important is to fill only weight histo and only have 1 process class

    float weight_pt = 1;
    if( currentProcessClass == "ttZ")
        weight_pt = reweightTop(myEvent.top_pt, myEvent.atop_pt);


    float nEvents =  myEvent.wNormalization.at(22);
    float ISRNJ = 1; 
    float ISRNJ_UP = 1; 
    float ISRNJ_DOWN = 1; 
    if( currentProcessClass == "ttZ")
    {
        ISRNJ = myEvent.weight_ISRnjets*( nEvents / myEvent.wNormalization.at(25));
        ISRNJ_UP = myEvent.weight_ISRnjets_UP*( nEvents / myEvent.wNormalization.at(26));
        ISRNJ_DOWN = myEvent.weight_ISRnjets_DN*( nEvents / myEvent.wNormalization.at(27));
    }
    float weightLumi = getWeight(currentProcessType, GetLumi())*ISRNJ; //@MJ@ TODO cross section form file?!
 
    vector<float> weightV;
    weightV.clear();
    //for number of SR
    for(uint32_t SR=0; SR<31; SR++) //@MJ@ TODO nr of sig regions changes
    {

        float w = 0;
        float btagmax = 0;
        //normal
        weightV.push_back(weightLumi);
        //LSFdown
        if(counter == 1) statnames << "lepSFDN" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_lepSF_down*( nEvents / myEvent.wNormalization.at(30)) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14));
        weightV.push_back(w);
        //LSFup
        if(counter == 1) statnames << "lepSFUP" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_lepSF_up*( nEvents / myEvent.wNormalization.at(29)) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) );
        weightV.push_back(w);
        //BTlightdown
        if(counter == 1) statnames << "btagLightDN" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_btagsf_light_DN*( nEvents / myEvent.wNormalization.at(18) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTlightup
        if(counter == 1) statnames << "btagLightUP" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_btagsf_light_UP*( nEvents / myEvent.wNormalization.at(16) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTheabydown
        if(counter == 1) statnames << "btagHeavyDN" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_btagsf_heavy_DN*( nEvents / myEvent.wNormalization.at(17) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTheavyup
        if(counter == 1) statnames << "btagHeavyUP" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_btagsf_heavy_UP*( nEvents / myEvent.wNormalization.at(15) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //PUdown
        if(counter == 1) statnames << "PUdown" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_PUdown; //@MJ@ TODO PU without any normalization?!
        weightV.push_back(w);
        //PUup
        if(counter == 1) statnames << "PUup"<< endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_PUup;
        weightV.push_back(w);
        //PDFdown
        if(counter == 1) statnames << "pdfDN" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs((myEvent.pdf_down_weight/myEvent.genweights->at(0)) * (  myEvent.wNormalization.at(1) / myEvent.wNormalization.at(11) ))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) );
        //cout << "pdfd  " << pdfd << " w " << w << endl;
        weightV.push_back(w);
        //PDFup
        if(counter == 1) statnames << "pdfUP" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * abs((myEvent.pdf_up_weight/myEvent.genweights->at(0)) * (  myEvent.wNormalization.at(1)/ myEvent.wNormalization.at(10) ))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) );
        weightV.push_back(w);
        //alphaSdown
        if(counter == 1) statnames << "alphaSDN" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * abs(myEvent.weight_alphas_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(13)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ; //TODO
        weightV.push_back(w);
        //alphaSup
        if(counter == 1) statnames << "alphaSUP" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs(myEvent.weight_alphas_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(12)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ; //TODO
        weightV.push_back(w);
        //Q2down
        if(counter == 1) statnames << "Q2DN" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs(myEvent.weight_q2_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(9)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )  ; //TODO
        weightV.push_back(w);
        //Q2up
        if(counter == 1) statnames << "Q2UP" << endl;
        w = GetLumi() *  myEvent.scale1fb *ISRNJ* myEvent.weight_PU *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28)) * abs(myEvent.weight_q2_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(5)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ); //TODO
        weightV.push_back(w);
        //ISRNjetsdown
        if(counter == 1) statnames << "ISRnjetsDown" << endl;
        w = GetLumi() *  myEvent.scale1fb* myEvent.weight_PU *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) *ISRNJ_DOWN * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ; //TODO
        weightV.push_back(w);
        //ISRnjetsup
        if(counter == 1) statnames << "ISRnjetsUp" << endl;
        w = GetLumi() *  myEvent.scale1fb* myEvent.weight_PU * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) *ISRNJ_UP* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )  ; //TODO
        weightV.push_back(w);
       //topptmodelling
        //if(counter == 1) statnames << "topPtModeling" << endl;
        //w = GetLumi() *   myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) * weight_pt;
        //weightV.push_back(w);
        //topptmodelling obsolete
        //if(counter == 1) statnames << "topPtmodeling2" << endl;
        //w = GetLumi() *   myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ); //not ready now
        //weightV.push_back(w);
    }

    //cout << "pass med b-tag " << myEvent.ak4pfjets_passMEDbtag->at(0) << endl;

    vector<string> theReg;
    GetRegionTagList(&theReg);
    if( weightV.size() != theReg.size())
        throw std::runtime_error("vector of weights does not have same size as regions, the weights will not be correctly assesed");

   //cout << "scale 1fb " << myEvent.scale1fb << " kfactor " << myEvent.kfactor << " xsec " << myEvent.crossSection << " 1/weight " << (myEvent.scale1fb)/(myEvent.kfactor*myEvent.crossSection) << " pdf down " << abs(myEvent.pdf_down_weight* ( nEvents / myEvent.wNormalization.at(11) )) << " pdf up " << abs(myEvent.pdf_up_weight*( nEvents / myEvent.wNormalization.at(10) )) << " gen[0]" << myEvent.genweights->at(0) << endl; 
 
    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields, weightV, true);


    if(counter % 10000 == 0)
    {
        cout << counter << endl;
    }

}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    statnames.close();

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plotsTest/");

    // ######################
    //  Tables and other stuff
    // ######################


vector<string> totYield = { "SR1l_A_250lessMETless350" , "SR1l_A_250lessMETless350LSFdown" , "SR1l_A_250lessMETless350LSFup" , "SR1l_A_250lessMETless350BTlightDown" , "SR1l_A_250lessMETless350BTlightUp" , "SR1l_A_250lessMETless350BTheavyDown" , "SR1l_A_250lessMETless350BTheavyUp" , "SR1l_A_250lessMETless350PUdown" , "SR1l_A_250lessMETless350PUup" , "SR1l_A_250lessMETless350PDFdown" , "SR1l_A_250lessMETless350PDFup" , "SR1l_A_250lessMETless350alphaSdown" , "SR1l_A_250lessMETless350alphaSup" , "SR1l_A_250lessMETless350Q2down" , "SR1l_A_250lessMETless350Q2up" , "SR1l_A_250lessMETless350ISRnjetsDown" , "SR1l_A_250lessMETless350ISRnjetsUp" , "SR1l_A_350lessMETless450" , "SR1l_A_350lessMETless450LSFdown" , "SR1l_A_350lessMETless450LSFup" , "SR1l_A_350lessMETless450BTlightDown" , "SR1l_A_350lessMETless450BTlightUp" , "SR1l_A_350lessMETless450BTheavyDown" , "SR1l_A_350lessMETless450BTheavyUp" , "SR1l_A_350lessMETless450PUdown" , "SR1l_A_350lessMETless450PUup" , "SR1l_A_350lessMETless450PDFdown" , "SR1l_A_350lessMETless450PDFup" , "SR1l_A_350lessMETless450alphaSdown" , "SR1l_A_350lessMETless450alphaSup" , "SR1l_A_350lessMETless450Q2down" , "SR1l_A_350lessMETless450Q2up" , "SR1l_A_350lessMETless450ISRnjetsDown" , "SR1l_A_350lessMETless450ISRnjetsUp" , "SR1l_A_450lessMETless600" , "SR1l_A_450lessMETless600LSFdown" , "SR1l_A_450lessMETless600LSFup" , "SR1l_A_450lessMETless600BTlightDown" , "SR1l_A_450lessMETless600BTlightUp" , "SR1l_A_450lessMETless600BTheavyDown" , "SR1l_A_450lessMETless600BTheavyUp" , "SR1l_A_450lessMETless600PUdown" , "SR1l_A_450lessMETless600PUup" , "SR1l_A_450lessMETless600PDFdown" , "SR1l_A_450lessMETless600PDFup" , "SR1l_A_450lessMETless600alphaSdown" , "SR1l_A_450lessMETless600alphaSup" , "SR1l_A_450lessMETless600Q2down" , "SR1l_A_450lessMETless600Q2up" , "SR1l_A_450lessMETless600ISRnjetsDown" , "SR1l_A_450lessMETless600ISRnjetsUp" , "SR1l_A_600lessMETlessInf" , "SR1l_A_600lessMETlessInfLSFdown" , "SR1l_A_600lessMETlessInfLSFup" , "SR1l_A_600lessMETlessInfBTlightDown" , "SR1l_A_600lessMETlessInfBTlightUp" , "SR1l_A_600lessMETlessInfBTheavyDown" , "SR1l_A_600lessMETlessInfBTheavyUp" , "SR1l_A_600lessMETlessInfPUdown" , "SR1l_A_600lessMETlessInfPUup" , "SR1l_A_600lessMETlessInfPDFdown" , "SR1l_A_600lessMETlessInfPDFup" , "SR1l_A_600lessMETlessInfalphaSdown" , "SR1l_A_600lessMETlessInfalphaSup" , "SR1l_A_600lessMETlessInfQ2down" , "SR1l_A_600lessMETlessInfQ2up" , "SR1l_A_600lessMETlessInfISRnjetsDown" , "SR1l_A_600lessMETlessInfISRnjetsUp" , "SR1l_B_250lessMETless450" , "SR1l_B_250lessMETless450LSFdown" , "SR1l_B_250lessMETless450LSFup" , "SR1l_B_250lessMETless450BTlightDown" , "SR1l_B_250lessMETless450BTlightUp" , "SR1l_B_250lessMETless450BTheavyDown" , "SR1l_B_250lessMETless450BTheavyUp" , "SR1l_B_250lessMETless450PUdown" , "SR1l_B_250lessMETless450PUup" , "SR1l_B_250lessMETless450PDFdown" , "SR1l_B_250lessMETless450PDFup" , "SR1l_B_250lessMETless450alphaSdown" , "SR1l_B_250lessMETless450alphaSup" , "SR1l_B_250lessMETless450Q2down" , "SR1l_B_250lessMETless450Q2up" , "SR1l_B_250lessMETless450ISRnjetsDown" , "SR1l_B_250lessMETless450ISRnjetsUp" , "SR1l_B_450lessMETless600" , "SR1l_B_450lessMETless600LSFdown" , "SR1l_B_450lessMETless600LSFup" , "SR1l_B_450lessMETless600BTlightDown" , "SR1l_B_450lessMETless600BTlightUp" , "SR1l_B_450lessMETless600BTheavyDown" , "SR1l_B_450lessMETless600BTheavyUp" , "SR1l_B_450lessMETless600PUdown" , "SR1l_B_450lessMETless600PUup" , "SR1l_B_450lessMETless600PDFdown" , "SR1l_B_450lessMETless600PDFup" , "SR1l_B_450lessMETless600alphaSdown" , "SR1l_B_450lessMETless600alphaSup" , "SR1l_B_450lessMETless600Q2down" , "SR1l_B_450lessMETless600Q2up" , "SR1l_B_450lessMETless600ISRnjetsDown" , "SR1l_B_450lessMETless600ISRnjetsUp" , "SR1l_B_600lessMETlessInf" , "SR1l_B_600lessMETlessInfLSFdown" , "SR1l_B_600lessMETlessInfLSFup" , "SR1l_B_600lessMETlessInfBTlightDown" , "SR1l_B_600lessMETlessInfBTlightUp" , "SR1l_B_600lessMETlessInfBTheavyDown" , "SR1l_B_600lessMETlessInfBTheavyUp" , "SR1l_B_600lessMETlessInfPUdown" , "SR1l_B_600lessMETlessInfPUup" , "SR1l_B_600lessMETlessInfPDFdown" , "SR1l_B_600lessMETlessInfPDFup" , "SR1l_B_600lessMETlessInfalphaSdown" , "SR1l_B_600lessMETlessInfalphaSup" , "SR1l_B_600lessMETlessInfQ2down" , "SR1l_B_600lessMETlessInfQ2up" , "SR1l_B_600lessMETlessInfISRnjetsDown" , "SR1l_B_600lessMETlessInfISRnjetsUp" , "SR1l_C_250lessMETless350" , "SR1l_C_250lessMETless350LSFdown" , "SR1l_C_250lessMETless350LSFup" , "SR1l_C_250lessMETless350BTlightDown" , "SR1l_C_250lessMETless350BTlightUp" , "SR1l_C_250lessMETless350BTheavyDown" , "SR1l_C_250lessMETless350BTheavyUp" , "SR1l_C_250lessMETless350PUdown" , "SR1l_C_250lessMETless350PUup" , "SR1l_C_250lessMETless350PDFdown" , "SR1l_C_250lessMETless350PDFup" , "SR1l_C_250lessMETless350alphaSdown" , "SR1l_C_250lessMETless350alphaSup" , "SR1l_C_250lessMETless350Q2down" , "SR1l_C_250lessMETless350Q2up" , "SR1l_C_250lessMETless350ISRnjetsDown" , "SR1l_C_250lessMETless350ISRnjetsUp" , "SR1l_C_350lessMETless450" , "SR1l_C_350lessMETless450LSFdown" , "SR1l_C_350lessMETless450LSFup" , "SR1l_C_350lessMETless450BTlightDown" , "SR1l_C_350lessMETless450BTlightUp" , "SR1l_C_350lessMETless450BTheavyDown" , "SR1l_C_350lessMETless450BTheavyUp" , "SR1l_C_350lessMETless450PUdown" , "SR1l_C_350lessMETless450PUup" , "SR1l_C_350lessMETless450PDFdown" , "SR1l_C_350lessMETless450PDFup" , "SR1l_C_350lessMETless450alphaSdown" , "SR1l_C_350lessMETless450alphaSup" , "SR1l_C_350lessMETless450Q2down" , "SR1l_C_350lessMETless450Q2up" , "SR1l_C_350lessMETless450ISRnjetsDown" , "SR1l_C_350lessMETless450ISRnjetsUp" , "SR1l_C_450lessMETless550" , "SR1l_C_450lessMETless550LSFdown" , "SR1l_C_450lessMETless550LSFup" , "SR1l_C_450lessMETless550BTlightDown" , "SR1l_C_450lessMETless550BTlightUp" , "SR1l_C_450lessMETless550BTheavyDown" , "SR1l_C_450lessMETless550BTheavyUp" , "SR1l_C_450lessMETless550PUdown" , "SR1l_C_450lessMETless550PUup" , "SR1l_C_450lessMETless550PDFdown" , "SR1l_C_450lessMETless550PDFup" , "SR1l_C_450lessMETless550alphaSdown" , "SR1l_C_450lessMETless550alphaSup" , "SR1l_C_450lessMETless550Q2down" , "SR1l_C_450lessMETless550Q2up" , "SR1l_C_450lessMETless550ISRnjetsDown" , "SR1l_C_450lessMETless550ISRnjetsUp" , "SR1l_C_550lessMETless650" , "SR1l_C_550lessMETless650LSFdown" , "SR1l_C_550lessMETless650LSFup" , "SR1l_C_550lessMETless650BTlightDown" , "SR1l_C_550lessMETless650BTlightUp" , "SR1l_C_550lessMETless650BTheavyDown" , "SR1l_C_550lessMETless650BTheavyUp" , "SR1l_C_550lessMETless650PUdown" , "SR1l_C_550lessMETless650PUup" , "SR1l_C_550lessMETless650PDFdown" , "SR1l_C_550lessMETless650PDFup" , "SR1l_C_550lessMETless650alphaSdown" , "SR1l_C_550lessMETless650alphaSup" , "SR1l_C_550lessMETless650Q2down" , "SR1l_C_550lessMETless650Q2up" , "SR1l_C_550lessMETless650ISRnjetsDown" , "SR1l_C_550lessMETless650ISRnjetsUp" , "SR1l_C_650lessMETlessInf" , "SR1l_C_650lessMETlessInfLSFdown" , "SR1l_C_650lessMETlessInfLSFup" , "SR1l_C_650lessMETlessInfBTlightDown" , "SR1l_C_650lessMETlessInfBTlightUp" , "SR1l_C_650lessMETlessInfBTheavyDown" , "SR1l_C_650lessMETlessInfBTheavyUp" , "SR1l_C_650lessMETlessInfPUdown" , "SR1l_C_650lessMETlessInfPUup" , "SR1l_C_650lessMETlessInfPDFdown" , "SR1l_C_650lessMETlessInfPDFup" , "SR1l_C_650lessMETlessInfalphaSdown" , "SR1l_C_650lessMETlessInfalphaSup" , "SR1l_C_650lessMETlessInfQ2down" , "SR1l_C_650lessMETlessInfQ2up" , "SR1l_C_650lessMETlessInfISRnjetsDown" , "SR1l_C_650lessMETlessInfISRnjetsUp" , "SR1l_D_250lessMETless350" , "SR1l_D_250lessMETless350LSFdown" , "SR1l_D_250lessMETless350LSFup" , "SR1l_D_250lessMETless350BTlightDown" , "SR1l_D_250lessMETless350BTlightUp" , "SR1l_D_250lessMETless350BTheavyDown" , "SR1l_D_250lessMETless350BTheavyUp" , "SR1l_D_250lessMETless350PUdown" , "SR1l_D_250lessMETless350PUup" , "SR1l_D_250lessMETless350PDFdown" , "SR1l_D_250lessMETless350PDFup" , "SR1l_D_250lessMETless350alphaSdown" , "SR1l_D_250lessMETless350alphaSup" , "SR1l_D_250lessMETless350Q2down" , "SR1l_D_250lessMETless350Q2up" , "SR1l_D_250lessMETless350ISRnjetsDown" , "SR1l_D_250lessMETless350ISRnjetsUp" , "SR1l_D_350lessMETless450" , "SR1l_D_350lessMETless450LSFdown" , "SR1l_D_350lessMETless450LSFup" , "SR1l_D_350lessMETless450BTlightDown" , "SR1l_D_350lessMETless450BTlightUp" , "SR1l_D_350lessMETless450BTheavyDown" , "SR1l_D_350lessMETless450BTheavyUp" , "SR1l_D_350lessMETless450PUdown" , "SR1l_D_350lessMETless450PUup" , "SR1l_D_350lessMETless450PDFdown" , "SR1l_D_350lessMETless450PDFup" , "SR1l_D_350lessMETless450alphaSdown" , "SR1l_D_350lessMETless450alphaSup" , "SR1l_D_350lessMETless450Q2down" , "SR1l_D_350lessMETless450Q2up" , "SR1l_D_350lessMETless450ISRnjetsDown" , "SR1l_D_350lessMETless450ISRnjetsUp" , "SR1l_D_450lessMETless550" , "SR1l_D_450lessMETless550LSFdown" , "SR1l_D_450lessMETless550LSFup" , "SR1l_D_450lessMETless550BTlightDown" , "SR1l_D_450lessMETless550BTlightUp" , "SR1l_D_450lessMETless550BTheavyDown" , "SR1l_D_450lessMETless550BTheavyUp" , "SR1l_D_450lessMETless550PUdown" , "SR1l_D_450lessMETless550PUup" , "SR1l_D_450lessMETless550PDFdown" , "SR1l_D_450lessMETless550PDFup" , "SR1l_D_450lessMETless550alphaSdown" , "SR1l_D_450lessMETless550alphaSup" , "SR1l_D_450lessMETless550Q2down" , "SR1l_D_450lessMETless550Q2up" , "SR1l_D_450lessMETless550ISRnjetsDown" , "SR1l_D_450lessMETless550ISRnjetsUp" , "SR1l_D_550lessMETlessInf" , "SR1l_D_550lessMETlessInfLSFdown" , "SR1l_D_550lessMETlessInfLSFup" , "SR1l_D_550lessMETlessInfBTlightDown" , "SR1l_D_550lessMETlessInfBTlightUp" , "SR1l_D_550lessMETlessInfBTheavyDown" , "SR1l_D_550lessMETlessInfBTheavyUp" , "SR1l_D_550lessMETlessInfPUdown" , "SR1l_D_550lessMETlessInfPUup" , "SR1l_D_550lessMETlessInfPDFdown" , "SR1l_D_550lessMETlessInfPDFup" , "SR1l_D_550lessMETlessInfalphaSdown" , "SR1l_D_550lessMETlessInfalphaSup" , "SR1l_D_550lessMETlessInfQ2down" , "SR1l_D_550lessMETlessInfQ2up" , "SR1l_D_550lessMETlessInfISRnjetsDown" , "SR1l_D_550lessMETlessInfISRnjetsUp" , "SR1l_E_250lessMETless350" , "SR1l_E_250lessMETless350LSFdown" , "SR1l_E_250lessMETless350LSFup" , "SR1l_E_250lessMETless350BTlightDown" , "SR1l_E_250lessMETless350BTlightUp" , "SR1l_E_250lessMETless350BTheavyDown" , "SR1l_E_250lessMETless350BTheavyUp" , "SR1l_E_250lessMETless350PUdown" , "SR1l_E_250lessMETless350PUup" , "SR1l_E_250lessMETless350PDFdown" , "SR1l_E_250lessMETless350PDFup" , "SR1l_E_250lessMETless350alphaSdown" , "SR1l_E_250lessMETless350alphaSup" , "SR1l_E_250lessMETless350Q2down" , "SR1l_E_250lessMETless350Q2up" , "SR1l_E_250lessMETless350ISRnjetsDown" , "SR1l_E_250lessMETless350ISRnjetsUp" , "SR1l_E_350lessMETless550" , "SR1l_E_350lessMETless550LSFdown" , "SR1l_E_350lessMETless550LSFup" , "SR1l_E_350lessMETless550BTlightDown" , "SR1l_E_350lessMETless550BTlightUp" , "SR1l_E_350lessMETless550BTheavyDown" , "SR1l_E_350lessMETless550BTheavyUp" , "SR1l_E_350lessMETless550PUdown" , "SR1l_E_350lessMETless550PUup" , "SR1l_E_350lessMETless550PDFdown" , "SR1l_E_350lessMETless550PDFup" , "SR1l_E_350lessMETless550alphaSdown" , "SR1l_E_350lessMETless550alphaSup" , "SR1l_E_350lessMETless550Q2down" , "SR1l_E_350lessMETless550Q2up" , "SR1l_E_350lessMETless550ISRnjetsDown" , "SR1l_E_350lessMETless550ISRnjetsUp" , "SR1l_E_550lessMETlessInf" , "SR1l_E_550lessMETlessInfLSFdown" , "SR1l_E_550lessMETlessInfLSFup" , "SR1l_E_550lessMETlessInfBTlightDown" , "SR1l_E_550lessMETlessInfBTlightUp" , "SR1l_E_550lessMETlessInfBTheavyDown" , "SR1l_E_550lessMETlessInfBTheavyUp" , "SR1l_E_550lessMETlessInfPUdown" , "SR1l_E_550lessMETlessInfPUup" , "SR1l_E_550lessMETlessInfPDFdown" , "SR1l_E_550lessMETlessInfPDFup" , "SR1l_E_550lessMETlessInfalphaSdown" , "SR1l_E_550lessMETlessInfalphaSup" , "SR1l_E_550lessMETlessInfQ2down" , "SR1l_E_550lessMETlessInfQ2up" , "SR1l_E_550lessMETlessInfISRnjetsDown" , "SR1l_E_550lessMETlessInfISRnjetsUp" , "SR1l_F_250lessMETless450" , "SR1l_F_250lessMETless450LSFdown" , "SR1l_F_250lessMETless450LSFup" , "SR1l_F_250lessMETless450BTlightDown" , "SR1l_F_250lessMETless450BTlightUp" , "SR1l_F_250lessMETless450BTheavyDown" , "SR1l_F_250lessMETless450BTheavyUp" , "SR1l_F_250lessMETless450PUdown" , "SR1l_F_250lessMETless450PUup" , "SR1l_F_250lessMETless450PDFdown" , "SR1l_F_250lessMETless450PDFup" , "SR1l_F_250lessMETless450alphaSdown" , "SR1l_F_250lessMETless450alphaSup" , "SR1l_F_250lessMETless450Q2down" , "SR1l_F_250lessMETless450Q2up" , "SR1l_F_250lessMETless450ISRnjetsDown" , "SR1l_F_250lessMETless450ISRnjetsUp" , "SR1l_F_450lessMETlessInf" , "SR1l_F_450lessMETlessInfLSFdown" , "SR1l_F_450lessMETlessInfLSFup" , "SR1l_F_450lessMETlessInfBTlightDown" , "SR1l_F_450lessMETlessInfBTlightUp" , "SR1l_F_450lessMETlessInfBTheavyDown" , "SR1l_F_450lessMETlessInfBTheavyUp" , "SR1l_F_450lessMETlessInfPUdown" , "SR1l_F_450lessMETlessInfPUup" , "SR1l_F_450lessMETlessInfPDFdown" , "SR1l_F_450lessMETlessInfPDFup" , "SR1l_F_450lessMETlessInfalphaSdown" , "SR1l_F_450lessMETlessInfalphaSup" , "SR1l_F_450lessMETlessInfQ2down" , "SR1l_F_450lessMETlessInfQ2up" , "SR1l_F_450lessMETlessInfISRnjetsDown" , "SR1l_F_450lessMETlessInfISRnjetsUp" , "SR1l_G_250lessMETless350" , "SR1l_G_250lessMETless350LSFdown" , "SR1l_G_250lessMETless350LSFup" , "SR1l_G_250lessMETless350BTlightDown" , "SR1l_G_250lessMETless350BTlightUp" , "SR1l_G_250lessMETless350BTheavyDown" , "SR1l_G_250lessMETless350BTheavyUp" , "SR1l_G_250lessMETless350PUdown" , "SR1l_G_250lessMETless350PUup" , "SR1l_G_250lessMETless350PDFdown" , "SR1l_G_250lessMETless350PDFup" , "SR1l_G_250lessMETless350alphaSdown" , "SR1l_G_250lessMETless350alphaSup" , "SR1l_G_250lessMETless350Q2down" , "SR1l_G_250lessMETless350Q2up" , "SR1l_G_250lessMETless350ISRnjetsDown" , "SR1l_G_250lessMETless350ISRnjetsUp" , "SR1l_G_350lessMETless450" , "SR1l_G_350lessMETless450LSFdown" , "SR1l_G_350lessMETless450LSFup" , "SR1l_G_350lessMETless450BTlightDown" , "SR1l_G_350lessMETless450BTlightUp" , "SR1l_G_350lessMETless450BTheavyDown" , "SR1l_G_350lessMETless450BTheavyUp" , "SR1l_G_350lessMETless450PUdown" , "SR1l_G_350lessMETless450PUup" , "SR1l_G_350lessMETless450PDFdown" , "SR1l_G_350lessMETless450PDFup" , "SR1l_G_350lessMETless450alphaSdown" , "SR1l_G_350lessMETless450alphaSup" , "SR1l_G_350lessMETless450Q2down" , "SR1l_G_350lessMETless450Q2up" , "SR1l_G_350lessMETless450ISRnjetsDown" , "SR1l_G_350lessMETless450ISRnjetsUp" , "SR1l_G_450lessMETless600" , "SR1l_G_450lessMETless600LSFdown" , "SR1l_G_450lessMETless600LSFup" , "SR1l_G_450lessMETless600BTlightDown" , "SR1l_G_450lessMETless600BTlightUp" , "SR1l_G_450lessMETless600BTheavyDown" , "SR1l_G_450lessMETless600BTheavyUp" , "SR1l_G_450lessMETless600PUdown" , "SR1l_G_450lessMETless600PUup" , "SR1l_G_450lessMETless600PDFdown" , "SR1l_G_450lessMETless600PDFup" , "SR1l_G_450lessMETless600alphaSdown" , "SR1l_G_450lessMETless600alphaSup" , "SR1l_G_450lessMETless600Q2down" , "SR1l_G_450lessMETless600Q2up" , "SR1l_G_450lessMETless600ISRnjetsDown" , "SR1l_G_450lessMETless600ISRnjetsUp" , "SR1l_G_600lessMETlessInf" , "SR1l_G_600lessMETlessInfLSFdown" , "SR1l_G_600lessMETlessInfLSFup" , "SR1l_G_600lessMETlessInfBTlightDown" , "SR1l_G_600lessMETlessInfBTlightUp" , "SR1l_G_600lessMETlessInfBTheavyDown" , "SR1l_G_600lessMETlessInfBTheavyUp" , "SR1l_G_600lessMETlessInfPUdown" , "SR1l_G_600lessMETlessInfPUup" , "SR1l_G_600lessMETlessInfPDFdown" , "SR1l_G_600lessMETlessInfPDFup" , "SR1l_G_600lessMETlessInfalphaSdown" , "SR1l_G_600lessMETlessInfalphaSup" , "SR1l_G_600lessMETlessInfQ2down" , "SR1l_G_600lessMETlessInfQ2up" , "SR1l_G_600lessMETlessInfISRnjetsDown" , "SR1l_G_600lessMETlessInfISRnjetsUp" , "SR1l_H_250lessMETless450" , "SR1l_H_250lessMETless450LSFdown" , "SR1l_H_250lessMETless450LSFup" , "SR1l_H_250lessMETless450BTlightDown" , "SR1l_H_250lessMETless450BTlightUp" , "SR1l_H_250lessMETless450BTheavyDown" , "SR1l_H_250lessMETless450BTheavyUp" , "SR1l_H_250lessMETless450PUdown" , "SR1l_H_250lessMETless450PUup" , "SR1l_H_250lessMETless450PDFdown" , "SR1l_H_250lessMETless450PDFup" , "SR1l_H_250lessMETless450alphaSdown" , "SR1l_H_250lessMETless450alphaSup" , "SR1l_H_250lessMETless450Q2down" , "SR1l_H_250lessMETless450Q2up" , "SR1l_H_250lessMETless450ISRnjetsDown" , "SR1l_H_250lessMETless450ISRnjetsUp" , "SR1l_H_450lessMETlessInf" , "SR1l_H_450lessMETlessInfLSFdown" , "SR1l_H_450lessMETlessInfLSFup" , "SR1l_H_450lessMETlessInfBTlightDown" , "SR1l_H_450lessMETlessInfBTlightUp" , "SR1l_H_450lessMETlessInfBTheavyDown" , "SR1l_H_450lessMETlessInfBTheavyUp" , "SR1l_H_450lessMETlessInfPUdown" , "SR1l_H_450lessMETlessInfPUup" , "SR1l_H_450lessMETlessInfPDFdown" , "SR1l_H_450lessMETlessInfPDFup" , "SR1l_H_450lessMETlessInfalphaSdown" , "SR1l_H_450lessMETlessInfalphaSup" , "SR1l_H_450lessMETlessInfQ2down" , "SR1l_H_450lessMETlessInfQ2up" , "SR1l_H_450lessMETlessInfISRnjetsDown" , "SR1l_H_450lessMETlessInfISRnjetsUp" , "SR1l_I_250lessMETless350" , "SR1l_I_250lessMETless350LSFdown" , "SR1l_I_250lessMETless350LSFup" , "SR1l_I_250lessMETless350BTlightDown" , "SR1l_I_250lessMETless350BTlightUp" , "SR1l_I_250lessMETless350BTheavyDown" , "SR1l_I_250lessMETless350BTheavyUp" , "SR1l_I_250lessMETless350PUdown" , "SR1l_I_250lessMETless350PUup" , "SR1l_I_250lessMETless350PDFdown" , "SR1l_I_250lessMETless350PDFup" , "SR1l_I_250lessMETless350alphaSdown" , "SR1l_I_250lessMETless350alphaSup" , "SR1l_I_250lessMETless350Q2down" , "SR1l_I_250lessMETless350Q2up" , "SR1l_I_250lessMETless350ISRnjetsDown" , "SR1l_I_250lessMETless350ISRnjetsUp" , "SR1l_I_350lessMETless450" , "SR1l_I_350lessMETless450LSFdown" , "SR1l_I_350lessMETless450LSFup" , "SR1l_I_350lessMETless450BTlightDown" , "SR1l_I_350lessMETless450BTlightUp" , "SR1l_I_350lessMETless450BTheavyDown" , "SR1l_I_350lessMETless450BTheavyUp" , "SR1l_I_350lessMETless450PUdown" , "SR1l_I_350lessMETless450PUup" , "SR1l_I_350lessMETless450PDFdown" , "SR1l_I_350lessMETless450PDFup" , "SR1l_I_350lessMETless450alphaSdown" , "SR1l_I_350lessMETless450alphaSup" , "SR1l_I_350lessMETless450Q2down" , "SR1l_I_350lessMETless450Q2up" , "SR1l_I_350lessMETless450ISRnjetsDown" , "SR1l_I_350lessMETless450ISRnjetsUp" , "SR1l_I_450lessMETless550" , "SR1l_I_450lessMETless550LSFdown" , "SR1l_I_450lessMETless550LSFup" , "SR1l_I_450lessMETless550BTlightDown" , "SR1l_I_450lessMETless550BTlightUp" , "SR1l_I_450lessMETless550BTheavyDown" , "SR1l_I_450lessMETless550BTheavyUp" , "SR1l_I_450lessMETless550PUdown" , "SR1l_I_450lessMETless550PUup" , "SR1l_I_450lessMETless550PDFdown" , "SR1l_I_450lessMETless550PDFup" , "SR1l_I_450lessMETless550alphaSdown" , "SR1l_I_450lessMETless550alphaSup" , "SR1l_I_450lessMETless550Q2down" , "SR1l_I_450lessMETless550Q2up" , "SR1l_I_450lessMETless550ISRnjetsDown" , "SR1l_I_450lessMETless550ISRnjetsUp" , "SR1l_I_550lessMETlessInf" , "SR1l_I_550lessMETlessInfLSFdown" , "SR1l_I_550lessMETlessInfLSFup" , "SR1l_I_550lessMETlessInfBTlightDown" , "SR1l_I_550lessMETlessInfBTlightUp" , "SR1l_I_550lessMETlessInfBTheavyDown" , "SR1l_I_550lessMETlessInfBTheavyUp" , "SR1l_I_550lessMETlessInfPUdown" , "SR1l_I_550lessMETlessInfPUup" , "SR1l_I_550lessMETlessInfPDFdown" , "SR1l_I_550lessMETlessInfPDFup" , "SR1l_I_550lessMETlessInfalphaSdown" , "SR1l_I_550lessMETlessInfalphaSup" , "SR1l_I_550lessMETlessInfQ2down" , "SR1l_I_550lessMETlessInfQ2up" , "SR1l_I_550lessMETlessInfISRnjetsDown" , "SR1l_I_550lessMETlessInfISRnjetsUp"};

    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+".tab", 6);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+".tex", 6);

    ofstream sigfile("signalRegMor.txt");
    if (sigfile.is_open())
    {
        for(uint32_t r=0; r<totYield.size(); r++)
        {
            sigfile << totYield.at(r) << endl;
        }
            sigfile.close();
    }


    cout << "end of processing" << endl;
 }


    float getWeight(string currentProcessType, float lumi)
    {
        float nEvents =  myEvent.wNormalization.at(22);
        float all_weights = lumi*  myEvent.scale1fb* myEvent.weight_PU * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) );
        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
        return all_weights;
    }

    
    float reweightTop(float topPt, float atopPt)
    {
        return sqrt(exp(0.0615-(0.0005*topPt))*exp(0.0615-(0.0005*atopPt))); 
    }
