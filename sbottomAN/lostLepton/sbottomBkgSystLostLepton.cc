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

#include "../../common/TFFactory.h"
#include "../../Selection/moriondLostLep.h"

using namespace std;


// ----------------------------------------------
// Should be called only here because many
// struct and fuctions have to be declare first
// ----------------------------------------------
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"

uint32_t counter = 0;
string empty = "";
string storedDataset = "";
TH2D *h2 = NULL;
TAxis *xaxis = NULL;
TAxis *yaxis = NULL;
bool checkNegativeYields = false;
uint32_t nthentry = 0;
string outputName = "";
float scale1fbS2 =1;

float getWeight(string currentProcessType, float lumi, float s1fb2=1);
void getscale1fb2(TString fleName, float* scale1fb2 );
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

    babyTuplePath = "/opt/sbg/data/data6/cms/mjansova/Stop1lSharedBabies/v22/skim/";
    totalNumberOfWorkers = 1;

    // ------------------
    // Datasets
    // ------------------

    AddProcessClass("bkgLostLepton", "bkgLostLepton", "background", kBlue); 
    	AddDataset("W1JetsToLNu_madgraph_pythia8_25ns","bkgLostLepton",0,0); //@MJ@ TODO which W+jets - check Indara's slides?
    	AddDataset("W2JetsToLNu_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("W3JetsToLNu_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("W4JetsToLNu_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns","bkgLostLepton",0,0); //@MJ@ TODO what about extentions, good normalization to total number of events?
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns","bkgLostLepton",0,0); 
    	AddDataset("ttbar_diLept_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ttbar_diLept_madgraph_pythia8_ext1_25ns","bkgLostLepton",0,0);
    	AddDataset("t_sch_4f_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	//AddDataset("t_tW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
    	//AddDataset("t_tbarW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
    	//AddDataset("t_tch_4f_powheg_pythia8_inclDecays_25ns","bkgLostLepton",0,0); //@MJ@ TODO no atop?!
    	AddDataset("ttWJets_13TeV_madgraphMLM","bkgLostLepton",0,0);
    	AddDataset("WWTo2l2Nu_powheg_25ns","bkgLostLepton",0,0);
    	AddDataset("WWToLNuQQ_powheg_25ns","bkgLostLepton",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM","bkgZnunu",0,0);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("WZTo1LNu2Q_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("WZTo3LNu_powheg_pythia8_25ns","bkgLostLepton",0,0);//@MJ@ TODO 3l missing
    	AddDataset("ZZTo2L2Nu_powheg_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ZZTo2L2Q_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ZZTo4L_powheg_pythia8_25ns","bkgLostLepton",0,0);
    


        //@MJ@ TODO regions
        //both signal and control -use both *.h and two different tables for CRs and SRs
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
//CRs
AddRegion("CR2l_A_250lessMETless350","CR2l_A_250lessMETless350",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350LSFdown","CR2l_A_250lessMETless350LSFdown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350LSFup","CR2l_A_250lessMETless350LSFup",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350BTlightDown","CR2l_A_250lessMETless350BTlightDown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350BTlightUp","CR2l_A_250lessMETless350BTlightUp",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350BTheavyDown","CR2l_A_250lessMETless350BTheavyDown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350BTheavyUp","CR2l_A_250lessMETless350BTheavyUp",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350PUdown","CR2l_A_250lessMETless350PUdown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350PUup","CR2l_A_250lessMETless350PUup",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350PDFdown","CR2l_A_250lessMETless350PDFdown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350PDFup","CR2l_A_250lessMETless350PDFup",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350alphaSdown","CR2l_A_250lessMETless350alphaSdown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350alphaSup","CR2l_A_250lessMETless350alphaSup",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350Q2down","CR2l_A_250lessMETless350Q2down",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350Q2up","CR2l_A_250lessMETless350Q2up",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350ISRnjetsDown","CR2l_A_250lessMETless350ISRnjetsDown",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_250lessMETless350ISRnjetsUp","CR2l_A_250lessMETless350ISRnjetsUp",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_350lessMETless450","CR2l_A_350lessMETless450",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450LSFdown","CR2l_A_350lessMETless450LSFdown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450LSFup","CR2l_A_350lessMETless450LSFup",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450BTlightDown","CR2l_A_350lessMETless450BTlightDown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450BTlightUp","CR2l_A_350lessMETless450BTlightUp",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450BTheavyDown","CR2l_A_350lessMETless450BTheavyDown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450BTheavyUp","CR2l_A_350lessMETless450BTheavyUp",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450PUdown","CR2l_A_350lessMETless450PUdown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450PUup","CR2l_A_350lessMETless450PUup",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450PDFdown","CR2l_A_350lessMETless450PDFdown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450PDFup","CR2l_A_350lessMETless450PDFup",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450alphaSdown","CR2l_A_350lessMETless450alphaSdown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450alphaSup","CR2l_A_350lessMETless450alphaSup",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450Q2down","CR2l_A_350lessMETless450Q2down",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450Q2up","CR2l_A_350lessMETless450Q2up",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450ISRnjetsDown","CR2l_A_350lessMETless450ISRnjetsDown",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_350lessMETless450ISRnjetsUp","CR2l_A_350lessMETless450ISRnjetsUp",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_450lessMETless600","CR2l_A_450lessMETless600",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600LSFdown","CR2l_A_450lessMETless600LSFdown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600LSFup","CR2l_A_450lessMETless600LSFup",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600BTlightDown","CR2l_A_450lessMETless600BTlightDown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600BTlightUp","CR2l_A_450lessMETless600BTlightUp",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600BTheavyDown","CR2l_A_450lessMETless600BTheavyDown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600BTheavyUp","CR2l_A_450lessMETless600BTheavyUp",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600PUdown","CR2l_A_450lessMETless600PUdown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600PUup","CR2l_A_450lessMETless600PUup",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600PDFdown","CR2l_A_450lessMETless600PDFdown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600PDFup","CR2l_A_450lessMETless600PDFup",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600alphaSdown","CR2l_A_450lessMETless600alphaSdown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600alphaSup","CR2l_A_450lessMETless600alphaSup",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600Q2down","CR2l_A_450lessMETless600Q2down",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600Q2up","CR2l_A_450lessMETless600Q2up",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600ISRnjetsDown","CR2l_A_450lessMETless600ISRnjetsDown",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_450lessMETless600ISRnjetsUp","CR2l_A_450lessMETless600ISRnjetsUp",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_600lessMETlessInf","CR2l_A_600lessMETlessInf",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfLSFdown","CR2l_A_600lessMETlessInfLSFdown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfLSFup","CR2l_A_600lessMETlessInfLSFup",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfBTlightDown","CR2l_A_600lessMETlessInfBTlightDown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfBTlightUp","CR2l_A_600lessMETlessInfBTlightUp",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfBTheavyDown","CR2l_A_600lessMETlessInfBTheavyDown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfBTheavyUp","CR2l_A_600lessMETlessInfBTheavyUp",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfPUdown","CR2l_A_600lessMETlessInfPUdown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfPUup","CR2l_A_600lessMETlessInfPUup",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfPDFdown","CR2l_A_600lessMETlessInfPDFdown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfPDFup","CR2l_A_600lessMETlessInfPDFup",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfalphaSdown","CR2l_A_600lessMETlessInfalphaSdown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfalphaSup","CR2l_A_600lessMETlessInfalphaSup",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfQ2down","CR2l_A_600lessMETlessInfQ2down",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfQ2up","CR2l_A_600lessMETlessInfQ2up",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfISRnjetsDown","CR2l_A_600lessMETlessInfISRnjetsDown",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_A_600lessMETlessInfISRnjetsUp","CR2l_A_600lessMETlessInfISRnjetsUp",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_B_250lessMETless450","CR2l_B_250lessMETless450",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450LSFdown","CR2l_B_250lessMETless450LSFdown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450LSFup","CR2l_B_250lessMETless450LSFup",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450BTlightDown","CR2l_B_250lessMETless450BTlightDown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450BTlightUp","CR2l_B_250lessMETless450BTlightUp",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450BTheavyDown","CR2l_B_250lessMETless450BTheavyDown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450BTheavyUp","CR2l_B_250lessMETless450BTheavyUp",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450PUdown","CR2l_B_250lessMETless450PUdown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450PUup","CR2l_B_250lessMETless450PUup",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450PDFdown","CR2l_B_250lessMETless450PDFdown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450PDFup","CR2l_B_250lessMETless450PDFup",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450alphaSdown","CR2l_B_250lessMETless450alphaSdown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450alphaSup","CR2l_B_250lessMETless450alphaSup",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450Q2down","CR2l_B_250lessMETless450Q2down",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450Q2up","CR2l_B_250lessMETless450Q2up",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450ISRnjetsDown","CR2l_B_250lessMETless450ISRnjetsDown",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_250lessMETless450ISRnjetsUp","CR2l_B_250lessMETless450ISRnjetsUp",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_450lessMETless600","CR2l_B_450lessMETless600",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600LSFdown","CR2l_B_450lessMETless600LSFdown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600LSFup","CR2l_B_450lessMETless600LSFup",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600BTlightDown","CR2l_B_450lessMETless600BTlightDown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600BTlightUp","CR2l_B_450lessMETless600BTlightUp",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600BTheavyDown","CR2l_B_450lessMETless600BTheavyDown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600BTheavyUp","CR2l_B_450lessMETless600BTheavyUp",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600PUdown","CR2l_B_450lessMETless600PUdown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600PUup","CR2l_B_450lessMETless600PUup",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600PDFdown","CR2l_B_450lessMETless600PDFdown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600PDFup","CR2l_B_450lessMETless600PDFup",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600alphaSdown","CR2l_B_450lessMETless600alphaSdown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600alphaSup","CR2l_B_450lessMETless600alphaSup",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600Q2down","CR2l_B_450lessMETless600Q2down",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600Q2up","CR2l_B_450lessMETless600Q2up",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600ISRnjetsDown","CR2l_B_450lessMETless600ISRnjetsDown",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_450lessMETless600ISRnjetsUp","CR2l_B_450lessMETless600ISRnjetsUp",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_600lessMETlessInf","CR2l_B_600lessMETlessInf",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfLSFdown","CR2l_B_600lessMETlessInfLSFdown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfLSFup","CR2l_B_600lessMETlessInfLSFup",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfBTlightDown","CR2l_B_600lessMETlessInfBTlightDown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfBTlightUp","CR2l_B_600lessMETlessInfBTlightUp",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfBTheavyDown","CR2l_B_600lessMETlessInfBTheavyDown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfBTheavyUp","CR2l_B_600lessMETlessInfBTheavyUp",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfPUdown","CR2l_B_600lessMETlessInfPUdown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfPUup","CR2l_B_600lessMETlessInfPUup",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfPDFdown","CR2l_B_600lessMETlessInfPDFdown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfPDFup","CR2l_B_600lessMETlessInfPDFup",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfalphaSdown","CR2l_B_600lessMETlessInfalphaSdown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfalphaSup","CR2l_B_600lessMETlessInfalphaSup",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfQ2down","CR2l_B_600lessMETlessInfQ2down",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfQ2up","CR2l_B_600lessMETlessInfQ2up",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfISRnjetsDown","CR2l_B_600lessMETlessInfISRnjetsDown",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_B_600lessMETlessInfISRnjetsUp","CR2l_B_600lessMETlessInfISRnjetsUp",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_C_250lessMETless350","CR2l_C_250lessMETless350",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350LSFdown","CR2l_C_250lessMETless350LSFdown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350LSFup","CR2l_C_250lessMETless350LSFup",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350BTlightDown","CR2l_C_250lessMETless350BTlightDown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350BTlightUp","CR2l_C_250lessMETless350BTlightUp",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350BTheavyDown","CR2l_C_250lessMETless350BTheavyDown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350BTheavyUp","CR2l_C_250lessMETless350BTheavyUp",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350PUdown","CR2l_C_250lessMETless350PUdown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350PUup","CR2l_C_250lessMETless350PUup",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350PDFdown","CR2l_C_250lessMETless350PDFdown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350PDFup","CR2l_C_250lessMETless350PDFup",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350alphaSdown","CR2l_C_250lessMETless350alphaSdown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350alphaSup","CR2l_C_250lessMETless350alphaSup",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350Q2down","CR2l_C_250lessMETless350Q2down",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350Q2up","CR2l_C_250lessMETless350Q2up",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350ISRnjetsDown","CR2l_C_250lessMETless350ISRnjetsDown",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_250lessMETless350ISRnjetsUp","CR2l_C_250lessMETless350ISRnjetsUp",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_350lessMETless450","CR2l_C_350lessMETless450",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450LSFdown","CR2l_C_350lessMETless450LSFdown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450LSFup","CR2l_C_350lessMETless450LSFup",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450BTlightDown","CR2l_C_350lessMETless450BTlightDown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450BTlightUp","CR2l_C_350lessMETless450BTlightUp",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450BTheavyDown","CR2l_C_350lessMETless450BTheavyDown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450BTheavyUp","CR2l_C_350lessMETless450BTheavyUp",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450PUdown","CR2l_C_350lessMETless450PUdown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450PUup","CR2l_C_350lessMETless450PUup",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450PDFdown","CR2l_C_350lessMETless450PDFdown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450PDFup","CR2l_C_350lessMETless450PDFup",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450alphaSdown","CR2l_C_350lessMETless450alphaSdown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450alphaSup","CR2l_C_350lessMETless450alphaSup",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450Q2down","CR2l_C_350lessMETless450Q2down",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450Q2up","CR2l_C_350lessMETless450Q2up",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450ISRnjetsDown","CR2l_C_350lessMETless450ISRnjetsDown",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_350lessMETless450ISRnjetsUp","CR2l_C_350lessMETless450ISRnjetsUp",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_450lessMETless550","CR2l_C_450lessMETless550",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550LSFdown","CR2l_C_450lessMETless550LSFdown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550LSFup","CR2l_C_450lessMETless550LSFup",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550BTlightDown","CR2l_C_450lessMETless550BTlightDown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550BTlightUp","CR2l_C_450lessMETless550BTlightUp",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550BTheavyDown","CR2l_C_450lessMETless550BTheavyDown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550BTheavyUp","CR2l_C_450lessMETless550BTheavyUp",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550PUdown","CR2l_C_450lessMETless550PUdown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550PUup","CR2l_C_450lessMETless550PUup",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550PDFdown","CR2l_C_450lessMETless550PDFdown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550PDFup","CR2l_C_450lessMETless550PDFup",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550alphaSdown","CR2l_C_450lessMETless550alphaSdown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550alphaSup","CR2l_C_450lessMETless550alphaSup",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550Q2down","CR2l_C_450lessMETless550Q2down",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550Q2up","CR2l_C_450lessMETless550Q2up",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550ISRnjetsDown","CR2l_C_450lessMETless550ISRnjetsDown",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_450lessMETless550ISRnjetsUp","CR2l_C_450lessMETless550ISRnjetsUp",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_550lessMETless650","CR2l_C_550lessMETless650",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650LSFdown","CR2l_C_550lessMETless650LSFdown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650LSFup","CR2l_C_550lessMETless650LSFup",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650BTlightDown","CR2l_C_550lessMETless650BTlightDown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650BTlightUp","CR2l_C_550lessMETless650BTlightUp",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650BTheavyDown","CR2l_C_550lessMETless650BTheavyDown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650BTheavyUp","CR2l_C_550lessMETless650BTheavyUp",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650PUdown","CR2l_C_550lessMETless650PUdown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650PUup","CR2l_C_550lessMETless650PUup",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650PDFdown","CR2l_C_550lessMETless650PDFdown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650PDFup","CR2l_C_550lessMETless650PDFup",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650alphaSdown","CR2l_C_550lessMETless650alphaSdown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650alphaSup","CR2l_C_550lessMETless650alphaSup",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650Q2down","CR2l_C_550lessMETless650Q2down",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650Q2up","CR2l_C_550lessMETless650Q2up",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650ISRnjetsDown","CR2l_C_550lessMETless650ISRnjetsDown",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_550lessMETless650ISRnjetsUp","CR2l_C_550lessMETless650ISRnjetsUp",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_650lessMETlessInf","CR2l_C_650lessMETlessInf",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfLSFdown","CR2l_C_650lessMETlessInfLSFdown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfLSFup","CR2l_C_650lessMETlessInfLSFup",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfBTlightDown","CR2l_C_650lessMETlessInfBTlightDown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfBTlightUp","CR2l_C_650lessMETlessInfBTlightUp",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfBTheavyDown","CR2l_C_650lessMETlessInfBTheavyDown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfBTheavyUp","CR2l_C_650lessMETlessInfBTheavyUp",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfPUdown","CR2l_C_650lessMETlessInfPUdown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfPUup","CR2l_C_650lessMETlessInfPUup",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfPDFdown","CR2l_C_650lessMETlessInfPDFdown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfPDFup","CR2l_C_650lessMETlessInfPDFup",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfalphaSdown","CR2l_C_650lessMETlessInfalphaSdown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfalphaSup","CR2l_C_650lessMETlessInfalphaSup",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfQ2down","CR2l_C_650lessMETlessInfQ2down",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfQ2up","CR2l_C_650lessMETlessInfQ2up",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfISRnjetsDown","CR2l_C_650lessMETlessInfISRnjetsDown",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_C_650lessMETlessInfISRnjetsUp","CR2l_C_650lessMETlessInfISRnjetsUp",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_D_250lessMETless350","CR2l_D_250lessMETless350",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350LSFdown","CR2l_D_250lessMETless350LSFdown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350LSFup","CR2l_D_250lessMETless350LSFup",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350BTlightDown","CR2l_D_250lessMETless350BTlightDown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350BTlightUp","CR2l_D_250lessMETless350BTlightUp",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350BTheavyDown","CR2l_D_250lessMETless350BTheavyDown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350BTheavyUp","CR2l_D_250lessMETless350BTheavyUp",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350PUdown","CR2l_D_250lessMETless350PUdown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350PUup","CR2l_D_250lessMETless350PUup",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350PDFdown","CR2l_D_250lessMETless350PDFdown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350PDFup","CR2l_D_250lessMETless350PDFup",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350alphaSdown","CR2l_D_250lessMETless350alphaSdown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350alphaSup","CR2l_D_250lessMETless350alphaSup",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350Q2down","CR2l_D_250lessMETless350Q2down",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350Q2up","CR2l_D_250lessMETless350Q2up",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350ISRnjetsDown","CR2l_D_250lessMETless350ISRnjetsDown",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_250lessMETless350ISRnjetsUp","CR2l_D_250lessMETless350ISRnjetsUp",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_350lessMETless450","CR2l_D_350lessMETless450",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450LSFdown","CR2l_D_350lessMETless450LSFdown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450LSFup","CR2l_D_350lessMETless450LSFup",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450BTlightDown","CR2l_D_350lessMETless450BTlightDown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450BTlightUp","CR2l_D_350lessMETless450BTlightUp",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450BTheavyDown","CR2l_D_350lessMETless450BTheavyDown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450BTheavyUp","CR2l_D_350lessMETless450BTheavyUp",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450PUdown","CR2l_D_350lessMETless450PUdown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450PUup","CR2l_D_350lessMETless450PUup",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450PDFdown","CR2l_D_350lessMETless450PDFdown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450PDFup","CR2l_D_350lessMETless450PDFup",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450alphaSdown","CR2l_D_350lessMETless450alphaSdown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450alphaSup","CR2l_D_350lessMETless450alphaSup",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450Q2down","CR2l_D_350lessMETless450Q2down",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450Q2up","CR2l_D_350lessMETless450Q2up",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450ISRnjetsDown","CR2l_D_350lessMETless450ISRnjetsDown",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_350lessMETless450ISRnjetsUp","CR2l_D_350lessMETless450ISRnjetsUp",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_450lessMETless550","CR2l_D_450lessMETless550",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550LSFdown","CR2l_D_450lessMETless550LSFdown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550LSFup","CR2l_D_450lessMETless550LSFup",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550BTlightDown","CR2l_D_450lessMETless550BTlightDown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550BTlightUp","CR2l_D_450lessMETless550BTlightUp",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550BTheavyDown","CR2l_D_450lessMETless550BTheavyDown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550BTheavyUp","CR2l_D_450lessMETless550BTheavyUp",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550PUdown","CR2l_D_450lessMETless550PUdown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550PUup","CR2l_D_450lessMETless550PUup",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550PDFdown","CR2l_D_450lessMETless550PDFdown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550PDFup","CR2l_D_450lessMETless550PDFup",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550alphaSdown","CR2l_D_450lessMETless550alphaSdown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550alphaSup","CR2l_D_450lessMETless550alphaSup",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550Q2down","CR2l_D_450lessMETless550Q2down",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550Q2up","CR2l_D_450lessMETless550Q2up",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550ISRnjetsDown","CR2l_D_450lessMETless550ISRnjetsDown",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_450lessMETless550ISRnjetsUp","CR2l_D_450lessMETless550ISRnjetsUp",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_550lessMETlessInf","CR2l_D_550lessMETlessInf",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfLSFdown","CR2l_D_550lessMETlessInfLSFdown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfLSFup","CR2l_D_550lessMETlessInfLSFup",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfBTlightDown","CR2l_D_550lessMETlessInfBTlightDown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfBTlightUp","CR2l_D_550lessMETlessInfBTlightUp",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfBTheavyDown","CR2l_D_550lessMETlessInfBTheavyDown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfBTheavyUp","CR2l_D_550lessMETlessInfBTheavyUp",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfPUdown","CR2l_D_550lessMETlessInfPUdown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfPUup","CR2l_D_550lessMETlessInfPUup",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfPDFdown","CR2l_D_550lessMETlessInfPDFdown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfPDFup","CR2l_D_550lessMETlessInfPDFup",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfalphaSdown","CR2l_D_550lessMETlessInfalphaSdown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfalphaSup","CR2l_D_550lessMETlessInfalphaSup",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfQ2down","CR2l_D_550lessMETlessInfQ2down",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfQ2up","CR2l_D_550lessMETlessInfQ2up",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfISRnjetsDown","CR2l_D_550lessMETlessInfISRnjetsDown",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_D_550lessMETlessInfISRnjetsUp","CR2l_D_550lessMETlessInfISRnjetsUp",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_E_250lessMETless350","CR2l_E_250lessMETless350",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350LSFdown","CR2l_E_250lessMETless350LSFdown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350LSFup","CR2l_E_250lessMETless350LSFup",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350BTlightDown","CR2l_E_250lessMETless350BTlightDown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350BTlightUp","CR2l_E_250lessMETless350BTlightUp",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350BTheavyDown","CR2l_E_250lessMETless350BTheavyDown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350BTheavyUp","CR2l_E_250lessMETless350BTheavyUp",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350PUdown","CR2l_E_250lessMETless350PUdown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350PUup","CR2l_E_250lessMETless350PUup",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350PDFdown","CR2l_E_250lessMETless350PDFdown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350PDFup","CR2l_E_250lessMETless350PDFup",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350alphaSdown","CR2l_E_250lessMETless350alphaSdown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350alphaSup","CR2l_E_250lessMETless350alphaSup",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350Q2down","CR2l_E_250lessMETless350Q2down",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350Q2up","CR2l_E_250lessMETless350Q2up",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350ISRnjetsDown","CR2l_E_250lessMETless350ISRnjetsDown",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_250lessMETless350ISRnjetsUp","CR2l_E_250lessMETless350ISRnjetsUp",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_350lessMETless550","CR2l_E_350lessMETless550",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550LSFdown","CR2l_E_350lessMETless550LSFdown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550LSFup","CR2l_E_350lessMETless550LSFup",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550BTlightDown","CR2l_E_350lessMETless550BTlightDown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550BTlightUp","CR2l_E_350lessMETless550BTlightUp",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550BTheavyDown","CR2l_E_350lessMETless550BTheavyDown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550BTheavyUp","CR2l_E_350lessMETless550BTheavyUp",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550PUdown","CR2l_E_350lessMETless550PUdown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550PUup","CR2l_E_350lessMETless550PUup",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550PDFdown","CR2l_E_350lessMETless550PDFdown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550PDFup","CR2l_E_350lessMETless550PDFup",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550alphaSdown","CR2l_E_350lessMETless550alphaSdown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550alphaSup","CR2l_E_350lessMETless550alphaSup",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550Q2down","CR2l_E_350lessMETless550Q2down",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550Q2up","CR2l_E_350lessMETless550Q2up",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550ISRnjetsDown","CR2l_E_350lessMETless550ISRnjetsDown",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_350lessMETless550ISRnjetsUp","CR2l_E_350lessMETless550ISRnjetsUp",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_550lessMETlessInf","CR2l_E_550lessMETlessInf",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfLSFdown","CR2l_E_550lessMETlessInfLSFdown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfLSFup","CR2l_E_550lessMETlessInfLSFup",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfBTlightDown","CR2l_E_550lessMETlessInfBTlightDown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfBTlightUp","CR2l_E_550lessMETlessInfBTlightUp",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfBTheavyDown","CR2l_E_550lessMETlessInfBTheavyDown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfBTheavyUp","CR2l_E_550lessMETlessInfBTheavyUp",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfPUdown","CR2l_E_550lessMETlessInfPUdown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfPUup","CR2l_E_550lessMETlessInfPUup",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfPDFdown","CR2l_E_550lessMETlessInfPDFdown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfPDFup","CR2l_E_550lessMETlessInfPDFup",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfalphaSdown","CR2l_E_550lessMETlessInfalphaSdown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfalphaSup","CR2l_E_550lessMETlessInfalphaSup",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfQ2down","CR2l_E_550lessMETlessInfQ2down",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfQ2up","CR2l_E_550lessMETlessInfQ2up",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfISRnjetsDown","CR2l_E_550lessMETlessInfISRnjetsDown",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_E_550lessMETlessInfISRnjetsUp","CR2l_E_550lessMETlessInfISRnjetsUp",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_F_250lessMETless450","CR2l_F_250lessMETless450",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450LSFdown","CR2l_F_250lessMETless450LSFdown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450LSFup","CR2l_F_250lessMETless450LSFup",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450BTlightDown","CR2l_F_250lessMETless450BTlightDown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450BTlightUp","CR2l_F_250lessMETless450BTlightUp",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450BTheavyDown","CR2l_F_250lessMETless450BTheavyDown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450BTheavyUp","CR2l_F_250lessMETless450BTheavyUp",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450PUdown","CR2l_F_250lessMETless450PUdown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450PUup","CR2l_F_250lessMETless450PUup",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450PDFdown","CR2l_F_250lessMETless450PDFdown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450PDFup","CR2l_F_250lessMETless450PDFup",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450alphaSdown","CR2l_F_250lessMETless450alphaSdown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450alphaSup","CR2l_F_250lessMETless450alphaSup",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450Q2down","CR2l_F_250lessMETless450Q2down",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450Q2up","CR2l_F_250lessMETless450Q2up",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450ISRnjetsDown","CR2l_F_250lessMETless450ISRnjetsDown",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_250lessMETless450ISRnjetsUp","CR2l_F_250lessMETless450ISRnjetsUp",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_450lessMETlessInf","CR2l_F_450lessMETlessInf",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfLSFdown","CR2l_F_450lessMETlessInfLSFdown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfLSFup","CR2l_F_450lessMETlessInfLSFup",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfBTlightDown","CR2l_F_450lessMETlessInfBTlightDown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfBTlightUp","CR2l_F_450lessMETlessInfBTlightUp",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfBTheavyDown","CR2l_F_450lessMETlessInfBTheavyDown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfBTheavyUp","CR2l_F_450lessMETlessInfBTheavyUp",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfPUdown","CR2l_F_450lessMETlessInfPUdown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfPUup","CR2l_F_450lessMETlessInfPUup",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfPDFdown","CR2l_F_450lessMETlessInfPDFdown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfPDFup","CR2l_F_450lessMETlessInfPDFup",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfalphaSdown","CR2l_F_450lessMETlessInfalphaSdown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfalphaSup","CR2l_F_450lessMETlessInfalphaSup",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfQ2down","CR2l_F_450lessMETlessInfQ2down",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfQ2up","CR2l_F_450lessMETlessInfQ2up",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfISRnjetsDown","CR2l_F_450lessMETlessInfISRnjetsDown",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_F_450lessMETlessInfISRnjetsUp","CR2l_F_450lessMETlessInfISRnjetsUp",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_G_250lessMETless350","CR2l_G_250lessMETless350",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350LSFdown","CR2l_G_250lessMETless350LSFdown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350LSFup","CR2l_G_250lessMETless350LSFup",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350BTlightDown","CR2l_G_250lessMETless350BTlightDown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350BTlightUp","CR2l_G_250lessMETless350BTlightUp",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350BTheavyDown","CR2l_G_250lessMETless350BTheavyDown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350BTheavyUp","CR2l_G_250lessMETless350BTheavyUp",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350PUdown","CR2l_G_250lessMETless350PUdown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350PUup","CR2l_G_250lessMETless350PUup",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350PDFdown","CR2l_G_250lessMETless350PDFdown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350PDFup","CR2l_G_250lessMETless350PDFup",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350alphaSdown","CR2l_G_250lessMETless350alphaSdown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350alphaSup","CR2l_G_250lessMETless350alphaSup",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350Q2down","CR2l_G_250lessMETless350Q2down",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350Q2up","CR2l_G_250lessMETless350Q2up",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350ISRnjetsDown","CR2l_G_250lessMETless350ISRnjetsDown",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_250lessMETless350ISRnjetsUp","CR2l_G_250lessMETless350ISRnjetsUp",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_350lessMETless450","CR2l_G_350lessMETless450",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450LSFdown","CR2l_G_350lessMETless450LSFdown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450LSFup","CR2l_G_350lessMETless450LSFup",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450BTlightDown","CR2l_G_350lessMETless450BTlightDown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450BTlightUp","CR2l_G_350lessMETless450BTlightUp",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450BTheavyDown","CR2l_G_350lessMETless450BTheavyDown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450BTheavyUp","CR2l_G_350lessMETless450BTheavyUp",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450PUdown","CR2l_G_350lessMETless450PUdown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450PUup","CR2l_G_350lessMETless450PUup",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450PDFdown","CR2l_G_350lessMETless450PDFdown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450PDFup","CR2l_G_350lessMETless450PDFup",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450alphaSdown","CR2l_G_350lessMETless450alphaSdown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450alphaSup","CR2l_G_350lessMETless450alphaSup",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450Q2down","CR2l_G_350lessMETless450Q2down",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450Q2up","CR2l_G_350lessMETless450Q2up",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450ISRnjetsDown","CR2l_G_350lessMETless450ISRnjetsDown",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_350lessMETless450ISRnjetsUp","CR2l_G_350lessMETless450ISRnjetsUp",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_450lessMETless600","CR2l_G_450lessMETless600",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600LSFdown","CR2l_G_450lessMETless600LSFdown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600LSFup","CR2l_G_450lessMETless600LSFup",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600BTlightDown","CR2l_G_450lessMETless600BTlightDown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600BTlightUp","CR2l_G_450lessMETless600BTlightUp",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600BTheavyDown","CR2l_G_450lessMETless600BTheavyDown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600BTheavyUp","CR2l_G_450lessMETless600BTheavyUp",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600PUdown","CR2l_G_450lessMETless600PUdown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600PUup","CR2l_G_450lessMETless600PUup",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600PDFdown","CR2l_G_450lessMETless600PDFdown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600PDFup","CR2l_G_450lessMETless600PDFup",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600alphaSdown","CR2l_G_450lessMETless600alphaSdown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600alphaSup","CR2l_G_450lessMETless600alphaSup",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600Q2down","CR2l_G_450lessMETless600Q2down",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600Q2up","CR2l_G_450lessMETless600Q2up",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600ISRnjetsDown","CR2l_G_450lessMETless600ISRnjetsDown",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_450lessMETless600ISRnjetsUp","CR2l_G_450lessMETless600ISRnjetsUp",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_600lessMETlessInf","CR2l_G_600lessMETlessInf",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfLSFdown","CR2l_G_600lessMETlessInfLSFdown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfLSFup","CR2l_G_600lessMETlessInfLSFup",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfBTlightDown","CR2l_G_600lessMETlessInfBTlightDown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfBTlightUp","CR2l_G_600lessMETlessInfBTlightUp",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfBTheavyDown","CR2l_G_600lessMETlessInfBTheavyDown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfBTheavyUp","CR2l_G_600lessMETlessInfBTheavyUp",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfPUdown","CR2l_G_600lessMETlessInfPUdown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfPUup","CR2l_G_600lessMETlessInfPUup",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfPDFdown","CR2l_G_600lessMETlessInfPDFdown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfPDFup","CR2l_G_600lessMETlessInfPDFup",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfalphaSdown","CR2l_G_600lessMETlessInfalphaSdown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfalphaSup","CR2l_G_600lessMETlessInfalphaSup",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfQ2down","CR2l_G_600lessMETlessInfQ2down",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfQ2up","CR2l_G_600lessMETlessInfQ2up",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfISRnjetsDown","CR2l_G_600lessMETlessInfISRnjetsDown",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_G_600lessMETlessInfISRnjetsUp","CR2l_G_600lessMETlessInfISRnjetsUp",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_H_250lessMETless450","CR2l_H_250lessMETless450",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450LSFdown","CR2l_H_250lessMETless450LSFdown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450LSFup","CR2l_H_250lessMETless450LSFup",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450BTlightDown","CR2l_H_250lessMETless450BTlightDown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450BTlightUp","CR2l_H_250lessMETless450BTlightUp",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450BTheavyDown","CR2l_H_250lessMETless450BTheavyDown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450BTheavyUp","CR2l_H_250lessMETless450BTheavyUp",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450PUdown","CR2l_H_250lessMETless450PUdown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450PUup","CR2l_H_250lessMETless450PUup",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450PDFdown","CR2l_H_250lessMETless450PDFdown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450PDFup","CR2l_H_250lessMETless450PDFup",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450alphaSdown","CR2l_H_250lessMETless450alphaSdown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450alphaSup","CR2l_H_250lessMETless450alphaSup",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450Q2down","CR2l_H_250lessMETless450Q2down",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450Q2up","CR2l_H_250lessMETless450Q2up",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450ISRnjetsDown","CR2l_H_250lessMETless450ISRnjetsDown",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_250lessMETless450ISRnjetsUp","CR2l_H_250lessMETless450ISRnjetsUp",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_450lessMETlessInf","CR2l_H_450lessMETlessInf",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfLSFdown","CR2l_H_450lessMETlessInfLSFdown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfLSFup","CR2l_H_450lessMETlessInfLSFup",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfBTlightDown","CR2l_H_450lessMETlessInfBTlightDown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfBTlightUp","CR2l_H_450lessMETlessInfBTlightUp",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfBTheavyDown","CR2l_H_450lessMETlessInfBTheavyDown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfBTheavyUp","CR2l_H_450lessMETlessInfBTheavyUp",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfPUdown","CR2l_H_450lessMETlessInfPUdown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfPUup","CR2l_H_450lessMETlessInfPUup",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfPDFdown","CR2l_H_450lessMETlessInfPDFdown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfPDFup","CR2l_H_450lessMETlessInfPDFup",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfalphaSdown","CR2l_H_450lessMETlessInfalphaSdown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfalphaSup","CR2l_H_450lessMETlessInfalphaSup",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfQ2down","CR2l_H_450lessMETlessInfQ2down",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfQ2up","CR2l_H_450lessMETlessInfQ2up",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfISRnjetsDown","CR2l_H_450lessMETlessInfISRnjetsDown",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_H_450lessMETlessInfISRnjetsUp","CR2l_H_450lessMETlessInfISRnjetsUp",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_I_250lessMETless350","CR2l_I_250lessMETless350",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350LSFdown","CR2l_I_250lessMETless350LSFdown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350LSFup","CR2l_I_250lessMETless350LSFup",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350BTlightDown","CR2l_I_250lessMETless350BTlightDown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350BTlightUp","CR2l_I_250lessMETless350BTlightUp",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350BTheavyDown","CR2l_I_250lessMETless350BTheavyDown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350BTheavyUp","CR2l_I_250lessMETless350BTheavyUp",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350PUdown","CR2l_I_250lessMETless350PUdown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350PUup","CR2l_I_250lessMETless350PUup",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350PDFdown","CR2l_I_250lessMETless350PDFdown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350PDFup","CR2l_I_250lessMETless350PDFup",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350alphaSdown","CR2l_I_250lessMETless350alphaSdown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350alphaSup","CR2l_I_250lessMETless350alphaSup",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350Q2down","CR2l_I_250lessMETless350Q2down",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350Q2up","CR2l_I_250lessMETless350Q2up",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350ISRnjetsDown","CR2l_I_250lessMETless350ISRnjetsDown",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_250lessMETless350ISRnjetsUp","CR2l_I_250lessMETless350ISRnjetsUp",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_350lessMETless450","CR2l_I_350lessMETless450",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450LSFdown","CR2l_I_350lessMETless450LSFdown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450LSFup","CR2l_I_350lessMETless450LSFup",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450BTlightDown","CR2l_I_350lessMETless450BTlightDown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450BTlightUp","CR2l_I_350lessMETless450BTlightUp",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450BTheavyDown","CR2l_I_350lessMETless450BTheavyDown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450BTheavyUp","CR2l_I_350lessMETless450BTheavyUp",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450PUdown","CR2l_I_350lessMETless450PUdown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450PUup","CR2l_I_350lessMETless450PUup",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450PDFdown","CR2l_I_350lessMETless450PDFdown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450PDFup","CR2l_I_350lessMETless450PDFup",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450alphaSdown","CR2l_I_350lessMETless450alphaSdown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450alphaSup","CR2l_I_350lessMETless450alphaSup",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450Q2down","CR2l_I_350lessMETless450Q2down",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450Q2up","CR2l_I_350lessMETless450Q2up",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450ISRnjetsDown","CR2l_I_350lessMETless450ISRnjetsDown",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_350lessMETless450ISRnjetsUp","CR2l_I_350lessMETless450ISRnjetsUp",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_450lessMETless550","CR2l_I_450lessMETless550",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550LSFdown","CR2l_I_450lessMETless550LSFdown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550LSFup","CR2l_I_450lessMETless550LSFup",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550BTlightDown","CR2l_I_450lessMETless550BTlightDown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550BTlightUp","CR2l_I_450lessMETless550BTlightUp",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550BTheavyDown","CR2l_I_450lessMETless550BTheavyDown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550BTheavyUp","CR2l_I_450lessMETless550BTheavyUp",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550PUdown","CR2l_I_450lessMETless550PUdown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550PUup","CR2l_I_450lessMETless550PUup",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550PDFdown","CR2l_I_450lessMETless550PDFdown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550PDFup","CR2l_I_450lessMETless550PDFup",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550alphaSdown","CR2l_I_450lessMETless550alphaSdown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550alphaSup","CR2l_I_450lessMETless550alphaSup",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550Q2down","CR2l_I_450lessMETless550Q2down",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550Q2up","CR2l_I_450lessMETless550Q2up",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550ISRnjetsDown","CR2l_I_450lessMETless550ISRnjetsDown",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_450lessMETless550ISRnjetsUp","CR2l_I_450lessMETless550ISRnjetsUp",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_550lessMETlessInf","CR2l_I_550lessMETlessInf",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfLSFdown","CR2l_I_550lessMETlessInfLSFdown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfLSFup","CR2l_I_550lessMETlessInfLSFup",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfBTlightDown","CR2l_I_550lessMETlessInfBTlightDown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfBTlightUp","CR2l_I_550lessMETlessInfBTlightUp",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfBTheavyDown","CR2l_I_550lessMETlessInfBTheavyDown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfBTheavyUp","CR2l_I_550lessMETlessInfBTheavyUp",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfPUdown","CR2l_I_550lessMETlessInfPUdown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfPUup","CR2l_I_550lessMETlessInfPUup",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfPDFdown","CR2l_I_550lessMETlessInfPDFdown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfPDFup","CR2l_I_550lessMETlessInfPDFup",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfalphaSdown","CR2l_I_550lessMETlessInfalphaSdown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfalphaSup","CR2l_I_550lessMETlessInfalphaSup",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfQ2down","CR2l_I_550lessMETlessInfQ2down",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfQ2up","CR2l_I_550lessMETlessInfQ2up",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfISRnjetsDown","CR2l_I_550lessMETlessInfISRnjetsDown",&CR2l_I_550lessMETlessInf);
AddRegion("CR2l_I_550lessMETlessInfISRnjetsUp","CR2l_I_550lessMETlessInfISRnjetsUp",&CR2l_I_550lessMETlessInf);         
    fillYieldsVector();
                                                                                                                       

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(35.867);

    Create1DHistos();

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
        cout << "checking histogram for negative values" << endl;
        nthentry =0;
        checkNegativeYields = true; //@MJ@ TODO be aware of this, I can use multiprocessing now but I can endup with incorrect MC yields!!!
    }

    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);
    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset); //@MJ@ TODO check this!!!

    if( (currentProcessClass == "bkgLostLepton")  && ( !(myEvent.is2lep) ))
    {
         currentProcessClass = "";
    }
    if(currentDataset == "ZZTo2L2Nu_powheg_pythia8_25ns")
    {
             currentProcessClass = "bkgLostLepton"; //@MJ@ TODO mistake in babies, maybe somewhere else too?
    }
    

    //in case of ext -recomputation
    if(currentDataset != storedDataset && currentProcessType == "background") //@MJ@ TODO this can work only with one signal dataset!!!
    {
        storedDataset = currentDataset;
        scale1fbS2 = 1;

        if(currentDataset == "ttbar_singleLeptFromTbar_madgraph_pythia8_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns.root"; 
            getscale1fb2(fBkgName, &scale1fbS2 );

        }
        else if(currentDataset == "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromTbar_madgraph_pythia8_25ns.root"; 
            getscale1fb2(fBkgName, &scale1fbS2 );

        }
        else if(currentDataset == "ttbar_singleLeptFromT_madgraph_pythia8_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns.root"; 
            getscale1fb2(fBkgName, &scale1fbS2 );
        }
        else if(currentDataset == "ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromT_madgraph_pythia8_25ns.root"; 
            getscale1fb2(fBkgName, &scale1fbS2 );
        }
        else if(currentDataset == "ttbar_diLept_madgraph_pythia8_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_diLept_madgraph_pythia8_ext1_25ns.root"; 
            getscale1fb2(fBkgName, &scale1fbS2 );
        }
        else if(currentDataset == "ttbar_diLept_madgraph_pythia8_ext1_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_diLept_madgraph_pythia8_25ns.root"; 
            getscale1fb2(fBkgName, &scale1fbS2 );
        }
        else
        {
            scale1fbS2 = 1;
        }
    }
    
    float weightLumi = getWeight(currentProcessType, GetLumi(), scale1fbS2); 
    vector<float> weightV;
    weightV.clear();

    //@MJ@ TODO could be done better
    float nEvents =  myEvent.wNormalization.at(22);
    float ISRNJ = myEvent.weight_ISRnjets*( nEvents / myEvent.wNormalization.at(25));
    float ISRNJ_UP = myEvent.weight_ISRnjets_UP*( nEvents / myEvent.wNormalization.at(26));
    float ISRNJ_DOWN = myEvent.weight_ISRnjets_DN*( nEvents / myEvent.wNormalization.at(27));


    for(uint32_t SR=0; SR<62; SR++) //@MJ@ both signal and control regions
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
        //ISRnjetsdown
        if(counter == 1) statnames << "ISRnjetsDown" << endl;
        w = GetLumi() *  myEvent.scale1fb* myEvent.weight_PU *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) *ISRNJ_DOWN * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ; //TODO
        weightV.push_back(w);
        //ISRnjetsup
        if(counter == 1) statnames << "ISRnjetsUp" << endl;
        w = GetLumi() *  myEvent.scale1fb* myEvent.weight_PU * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) *ISRNJ_UP* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )  ; //TODO
        weightV.push_back(w);
//@MJ@ TODO misssing tau efficiency, CR2l trigger,JES, MET resolution, MET SF
    }    
      
    vector<string> theReg;
    GetRegionTagList(&theReg);
    if( weightV.size() != theReg.size())
        throw std::runtime_error("vector of weights does not have same size as regions, the weights will not be correctly assesed");

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
    // ######################
    //  Plot configuration and production
    // ######################

    // Schedule plots
    //

    //SchedulePlots("1DSuperimposed");

    // Config plots

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

vector<string> yieldSR = { "SR1l_A_250lessMETless350" , "SR1l_A_250lessMETless350LSFdown" , "SR1l_A_250lessMETless350LSFup" , "SR1l_A_250lessMETless350BTlightDown" , "SR1l_A_250lessMETless350BTlightUp" , "SR1l_A_250lessMETless350BTheavyDown" , "SR1l_A_250lessMETless350BTheavyUp" , "SR1l_A_250lessMETless350PUdown" , "SR1l_A_250lessMETless350PUup" , "SR1l_A_250lessMETless350PDFdown" , "SR1l_A_250lessMETless350PDFup" , "SR1l_A_250lessMETless350alphaSdown" , "SR1l_A_250lessMETless350alphaSup" , "SR1l_A_250lessMETless350Q2down" , "SR1l_A_250lessMETless350Q2up" , "SR1l_A_250lessMETless350ISRnjetsDown" , "SR1l_A_250lessMETless350ISRnjetsUp" , "SR1l_A_350lessMETless450" , "SR1l_A_350lessMETless450LSFdown" , "SR1l_A_350lessMETless450LSFup" , "SR1l_A_350lessMETless450BTlightDown" , "SR1l_A_350lessMETless450BTlightUp" , "SR1l_A_350lessMETless450BTheavyDown" , "SR1l_A_350lessMETless450BTheavyUp" , "SR1l_A_350lessMETless450PUdown" , "SR1l_A_350lessMETless450PUup" , "SR1l_A_350lessMETless450PDFdown" , "SR1l_A_350lessMETless450PDFup" , "SR1l_A_350lessMETless450alphaSdown" , "SR1l_A_350lessMETless450alphaSup" , "SR1l_A_350lessMETless450Q2down" , "SR1l_A_350lessMETless450Q2up" , "SR1l_A_350lessMETless450ISRnjetsDown" , "SR1l_A_350lessMETless450ISRnjetsUp" , "SR1l_A_450lessMETless600" , "SR1l_A_450lessMETless600LSFdown" , "SR1l_A_450lessMETless600LSFup" , "SR1l_A_450lessMETless600BTlightDown" , "SR1l_A_450lessMETless600BTlightUp" , "SR1l_A_450lessMETless600BTheavyDown" , "SR1l_A_450lessMETless600BTheavyUp" , "SR1l_A_450lessMETless600PUdown" , "SR1l_A_450lessMETless600PUup" , "SR1l_A_450lessMETless600PDFdown" , "SR1l_A_450lessMETless600PDFup" , "SR1l_A_450lessMETless600alphaSdown" , "SR1l_A_450lessMETless600alphaSup" , "SR1l_A_450lessMETless600Q2down" , "SR1l_A_450lessMETless600Q2up" , "SR1l_A_450lessMETless600ISRnjetsDown" , "SR1l_A_450lessMETless600ISRnjetsUp" , "SR1l_A_600lessMETlessInf" , "SR1l_A_600lessMETlessInfLSFdown" , "SR1l_A_600lessMETlessInfLSFup" , "SR1l_A_600lessMETlessInfBTlightDown" , "SR1l_A_600lessMETlessInfBTlightUp" , "SR1l_A_600lessMETlessInfBTheavyDown" , "SR1l_A_600lessMETlessInfBTheavyUp" , "SR1l_A_600lessMETlessInfPUdown" , "SR1l_A_600lessMETlessInfPUup" , "SR1l_A_600lessMETlessInfPDFdown" , "SR1l_A_600lessMETlessInfPDFup" , "SR1l_A_600lessMETlessInfalphaSdown" , "SR1l_A_600lessMETlessInfalphaSup" , "SR1l_A_600lessMETlessInfQ2down" , "SR1l_A_600lessMETlessInfQ2up" , "SR1l_A_600lessMETlessInfISRnjetsDown" , "SR1l_A_600lessMETlessInfISRnjetsUp" , "SR1l_B_250lessMETless450" , "SR1l_B_250lessMETless450LSFdown" , "SR1l_B_250lessMETless450LSFup" , "SR1l_B_250lessMETless450BTlightDown" , "SR1l_B_250lessMETless450BTlightUp" , "SR1l_B_250lessMETless450BTheavyDown" , "SR1l_B_250lessMETless450BTheavyUp" , "SR1l_B_250lessMETless450PUdown" , "SR1l_B_250lessMETless450PUup" , "SR1l_B_250lessMETless450PDFdown" , "SR1l_B_250lessMETless450PDFup" , "SR1l_B_250lessMETless450alphaSdown" , "SR1l_B_250lessMETless450alphaSup" , "SR1l_B_250lessMETless450Q2down" , "SR1l_B_250lessMETless450Q2up" , "SR1l_B_250lessMETless450ISRnjetsDown" , "SR1l_B_250lessMETless450ISRnjetsUp" , "SR1l_B_450lessMETless600" , "SR1l_B_450lessMETless600LSFdown" , "SR1l_B_450lessMETless600LSFup" , "SR1l_B_450lessMETless600BTlightDown" , "SR1l_B_450lessMETless600BTlightUp" , "SR1l_B_450lessMETless600BTheavyDown" , "SR1l_B_450lessMETless600BTheavyUp" , "SR1l_B_450lessMETless600PUdown" , "SR1l_B_450lessMETless600PUup" , "SR1l_B_450lessMETless600PDFdown" , "SR1l_B_450lessMETless600PDFup" , "SR1l_B_450lessMETless600alphaSdown" , "SR1l_B_450lessMETless600alphaSup" , "SR1l_B_450lessMETless600Q2down" , "SR1l_B_450lessMETless600Q2up" , "SR1l_B_450lessMETless600ISRnjetsDown" , "SR1l_B_450lessMETless600ISRnjetsUp" , "SR1l_B_600lessMETlessInf" , "SR1l_B_600lessMETlessInfLSFdown" , "SR1l_B_600lessMETlessInfLSFup" , "SR1l_B_600lessMETlessInfBTlightDown" , "SR1l_B_600lessMETlessInfBTlightUp" , "SR1l_B_600lessMETlessInfBTheavyDown" , "SR1l_B_600lessMETlessInfBTheavyUp" , "SR1l_B_600lessMETlessInfPUdown" , "SR1l_B_600lessMETlessInfPUup" , "SR1l_B_600lessMETlessInfPDFdown" , "SR1l_B_600lessMETlessInfPDFup" , "SR1l_B_600lessMETlessInfalphaSdown" , "SR1l_B_600lessMETlessInfalphaSup" , "SR1l_B_600lessMETlessInfQ2down" , "SR1l_B_600lessMETlessInfQ2up" , "SR1l_B_600lessMETlessInfISRnjetsDown" , "SR1l_B_600lessMETlessInfISRnjetsUp" , "SR1l_C_250lessMETless350" , "SR1l_C_250lessMETless350LSFdown" , "SR1l_C_250lessMETless350LSFup" , "SR1l_C_250lessMETless350BTlightDown" , "SR1l_C_250lessMETless350BTlightUp" , "SR1l_C_250lessMETless350BTheavyDown" , "SR1l_C_250lessMETless350BTheavyUp" , "SR1l_C_250lessMETless350PUdown" , "SR1l_C_250lessMETless350PUup" , "SR1l_C_250lessMETless350PDFdown" , "SR1l_C_250lessMETless350PDFup" , "SR1l_C_250lessMETless350alphaSdown" , "SR1l_C_250lessMETless350alphaSup" , "SR1l_C_250lessMETless350Q2down" , "SR1l_C_250lessMETless350Q2up" , "SR1l_C_250lessMETless350ISRnjetsDown" , "SR1l_C_250lessMETless350ISRnjetsUp" , "SR1l_C_350lessMETless450" , "SR1l_C_350lessMETless450LSFdown" , "SR1l_C_350lessMETless450LSFup" , "SR1l_C_350lessMETless450BTlightDown" , "SR1l_C_350lessMETless450BTlightUp" , "SR1l_C_350lessMETless450BTheavyDown" , "SR1l_C_350lessMETless450BTheavyUp" , "SR1l_C_350lessMETless450PUdown" , "SR1l_C_350lessMETless450PUup" , "SR1l_C_350lessMETless450PDFdown" , "SR1l_C_350lessMETless450PDFup" , "SR1l_C_350lessMETless450alphaSdown" , "SR1l_C_350lessMETless450alphaSup" , "SR1l_C_350lessMETless450Q2down" , "SR1l_C_350lessMETless450Q2up" , "SR1l_C_350lessMETless450ISRnjetsDown" , "SR1l_C_350lessMETless450ISRnjetsUp" , "SR1l_C_450lessMETless550" , "SR1l_C_450lessMETless550LSFdown" , "SR1l_C_450lessMETless550LSFup" , "SR1l_C_450lessMETless550BTlightDown" , "SR1l_C_450lessMETless550BTlightUp" , "SR1l_C_450lessMETless550BTheavyDown" , "SR1l_C_450lessMETless550BTheavyUp" , "SR1l_C_450lessMETless550PUdown" , "SR1l_C_450lessMETless550PUup" , "SR1l_C_450lessMETless550PDFdown" , "SR1l_C_450lessMETless550PDFup" , "SR1l_C_450lessMETless550alphaSdown" , "SR1l_C_450lessMETless550alphaSup" , "SR1l_C_450lessMETless550Q2down" , "SR1l_C_450lessMETless550Q2up" , "SR1l_C_450lessMETless550ISRnjetsDown" , "SR1l_C_450lessMETless550ISRnjetsUp" , "SR1l_C_550lessMETless650" , "SR1l_C_550lessMETless650LSFdown" , "SR1l_C_550lessMETless650LSFup" , "SR1l_C_550lessMETless650BTlightDown" , "SR1l_C_550lessMETless650BTlightUp" , "SR1l_C_550lessMETless650BTheavyDown" , "SR1l_C_550lessMETless650BTheavyUp" , "SR1l_C_550lessMETless650PUdown" , "SR1l_C_550lessMETless650PUup" , "SR1l_C_550lessMETless650PDFdown" , "SR1l_C_550lessMETless650PDFup" , "SR1l_C_550lessMETless650alphaSdown" , "SR1l_C_550lessMETless650alphaSup" , "SR1l_C_550lessMETless650Q2down" , "SR1l_C_550lessMETless650Q2up" , "SR1l_C_550lessMETless650ISRnjetsDown" , "SR1l_C_550lessMETless650ISRnjetsUp" , "SR1l_C_650lessMETlessInf" , "SR1l_C_650lessMETlessInfLSFdown" , "SR1l_C_650lessMETlessInfLSFup" , "SR1l_C_650lessMETlessInfBTlightDown" , "SR1l_C_650lessMETlessInfBTlightUp" , "SR1l_C_650lessMETlessInfBTheavyDown" , "SR1l_C_650lessMETlessInfBTheavyUp" , "SR1l_C_650lessMETlessInfPUdown" , "SR1l_C_650lessMETlessInfPUup" , "SR1l_C_650lessMETlessInfPDFdown" , "SR1l_C_650lessMETlessInfPDFup" , "SR1l_C_650lessMETlessInfalphaSdown" , "SR1l_C_650lessMETlessInfalphaSup" , "SR1l_C_650lessMETlessInfQ2down" , "SR1l_C_650lessMETlessInfQ2up" , "SR1l_C_650lessMETlessInfISRnjetsDown" , "SR1l_C_650lessMETlessInfISRnjetsUp" , "SR1l_D_250lessMETless350" , "SR1l_D_250lessMETless350LSFdown" , "SR1l_D_250lessMETless350LSFup" , "SR1l_D_250lessMETless350BTlightDown" , "SR1l_D_250lessMETless350BTlightUp" , "SR1l_D_250lessMETless350BTheavyDown" , "SR1l_D_250lessMETless350BTheavyUp" , "SR1l_D_250lessMETless350PUdown" , "SR1l_D_250lessMETless350PUup" , "SR1l_D_250lessMETless350PDFdown" , "SR1l_D_250lessMETless350PDFup" , "SR1l_D_250lessMETless350alphaSdown" , "SR1l_D_250lessMETless350alphaSup" , "SR1l_D_250lessMETless350Q2down" , "SR1l_D_250lessMETless350Q2up" , "SR1l_D_250lessMETless350ISRnjetsDown" , "SR1l_D_250lessMETless350ISRnjetsUp" , "SR1l_D_350lessMETless450" , "SR1l_D_350lessMETless450LSFdown" , "SR1l_D_350lessMETless450LSFup" , "SR1l_D_350lessMETless450BTlightDown" , "SR1l_D_350lessMETless450BTlightUp" , "SR1l_D_350lessMETless450BTheavyDown" , "SR1l_D_350lessMETless450BTheavyUp" , "SR1l_D_350lessMETless450PUdown" , "SR1l_D_350lessMETless450PUup" , "SR1l_D_350lessMETless450PDFdown" , "SR1l_D_350lessMETless450PDFup" , "SR1l_D_350lessMETless450alphaSdown" , "SR1l_D_350lessMETless450alphaSup" , "SR1l_D_350lessMETless450Q2down" , "SR1l_D_350lessMETless450Q2up" , "SR1l_D_350lessMETless450ISRnjetsDown" , "SR1l_D_350lessMETless450ISRnjetsUp" , "SR1l_D_450lessMETless550" , "SR1l_D_450lessMETless550LSFdown" , "SR1l_D_450lessMETless550LSFup" , "SR1l_D_450lessMETless550BTlightDown" , "SR1l_D_450lessMETless550BTlightUp" , "SR1l_D_450lessMETless550BTheavyDown" , "SR1l_D_450lessMETless550BTheavyUp" , "SR1l_D_450lessMETless550PUdown" , "SR1l_D_450lessMETless550PUup" , "SR1l_D_450lessMETless550PDFdown" , "SR1l_D_450lessMETless550PDFup" , "SR1l_D_450lessMETless550alphaSdown" , "SR1l_D_450lessMETless550alphaSup" , "SR1l_D_450lessMETless550Q2down" , "SR1l_D_450lessMETless550Q2up" , "SR1l_D_450lessMETless550ISRnjetsDown" , "SR1l_D_450lessMETless550ISRnjetsUp" , "SR1l_D_550lessMETlessInf" , "SR1l_D_550lessMETlessInfLSFdown" , "SR1l_D_550lessMETlessInfLSFup" , "SR1l_D_550lessMETlessInfBTlightDown" , "SR1l_D_550lessMETlessInfBTlightUp" , "SR1l_D_550lessMETlessInfBTheavyDown" , "SR1l_D_550lessMETlessInfBTheavyUp" , "SR1l_D_550lessMETlessInfPUdown" , "SR1l_D_550lessMETlessInfPUup" , "SR1l_D_550lessMETlessInfPDFdown" , "SR1l_D_550lessMETlessInfPDFup" , "SR1l_D_550lessMETlessInfalphaSdown" , "SR1l_D_550lessMETlessInfalphaSup" , "SR1l_D_550lessMETlessInfQ2down" , "SR1l_D_550lessMETlessInfQ2up" , "SR1l_D_550lessMETlessInfISRnjetsDown" , "SR1l_D_550lessMETlessInfISRnjetsUp" , "SR1l_E_250lessMETless350" , "SR1l_E_250lessMETless350LSFdown" , "SR1l_E_250lessMETless350LSFup" , "SR1l_E_250lessMETless350BTlightDown" , "SR1l_E_250lessMETless350BTlightUp" , "SR1l_E_250lessMETless350BTheavyDown" , "SR1l_E_250lessMETless350BTheavyUp" , "SR1l_E_250lessMETless350PUdown" , "SR1l_E_250lessMETless350PUup" , "SR1l_E_250lessMETless350PDFdown" , "SR1l_E_250lessMETless350PDFup" , "SR1l_E_250lessMETless350alphaSdown" , "SR1l_E_250lessMETless350alphaSup" , "SR1l_E_250lessMETless350Q2down" , "SR1l_E_250lessMETless350Q2up" , "SR1l_E_250lessMETless350ISRnjetsDown" , "SR1l_E_250lessMETless350ISRnjetsUp" , "SR1l_E_350lessMETless550" , "SR1l_E_350lessMETless550LSFdown" , "SR1l_E_350lessMETless550LSFup" , "SR1l_E_350lessMETless550BTlightDown" , "SR1l_E_350lessMETless550BTlightUp" , "SR1l_E_350lessMETless550BTheavyDown" , "SR1l_E_350lessMETless550BTheavyUp" , "SR1l_E_350lessMETless550PUdown" , "SR1l_E_350lessMETless550PUup" , "SR1l_E_350lessMETless550PDFdown" , "SR1l_E_350lessMETless550PDFup" , "SR1l_E_350lessMETless550alphaSdown" , "SR1l_E_350lessMETless550alphaSup" , "SR1l_E_350lessMETless550Q2down" , "SR1l_E_350lessMETless550Q2up" , "SR1l_E_350lessMETless550ISRnjetsDown" , "SR1l_E_350lessMETless550ISRnjetsUp" , "SR1l_E_550lessMETlessInf" , "SR1l_E_550lessMETlessInfLSFdown" , "SR1l_E_550lessMETlessInfLSFup" , "SR1l_E_550lessMETlessInfBTlightDown" , "SR1l_E_550lessMETlessInfBTlightUp" , "SR1l_E_550lessMETlessInfBTheavyDown" , "SR1l_E_550lessMETlessInfBTheavyUp" , "SR1l_E_550lessMETlessInfPUdown" , "SR1l_E_550lessMETlessInfPUup" , "SR1l_E_550lessMETlessInfPDFdown" , "SR1l_E_550lessMETlessInfPDFup" , "SR1l_E_550lessMETlessInfalphaSdown" , "SR1l_E_550lessMETlessInfalphaSup" , "SR1l_E_550lessMETlessInfQ2down" , "SR1l_E_550lessMETlessInfQ2up" , "SR1l_E_550lessMETlessInfISRnjetsDown" , "SR1l_E_550lessMETlessInfISRnjetsUp" , "SR1l_F_250lessMETless450" , "SR1l_F_250lessMETless450LSFdown" , "SR1l_F_250lessMETless450LSFup" , "SR1l_F_250lessMETless450BTlightDown" , "SR1l_F_250lessMETless450BTlightUp" , "SR1l_F_250lessMETless450BTheavyDown" , "SR1l_F_250lessMETless450BTheavyUp" , "SR1l_F_250lessMETless450PUdown" , "SR1l_F_250lessMETless450PUup" , "SR1l_F_250lessMETless450PDFdown" , "SR1l_F_250lessMETless450PDFup" , "SR1l_F_250lessMETless450alphaSdown" , "SR1l_F_250lessMETless450alphaSup" , "SR1l_F_250lessMETless450Q2down" , "SR1l_F_250lessMETless450Q2up" , "SR1l_F_250lessMETless450ISRnjetsDown" , "SR1l_F_250lessMETless450ISRnjetsUp" , "SR1l_F_450lessMETlessInf" , "SR1l_F_450lessMETlessInfLSFdown" , "SR1l_F_450lessMETlessInfLSFup" , "SR1l_F_450lessMETlessInfBTlightDown" , "SR1l_F_450lessMETlessInfBTlightUp" , "SR1l_F_450lessMETlessInfBTheavyDown" , "SR1l_F_450lessMETlessInfBTheavyUp" , "SR1l_F_450lessMETlessInfPUdown" , "SR1l_F_450lessMETlessInfPUup" , "SR1l_F_450lessMETlessInfPDFdown" , "SR1l_F_450lessMETlessInfPDFup" , "SR1l_F_450lessMETlessInfalphaSdown" , "SR1l_F_450lessMETlessInfalphaSup" , "SR1l_F_450lessMETlessInfQ2down" , "SR1l_F_450lessMETlessInfQ2up" , "SR1l_F_450lessMETlessInfISRnjetsDown" , "SR1l_F_450lessMETlessInfISRnjetsUp" , "SR1l_G_250lessMETless350" , "SR1l_G_250lessMETless350LSFdown" , "SR1l_G_250lessMETless350LSFup" , "SR1l_G_250lessMETless350BTlightDown" , "SR1l_G_250lessMETless350BTlightUp" , "SR1l_G_250lessMETless350BTheavyDown" , "SR1l_G_250lessMETless350BTheavyUp" , "SR1l_G_250lessMETless350PUdown" , "SR1l_G_250lessMETless350PUup" , "SR1l_G_250lessMETless350PDFdown" , "SR1l_G_250lessMETless350PDFup" , "SR1l_G_250lessMETless350alphaSdown" , "SR1l_G_250lessMETless350alphaSup" , "SR1l_G_250lessMETless350Q2down" , "SR1l_G_250lessMETless350Q2up" , "SR1l_G_250lessMETless350ISRnjetsDown" , "SR1l_G_250lessMETless350ISRnjetsUp" , "SR1l_G_350lessMETless450" , "SR1l_G_350lessMETless450LSFdown" , "SR1l_G_350lessMETless450LSFup" , "SR1l_G_350lessMETless450BTlightDown" , "SR1l_G_350lessMETless450BTlightUp" , "SR1l_G_350lessMETless450BTheavyDown" , "SR1l_G_350lessMETless450BTheavyUp" , "SR1l_G_350lessMETless450PUdown" , "SR1l_G_350lessMETless450PUup" , "SR1l_G_350lessMETless450PDFdown" , "SR1l_G_350lessMETless450PDFup" , "SR1l_G_350lessMETless450alphaSdown" , "SR1l_G_350lessMETless450alphaSup" , "SR1l_G_350lessMETless450Q2down" , "SR1l_G_350lessMETless450Q2up" , "SR1l_G_350lessMETless450ISRnjetsDown" , "SR1l_G_350lessMETless450ISRnjetsUp" , "SR1l_G_450lessMETless600" , "SR1l_G_450lessMETless600LSFdown" , "SR1l_G_450lessMETless600LSFup" , "SR1l_G_450lessMETless600BTlightDown" , "SR1l_G_450lessMETless600BTlightUp" , "SR1l_G_450lessMETless600BTheavyDown" , "SR1l_G_450lessMETless600BTheavyUp" , "SR1l_G_450lessMETless600PUdown" , "SR1l_G_450lessMETless600PUup" , "SR1l_G_450lessMETless600PDFdown" , "SR1l_G_450lessMETless600PDFup" , "SR1l_G_450lessMETless600alphaSdown" , "SR1l_G_450lessMETless600alphaSup" , "SR1l_G_450lessMETless600Q2down" , "SR1l_G_450lessMETless600Q2up" , "SR1l_G_450lessMETless600ISRnjetsDown" , "SR1l_G_450lessMETless600ISRnjetsUp" , "SR1l_G_600lessMETlessInf" , "SR1l_G_600lessMETlessInfLSFdown" , "SR1l_G_600lessMETlessInfLSFup" , "SR1l_G_600lessMETlessInfBTlightDown" , "SR1l_G_600lessMETlessInfBTlightUp" , "SR1l_G_600lessMETlessInfBTheavyDown" , "SR1l_G_600lessMETlessInfBTheavyUp" , "SR1l_G_600lessMETlessInfPUdown" , "SR1l_G_600lessMETlessInfPUup" , "SR1l_G_600lessMETlessInfPDFdown" , "SR1l_G_600lessMETlessInfPDFup" , "SR1l_G_600lessMETlessInfalphaSdown" , "SR1l_G_600lessMETlessInfalphaSup" , "SR1l_G_600lessMETlessInfQ2down" , "SR1l_G_600lessMETlessInfQ2up" , "SR1l_G_600lessMETlessInfISRnjetsDown" , "SR1l_G_600lessMETlessInfISRnjetsUp" , "SR1l_H_250lessMETless450" , "SR1l_H_250lessMETless450LSFdown" , "SR1l_H_250lessMETless450LSFup" , "SR1l_H_250lessMETless450BTlightDown" , "SR1l_H_250lessMETless450BTlightUp" , "SR1l_H_250lessMETless450BTheavyDown" , "SR1l_H_250lessMETless450BTheavyUp" , "SR1l_H_250lessMETless450PUdown" , "SR1l_H_250lessMETless450PUup" , "SR1l_H_250lessMETless450PDFdown" , "SR1l_H_250lessMETless450PDFup" , "SR1l_H_250lessMETless450alphaSdown" , "SR1l_H_250lessMETless450alphaSup" , "SR1l_H_250lessMETless450Q2down" , "SR1l_H_250lessMETless450Q2up" , "SR1l_H_250lessMETless450ISRnjetsDown" , "SR1l_H_250lessMETless450ISRnjetsUp" , "SR1l_H_450lessMETlessInf" , "SR1l_H_450lessMETlessInfLSFdown" , "SR1l_H_450lessMETlessInfLSFup" , "SR1l_H_450lessMETlessInfBTlightDown" , "SR1l_H_450lessMETlessInfBTlightUp" , "SR1l_H_450lessMETlessInfBTheavyDown" , "SR1l_H_450lessMETlessInfBTheavyUp" , "SR1l_H_450lessMETlessInfPUdown" , "SR1l_H_450lessMETlessInfPUup" , "SR1l_H_450lessMETlessInfPDFdown" , "SR1l_H_450lessMETlessInfPDFup" , "SR1l_H_450lessMETlessInfalphaSdown" , "SR1l_H_450lessMETlessInfalphaSup" , "SR1l_H_450lessMETlessInfQ2down" , "SR1l_H_450lessMETlessInfQ2up" , "SR1l_H_450lessMETlessInfISRnjetsDown" , "SR1l_H_450lessMETlessInfISRnjetsUp" , "SR1l_I_250lessMETless350" , "SR1l_I_250lessMETless350LSFdown" , "SR1l_I_250lessMETless350LSFup" , "SR1l_I_250lessMETless350BTlightDown" , "SR1l_I_250lessMETless350BTlightUp" , "SR1l_I_250lessMETless350BTheavyDown" , "SR1l_I_250lessMETless350BTheavyUp" , "SR1l_I_250lessMETless350PUdown" , "SR1l_I_250lessMETless350PUup" , "SR1l_I_250lessMETless350PDFdown" , "SR1l_I_250lessMETless350PDFup" , "SR1l_I_250lessMETless350alphaSdown" , "SR1l_I_250lessMETless350alphaSup" , "SR1l_I_250lessMETless350Q2down" , "SR1l_I_250lessMETless350Q2up" , "SR1l_I_250lessMETless350ISRnjetsDown" , "SR1l_I_250lessMETless350ISRnjetsUp" , "SR1l_I_350lessMETless450" , "SR1l_I_350lessMETless450LSFdown" , "SR1l_I_350lessMETless450LSFup" , "SR1l_I_350lessMETless450BTlightDown" , "SR1l_I_350lessMETless450BTlightUp" , "SR1l_I_350lessMETless450BTheavyDown" , "SR1l_I_350lessMETless450BTheavyUp" , "SR1l_I_350lessMETless450PUdown" , "SR1l_I_350lessMETless450PUup" , "SR1l_I_350lessMETless450PDFdown" , "SR1l_I_350lessMETless450PDFup" , "SR1l_I_350lessMETless450alphaSdown" , "SR1l_I_350lessMETless450alphaSup" , "SR1l_I_350lessMETless450Q2down" , "SR1l_I_350lessMETless450Q2up" , "SR1l_I_350lessMETless450ISRnjetsDown" , "SR1l_I_350lessMETless450ISRnjetsUp" , "SR1l_I_450lessMETless550" , "SR1l_I_450lessMETless550LSFdown" , "SR1l_I_450lessMETless550LSFup" , "SR1l_I_450lessMETless550BTlightDown" , "SR1l_I_450lessMETless550BTlightUp" , "SR1l_I_450lessMETless550BTheavyDown" , "SR1l_I_450lessMETless550BTheavyUp" , "SR1l_I_450lessMETless550PUdown" , "SR1l_I_450lessMETless550PUup" , "SR1l_I_450lessMETless550PDFdown" , "SR1l_I_450lessMETless550PDFup" , "SR1l_I_450lessMETless550alphaSdown" , "SR1l_I_450lessMETless550alphaSup" , "SR1l_I_450lessMETless550Q2down" , "SR1l_I_450lessMETless550Q2up" , "SR1l_I_450lessMETless550ISRnjetsDown" , "SR1l_I_450lessMETless550ISRnjetsUp" , "SR1l_I_550lessMETlessInf" , "SR1l_I_550lessMETlessInfLSFdown" , "SR1l_I_550lessMETlessInfLSFup" , "SR1l_I_550lessMETlessInfBTlightDown" , "SR1l_I_550lessMETlessInfBTlightUp" , "SR1l_I_550lessMETlessInfBTheavyDown" , "SR1l_I_550lessMETlessInfBTheavyUp" , "SR1l_I_550lessMETlessInfPUdown" , "SR1l_I_550lessMETlessInfPUup" , "SR1l_I_550lessMETlessInfPDFdown" , "SR1l_I_550lessMETlessInfPDFup" , "SR1l_I_550lessMETlessInfalphaSdown" , "SR1l_I_550lessMETlessInfalphaSup" , "SR1l_I_550lessMETlessInfQ2down" , "SR1l_I_550lessMETlessInfQ2up" , "SR1l_I_550lessMETlessInfISRnjetsDown" , "SR1l_I_550lessMETlessInfISRnjetsUp"   };

    TableDataMC(this, yieldSR,"lepChannel",  "keepNegative" ).Print(outputName+ "yieldsSRsyst.tab", 4);
    TableDataMC(this, yieldSR,"lepChannel", "keepNegative" ).PrintLatex(outputName+ "yieldsSRsyst.tex", 4);

vector<string> yieldCR = { "CR2l_A_250lessMETless350" , "CR2l_A_250lessMETless350LSFdown" , "CR2l_A_250lessMETless350LSFup" , "CR2l_A_250lessMETless350BTlightDown" , "CR2l_A_250lessMETless350BTlightUp" , "CR2l_A_250lessMETless350BTheavyDown" , "CR2l_A_250lessMETless350BTheavyUp" , "CR2l_A_250lessMETless350PUdown" , "CR2l_A_250lessMETless350PUup" , "CR2l_A_250lessMETless350PDFdown" , "CR2l_A_250lessMETless350PDFup" , "CR2l_A_250lessMETless350alphaSdown" , "CR2l_A_250lessMETless350alphaSup" , "CR2l_A_250lessMETless350Q2down" , "CR2l_A_250lessMETless350Q2up" , "CR2l_A_250lessMETless350ISRnjetsDown" , "CR2l_A_250lessMETless350ISRnjetsUp" , "CR2l_A_350lessMETless450" , "CR2l_A_350lessMETless450LSFdown" , "CR2l_A_350lessMETless450LSFup" , "CR2l_A_350lessMETless450BTlightDown" , "CR2l_A_350lessMETless450BTlightUp" , "CR2l_A_350lessMETless450BTheavyDown" , "CR2l_A_350lessMETless450BTheavyUp" , "CR2l_A_350lessMETless450PUdown" , "CR2l_A_350lessMETless450PUup" , "CR2l_A_350lessMETless450PDFdown" , "CR2l_A_350lessMETless450PDFup" , "CR2l_A_350lessMETless450alphaSdown" , "CR2l_A_350lessMETless450alphaSup" , "CR2l_A_350lessMETless450Q2down" , "CR2l_A_350lessMETless450Q2up" , "CR2l_A_350lessMETless450ISRnjetsDown" , "CR2l_A_350lessMETless450ISRnjetsUp" , "CR2l_A_450lessMETless600" , "CR2l_A_450lessMETless600LSFdown" , "CR2l_A_450lessMETless600LSFup" , "CR2l_A_450lessMETless600BTlightDown" , "CR2l_A_450lessMETless600BTlightUp" , "CR2l_A_450lessMETless600BTheavyDown" , "CR2l_A_450lessMETless600BTheavyUp" , "CR2l_A_450lessMETless600PUdown" , "CR2l_A_450lessMETless600PUup" , "CR2l_A_450lessMETless600PDFdown" , "CR2l_A_450lessMETless600PDFup" , "CR2l_A_450lessMETless600alphaSdown" , "CR2l_A_450lessMETless600alphaSup" , "CR2l_A_450lessMETless600Q2down" , "CR2l_A_450lessMETless600Q2up" , "CR2l_A_450lessMETless600ISRnjetsDown" , "CR2l_A_450lessMETless600ISRnjetsUp" , "CR2l_A_600lessMETlessInf" , "CR2l_A_600lessMETlessInfLSFdown" , "CR2l_A_600lessMETlessInfLSFup" , "CR2l_A_600lessMETlessInfBTlightDown" , "CR2l_A_600lessMETlessInfBTlightUp" , "CR2l_A_600lessMETlessInfBTheavyDown" , "CR2l_A_600lessMETlessInfBTheavyUp" , "CR2l_A_600lessMETlessInfPUdown" , "CR2l_A_600lessMETlessInfPUup" , "CR2l_A_600lessMETlessInfPDFdown" , "CR2l_A_600lessMETlessInfPDFup" , "CR2l_A_600lessMETlessInfalphaSdown" , "CR2l_A_600lessMETlessInfalphaSup" , "CR2l_A_600lessMETlessInfQ2down" , "CR2l_A_600lessMETlessInfQ2up" , "CR2l_A_600lessMETlessInfISRnjetsDown" , "CR2l_A_600lessMETlessInfISRnjetsUp" , "CR2l_B_250lessMETless450" , "CR2l_B_250lessMETless450LSFdown" , "CR2l_B_250lessMETless450LSFup" , "CR2l_B_250lessMETless450BTlightDown" , "CR2l_B_250lessMETless450BTlightUp" , "CR2l_B_250lessMETless450BTheavyDown" , "CR2l_B_250lessMETless450BTheavyUp" , "CR2l_B_250lessMETless450PUdown" , "CR2l_B_250lessMETless450PUup" , "CR2l_B_250lessMETless450PDFdown" , "CR2l_B_250lessMETless450PDFup" , "CR2l_B_250lessMETless450alphaSdown" , "CR2l_B_250lessMETless450alphaSup" , "CR2l_B_250lessMETless450Q2down" , "CR2l_B_250lessMETless450Q2up" , "CR2l_B_250lessMETless450ISRnjetsDown" , "CR2l_B_250lessMETless450ISRnjetsUp" , "CR2l_B_450lessMETless600" , "CR2l_B_450lessMETless600LSFdown" , "CR2l_B_450lessMETless600LSFup" , "CR2l_B_450lessMETless600BTlightDown" , "CR2l_B_450lessMETless600BTlightUp" , "CR2l_B_450lessMETless600BTheavyDown" , "CR2l_B_450lessMETless600BTheavyUp" , "CR2l_B_450lessMETless600PUdown" , "CR2l_B_450lessMETless600PUup" , "CR2l_B_450lessMETless600PDFdown" , "CR2l_B_450lessMETless600PDFup" , "CR2l_B_450lessMETless600alphaSdown" , "CR2l_B_450lessMETless600alphaSup" , "CR2l_B_450lessMETless600Q2down" , "CR2l_B_450lessMETless600Q2up" , "CR2l_B_450lessMETless600ISRnjetsDown" , "CR2l_B_450lessMETless600ISRnjetsUp" , "CR2l_B_600lessMETlessInf" , "CR2l_B_600lessMETlessInfLSFdown" , "CR2l_B_600lessMETlessInfLSFup" , "CR2l_B_600lessMETlessInfBTlightDown" , "CR2l_B_600lessMETlessInfBTlightUp" , "CR2l_B_600lessMETlessInfBTheavyDown" , "CR2l_B_600lessMETlessInfBTheavyUp" , "CR2l_B_600lessMETlessInfPUdown" , "CR2l_B_600lessMETlessInfPUup" , "CR2l_B_600lessMETlessInfPDFdown" , "CR2l_B_600lessMETlessInfPDFup" , "CR2l_B_600lessMETlessInfalphaSdown" , "CR2l_B_600lessMETlessInfalphaSup" , "CR2l_B_600lessMETlessInfQ2down" , "CR2l_B_600lessMETlessInfQ2up" , "CR2l_B_600lessMETlessInfISRnjetsDown" , "CR2l_B_600lessMETlessInfISRnjetsUp" , "CR2l_C_250lessMETless350" , "CR2l_C_250lessMETless350LSFdown" , "CR2l_C_250lessMETless350LSFup" , "CR2l_C_250lessMETless350BTlightDown" , "CR2l_C_250lessMETless350BTlightUp" , "CR2l_C_250lessMETless350BTheavyDown" , "CR2l_C_250lessMETless350BTheavyUp" , "CR2l_C_250lessMETless350PUdown" , "CR2l_C_250lessMETless350PUup" , "CR2l_C_250lessMETless350PDFdown" , "CR2l_C_250lessMETless350PDFup" , "CR2l_C_250lessMETless350alphaSdown" , "CR2l_C_250lessMETless350alphaSup" , "CR2l_C_250lessMETless350Q2down" , "CR2l_C_250lessMETless350Q2up" , "CR2l_C_250lessMETless350ISRnjetsDown" , "CR2l_C_250lessMETless350ISRnjetsUp" , "CR2l_C_350lessMETless450" , "CR2l_C_350lessMETless450LSFdown" , "CR2l_C_350lessMETless450LSFup" , "CR2l_C_350lessMETless450BTlightDown" , "CR2l_C_350lessMETless450BTlightUp" , "CR2l_C_350lessMETless450BTheavyDown" , "CR2l_C_350lessMETless450BTheavyUp" , "CR2l_C_350lessMETless450PUdown" , "CR2l_C_350lessMETless450PUup" , "CR2l_C_350lessMETless450PDFdown" , "CR2l_C_350lessMETless450PDFup" , "CR2l_C_350lessMETless450alphaSdown" , "CR2l_C_350lessMETless450alphaSup" , "CR2l_C_350lessMETless450Q2down" , "CR2l_C_350lessMETless450Q2up" , "CR2l_C_350lessMETless450ISRnjetsDown" , "CR2l_C_350lessMETless450ISRnjetsUp" , "CR2l_C_450lessMETless550" , "CR2l_C_450lessMETless550LSFdown" , "CR2l_C_450lessMETless550LSFup" , "CR2l_C_450lessMETless550BTlightDown" , "CR2l_C_450lessMETless550BTlightUp" , "CR2l_C_450lessMETless550BTheavyDown" , "CR2l_C_450lessMETless550BTheavyUp" , "CR2l_C_450lessMETless550PUdown" , "CR2l_C_450lessMETless550PUup" , "CR2l_C_450lessMETless550PDFdown" , "CR2l_C_450lessMETless550PDFup" , "CR2l_C_450lessMETless550alphaSdown" , "CR2l_C_450lessMETless550alphaSup" , "CR2l_C_450lessMETless550Q2down" , "CR2l_C_450lessMETless550Q2up" , "CR2l_C_450lessMETless550ISRnjetsDown" , "CR2l_C_450lessMETless550ISRnjetsUp" , "CR2l_C_550lessMETless650" , "CR2l_C_550lessMETless650LSFdown" , "CR2l_C_550lessMETless650LSFup" , "CR2l_C_550lessMETless650BTlightDown" , "CR2l_C_550lessMETless650BTlightUp" , "CR2l_C_550lessMETless650BTheavyDown" , "CR2l_C_550lessMETless650BTheavyUp" , "CR2l_C_550lessMETless650PUdown" , "CR2l_C_550lessMETless650PUup" , "CR2l_C_550lessMETless650PDFdown" , "CR2l_C_550lessMETless650PDFup" , "CR2l_C_550lessMETless650alphaSdown" , "CR2l_C_550lessMETless650alphaSup" , "CR2l_C_550lessMETless650Q2down" , "CR2l_C_550lessMETless650Q2up" , "CR2l_C_550lessMETless650ISRnjetsDown" , "CR2l_C_550lessMETless650ISRnjetsUp" , "CR2l_C_650lessMETlessInf" , "CR2l_C_650lessMETlessInfLSFdown" , "CR2l_C_650lessMETlessInfLSFup" , "CR2l_C_650lessMETlessInfBTlightDown" , "CR2l_C_650lessMETlessInfBTlightUp" , "CR2l_C_650lessMETlessInfBTheavyDown" , "CR2l_C_650lessMETlessInfBTheavyUp" , "CR2l_C_650lessMETlessInfPUdown" , "CR2l_C_650lessMETlessInfPUup" , "CR2l_C_650lessMETlessInfPDFdown" , "CR2l_C_650lessMETlessInfPDFup" , "CR2l_C_650lessMETlessInfalphaSdown" , "CR2l_C_650lessMETlessInfalphaSup" , "CR2l_C_650lessMETlessInfQ2down" , "CR2l_C_650lessMETlessInfQ2up" , "CR2l_C_650lessMETlessInfISRnjetsDown" , "CR2l_C_650lessMETlessInfISRnjetsUp" , "CR2l_D_250lessMETless350" , "CR2l_D_250lessMETless350LSFdown" , "CR2l_D_250lessMETless350LSFup" , "CR2l_D_250lessMETless350BTlightDown" , "CR2l_D_250lessMETless350BTlightUp" , "CR2l_D_250lessMETless350BTheavyDown" , "CR2l_D_250lessMETless350BTheavyUp" , "CR2l_D_250lessMETless350PUdown" , "CR2l_D_250lessMETless350PUup" , "CR2l_D_250lessMETless350PDFdown" , "CR2l_D_250lessMETless350PDFup" , "CR2l_D_250lessMETless350alphaSdown" , "CR2l_D_250lessMETless350alphaSup" , "CR2l_D_250lessMETless350Q2down" , "CR2l_D_250lessMETless350Q2up" , "CR2l_D_250lessMETless350ISRnjetsDown" , "CR2l_D_250lessMETless350ISRnjetsUp" , "CR2l_D_350lessMETless450" , "CR2l_D_350lessMETless450LSFdown" , "CR2l_D_350lessMETless450LSFup" , "CR2l_D_350lessMETless450BTlightDown" , "CR2l_D_350lessMETless450BTlightUp" , "CR2l_D_350lessMETless450BTheavyDown" , "CR2l_D_350lessMETless450BTheavyUp" , "CR2l_D_350lessMETless450PUdown" , "CR2l_D_350lessMETless450PUup" , "CR2l_D_350lessMETless450PDFdown" , "CR2l_D_350lessMETless450PDFup" , "CR2l_D_350lessMETless450alphaSdown" , "CR2l_D_350lessMETless450alphaSup" , "CR2l_D_350lessMETless450Q2down" , "CR2l_D_350lessMETless450Q2up" , "CR2l_D_350lessMETless450ISRnjetsDown" , "CR2l_D_350lessMETless450ISRnjetsUp" , "CR2l_D_450lessMETless550" , "CR2l_D_450lessMETless550LSFdown" , "CR2l_D_450lessMETless550LSFup" , "CR2l_D_450lessMETless550BTlightDown" , "CR2l_D_450lessMETless550BTlightUp" , "CR2l_D_450lessMETless550BTheavyDown" , "CR2l_D_450lessMETless550BTheavyUp" , "CR2l_D_450lessMETless550PUdown" , "CR2l_D_450lessMETless550PUup" , "CR2l_D_450lessMETless550PDFdown" , "CR2l_D_450lessMETless550PDFup" , "CR2l_D_450lessMETless550alphaSdown" , "CR2l_D_450lessMETless550alphaSup" , "CR2l_D_450lessMETless550Q2down" , "CR2l_D_450lessMETless550Q2up" , "CR2l_D_450lessMETless550ISRnjetsDown" , "CR2l_D_450lessMETless550ISRnjetsUp" , "CR2l_D_550lessMETlessInf" , "CR2l_D_550lessMETlessInfLSFdown" , "CR2l_D_550lessMETlessInfLSFup" , "CR2l_D_550lessMETlessInfBTlightDown" , "CR2l_D_550lessMETlessInfBTlightUp" , "CR2l_D_550lessMETlessInfBTheavyDown" , "CR2l_D_550lessMETlessInfBTheavyUp" , "CR2l_D_550lessMETlessInfPUdown" , "CR2l_D_550lessMETlessInfPUup" , "CR2l_D_550lessMETlessInfPDFdown" , "CR2l_D_550lessMETlessInfPDFup" , "CR2l_D_550lessMETlessInfalphaSdown" , "CR2l_D_550lessMETlessInfalphaSup" , "CR2l_D_550lessMETlessInfQ2down" , "CR2l_D_550lessMETlessInfQ2up" , "CR2l_D_550lessMETlessInfISRnjetsDown" , "CR2l_D_550lessMETlessInfISRnjetsUp" , "CR2l_E_250lessMETless350" , "CR2l_E_250lessMETless350LSFdown" , "CR2l_E_250lessMETless350LSFup" , "CR2l_E_250lessMETless350BTlightDown" , "CR2l_E_250lessMETless350BTlightUp" , "CR2l_E_250lessMETless350BTheavyDown" , "CR2l_E_250lessMETless350BTheavyUp" , "CR2l_E_250lessMETless350PUdown" , "CR2l_E_250lessMETless350PUup" , "CR2l_E_250lessMETless350PDFdown" , "CR2l_E_250lessMETless350PDFup" , "CR2l_E_250lessMETless350alphaSdown" , "CR2l_E_250lessMETless350alphaSup" , "CR2l_E_250lessMETless350Q2down" , "CR2l_E_250lessMETless350Q2up" , "CR2l_E_250lessMETless350ISRnjetsDown" , "CR2l_E_250lessMETless350ISRnjetsUp" , "CR2l_E_350lessMETless550" , "CR2l_E_350lessMETless550LSFdown" , "CR2l_E_350lessMETless550LSFup" , "CR2l_E_350lessMETless550BTlightDown" , "CR2l_E_350lessMETless550BTlightUp" , "CR2l_E_350lessMETless550BTheavyDown" , "CR2l_E_350lessMETless550BTheavyUp" , "CR2l_E_350lessMETless550PUdown" , "CR2l_E_350lessMETless550PUup" , "CR2l_E_350lessMETless550PDFdown" , "CR2l_E_350lessMETless550PDFup" , "CR2l_E_350lessMETless550alphaSdown" , "CR2l_E_350lessMETless550alphaSup" , "CR2l_E_350lessMETless550Q2down" , "CR2l_E_350lessMETless550Q2up" , "CR2l_E_350lessMETless550ISRnjetsDown" , "CR2l_E_350lessMETless550ISRnjetsUp" , "CR2l_E_550lessMETlessInf" , "CR2l_E_550lessMETlessInfLSFdown" , "CR2l_E_550lessMETlessInfLSFup" , "CR2l_E_550lessMETlessInfBTlightDown" , "CR2l_E_550lessMETlessInfBTlightUp" , "CR2l_E_550lessMETlessInfBTheavyDown" , "CR2l_E_550lessMETlessInfBTheavyUp" , "CR2l_E_550lessMETlessInfPUdown" , "CR2l_E_550lessMETlessInfPUup" , "CR2l_E_550lessMETlessInfPDFdown" , "CR2l_E_550lessMETlessInfPDFup" , "CR2l_E_550lessMETlessInfalphaSdown" , "CR2l_E_550lessMETlessInfalphaSup" , "CR2l_E_550lessMETlessInfQ2down" , "CR2l_E_550lessMETlessInfQ2up" , "CR2l_E_550lessMETlessInfISRnjetsDown" , "CR2l_E_550lessMETlessInfISRnjetsUp" , "CR2l_F_250lessMETless450" , "CR2l_F_250lessMETless450LSFdown" , "CR2l_F_250lessMETless450LSFup" , "CR2l_F_250lessMETless450BTlightDown" , "CR2l_F_250lessMETless450BTlightUp" , "CR2l_F_250lessMETless450BTheavyDown" , "CR2l_F_250lessMETless450BTheavyUp" , "CR2l_F_250lessMETless450PUdown" , "CR2l_F_250lessMETless450PUup" , "CR2l_F_250lessMETless450PDFdown" , "CR2l_F_250lessMETless450PDFup" , "CR2l_F_250lessMETless450alphaSdown" , "CR2l_F_250lessMETless450alphaSup" , "CR2l_F_250lessMETless450Q2down" , "CR2l_F_250lessMETless450Q2up" , "CR2l_F_250lessMETless450ISRnjetsDown" , "CR2l_F_250lessMETless450ISRnjetsUp" , "CR2l_F_450lessMETlessInf" , "CR2l_F_450lessMETlessInfLSFdown" , "CR2l_F_450lessMETlessInfLSFup" , "CR2l_F_450lessMETlessInfBTlightDown" , "CR2l_F_450lessMETlessInfBTlightUp" , "CR2l_F_450lessMETlessInfBTheavyDown" , "CR2l_F_450lessMETlessInfBTheavyUp" , "CR2l_F_450lessMETlessInfPUdown" , "CR2l_F_450lessMETlessInfPUup" , "CR2l_F_450lessMETlessInfPDFdown" , "CR2l_F_450lessMETlessInfPDFup" , "CR2l_F_450lessMETlessInfalphaSdown" , "CR2l_F_450lessMETlessInfalphaSup" , "CR2l_F_450lessMETlessInfQ2down" , "CR2l_F_450lessMETlessInfQ2up" , "CR2l_F_450lessMETlessInfISRnjetsDown" , "CR2l_F_450lessMETlessInfISRnjetsUp" , "CR2l_G_250lessMETless350" , "CR2l_G_250lessMETless350LSFdown" , "CR2l_G_250lessMETless350LSFup" , "CR2l_G_250lessMETless350BTlightDown" , "CR2l_G_250lessMETless350BTlightUp" , "CR2l_G_250lessMETless350BTheavyDown" , "CR2l_G_250lessMETless350BTheavyUp" , "CR2l_G_250lessMETless350PUdown" , "CR2l_G_250lessMETless350PUup" , "CR2l_G_250lessMETless350PDFdown" , "CR2l_G_250lessMETless350PDFup" , "CR2l_G_250lessMETless350alphaSdown" , "CR2l_G_250lessMETless350alphaSup" , "CR2l_G_250lessMETless350Q2down" , "CR2l_G_250lessMETless350Q2up" , "CR2l_G_250lessMETless350ISRnjetsDown" , "CR2l_G_250lessMETless350ISRnjetsUp" , "CR2l_G_350lessMETless450" , "CR2l_G_350lessMETless450LSFdown" , "CR2l_G_350lessMETless450LSFup" , "CR2l_G_350lessMETless450BTlightDown" , "CR2l_G_350lessMETless450BTlightUp" , "CR2l_G_350lessMETless450BTheavyDown" , "CR2l_G_350lessMETless450BTheavyUp" , "CR2l_G_350lessMETless450PUdown" , "CR2l_G_350lessMETless450PUup" , "CR2l_G_350lessMETless450PDFdown" , "CR2l_G_350lessMETless450PDFup" , "CR2l_G_350lessMETless450alphaSdown" , "CR2l_G_350lessMETless450alphaSup" , "CR2l_G_350lessMETless450Q2down" , "CR2l_G_350lessMETless450Q2up" , "CR2l_G_350lessMETless450ISRnjetsDown" , "CR2l_G_350lessMETless450ISRnjetsUp" , "CR2l_G_450lessMETless600" , "CR2l_G_450lessMETless600LSFdown" , "CR2l_G_450lessMETless600LSFup" , "CR2l_G_450lessMETless600BTlightDown" , "CR2l_G_450lessMETless600BTlightUp" , "CR2l_G_450lessMETless600BTheavyDown" , "CR2l_G_450lessMETless600BTheavyUp" , "CR2l_G_450lessMETless600PUdown" , "CR2l_G_450lessMETless600PUup" , "CR2l_G_450lessMETless600PDFdown" , "CR2l_G_450lessMETless600PDFup" , "CR2l_G_450lessMETless600alphaSdown" , "CR2l_G_450lessMETless600alphaSup" , "CR2l_G_450lessMETless600Q2down" , "CR2l_G_450lessMETless600Q2up" , "CR2l_G_450lessMETless600ISRnjetsDown" , "CR2l_G_450lessMETless600ISRnjetsUp" , "CR2l_G_600lessMETlessInf" , "CR2l_G_600lessMETlessInfLSFdown" , "CR2l_G_600lessMETlessInfLSFup" , "CR2l_G_600lessMETlessInfBTlightDown" , "CR2l_G_600lessMETlessInfBTlightUp" , "CR2l_G_600lessMETlessInfBTheavyDown" , "CR2l_G_600lessMETlessInfBTheavyUp" , "CR2l_G_600lessMETlessInfPUdown" , "CR2l_G_600lessMETlessInfPUup" , "CR2l_G_600lessMETlessInfPDFdown" , "CR2l_G_600lessMETlessInfPDFup" , "CR2l_G_600lessMETlessInfalphaSdown" , "CR2l_G_600lessMETlessInfalphaSup" , "CR2l_G_600lessMETlessInfQ2down" , "CR2l_G_600lessMETlessInfQ2up" , "CR2l_G_600lessMETlessInfISRnjetsDown" , "CR2l_G_600lessMETlessInfISRnjetsUp" , "CR2l_H_250lessMETless450" , "CR2l_H_250lessMETless450LSFdown" , "CR2l_H_250lessMETless450LSFup" , "CR2l_H_250lessMETless450BTlightDown" , "CR2l_H_250lessMETless450BTlightUp" , "CR2l_H_250lessMETless450BTheavyDown" , "CR2l_H_250lessMETless450BTheavyUp" , "CR2l_H_250lessMETless450PUdown" , "CR2l_H_250lessMETless450PUup" , "CR2l_H_250lessMETless450PDFdown" , "CR2l_H_250lessMETless450PDFup" , "CR2l_H_250lessMETless450alphaSdown" , "CR2l_H_250lessMETless450alphaSup" , "CR2l_H_250lessMETless450Q2down" , "CR2l_H_250lessMETless450Q2up" , "CR2l_H_250lessMETless450ISRnjetsDown" , "CR2l_H_250lessMETless450ISRnjetsUp" , "CR2l_H_450lessMETlessInf" , "CR2l_H_450lessMETlessInfLSFdown" , "CR2l_H_450lessMETlessInfLSFup" , "CR2l_H_450lessMETlessInfBTlightDown" , "CR2l_H_450lessMETlessInfBTlightUp" , "CR2l_H_450lessMETlessInfBTheavyDown" , "CR2l_H_450lessMETlessInfBTheavyUp" , "CR2l_H_450lessMETlessInfPUdown" , "CR2l_H_450lessMETlessInfPUup" , "CR2l_H_450lessMETlessInfPDFdown" , "CR2l_H_450lessMETlessInfPDFup" , "CR2l_H_450lessMETlessInfalphaSdown" , "CR2l_H_450lessMETlessInfalphaSup" , "CR2l_H_450lessMETlessInfQ2down" , "CR2l_H_450lessMETlessInfQ2up" , "CR2l_H_450lessMETlessInfISRnjetsDown" , "CR2l_H_450lessMETlessInfISRnjetsUp" , "CR2l_I_250lessMETless350" , "CR2l_I_250lessMETless350LSFdown" , "CR2l_I_250lessMETless350LSFup" , "CR2l_I_250lessMETless350BTlightDown" , "CR2l_I_250lessMETless350BTlightUp" , "CR2l_I_250lessMETless350BTheavyDown" , "CR2l_I_250lessMETless350BTheavyUp" , "CR2l_I_250lessMETless350PUdown" , "CR2l_I_250lessMETless350PUup" , "CR2l_I_250lessMETless350PDFdown" , "CR2l_I_250lessMETless350PDFup" , "CR2l_I_250lessMETless350alphaSdown" , "CR2l_I_250lessMETless350alphaSup" , "CR2l_I_250lessMETless350Q2down" , "CR2l_I_250lessMETless350Q2up" , "CR2l_I_250lessMETless350ISRnjetsDown" , "CR2l_I_250lessMETless350ISRnjetsUp" , "CR2l_I_350lessMETless450" , "CR2l_I_350lessMETless450LSFdown" , "CR2l_I_350lessMETless450LSFup" , "CR2l_I_350lessMETless450BTlightDown" , "CR2l_I_350lessMETless450BTlightUp" , "CR2l_I_350lessMETless450BTheavyDown" , "CR2l_I_350lessMETless450BTheavyUp" , "CR2l_I_350lessMETless450PUdown" , "CR2l_I_350lessMETless450PUup" , "CR2l_I_350lessMETless450PDFdown" , "CR2l_I_350lessMETless450PDFup" , "CR2l_I_350lessMETless450alphaSdown" , "CR2l_I_350lessMETless450alphaSup" , "CR2l_I_350lessMETless450Q2down" , "CR2l_I_350lessMETless450Q2up" , "CR2l_I_350lessMETless450ISRnjetsDown" , "CR2l_I_350lessMETless450ISRnjetsUp" , "CR2l_I_450lessMETless550" , "CR2l_I_450lessMETless550LSFdown" , "CR2l_I_450lessMETless550LSFup" , "CR2l_I_450lessMETless550BTlightDown" , "CR2l_I_450lessMETless550BTlightUp" , "CR2l_I_450lessMETless550BTheavyDown" , "CR2l_I_450lessMETless550BTheavyUp" , "CR2l_I_450lessMETless550PUdown" , "CR2l_I_450lessMETless550PUup" , "CR2l_I_450lessMETless550PDFdown" , "CR2l_I_450lessMETless550PDFup" , "CR2l_I_450lessMETless550alphaSdown" , "CR2l_I_450lessMETless550alphaSup" , "CR2l_I_450lessMETless550Q2down" , "CR2l_I_450lessMETless550Q2up" , "CR2l_I_450lessMETless550ISRnjetsDown" , "CR2l_I_450lessMETless550ISRnjetsUp" , "CR2l_I_550lessMETlessInf" , "CR2l_I_550lessMETlessInfLSFdown" , "CR2l_I_550lessMETlessInfLSFup" , "CR2l_I_550lessMETlessInfBTlightDown" , "CR2l_I_550lessMETlessInfBTlightUp" , "CR2l_I_550lessMETlessInfBTheavyDown" , "CR2l_I_550lessMETlessInfBTheavyUp" , "CR2l_I_550lessMETlessInfPUdown" , "CR2l_I_550lessMETlessInfPUup" , "CR2l_I_550lessMETlessInfPDFdown" , "CR2l_I_550lessMETlessInfPDFup" , "CR2l_I_550lessMETlessInfalphaSdown" , "CR2l_I_550lessMETlessInfalphaSup" , "CR2l_I_550lessMETlessInfQ2down" , "CR2l_I_550lessMETlessInfQ2up" , "CR2l_I_550lessMETlessInfISRnjetsDown" , "CR2l_I_550lessMETlessInfISRnjetsUp"   };

    TableDataMC(this, yieldCR,"lepChannel",  "keepNegative" ).Print(outputName+ "yieldsCRsyst.tab", 4);
    TableDataMC(this, yieldCR,"lepChannel", "keepNegative" ).PrintLatex(outputName+ "yieldsCRsyst.tex", 4);
/*    ofstream sigfile("controlRegMor.txt");
    if (sigfile.is_open())
    {
        for(uint32_t r=0; r<totYield.size(); r++)
        {
            sigfile << totYield.at(r) << endl;
        }
            sigfile.close();
    } */ 

    cout << "end of processing" << endl;
 }


    float getWeight(string currentProcessType, float lumi, float s1fb2)
    {
        float nEvents =  myEvent.wNormalization.at(22);
        float all_weights = lumi*  myEvent.scale1fb * myEvent.weight_PU  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31))* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ;

        //cout << "weight normal " <<all_weights << endl;
        if(s1fb2 != 1)
        {
            float scaleUpdate = (myEvent.scale1fb/(myEvent.scale1fb + s1fb2 )  );
            //cout << "scaleUpdate " << scaleUpdate << "sclae 1fb " << myEvent.scale1fb << "scale 1fb2 " << s1fb2 << endl;
            all_weights *= scaleUpdate;
            //cout << "weight updated " <<all_weights << endl;
        }

        if(currentProcessType == "signal")
        {
            all_weights = lumi; 
        }
            
        return all_weights;
    }

    void getscale1fb2(TString fleName, float* scale1fb2 )
    {
	TFile *fbkg = NULL;
	fbkg = TFile::Open(fleName);
	TTree* tBkg = NULL;
	tBkg =  (TTree*) fbkg->Get("t");
	if(tBkg->GetListOfBranches()->FindObject("scale1fb"))
        {
            cout << "filling nevt " << endl;
	    tBkg->SetBranchAddress("scale1fb",      scale1fb2);
        }
        tBkg->GetEntry(1);
	fbkg->Close();

        cout << "scale1fb2 after filling " << *scale1fb2 << " file " << fleName << endl;
    }
