#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

#define USE_TTZ
//#define USE_WZ
//#define USE_ZZ

#define USE_VAR_BASELINE
#define USE_LEP1
#define USE_LEP2
#define USE_JETS
#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_GLOBAL_VAR

//#include "../../common/common.h"
#include "../../common/TFFactory.h"
#include "../../Selection/moriondPU.h"
#define _METCUT_ 50
#define _LEPTONPTCUT_ 40

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
map< pair<uint32_t,uint32_t>, string > scanMap;

bool lepChannel() 
{ 
    return true; 
}
    
ofstream statnames("statNamesPU.txt");
//Add this as a global variable 

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11/NotSkimmed";
    totalNumberOfWorkers = 1;



    #ifdef USE_WZ
    AddProcessClass("WZ", "WZ", "background", kBlue);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","WZ",0,0);
    outputName = "yieldMorWZ";
    #endif

    #ifdef USE_ZZ
    AddProcessClass("ZZ", "ZZ", "background", kBlue);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZ",0,0);
    outputName = "yieldMorWZ";
    #endif

    #ifdef USE_TTZ
    AddProcessClass("ttZ", "ttZ", "background", kBlue);
        AddDataset("ttZJets_13TeV_madgraphMLM","ttZ",0,0);
        AddDataset("ttZJets_13TeV_madgraphMLM_1","ttZ",0,0);
        AddDataset("ttZJets_13TeV_madgraphMLM_2","ttZ",0,0);
        AddDataset("ttZJets_13TeV_madgraphMLM_3","ttZ",0,0);
        AddDataset("ttZJets_13TeV_madgraphMLM_4","ttZ",0,0);
        AddDataset("ttZJets_13TeV_madgraphMLM_5","ttZ",0,0);
        AddDataset("ttZJets_13TeV_madgraphMLM_6","ttZ",0,0);
    outputName = "yieldMorttZ";
    #endif



//regions //@MJ@ TODO add I
AddRegion("SR1l_A_250lessMETlessInf","SR1l_A_250lessMETlessInf",&SR1l_A_250lessMETlessInf);
AddRegion("SR1l_A_250lessMETlessInfPUdown","SR1l_A_250lessMETlessInfPUdown",&SR1l_A_250lessMETlessInf);
AddRegion("SR1l_A_250lessMETlessInfPUup","SR1l_A_250lessMETlessInfPUup",&SR1l_A_250lessMETlessInf);
AddRegion("SR1l_B_250lessMETlessInf","SR1l_B_250lessMETlessInf",&SR1l_B_250lessMETlessInf);
AddRegion("SR1l_B_250lessMETlessInfPUdown","SR1l_B_250lessMETlessInfPUdown",&SR1l_B_250lessMETlessInf);
AddRegion("SR1l_B_250lessMETlessInfPUup","SR1l_B_250lessMETlessInfPUup",&SR1l_B_250lessMETlessInf);
AddRegion("SR1l_C_250lessMETlessInf","SR1l_C_250lessMETlessInf",&SR1l_C_250lessMETlessInf);
AddRegion("SR1l_C_250lessMETlessInfPUdown","SR1l_C_250lessMETlessInfPUdown",&SR1l_C_250lessMETlessInf);
AddRegion("SR1l_C_250lessMETlessInfPUup","SR1l_C_250lessMETlessInfPUup",&SR1l_C_250lessMETlessInf);
AddRegion("SR1l_D_250lessMETlessInf","SR1l_D_250lessMETlessInf",&SR1l_D_250lessMETlessInf);
AddRegion("SR1l_D_250lessMETlessInfPUdown","SR1l_D_250lessMETlessInfPUdown",&SR1l_D_250lessMETlessInf);
AddRegion("SR1l_D_250lessMETlessInfPUup","SR1l_D_250lessMETlessInfPUup",&SR1l_D_250lessMETlessInf);
AddRegion("SR1l_E_250lessMETlessInf","SR1l_E_250lessMETlessInf",&SR1l_E_250lessMETlessInf);
AddRegion("SR1l_E_250lessMETlessInfPUdown","SR1l_E_250lessMETlessInfPUdown",&SR1l_E_250lessMETlessInf);
AddRegion("SR1l_E_250lessMETlessInfPUup","SR1l_E_250lessMETlessInfPUup",&SR1l_E_250lessMETlessInf);
AddRegion("SR1l_F_250lessMETlessInf","SR1l_F_250lessMETlessInf",&SR1l_F_250lessMETlessInf);
AddRegion("SR1l_F_250lessMETlessInfPUdown","SR1l_F_250lessMETlessInfPUdown",&SR1l_F_250lessMETlessInf);
AddRegion("SR1l_F_250lessMETlessInfPUup","SR1l_F_250lessMETlessInfPUup",&SR1l_F_250lessMETlessInf);
AddRegion("SR1l_G_250lessMETlessInf","SR1l_G_250lessMETlessInf",&SR1l_G_250lessMETlessInf);
AddRegion("SR1l_G_250lessMETlessInfPUdown","SR1l_G_250lessMETlessInfPUdown",&SR1l_G_250lessMETlessInf);
AddRegion("SR1l_G_250lessMETlessInfPUup","SR1l_G_250lessMETlessInfPUup",&SR1l_G_250lessMETlessInf);
AddRegion("SR1l_H_250lessMETlessInf","SR1l_H_250lessMETlessInf",&SR1l_H_250lessMETlessInf);
AddRegion("SR1l_H_250lessMETlessInfPUdown","SR1l_H_250lessMETlessInfPUdown",&SR1l_H_250lessMETlessInf);
AddRegion("SR1l_H_250lessMETlessInfPUup","SR1l_H_250lessMETlessInfPUup",&SR1l_H_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInf","SR1l_I_250lessMETlessInf",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfPUdown","SR1l_I_250lessMETlessInfPUdown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfPUup","SR1l_I_250lessMETlessInfPUup",&SR1l_I_250lessMETlessInf);

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(35.867);

fillYieldsVector();
    Create1DHistos();

    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);
    bool useTriggerInfo = currentProcessType == "data" ? true: false; //@MJ@ TODO fix the trigger

    nthentry++;


    checkNegativeYields = false;
    if(nthentry == myEvent.nentries)
    {
        cout << " checking yields for negative values " << endl;
        nthentry =0;
        checkNegativeYields = true;
    }

    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);

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
    for(uint32_t SR=0; SR<8; SR++) //@MJ@ TODO nr of sig regions changes
    {

        float w = 0;
        float btagmax = 0;
        //normal
        weightV.push_back(weightLumi);
        //PUdown
        if(counter == 1) statnames << "PUdown" << endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )*myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )* myEvent.weight_PUdown; //@MJ@ TODO PU without any normalization?!
        weightV.push_back(w);
        //PUup
        if(counter == 1) statnames << "PUup"<< endl;
        w = GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )*myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )* myEvent.weight_PUup;
        weightV.push_back(w);
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

vector<string> totYield = { "SR1l_A_250lessMETlessInf" , "SR1l_A_250lessMETlessInfPUdown" , "SR1l_A_250lessMETlessInfPUup" , "SR1l_B_250lessMETlessInf" , "SR1l_B_250lessMETlessInfPUdown" , "SR1l_B_250lessMETlessInfPUup" , "SR1l_C_250lessMETlessInf" , "SR1l_C_250lessMETlessInfPUdown" , "SR1l_C_250lessMETlessInfPUup" , "SR1l_D_250lessMETlessInf" , "SR1l_D_250lessMETlessInfPUdown" , "SR1l_D_250lessMETlessInfPUup" , "SR1l_E_250lessMETlessInf" , "SR1l_E_250lessMETlessInfPUdown" , "SR1l_E_250lessMETlessInfPUup" , "SR1l_F_250lessMETlessInf" , "SR1l_F_250lessMETlessInfPUdown" , "SR1l_F_250lessMETlessInfPUup" , "SR1l_G_250lessMETlessInf" , "SR1l_G_250lessMETlessInfPUdown" , "SR1l_G_250lessMETlessInfPUup" , "SR1l_H_250lessMETlessInf" , "SR1l_H_250lessMETlessInfPUdown" , "SR1l_H_250lessMETlessInfPUup", "SR1l_I_250lessMETlessInf" , "SR1l_I_250lessMETlessInfPUdown" , "SR1l_I_250lessMETlessInfPUup" };


    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+"PU.tab", 6);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+"yieldMorPU.tex", 6);

    ofstream sigfile("signalRegMorPU.txt");
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
        float all_weights = lumi*  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) )* myEvent.weight_PU;
        //cout << "cs: " << myEvent.crossSection << "scale 1 fb " << myEvent.scale1fb << " weight total: " << all_weights << endl;
        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
        return all_weights;
    }


