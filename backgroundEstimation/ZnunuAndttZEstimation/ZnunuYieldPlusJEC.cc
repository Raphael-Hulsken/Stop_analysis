#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

#define USE_TTZ
//#define USE_TTZNLO
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
string storedDataset = "";
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
    
//Add this as a global variable 

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11/NotSkimmed";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11/";
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/haweber_v22/";
    totalNumberOfWorkers = 1;



    int trackveto = (int) myEvent.PassTrackVeto; 
    int tauveto = (int) myEvent.PassTauVeto; 
    AddVariable("MET", "MET",  "MET", 100 ,200,1000,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "MT2W", 100 ,0,500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("MT", "MT",  "MT", 100 ,100,1000,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("nJets","nJets","nJets",5,1,5,&(myEvent.ngoodjets),"noUnderflowInFirstBin");
    AddVariable("nBJets","nBJets","nBJets",5,1,5,&(myEvent.ngoodbtags),"noUnderflowInFirstBin");
    AddVariable("topnessMod","topnessMod","topnessMod",20,-20,20,&(myEvent.topnessMod),"noUnderflowInFirstBin");
    AddVariable("dphi","dphi","dphi", 100,0,3.5,&(myEvent.dphi_ak4pfjets_met),"noUnderflowInFirstBin");
    AddVariable("nvertex","nvertex","nvertex",50,0,50,&(myEvent.nvertex),"noUnderflowInFirstBin");
    AddVariable("trackveto","trackveto","trackveto",3,0,3,&(trackveto),"noUnderflowInFirstBin");
    AddVariable("tauveto","tauveto","tauveto",3,0,3,&(tauveto),"noUnderflowInFirstBin");


    // ------------------
    // Datasets
    // ------------------
    #ifdef USE_TTZ
    AddProcessClass("ttZ", "ttZ", "background", kBlue);
    	AddDataset("ttZJets_13TeV_madgraphMLM","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_1","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_2","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_3","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_4","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_5","ttZ",0,0);
    	//AddDataset("ttZJets_13TeV_madgraphMLM_6","ttZ",0,0);
    outputName = "yieldZnunuMorTTZ";
    #endif

    #ifdef USE_TTZNLO
    AddProcessClass("ttZNLO", "ttZNLO", "background", kBlue);
    	AddDataset("TTZToLLNuNu_M-10_amcnlo_pythia8_25ns","ttZNLO",0,0);
    outputName = "yieldZnunuMorTTZNLO";
    #endif


    //AddProcessClass("ZZ", "ZZ", "background", kBlue);
    	//AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZ",0,0);
    AddProcessClass("WZ", "WZ", "background", kBlue);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","WZ",0,0);
    //AddProcessClass("tZq", "tZq", "background", kBlue);
    	//AddDataset("tZq_ll_4f_amcatnlo-pythia8_25ns","tZq",0,0);

AddRegion("SR1l_A_250lessMETless350","SR1l_A_250lessMETless350",&SR1l_A_250lessMETless350);
AddRegion("SR1l_A_350lessMETless450","SR1l_A_350lessMETless450",&SR1l_A_350lessMETless450);
AddRegion("SR1l_A_450lessMETless600","SR1l_A_450lessMETless600",&SR1l_A_450lessMETless600);
AddRegion("SR1l_A_600lessMETlessInf","SR1l_A_600lessMETlessInf",&SR1l_A_600lessMETlessInf);
AddRegion("SR1l_B_250lessMETless450","SR1l_B_250lessMETless450",&SR1l_B_250lessMETless450);
AddRegion("SR1l_B_450lessMETless600","SR1l_B_450lessMETless600",&SR1l_B_450lessMETless600);
AddRegion("SR1l_B_600lessMETlessInf","SR1l_B_600lessMETlessInf",&SR1l_B_600lessMETlessInf);
AddRegion("SR1l_C_250lessMETless350","SR1l_C_250lessMETless350",&SR1l_C_250lessMETless350);
AddRegion("SR1l_C_350lessMETless450","SR1l_C_350lessMETless450",&SR1l_C_350lessMETless450);
AddRegion("SR1l_C_450lessMETless550","SR1l_C_450lessMETless550",&SR1l_C_450lessMETless550);
AddRegion("SR1l_C_550lessMETless650","SR1l_C_550lessMETless650",&SR1l_C_550lessMETless650);
AddRegion("SR1l_C_650lessMETlessInf","SR1l_C_650lessMETlessInf",&SR1l_C_650lessMETlessInf);
AddRegion("SR1l_D_250lessMETless350","SR1l_D_250lessMETless350",&SR1l_D_250lessMETless350);
AddRegion("SR1l_D_350lessMETless450","SR1l_D_350lessMETless450",&SR1l_D_350lessMETless450);
AddRegion("SR1l_D_450lessMETless550","SR1l_D_450lessMETless550",&SR1l_D_450lessMETless550);
AddRegion("SR1l_D_550lessMETlessInf","SR1l_D_550lessMETlessInf",&SR1l_D_550lessMETlessInf);
AddRegion("SR1l_E_250lessMETless350","SR1l_E_250lessMETless350",&SR1l_E_250lessMETless350);
AddRegion("SR1l_E_350lessMETless550","SR1l_E_350lessMETless550",&SR1l_E_350lessMETless550);
AddRegion("SR1l_E_550lessMETlessInf","SR1l_E_550lessMETlessInf",&SR1l_E_550lessMETlessInf);
AddRegion("SR1l_F_250lessMETless450","SR1l_F_250lessMETless450",&SR1l_F_250lessMETless450);
AddRegion("SR1l_F_450lessMETlessInf","SR1l_F_450lessMETlessInf",&SR1l_F_450lessMETlessInf);
AddRegion("SR1l_G_250lessMETless350","SR1l_G_250lessMETless350",&SR1l_G_250lessMETless350);
AddRegion("SR1l_G_350lessMETless450","SR1l_G_350lessMETless450",&SR1l_G_350lessMETless450);
AddRegion("SR1l_G_450lessMETless600","SR1l_G_450lessMETless600",&SR1l_G_450lessMETless600);
AddRegion("SR1l_G_600lessMETlessInf","SR1l_G_600lessMETlessInf",&SR1l_G_600lessMETlessInf);
AddRegion("SR1l_H_250lessMETless450","SR1l_H_250lessMETless450",&SR1l_H_250lessMETless450);
AddRegion("SR1l_H_450lessMETlessInf","SR1l_H_450lessMETlessInf",&SR1l_H_450lessMETlessInf);
AddRegion("SR1l_I_250lessMETless350","SR1l_I_250lessMETless350",&SR1l_I_250lessMETless350);
AddRegion("SR1l_I_350lessMETless450","SR1l_I_350lessMETless450",&SR1l_I_350lessMETless450);
AddRegion("SR1l_I_450lessMETless550","SR1l_I_450lessMETless550",&SR1l_I_450lessMETless550);
AddRegion("SR1l_I_550lessMETlessInf","SR1l_I_550lessMETlessInf",&SR1l_I_550lessMETlessInf);


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


    checkNegativeYields = false;
    if(nthentry == myEvent.nentries)
    {
        cout << "checking histogram for negative values" << endl;
        nthentry =0;
        checkNegativeYields = true;
    }


    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);


    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);


    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
    if( (currentProcessClass == "Znunu" || currentProcessClass =="ZZ" || currentProcessClass == "ttZ" || currentProcessClass == "WZ" ||  currentProcessClass == "ttZNLO"||  currentProcessClass == "tZq")  && ( !(myEvent.isZtoNuNu) ))
    {
         currentProcessClass = "";
    }
    
    float nEventsN =  myEvent.wNormalization.at(22);
    float ISRNJ = 1;
    if( currentProcessClass == "ttZ")
    {
        ISRNJ = myEvent.weight_ISRnjets*( nEventsN / myEvent.wNormalization.at(25));
        //ISRNJ=1;
    }
    cout <<"class " << currentProcessClass  << " ngoodjets" << myEvent.ngoodjets << "ISR njets weight "  << myEvent.weight_ISRnjets*( nEventsN / myEvent.wNormalization.at(25)) << endl;
    float weightLumi = getWeight(currentProcessType, GetLumi())*ISRNJ; 
    float weight     = weightLumi;

    if (currentProcessType == "data") weight = 1.0;

    AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields);//, dummy, false);



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
vector<string> totYield = { "SR1l_A_250lessMETless350" , "SR1l_A_350lessMETless450" , "SR1l_A_450lessMETless600" , "SR1l_A_600lessMETlessInf" , "SR1l_B_250lessMETless450" , "SR1l_B_450lessMETless600" , "SR1l_B_600lessMETlessInf" , "SR1l_C_250lessMETless350" , "SR1l_C_350lessMETless450" , "SR1l_C_450lessMETless550" , "SR1l_C_550lessMETless650" , "SR1l_C_650lessMETlessInf" , "SR1l_D_250lessMETless350" , "SR1l_D_350lessMETless450" , "SR1l_D_450lessMETless550" , "SR1l_D_550lessMETlessInf" , "SR1l_E_250lessMETless350" , "SR1l_E_350lessMETless550" , "SR1l_E_550lessMETlessInf" , "SR1l_F_250lessMETless450" , "SR1l_F_450lessMETlessInf" , "SR1l_G_250lessMETless350" , "SR1l_G_350lessMETless450" , "SR1l_G_450lessMETless600" , "SR1l_G_600lessMETlessInf" , "SR1l_H_250lessMETless450" , "SR1l_H_450lessMETlessInf", "SR1l_I_250lessMETless350" , "SR1l_I_350lessMETless450" , "SR1l_I_450lessMETless550" , "SR1l_I_550lessMETlessInf"};


    #ifdef USE_VAR_BASELINE_DOWN
    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+ "JECDown.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+ "JECdown.tex", 4);
    #endif
    #ifdef USE_VAR_BASELINE_UP
    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+ "JECUp.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+ "JECUp.tex", 4);
    #endif
    #ifdef USE_VAR_BASELINE
    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+ "SR.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+ "SR.tex", 4);
    #endif

    cout << "end of processing" << endl;
 }


    float getWeight(string currentProcessType, float lumi)
    {
        float nEvents =  myEvent.wNormalization.at(22);
        float all_weights = lumi*  myEvent.scale1fb /** myEvent.weight_PU*/  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31))* myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ;

        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
            //all_weights = lumi * myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) * myEvent.weight_ISR*( nEvents / myEvent.wNormalization.at(19) ) ;
            
        return all_weights;
    }

