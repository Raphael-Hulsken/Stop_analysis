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

using namespace std;

#include "signalCS.h"

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

    string distingusihClassBkg(string currentProcessClass, vector<string> classLabels); // class, type, list of labels, ...
    void linkTheNameWithProcessClass(string PC);
    void distingusihcClassSig( string currentprocessclass, string currentprocesstype, vector<string> classlabels); // class, type, list of labels, ...
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
    //AddProcessClass("rare", "rare", "background", kBlue);//@MJ@ TODO K-factor?
    //AddProcessClass("Znunu", "Znunu", "background", kBlue);//@MJ@ TODO K-factor?
    AddProcessClass("ttZ", "ttZ", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("ttZJets_13TeV_madgraphMLM","ttZ",0,0);
    AddProcessClass("ZZTo2Q2nu", "ZZTo2Q2nu", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZTo2Q2nu",0,0);
    AddProcessClass("WZ", "WZ", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","WZ",0,0);


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

fillYieldsVector();
                                                                                                                       

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(36.46);
    //SetLumi(18.1);

    Create1DHistos();
    //Add2DHisto("nJets","MET");

    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;
    nthentry++;

    cout << "nentries " << myEvent.nentries << " nthentry "<< nthentry << endl;

    checkNegativeYields = false;
    if(nthentry == myEvent.nentries)
    {
        nthentry =0;
        checkNegativeYields = true;
    }
       

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);


    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);


    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
    if( (currentProcessClass == "Znunu" || currentProcessClass =="ZZTo2Q2nu" || currentProcessClass == "ttZ" || currentProcessClass == "WZ")  && ( !(myEvent.isZtoNuNu) ))
    {
         cout << "not total for " << currentProcessClass << endl;
         currentProcessClass = "";
    }
  

    float weightLumi = getWeight(currentProcessType, GetLumi()); //@MJ@ TODO cross section form file?!
    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    //if (currentProcessType == "signal") weight = 1;
    vector<float> dummy;
    AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields);//, dummy, false);


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
vector<string> totYield = { "SR1l_A_250lessMETless350" , "SR1l_A_350lessMETless450" , "SR1l_A_450lessMETless600" , "SR1l_A_600lessMETlessInf" , "SR1l_B_250lessMETless450" , "SR1l_B_450lessMETless600" , "SR1l_B_600lessMETlessInf" , "SR1l_C_250lessMETless350" , "SR1l_C_350lessMETless450" , "SR1l_C_450lessMETless550" , "SR1l_C_550lessMETless650" , "SR1l_C_650lessMETlessInf" , "SR1l_D_250lessMETless350" , "SR1l_D_350lessMETless450" , "SR1l_D_450lessMETless550" , "SR1l_D_550lessMETlessInf" , "SR1l_E_250lessMETless350" , "SR1l_E_350lessMETless550" , "SR1l_E_550lessMETlessInf" , "SR1l_F_250lessMETless450" , "SR1l_F_450lessMETlessInf" , "SR1l_G_250lessMETless350" , "SR1l_G_350lessMETless450" , "SR1l_G_450lessMETless600" , "SR1l_G_600lessMETlessInf" , "SR1l_H_250lessMETless450" , "SR1l_H_450lessMETlessInf"};


    #ifdef USE_VAR_BASELINE_DOWN
    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print("yieldJECDownMor.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex("yieldJECdownMor.tex", 4);
    #endif
    #ifdef USE_VAR_BASELINE_UP
    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print("yieldJECUpMor.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex("yieldJECUpMor.tex", 4);
    #endif
    #ifdef USE_VAR_BASELINE
    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print("yieldSRMor.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex("yieldSRMor.tex", 4);
    #endif



    vector <string> tfreg{};


    TFProducer prod(tfreg, "lostLepton", "CR2l");
    prod.produceTFTable("yieldICHEP.tab", "lostLeptonTF");
    TFProducer prod2(tfreg, "singleLepton", "CR1l");
    prod2.produceTFTable("yieldICHEP.tab", "singleLeptonTF");
    TFProducer prod3(tfreg, "singleLeptonFromT", "CR1l");
    prod3.produceTFTable("yieldICHEP.tab", "singleLeptonFromTTF");

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
    //cout << "gen weight: " << myEvent.genweights->at(0) << " gen weight tot " << myEvent.wNormalization.at(1) << " nEvents " << nEvents << " total weight " << myEvent.genweights->at(0)*(nEvents/ myEvent.wNormalization.at(1)) << endl;
        //float all_weights = lumi*  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) ;
        float all_weights = lumi*  myEvent.scale1fb  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        //cout << "cs: " << myEvent.crossSection << "scale 1 fb " << myEvent.scale1fb << " weight total: " << all_weights << endl;
        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
            //all_weights = lumi * myEvent.scale1fb * myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) * myEvent.weight_ISR*( nEvents / myEvent.wNormalization.at(19) ) ;
        return all_weights;
    }

