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
#include "../../Selection/moriondCRs.h"

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

TFile *fileX = new TFile("../../common/xsec_stop_13TeV.root");
TH1D* stopXSEC = (TH1D*)fileX->Get("stop")->Clone();

bool lepChannel() 
{ 
    return true; 
}
    
//Add this as a global variable 

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/data/data6/cms/mjansova/Stop1lSharedBabies/v22/skim/";
    totalNumberOfWorkers = 1;

    AddVariable("MET", "MET",  "MET", 100 ,200,1000,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "MT2W", 100 ,0,500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("MT", "MT",  "MT", 100 ,100,1000,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("nJets","nJets","nJets",10,1,10,&(myEvent.ngoodjets),"noUnderflowInFirstBin");
    AddVariable("nBJets","nBJets","nBJets",5,1,5,&(myEvent.ngoodbtags),"noUnderflowInFirstBin");
    AddVariable("topnessMod","topnessMod","topnessMod",20,-20,20,&(myEvent.topnessMod),"noUnderflowInFirstBin");
    AddVariable("dphi","dphi","dphi", 100,0,3.5,&(myEvent.dphi_ak4pfjets_met),"noUnderflowInFirstBin");
    AddVariable("Mlb","Mlb","Mlb", 100,0, 500,&(myEvent.Mlb),"noUnderflowInFirstBin");

    // ------------------
    // Datasets
    // ------------------
    AddProcessClass("data1", "data1", "data", kBlack); ///@MJ@ TODO discared this!!
        AddDataset("data_double_eg_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1");
        AddDataset("data_double_eg_Run2016C_MINIAOD_03Feb2017-v1","data1");
        AddDataset("data_double_eg_Run2016D_MINIAOD_03Feb2017-v1","data1");
        AddDataset("data_double_eg_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_double_eg_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_double_eg_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_double_eg_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);
        AddDataset("data_double_mu_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_double_mu_Run2016C_MINIAOD_03Feb2017-v1","data1");
        AddDataset("data_double_mu_Run2016D_MINIAOD_03Feb2017-v1","data1");
        AddDataset("data_double_mu_Run2016E_MINIAOD_03Feb2017-v1","data1");
        AddDataset("data_double_mu_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_double_mu_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_double_mu_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);
        AddDataset("data_muon_eg_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_muon_eg_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_muon_eg_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_muon_eg_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_muon_eg_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_muon_eg_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_muon_eg_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);
        AddDataset("data_met_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_met_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_single_muon_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_single_electron_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);

    AddProcessClass("bkgOneLepFromW", "bkgOneLepFromW", "background", kBlue); ///@MJ@ TODO discared this!!
    	AddDataset("W1JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0); //@MJ@ TODO which W+jets - check Indara's slides?
    	AddDataset("W2JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
    	AddDataset("W3JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
    	AddDataset("W4JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
    AddProcessClass("bkgOneLepFromTop", "bkgOneLepFromTop", "background", kBlue); ///@MJ@ TODO discared this!!
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); //@MJ@ TODO what about extentions, good normalization to total number of events?
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); 
    AddProcessClass("bkgLostLepton", "bkgLostLepton", "background", kBlue); 
    	AddDataset("ttbar_diLept_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ttbar_diLept_madgraph_pythia8_ext1_25ns","bkgLostLepton",0,0);
    	AddDataset("t_sch_4f_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("t_tW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
    	AddDataset("t_tbarW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
    	AddDataset("t_tch_4f_powheg_pythia8_inclDecays_25ns","bkgLostLepton",0,0); //@MJ@ TODO no atop?!
    	AddDataset("ttWJets_13TeV_madgraphMLM","bkgLostLepton",0,0);
    	AddDataset("WWTo2l2Nu_powheg_25ns","bkgLostLepton",0,0);
    	AddDataset("WWToLNuQQ_powheg_25ns","bkgLostLepton",0,0);
    AddProcessClass("bkgZnunu", "bkgZnunu", "background", kBlue);
    	AddDataset("ttZJets_13TeV_madgraphMLM","bkgZnunu",0,0);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","bkgZnunu",0,0);
    	AddDataset("WZTo1LNu2Q_amcnlo_pythia8_25ns","bkgZnunu",0,0);
    	AddDataset("WZTo3LNu_powheg_pythia8_25ns","bkgZnunu",0,0);//@MJ@ TODO 3l missing
    	AddDataset("ZZTo2L2Nu_powheg_pythia8_25ns","bkgZnunu",0,0);
    	AddDataset("ZZTo2L2Q_amcnlo_pythia8_25ns","bkgZnunu",0,0);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","bkgZnunu",0,0);
    	AddDataset("ZZTo4L_powheg_pythia8_25ns","bkgZnunu",0,0);
    
AddRegion("CR2l_A_250lessMETless350","CR2l_A_250lessMETless350",&CR2l_A_250lessMETless350);
AddRegion("CR2l_A_350lessMETless450","CR2l_A_350lessMETless450",&CR2l_A_350lessMETless450);
AddRegion("CR2l_A_450lessMETless600","CR2l_A_450lessMETless600",&CR2l_A_450lessMETless600);
AddRegion("CR2l_A_600lessMETlessInf","CR2l_A_600lessMETlessInf",&CR2l_A_600lessMETlessInf);
AddRegion("CR2l_B_250lessMETless450","CR2l_B_250lessMETless450",&CR2l_B_250lessMETless450);
AddRegion("CR2l_B_450lessMETless600","CR2l_B_450lessMETless600",&CR2l_B_450lessMETless600);
AddRegion("CR2l_B_600lessMETlessInf","CR2l_B_600lessMETlessInf",&CR2l_B_600lessMETlessInf);
AddRegion("CR2l_C_250lessMETless350","CR2l_C_250lessMETless350",&CR2l_C_250lessMETless350);
AddRegion("CR2l_C_350lessMETless450","CR2l_C_350lessMETless450",&CR2l_C_350lessMETless450);
AddRegion("CR2l_C_450lessMETless550","CR2l_C_450lessMETless550",&CR2l_C_450lessMETless550);
AddRegion("CR2l_C_550lessMETless650","CR2l_C_550lessMETless650",&CR2l_C_550lessMETless650);
AddRegion("CR2l_C_650lessMETlessInf","CR2l_C_650lessMETlessInf",&CR2l_C_650lessMETlessInf);
AddRegion("CR2l_D_250lessMETless350","CR2l_D_250lessMETless350",&CR2l_D_250lessMETless350);
AddRegion("CR2l_D_350lessMETless450","CR2l_D_350lessMETless450",&CR2l_D_350lessMETless450);
AddRegion("CR2l_D_450lessMETless550","CR2l_D_450lessMETless550",&CR2l_D_450lessMETless550);
AddRegion("CR2l_D_550lessMETlessInf","CR2l_D_550lessMETlessInf",&CR2l_D_550lessMETlessInf);
AddRegion("CR2l_E_250lessMETless350","CR2l_E_250lessMETless350",&CR2l_E_250lessMETless350);
AddRegion("CR2l_E_350lessMETless550","CR2l_E_350lessMETless550",&CR2l_E_350lessMETless550);
AddRegion("CR2l_E_550lessMETlessInf","CR2l_E_550lessMETlessInf",&CR2l_E_550lessMETlessInf);
AddRegion("CR2l_F_250lessMETless450","CR2l_F_250lessMETless450",&CR2l_F_250lessMETless450);
AddRegion("CR2l_F_450lessMETlessInf","CR2l_F_450lessMETlessInf",&CR2l_F_450lessMETlessInf);
AddRegion("CR2l_G_250lessMETless350","CR2l_G_250lessMETless350",&CR2l_G_250lessMETless350);
AddRegion("CR2l_G_350lessMETless450","CR2l_G_350lessMETless450",&CR2l_G_350lessMETless450);
AddRegion("CR2l_G_450lessMETless600","CR2l_G_450lessMETless600",&CR2l_G_450lessMETless600);
AddRegion("CR2l_G_600lessMETlessInf","CR2l_G_600lessMETlessInf",&CR2l_G_600lessMETlessInf);
AddRegion("CR2l_H_250lessMETless450","CR2l_H_250lessMETless450",&CR2l_H_250lessMETless450);
AddRegion("CR2l_H_450lessMETlessInf","CR2l_H_450lessMETlessInf",&CR2l_H_450lessMETlessInf);
AddRegion("CR2l_I_250lessMETless350","CR2l_I_250lessMETless350",&CR2l_I_250lessMETless350);
AddRegion("CR2l_I_350lessMETless450","CR2l_I_350lessMETless450",&CR2l_I_350lessMETless450);
AddRegion("CR2l_I_450lessMETless550","CR2l_I_450lessMETless550",&CR2l_I_450lessMETless550);
AddRegion("CR2l_I_550lessMETlessInf","CR2l_I_550lessMETlessInf",&CR2l_I_550lessMETlessInf);

AddRegion("CR0b_A_250lessMETless350","CR0b_A_250lessMETless350",&CR0b_A_250lessMETless350);
AddRegion("CR0b_A_350lessMETless450","CR0b_A_350lessMETless450",&CR0b_A_350lessMETless450);
AddRegion("CR0b_A_450lessMETless600","CR0b_A_450lessMETless600",&CR0b_A_450lessMETless600);
AddRegion("CR0b_A_600lessMETlessInf","CR0b_A_600lessMETlessInf",&CR0b_A_600lessMETlessInf);
AddRegion("CR0b_B_250lessMETless450","CR0b_B_250lessMETless450",&CR0b_B_250lessMETless450);
AddRegion("CR0b_B_450lessMETless600","CR0b_B_450lessMETless600",&CR0b_B_450lessMETless600);
AddRegion("CR0b_B_600lessMETlessInf","CR0b_B_600lessMETlessInf",&CR0b_B_600lessMETlessInf);
AddRegion("CR0b_C_250lessMETless350","CR0b_C_250lessMETless350",&CR0b_C_250lessMETless350);
AddRegion("CR0b_C_350lessMETless450","CR0b_C_350lessMETless450",&CR0b_C_350lessMETless450);
AddRegion("CR0b_C_450lessMETless550","CR0b_C_450lessMETless550",&CR0b_C_450lessMETless550);
AddRegion("CR0b_C_550lessMETless650","CR0b_C_550lessMETless650",&CR0b_C_550lessMETless650);
AddRegion("CR0b_C_650lessMETlessInf","CR0b_C_650lessMETlessInf",&CR0b_C_650lessMETlessInf);
AddRegion("CR0b_D_250lessMETless350","CR0b_D_250lessMETless350",&CR0b_D_250lessMETless350);
AddRegion("CR0b_D_350lessMETless450","CR0b_D_350lessMETless450",&CR0b_D_350lessMETless450);
AddRegion("CR0b_D_450lessMETless550","CR0b_D_450lessMETless550",&CR0b_D_450lessMETless550);
AddRegion("CR0b_D_550lessMETlessInf","CR0b_D_550lessMETlessInf",&CR0b_D_550lessMETlessInf);
AddRegion("CR0b_E_250lessMETless350","CR0b_E_250lessMETless350",&CR0b_E_250lessMETless350);
AddRegion("CR0b_E_350lessMETless550","CR0b_E_350lessMETless550",&CR0b_E_350lessMETless550);
AddRegion("CR0b_E_550lessMETlessInf","CR0b_E_550lessMETlessInf",&CR0b_E_550lessMETlessInf);
AddRegion("CR0b_F_250lessMETless450","CR0b_F_250lessMETless450",&CR0b_F_250lessMETless450);
AddRegion("CR0b_F_450lessMETlessInf","CR0b_F_450lessMETlessInf",&CR0b_F_450lessMETlessInf);
AddRegion("CR0b_G_250lessMETless350","CR0b_G_250lessMETless350",&CR0b_G_250lessMETless350);
AddRegion("CR0b_G_350lessMETless450","CR0b_G_350lessMETless450",&CR0b_G_350lessMETless450);
AddRegion("CR0b_G_450lessMETless600","CR0b_G_450lessMETless600",&CR0b_G_450lessMETless600);
AddRegion("CR0b_G_600lessMETlessInf","CR0b_G_600lessMETlessInf",&CR0b_G_600lessMETlessInf);
AddRegion("CR0b_H_250lessMETless450","CR0b_H_250lessMETless450",&CR0b_H_250lessMETless450);
AddRegion("CR0b_H_450lessMETlessInf","CR0b_H_450lessMETlessInf",&CR0b_H_450lessMETlessInf);
AddRegion("CR0b_I_250lessMETless350","CR0b_I_250lessMETless350",&CR0b_I_250lessMETless350);
AddRegion("CR0b_I_350lessMETless450","CR0b_I_350lessMETless450",&CR0b_I_350lessMETless450);
AddRegion("CR0b_I_450lessMETless550","CR0b_I_450lessMETless550",&CR0b_I_450lessMETless550);
AddRegion("CR0b_I_550lessMETlessInf","CR0b_I_550lessMETlessInf",&CR0b_I_550lessMETlessInf);


    //fillYieldsVector(); @MJ@ TODO probably not needed when I do not need to zero negative
                                                                                                                       

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
        //checkNegativeYields = true; //@MJ@ TODO be aware of this, I can use multiprocessing now but I can endup with incorrect MC yields!!!
    }


    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);


    //vector<string> classLabels;
    //GetProcessClassLabelList(&classLabels);
    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset); //@MJ@ TODO check this!!!

    if( (currentProcessClass == "bkgZnunu")  && ( !(myEvent.isZtoNuNu) ))
    {
         if(myEvent.is2lep)
             currentProcessClass = "bkgLostLepton";
         else if(myEvent.is1lepFromTop)
             currentProcessClass = "bkgOneLepFromTop";
         else if(myEvent.is1lepFromW)
             currentProcessClass = "bkgOneLepFromW";
         else
             cout<< "there are events which are not lost lepton and or 1lep from top/W" << endl;
    }
    if( (currentProcessClass == "bkgLostLepton")  && ( !(myEvent.is2lep) ))
    {
         if(myEvent.is1lepFromTop)
             currentProcessClass = "bkgOneLepFromTop";
         else if(myEvent.is1lepFromW)
             currentProcessClass = "bkgOneLepFromW";
         else if(myEvent.isZtoNuNu)
             currentProcessClass = "bkgZnunu";
         else
             cout<< "there are events which are not Znunu and or 1lep from top/W" << endl;
    }
    if(currentDataset == "ZZTo2L2Nu_powheg_pythia8_25ns")
    {
             currentProcessClass = "bkgLostLepton"; //@MJ@ TODO mistake in babies, maybe somewhere else too?
    }
    

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
vector<string> totYield = { "CR2l_A_250lessMETless350" , "CR2l_A_350lessMETless450" , "CR2l_A_450lessMETless600" , "CR2l_A_600lessMETlessInf" , "CR2l_B_250lessMETless450" , "CR2l_B_450lessMETless600" , "CR2l_B_600lessMETlessInf" , "CR2l_C_250lessMETless350" , "CR2l_C_350lessMETless450" , "CR2l_C_450lessMETless550" , "CR2l_C_550lessMETless650" , "CR2l_C_650lessMETlessInf" , "CR2l_D_250lessMETless350" , "CR2l_D_350lessMETless450" , "CR2l_D_450lessMETless550" , "CR2l_D_550lessMETlessInf" , "CR2l_E_250lessMETless350" , "CR2l_E_350lessMETless550" , "CR2l_E_550lessMETlessInf" , "CR2l_F_250lessMETless450" , "CR2l_F_450lessMETlessInf" , "CR2l_G_250lessMETless350" , "CR2l_G_350lessMETless450" , "CR2l_G_450lessMETless600" , "CR2l_G_600lessMETlessInf" , "CR2l_H_250lessMETless450" , "CR2l_H_450lessMETlessInf" , "CR2l_I_250lessMETless350" , "CR2l_I_350lessMETless450" , "CR2l_I_450lessMETless550" , "CR2l_I_550lessMETlessInf" , "CR0b_A_250lessMETless350" , "CR0b_A_350lessMETless450" , "CR0b_A_450lessMETless600" , "CR0b_A_600lessMETlessInf" , "CR0b_B_250lessMETless450" , "CR0b_B_450lessMETless600" , "CR0b_B_600lessMETlessInf" , "CR0b_C_250lessMETless350" , "CR0b_C_350lessMETless450" , "CR0b_C_450lessMETless550" , "CR0b_C_550lessMETless650" , "CR0b_C_650lessMETlessInf" , "CR0b_D_250lessMETless350" , "CR0b_D_350lessMETless450" , "CR0b_D_450lessMETless550" , "CR0b_D_550lessMETlessInf" , "CR0b_E_250lessMETless350" , "CR0b_E_350lessMETless550" , "CR0b_E_550lessMETlessInf" , "CR0b_F_250lessMETless450" , "CR0b_F_450lessMETlessInf" , "CR0b_G_250lessMETless350" , "CR0b_G_350lessMETless450" , "CR0b_G_450lessMETless600" , "CR0b_G_600lessMETlessInf" , "CR0b_H_250lessMETless450" , "CR0b_H_450lessMETlessInf" , "CR0b_I_250lessMETless350" , "CR0b_I_350lessMETless450" , "CR0b_I_450lessMETless550" , "CR0b_I_550lessMETlessInf" };


    TableDataMC(this, totYield,"lepChannel",  "keepNegative" ).Print(outputName+ "yieldsCRs.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "keepNegative" ).PrintLatex(outputName+ "yieldsCRs.tex", 4);


    ofstream sigfile("controlRegMor.txt");
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
