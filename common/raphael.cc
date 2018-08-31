
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
#include "goodrun.cc"
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


    // ------------------
    // Datasets
    // ------------------
 
    AddProcessClass("data1", "data1", "data", 1); ///@MJ@ TODO discared this!!
        AddDataset("data_met_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_met_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_met_Run2016H_MINIAOD_03Feb2017_ver2-v1","data1",0,0);
        AddDataset("data_met_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);

        AddDataset("data_single_muon_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_single_muon_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016H_MINIAOD_03Feb2017_ver2-v1","data1",0,0);
        AddDataset("data_single_muon_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);

        AddDataset("data_single_electron_Run2016B_MINIAOD_03Feb2017_ver2-v2","data1",0,0);
        AddDataset("data_single_electron_Run2016C_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016D_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016E_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016F_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016G_MINIAOD_03Feb2017-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016H_MINIAOD_03Feb2017_ver2-v1","data1",0,0);
        AddDataset("data_single_electron_Run2016H_MINIAOD_03Feb2017_ver3-v1","data1",0,0);


    AddProcessClass("bkgOneLepFromW", "bkgOneLepFromW", "background", kBlue); 
    	AddDataset("W1JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0); 
    	AddDataset("W2JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
    	AddDataset("W3JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
    	AddDataset("W4JetsToLNu_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);

        AddDataset("W1JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
        AddDataset("W2JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
        AddDataset("W3JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
        AddDataset("W4JetsToLNu_nupT200_madgraph_pythia8_25ns","bkgOneLepFromW",0,0);
    	AddDataset("WWToLNuQQ_powheg_25ns","bkgOneLepFromW",0,0);

    AddProcessClass("bkgOneLepFromTop", "bkgOneLepFromTop", "background", kRed);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); 
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); 


    AddProcessClass("bkgLostLepton", "bkgLostLepton", "background", kGreen+2); 
    	AddDataset("ttbar_diLept_madgraph_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("ttbar_diLept_madgraph_pythia8_ext1_25ns","bkgLostLepton",0,0);
        AddDataset("tbar_tch_4f_powheg_pythia8_inclDecays_25ns","bkgLostLepton",0,0);
    	AddDataset("t_sch_4f_amcnlo_pythia8_25ns","bkgLostLepton",0,0);
    	AddDataset("t_tW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
    	AddDataset("t_tbarW_5f_powheg_pythia8_noHadDecays_25ns","bkgLostLepton",0,0);
    	AddDataset("WWTo2l2Nu_powheg_25ns","bkgLostLepton",0,0);
    	AddDataset("ttWJets_13TeV_madgraphMLM","bkgLostLepton",0,0);
    AddProcessClass("bkgZnunu", "bkgZnunu", "background", kViolet);
    	AddDataset("ZZTo2L2Nu_powheg_pythia8_25ns","bkgZnunu",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM","bkgZnunu",0,0);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","bkgZnunu",0,0);
  

   //stop 
   AddProcessClass("throw", "throw", "signal", kOrange);
     	AddDataset("Signal_T2tt_mStop_400to1200_1", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_400to1200_2", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_400to1200", "throw", 0, 0 );


       AddProcessClass("signal_500_325",  "T2tt (500/325)",             "signal",kViolet);
       AddProcessClass("signal_800_400",  "T2tt (800/400)",             "signal",kBlue+3);
       AddProcessClass("signal_1000_50",  "T2tt (1000/50)",             "signal",kBlue+6);

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

    AddRegion("SR1l","SR1l",&SR1l);



    //fillYieldsVector(); @MJ@ TODO probably not needed when I do not need to zero negative
                                                                                                                       

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(35.9);

    Create1DHistos();

    Add2DHisto("StopMass","NeutralinoMass");
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

    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);

    TFile *file = NULL;
    double weightSignal = -13;
    double weightSignal_Tight = -13;
    if(currentProcessType == "signal")     {
        if(currentDataset != storedDataset && h2 == NULL) //@MJ@ TODO this can work only with one signal dataset!!!
        {
            storedDataset = currentDataset;
            TString fName =  babyTuplePath+currentDataset+".root";
            file = new TFile(fName);
            h2 = (TH2D*)file->Get("histNEvts")->Clone();
            h3 = (TH3D*)file->Get("h_counterSMS")->Clone();
        }
        if (h2 == NULL) throw std::runtime_error("The histogram used for CS was not filled!");
        if ((myEvent.mass_stop == 1000 && myEvent.mass_lsp == 50) ||
            (myEvent.mass_stop == 800 && myEvent.mass_lsp == 400)  ||
            (myEvent.mass_stop == 500 && myEvent.mass_lsp == 325)) {

            int totalNrOfEvents = (int)h2->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp));
            float xsec = stopXSEC->GetBinContent(stopXSEC->FindBin(myEvent.mass_stop));
            weightSignal =  1000* xsec * GetLumi() /totalNrOfEvents; //@MJ@ TODO improve my method
            weightSignal *=  myEvent.weight_lepSF*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,27))  ) ;
//            weightSignal *=  myEvent.weight_vetoLepSF*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,30))  ) ;
            weightSignal_Tight = weightSignal;
            weightSignal *=  myEvent.weight_btagsf*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,14))  ) ;
            weightSignal_Tight *=  myEvent.weight_tightbtagsf*( totalNrOfEvents / h3->GetBinContent(h3->FindBin(myEvent.mass_stop,myEvent.mass_lsp,37))  ) ;
        }
        else {
           weightSignal=0;
           weightSignal_Tight=0;
        }
    }


    // test 
    if( (currentProcessClass == "bkgLostLepton") || (currentProcessClass == "bkgOneLepFromTop")) {
      if (myEvent.is1lepFromTop) currentProcessClass = "bkgOneLepFromTop";
      else currentProcessClass = "bkgLostLepton";
    }


    if(currentDataset != storedDataset && currentProcessType == "background") //@MJ@ TODO this can work only with one signal dataset!!!
    {
        storedDataset = currentDataset;
        scale1fbS2 = 1;

        if(currentDataset == "ttbar_singleLeptFromTbar_madgraph_pythia8_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns.root"; 
              scale1fbS2= 1.13617e+07/(1.13617e+07 + 4.63189e+07);

        }
        else if(currentDataset == "ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromTbar_madgraph_pythia8_25ns.root"; 
            scale1fbS2= 4.63189e+07/(1.13617e+07 + 4.63189e+07);

        }
        else if(currentDataset == "ttbar_singleLeptFromT_madgraph_pythia8_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns.root"; 
            scale1fbS2=1.16509e+07/(1.16509e+07 + 4.08199e+07);

        }
        else if(currentDataset == "ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_singleLeptFromT_madgraph_pythia8_25ns.root"; 
            scale1fbS2=4.08199e+07/(1.16509e+07 + 4.08199e+07);
        }
        else if(currentDataset == "ttbar_diLept_madgraph_pythia8_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_diLept_madgraph_pythia8_ext1_25ns.root"; 
            scale1fbS2=5.77109e+06/(5.77109e+06 + 2.34556e+07);
        }
        else if(currentDataset == "ttbar_diLept_madgraph_pythia8_ext1_25ns")
        {
            TString fBkgName =  babyTuplePath+"ttbar_diLept_madgraph_pythia8_25ns.root"; 
            scale1fbS2=2.34556e+07/(5.77109e+06 + 2.34556e+07);
        }
        else
        {
            scale1fbS2 = 1;
        }
    }
   
 
    double weight = 0; 

    if(currentProcessType == "signal")
    {
        weight = weightSignal;
    }
    else if (currentProcessType == "data") 
    {
        weight = 1.0;
    }
    else { 
        weight = getWeight(currentProcessType, GetLumi(), scale1fbS2,true); 
    }


// Check WNJets genPt
   if( currentDataset == "W1JetsToLNu_madgraph_pythia8_25ns" || currentDataset == "W2JetsToLNu_madgraph_pythia8_25ns" 
     || currentDataset == "W3JetsToLNu_madgraph_pythia8_25ns" || currentDataset == "W4JetsToLNu_madgraph_pythia8_25ns") {
     if( myEvent.nupt>200.0 ) { weight = 0.;}
   }
   if (currentDataset == "W1JetsToLNu_nupT200_madgraph_pythia8_25ns" || currentDataset == "W2JetsToLNu_nupT200_madgraph_pythia8_25ns"
     || currentDataset == "W3JetsToLNu_nupT200_madgraph_pythia8_25ns" || currentDataset == "W4JetsToLNu_nupT200_madgraph_pythia8_25ns") {
     if (myEvent.nupt<200.0) { weight = 0.; }
   }

    AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields);

    if (myEvent.mass_stop == 1000 && myEvent.mass_lsp == 50) AutoFillProcessClass("signal_1000_50",weight,checkNegativeYields);
    if (myEvent.mass_stop ==  800 && myEvent.mass_lsp == 400) AutoFillProcessClass("signal_800_400",weight,checkNegativeYields);
    if (myEvent.mass_stop ==  500 && myEvent.mass_lsp == 325) AutoFillProcessClass("signal_500_325",weight,checkNegativeYields);

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


    SchedulePlots("1DSuperimposed");
    SchedulePlots("1DStack");
    SchedulePlots("1DDataMCComparison");

    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", false);
    SetGlobalBoolOption("Plot", "exportEps", false);
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


    vector<string> totYield = {"SR1l"};

    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+ "table.tab", 4);

    cout << "end of processing" << endl;
 }


    double getWeight(string currentProcessType, float lumi, double s1fb2, bool mediumbtag)
    {
        double nEvents =  myEvent.wNormalization.at(22);
        double all_weights = lumi*  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) );
        //    all_weights*= myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31));
        if (mediumbtag) all_weights *= myEvent.weight_btagsf*( nEvents / myEvent.wNormalization.at(14) ) ;
        else all_weights *= myEvent.weight_tightbtagsf*( nEvents / myEvent.wNormalization.at(37) ) ;
		// --> with the sytematics, the weights change but also the Normalization bins !!!! To complexify in future...

        if(s1fb2 != 1) all_weights *= s1fb2;

        if(currentProcessType == "signal")
        {
            all_weights = lumi; 
        }
            
        return all_weights;
    }

