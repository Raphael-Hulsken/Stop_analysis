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
    totalNumberOfWorkers = 5;

    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    TString fNameTmp =  babyTuplePath+"Signal_T6ttWW.root"; //sbottom
    ///TString fNameTmp =  babyTuplePath+"signal_T2tt_forBinning.root"; //stop
    ftmp = new TFile(fNameTmp);
    htmp = (TH2D*)ftmp->Get("histNEvts")->Clone();
    TH1D* pX = htmp->ProjectionX(); 
    TH1D* pY = htmp->ProjectionY();
    
   

    AddVariable("MET", "MET",  "MET", 100 ,200,1000,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "MT2W", 100 ,0,500,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("MT", "MT",  "MT", 100 ,100,1000,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
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
    AddProcessClass("data1", "data1", "data", kBlue); ///@MJ@ TODO discared this!!
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
    	//AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns_full","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromTbar_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); //@MJ@ TODO what about extentions, good normalization to total number of events?
    	//AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns_full","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_25ns","bkgOneLepFromTop",0,0);
    	AddDataset("ttbar_singleLeptFromT_madgraph_pythia8_ext1_25ns","bkgOneLepFromTop",0,0); 
    AddProcessClass("bkgLostLepton", "bkgLostLepton", "background", kBlue); 
    	//AddDataset("ttbar_diLept_madgraph_pythia8_25ns_full","bkgLostLepton",0,0);
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
    
    

   //sbottom
   AddProcessClass("throw", "throw", "signal", kBlue);
     	AddDataset("Signal_T6ttWW", "throw", 0, 0 );
     	AddDataset("Signal_T6ttWW_1", "throw", 0, 0 );
     	AddDataset("Signal_T6ttWW_2", "throw", 0, 0 );
     	AddDataset("Signal_T6ttWW_3", "throw", 0, 0 );
     	AddDataset("Signal_T6ttWW_4", "throw", 0, 0 );
     	AddDataset("Signal_T6ttWW_5", "throw", 0, 0 );
     	AddDataset("Signal_T6ttWW_6", "throw", 0, 0 );
   
   //stop 
   /*AddProcessClass("throw", "throw", "signal", kBlue);
     	AddDataset("Signal_T2tt_mStop_150to250", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_250to350", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_350to400_1", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_350to400", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_400to1200_1", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_400to1200_2", "throw", 0, 0 );
     	AddDataset("Signal_T2tt_mStop_400to1200", "throw", 0, 0 );
*/
    int citer = 0;
    for(uint32_t bx = 0; bx < htmp->GetNbinsX(); bx++)
    {
        for(uint32_t by = 0; by < htmp->GetNbinsY(); by++)
        {
            if(htmp->GetBinContent(bx+1,by+1))
            {                
                    if( (htmp->GetXaxis()->GetBinCenter(bx+1) == 1000 && htmp->GetYaxis()->GetBinCenter(by+1) == 50) ||
                    (htmp->GetXaxis()->GetBinCenter(bx+1) == 800 && htmp->GetYaxis()->GetBinCenter(by+1) == 400)  ||
		    (htmp->GetXaxis()->GetBinCenter(bx+1) == 500 && htmp->GetYaxis()->GetBinCenter(by+1) == 325)) //@MJ@ TODO just slimming in here
		    {
                    //std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinCenter(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinCenter(by+1) << std::endl;
                    std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinCenter(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinCenter(by+1) << std::endl;
                    pair<uint32_t, uint32_t> key = make_pair( htmp->GetXaxis()->GetBinCenter(bx+1), htmp->GetYaxis()->GetBinCenter(by+1));
                    string sbottoms = to_string(htmp->GetXaxis()->GetBinCenter(bx+1));
                    string neutrs = to_string( htmp->GetYaxis()->GetBinCenter(by+1));
                    scanMap[key] = sbottoms+"_"+neutrs;
                    std::string bval = sbottoms.substr(0,sbottoms.length()-7); //6 zeroes usually
                    std::string lspval = neutrs.substr(0,neutrs.length()-7); //6 zeroes usually
                    cout << "creating mass point: " <<"("+bval+","+lspval+")" << endl;
                    //AddProcessClass( sbottoms+"_"+neutrs, "("+bval+","+lspval+")", "signal", kViolet+citer++); //@MJ@ TODO normally should be signal, but it would fail
                    cout << "process class added" << endl;
		    }
            }
        }

    }

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

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




    //fillYieldsVector(); @MJ@ TODO probably not needed when I do not need to zero negative
                                                                                                                       

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(35.867);

    Create1DHistos();
    Add2DHisto("StopMass","NeutralinoMass");

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

    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
    cout << "returned trigger " << myEvent.trigger << endl;

    TFile *file = NULL;
    float weightSignal = -13;
    if(currentProcessType == "signal")     {
        if(currentDataset != storedDataset && h2 == NULL) //@MJ@ TODO this can work only with one signal dataset!!!
        {
            storedDataset = currentDataset;
            TString fName =  babyTuplePath+currentDataset+".root";
            file = new TFile(fName);
            h2 = (TH2D*)file->Get("histNEvts")->Clone();
            xaxis = h2->GetXaxis();
            yaxis = h2->GetYaxis();
        }
        if (h2 == NULL) throw std::runtime_error("The histogram used for CS was not filled!");
        float neutralinoMass = myEvent.mass_lsp;
        float sbottomMass = myEvent.mass_stop;
        //pair<uint32_t, uint32_t> yek = make_pair( sbottomMass,neutralinoMass );
        //auto it = scanMap.find( yek );
        //&& ( (myEvent.mass_stop == 1000 && myEvent.mass_lsp == 50) ||
                    //(myEvent.mass_stop == 800 && myEvent.mass_lsp == 400)  ||
	/*{	    //(myEvent.mass_stop == 500 && myEvent.mass_lsp == 325))) //@MJ@ TODO just slimming in here

        if ( it != scanMap.end() )
        {
            //cout << "class found" << endl;
            //currentProcessClass = it->second;
            cout << "process class " << currentProcessClass << endl;
        }
        else
        {
                cout << "signal point was not found" << endl;
        	//return;
	}
        }*/
        Int_t binx = xaxis->FindBin(sbottomMass);
        Int_t biny = yaxis->FindBin(neutralinoMass);
        uint32_t totalNrOfEvents = h2->GetBinContent(binx, biny);
        uint32_t binXSEC = stopXSEC->FindBin(sbottomMass);
        float xsec = stopXSEC->GetBinContent(binXSEC);
        //myEvent.totalNumberOfInitialEvent = h2->GetEntries();
        weightSignal =  1000* xsec *getWeight("signal", GetLumi()) /totalNrOfEvents; //@MJ@ TODO improve my method
        //weightSignal *=  myEvent.weight_lepSF*( totalNrOfEvents / myEvent.wNormalization.at(27) ) * myEvent.weight_vetoLepSF*( totalNrOfEvents / myEvent.wNormalization.at(30))* myEvent.weight_btagsf*( totalNrOfEvents / myEvent.wNormalization.at(14) ) ;
     }


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
    //if( (currentProcessClass == "throw"))
    //{
    //     currentProcessClass = "";
    //}
    
    //float nEventsN =  myEvent.wNormalization.at(22);
    

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
   
    //cout << currentDataset << " nevt2 " << scale1fbS2 << endl;
 
    float weightLumi = getWeight(currentProcessType, GetLumi(), scale1fbS2); 
    float weight     = weightLumi;
    if(currentProcessType == "signal")
    {
        weight = weightSignal;
    //cout << "weight signal latest " << weight << endl;
    }

    //cout << "weight " << weight << endl;

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
vector<string> totYield = {"SR1l_A_250lessMETless350" , "SR1l_A_350lessMETless450" , "SR1l_A_450lessMETless600" , "SR1l_A_600lessMETlessInf" , "SR1l_B_250lessMETless450" , "SR1l_B_450lessMETless600" , "SR1l_B_600lessMETlessInf" , "SR1l_C_250lessMETless350" , "SR1l_C_350lessMETless450" , "SR1l_C_450lessMETless550" , "SR1l_C_550lessMETless650" , "SR1l_C_650lessMETlessInf" , "SR1l_D_250lessMETless350" , "SR1l_D_350lessMETless450" , "SR1l_D_450lessMETless550" , "SR1l_D_550lessMETlessInf" , "SR1l_E_250lessMETless350" , "SR1l_E_350lessMETless550" , "SR1l_E_550lessMETlessInf" , "SR1l_F_250lessMETless450" , "SR1l_F_450lessMETlessInf" , "SR1l_G_250lessMETless350" , "SR1l_G_350lessMETless450" , "SR1l_G_450lessMETless600" , "SR1l_G_600lessMETlessInf" , "SR1l_H_250lessMETless450" , "SR1l_H_450lessMETlessInf", "SR1l_I_250lessMETless350" , "SR1l_I_350lessMETless450" , "SR1l_I_450lessMETless550" , "SR1l_I_550lessMETlessInf"};

    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+ "sbottomSignal.tab", 4);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+ "sbottomSignal.tex", 4);

    CombineCardMaker card(this, totYield , "lepChannel", "throw", "StopMass", "NeutralinoMass");
    card.Print("card.tab", 4);
    card.ProduceCard("card.tab" ,"../cards");
    //version for updating cards
    //CombineCardMaker card(this, totYield , "lepChannel", "none", "StopMass", "NeutralinoMass");
    //card.UpdateCardTable("card.tab");


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
