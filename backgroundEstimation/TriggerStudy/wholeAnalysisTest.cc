#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

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
#include "../../common/selection2016.h"
#include "../../Tools/Weighting/WeightFactory.h"
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
vector<string> process_name;
vector<string> process_table_name;
vector<double> process_syst;
vector<uint16_t> process_stat;
vector<string> process_signal;
vector<double> process_signal_syst;
vector<uint16_t> process_signal_stat;
string* process = new string();
//*process = empty.to_string();
string storedDataset;
TH2D *h2 = NULL;
TAxis *xaxis = NULL;
TAxis *yaxis = NULL;

map< pair<uint32_t,uint32_t>, string > scanMap;

bool lepChannel() 
{ 
    return true; 
}
    
//Add this as a global variable 
WeightFactory WFact;

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/0909/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/Synchro/CMSSW_8_0_5/src/store/babyTuples/TriggerStudy/";
    totalNumberOfWorkers = 10;


    //@EC@: to be done: put default files in the code itself or load info from a external config file
    string WPath = "../../Tools/Weighting/files/";
    //string WPath = "files/";
    //string CSVfileFullSim = WPath+"CSVv2_ichep_slimmed.csv";
    /*string CSVfileFullSim = WPath+"CSVv2_ichep.csv";
    string CSVfileFastSim = WPath+"CSV_13TEV_Combined_14_7_2016.csv";
    string BtagEffFullSimFilename = WPath+"btageff__ttbar_powheg_pythia8_25ns.root";
    string BtagEffFastSimFilename = WPath+"btageff__SMS-T1bbbb-T1qqqq_fastsim.root";
    WFact.InitializeBtagSFTool( CSVfileFullSim,  CSVfileFastSim,  BtagEffFullSimFilename,  BtagEffFastSimFilename);
    WFact.InitializeLeptonSFTool();
    PrintBoxedMessage("Weighting:Initialization ended");*/

    // ------------------
    // Histograms 
    // ------------------
    //vector<float> METBins = {150,175,200,225,250,275,300,350,400,500,800};
    vector<float> METBins = {150,200,250,350,450,600,800};
    vector<float> LepPtBins = {20,25,30,35,40,50,80,100,125,150,200,300};
    //vector<float> METBins = {250,350,450,800};
    //vector<float> LepPtBins = {30,40,50,100,200,300};

    // ------------------
    // Variables
    // ------------------
    // fixed binning
    //AddVariable("MET", "MET",  "MET", 20,   20, 500,  &(myEvent.MET), "noUnderflowInFirstBin,noOverflowInLastBin");
    //AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", 20,   0, 200,  &(myEvent.leadingLeptonPt), "noUnderflowInFirstBin,noOverflowInLastBin");
    //AddVariable("MET", "MET",  "MET", 20,   20, 500,  &(myEvent.MET), "");
    //AddVariable("LeptonPT", "LeptonPT",  "Lepton PT", 20,   0, 200,  &(myEvent.leadingLeptonPt), "");
    // defined binning
    //AddVariable("MET", "MET",  "MET", (int) (METBins.size()-1), METBins.data(),  &(myEvent.MET), "");
    AddVariable("MET", "MET",  "MET", 100 ,200,1000,  &(myEvent.pfmet), "");
    AddVariable("MT2W", "MT2W",  "MT2W", 100 ,0,500,  &(myEvent.MT2W), "");
    AddVariable("MT", "MT",  "MT", 100 ,100,1000,  &(myEvent.mt_met_lep), "");
    AddVariable("topness","topness","topness",100,-20,20,&(myEvent.topness),"");
    AddVariable("nJets","nJets","nJets",5,1,5,&(myEvent.ngoodjets),"");
    AddVariable("dphi","dphi","dphi", 100,0,3.5,&(myEvent.dphi_ak4pfjets_met),"");
    AddVariable("nvertex","nvertex","nvertex",50,0,50,&(myEvent.nvertex),"");

    // ------------------
    // Datasets
    // ------------------
    AddProcessClass("rare", "rare", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("ttZ","rare",0,0.7826);
    //	AddDataset("tZq","rare",0,0.0758);
    //	AddDataset("ZZ","rare",0,0.564);
    //	AddDataset("WZ","rare",0,3.06);

    AddProcessClass("throw", "throw", "signal", kBlue);
     	//AddDataset("T2tt_400to1200", "throw", 0, 0 );
     	AddDataset("signalMerged", "throw", 0, 0 );
    //
    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    //TString fNameTmp =  babyTuplePath+"T2tt_400to1200.root";
    TString fNameTmp =  babyTuplePath+"signalMerged.root";
    ftmp = new TFile(fNameTmp);
    htmp = (TH2D*)ftmp->Get("hStopNeutralino")->Clone();
   
    
    for(uint32_t bx = 0; bx < htmp->GetNbinsX(); bx++)
    {
        for(uint32_t by = 0; by < htmp->GetNbinsY(); by++)
        {
            if(htmp->GetBinContent(bx+1,by+1))
            {
                if( htmp->GetXaxis()->GetBinLowEdge(bx+1) == 600 || htmp->GetXaxis()->GetBinLowEdge(bx+1) == 800) //@MJ@ TODO avoid too many regions
                {
                    std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinLowEdge(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinLowEdge(by+1) << " content: " << htmp->GetBinContent(bx+1,by+1) << std::endl;
                    pair<uint32_t, uint32_t> key = make_pair( htmp->GetXaxis()->GetBinLowEdge(bx+1), htmp->GetYaxis()->GetBinLowEdge(by+1));
                    string stops = to_string(htmp->GetXaxis()->GetBinLowEdge(bx+1));
                    string neutrs = to_string( htmp->GetYaxis()->GetBinLowEdge(by+1));
                    scanMap[key] = stops+"_"+neutrs;
                    AddProcessClass( stops+"_"+neutrs, stops+"_"+neutrs, "signal", kBlue);
                }
            }
        }

    }

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

    AddProcessClass( "grouped", "grouped", "signal", kBlue);
    AddProcessClass("data", "data", "data", kViolet);
    	AddDataset("SE_0", "data", 0, 0 );
    	//AddDataset("SE_1", "data", 0, 0 );
       // AddDataset("SM_0", "data", 0, 0 );
       // AddDataset("SM_1", "data", 0, 0 );
        //AddDataset("MET_0", "data", 0, 0 );
        //AddDataset("MET_1", "data", 0, 0 );
    
    AddProcessClass("test", "test", "background", kRed);
    	AddDataset("ST_s","test",0,10.11*0.364176);
    	//AddDataset("ST_tW_top","test",0,38.09*0.5135);
    	//AddDataset("ST_tW_atop","test",0,38.09*0.5135);
    	//AddDataset("ST_t","test",0,80.95*0.324);
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

    
    AddProcessClass("lostLepton", "lostLepton", "background", kPink);
    AddProcessClass("singleLepton", "singleLepton", "background", kGreen);
    AddProcessClass("singleLeptonFromT", "singleLeptonFromT", "background", kGreen);
    
    // ------------------
    // Regions
    // ------------------
    
    AddRegion("SR1l","SR1l",&SR1l);
    AddRegion("CR1l","CR1l",&CR1l);
    AddRegion("CR2l","CR2l",&CR2l);

    AddRegion("SR1l2jMET250to350","SR1l2jMET250to350",&SR1l2jMET250to350); //@MJ@ TODO thing in time about more flexible naming/regions
    AddRegion("SR1l2jMET350to450","SR1l2jMET350to450",&SR1l2jMET350to450);
    AddRegion("SR1l2jMET450toInf","SR1l2jMET450toInf",&SR1l2jMET450toInf);
    
    AddRegion("SR1l3jMET250to350","SR1l3jMET250to350", &SR1l3jMET250to350);
    AddRegion("SR1l3jMET350to450","SR1l3jMET350to450", &SR1l3jMET350to450);
    AddRegion("SR1l3jMET450to550","SR1l3jMET450to550", &SR1l3jMET450to550);
    AddRegion("SR1l3jMET550toInf","SR1l3jMET550toInf", &SR1l3jMET550toInf);

    AddRegion("SR1l4jMET250to350lowMT2W","SR1l4jMET250to350lowMT2W", &SR1l4jMET250to350lowMT2W);
    AddRegion("SR1l4jMET350to450lowMT2W","SR1l4jMET350to450lowMT2W", &SR1l4jMET350to450lowMT2W);
    AddRegion("SR1l4jMET450toInflowMT2W","SR1l4jMET450toInflowMT2W", &SR1l4jMET450toInflowMT2W);

    AddRegion("SR1l4jMET250to350highMT2W","SR1l4jMET250to350highMT2W", &SR1l4jMET250to350highMT2W);
    AddRegion("SR1l4jMET350to450highMT2W","SR1l4jMET350to450highMT2W", &SR1l4jMET350to450highMT2W);
    AddRegion("SR1l4jMET450to550highMT2W","SR1l4jMET450to550highMT2W", &SR1l4jMET450to550highMT2W);
    AddRegion("SR1l4jMET550to650highMT2W","SR1l4jMET550to650highMT2W", &SR1l4jMET550to650highMT2W);
    AddRegion("SR1l4jMET650toInfhighMT2W","SR1l4jMET650toInfhighMT2W", &SR1l4jMET650toInfhighMT2W);


    AddRegion("CR2l2jMET250to350","CR2l2jMET250to350",&CR2l2jMET250to350);
    AddRegion("CR2l2jMET350to450","CR2l2jMET350to450",&CR2l2jMET350to450);
    AddRegion("CR2l2jMET450toInf","CR2l2jMET450toInf",&CR2l2jMET450toInf);
    
    AddRegion("CR2l3jMET250to350","CR2l3jMET250to350", &CR2l3jMET250to350);
    AddRegion("CR2l3jMET350to450","CR2l3jMET350to450", &CR2l3jMET350to450);
    AddRegion("CR2l3jMET450to550","CR2l3jMET450to550", &CR2l3jMET450to550);
    AddRegion("CR2l3jMET550toInf","CR2l3jMET550toInf", &CR2l3jMET550toInf);

    AddRegion("CR2l4jMET250to350lowMT2W","CR2l4jMET250to350lowMT2W", &CR2l4jMET250to350lowMT2W);
    AddRegion("CR2l4jMET350to450lowMT2W","CR2l4jMET350to450lowMT2W", &CR2l4jMET350to450lowMT2W);
    AddRegion("CR2l4jMET450toInflowMT2W","CR2l4jMET450toInflowMT2W", &CR2l4jMET450toInflowMT2W);

    AddRegion("CR2l4jMET250to350highMT2W","CR2l4jMET250to350highMT2W", &CR2l4jMET250to350highMT2W);
    AddRegion("CR2l4jMET350to450highMT2W","CR2l4jMET350to450highMT2W", &CR2l4jMET350to450highMT2W);
    AddRegion("CR2l4jMET450to550highMT2W","CR2l4jMET450to550highMT2W", &CR2l4jMET450to550highMT2W);
    AddRegion("CR2l4jMET550to650highMT2W","CR2l4jMET550to650highMT2W", &CR2l4jMET550to650highMT2W);
    AddRegion("CR2l4jMET650toInfhighMT2W","CR2l4jMET650toInfhighMT2W", &CR2l4jMET650toInfhighMT2W);


    AddRegion("CR1l2jMET250to350","CR1l2jMET250to350",&CR1l2jMET250to350);
    AddRegion("CR1l2jMET350to450","CR1l2jMET350to450",&CR1l2jMET350to450);
    AddRegion("CR1l2jMET450toInf","CR1l2jMET450toInf",&CR1l2jMET450toInf);
    
    AddRegion("CR1l3jMET250to350","CR1l3jMET250to350", &CR1l3jMET250to350);
    AddRegion("CR1l3jMET350to450","CR1l3jMET350to450", &CR1l3jMET350to450);
    AddRegion("CR1l3jMET450to550","CR1l3jMET450to550", &CR1l3jMET450to550);
    AddRegion("CR1l3jMET550toInf","CR1l3jMET550toInf", &CR1l3jMET550toInf);

    AddRegion("CR1l4jMET250to350lowMT2W","CR1l4jMET250to350lowMT2W", &CR1l4jMET250to350lowMT2W);
    AddRegion("CR1l4jMET350to450lowMT2W","CR1l4jMET350to450lowMT2W", &CR1l4jMET350to450lowMT2W);
    AddRegion("CR1l4jMET450toInflowMT2W","CR1l4jMET450toInflowMT2W", &CR1l4jMET450toInflowMT2W);

    AddRegion("CR1l4jMET250to350highMT2W","CR1l4jMET250to350highMT2W", &CR1l4jMET250to350highMT2W);
    AddRegion("CR1l4jMET350to450highMT2W","CR1l4jMET350to450highMT2W", &CR1l4jMET350to450highMT2W);
    AddRegion("CR1l4jMET450to550highMT2W","CR1l4jMET450to550highMT2W", &CR1l4jMET450to550highMT2W);
    AddRegion("CR1l4jMET550to650highMT2W","CR1l4jMET550to650highMT2W", &CR1l4jMET550to650highMT2W);
    AddRegion("CR1l4jMET650toInfhighMT2W","CR1l4jMET650toInfhighMT2W", &CR1l4jMET650toInfhighMT2W);



    vector<Cut> lowPU;
    lowPU.push_back(Cut("nvertex", '<', 20));
    vector<Cut> highPU;
    highPU.push_back(Cut("nvertex", '>', 19));

    AddRegion("allRegions","allRegions",&allRegions);
    AddRegion("allRegionsLowPU","allRegionsLowPU","allRegions", lowPU);
    AddRegion("allRegionsHighPU","allRegionsHighPU","allRegions", highPU);

    //low PU //@MJ@ TODO other cuts

    AddRegion("SR1l2jMET250to350LowPU","SR1l2jMET250to350LowPU","SR1l2jMET250to350",lowPU); //@MJ@ TODO thing in time about more flexible naming/regions
    AddRegion("SR1l2jMET350to450LowPU","SR1l2jMET350to450LowPU","SR1l2jMET350to450",lowPU);
    AddRegion("SR1l2jMET450toInfLowPU","SR1l2jMET450toInfLowPU","SR1l2jMET450toInf",lowPU);
    
    AddRegion("SR1l3jMET250to350LowPU","SR1l3jMET250to350LowPU", "SR1l3jMET250to350",lowPU);
    AddRegion("SR1l3jMET350to450LowPU","SR1l3jMET350to450LowPU", "SR1l3jMET350to450",lowPU);
    AddRegion("SR1l3jMET450to550LowPU","SR1l3jMET450to550LowPU", "SR1l3jMET450to550",lowPU);
    AddRegion("SR1l3jMET550toInfLowPU","SR1l3jMET550toInfLowPU", "SR1l3jMET550toInf",lowPU);

    AddRegion("SR1l4jMET250to350lowMT2WLowPU","SR1l4jMET250to350lowMT2WLowPU", "SR1l4jMET250to350lowMT2W",lowPU);
    AddRegion("SR1l4jMET350to450lowMT2WLowPU","SR1l4jMET350to450lowMT2WLowPU", "SR1l4jMET350to450lowMT2W",lowPU);
    AddRegion("SR1l4jMET450toInflowMT2WLowPU","SR1l4jMET450toInflowMT2WLowPU", "SR1l4jMET450toInflowMT2W",lowPU);

    AddRegion("SR1l4jMET250to350highMT2WLowPU","SR1l4jMET250to350highMT2WLowPU", "SR1l4jMET250to350highMT2W",lowPU);
    AddRegion("SR1l4jMET350to450highMT2WLowPU","SR1l4jMET350to450highMT2WLowPU", "SR1l4jMET350to450highMT2W",lowPU);
    AddRegion("SR1l4jMET450to550highMT2WLowPU","SR1l4jMET450to550highMT2WLowPU", "SR1l4jMET450to550highMT2W",lowPU);
    AddRegion("SR1l4jMET550to650highMT2WLowPU","SR1l4jMET550to650highMT2WLowPU", "SR1l4jMET550to650highMT2W",lowPU);
    AddRegion("SR1l4jMET650toInfhighMT2WLowPU","SR1l4jMET650toInfhighMT2WLowPU", "SR1l4jMET650toInfhighMT2W",lowPU);

    //high PU

    AddRegion("SR1l2jMET250to350HighPU","SR1l2jMET250to350HighPU","SR1l2jMET250to350",highPU); //@MJ@ TODO thing in time about more flexible naming/regions
    AddRegion("SR1l2jMET350to450HighPU","SR1l2jMET350to450HighPU","SR1l2jMET350to450",highPU);
    AddRegion("SR1l2jMET450toInfHighPU","SR1l2jMET450toInfHighPU","SR1l2jMET450toInf",highPU);
    
    AddRegion("SR1l3jMET250to350HighPU","SR1l3jMET250to350HighPU", "SR1l3jMET250to350",highPU);
    AddRegion("SR1l3jMET350to450HighPU","SR1l3jMET350to450HighPU", "SR1l3jMET350to450",highPU);
    AddRegion("SR1l3jMET450to550HighPU","SR1l3jMET450to550HighPU", "SR1l3jMET450to550",highPU);
    AddRegion("SR1l3jMET550toInfHighPU","SR1l3jMET550toInfHighPU", "SR1l3jMET550toInf",highPU);

    AddRegion("SR1l4jMET250to350lowMT2WHighPU","SR1l4jMET250to350lowMT2WHighPU", "SR1l4jMET250to350lowMT2W",highPU);
    AddRegion("SR1l4jMET350to450lowMT2WHighPU","SR1l4jMET350to450lowMT2WHighPU", "SR1l4jMET350to450lowMT2W",highPU);
    AddRegion("SR1l4jMET450toInflowMT2WHighPU","SR1l4jMET450toInflowMT2WHighPU", "SR1l4jMET450toInflowMT2W",highPU);

    AddRegion("SR1l4jMET250to350highMT2WHighPU","SR1l4jMET250to350highMT2WHighPU", "SR1l4jMET250to350highMT2W",highPU);
    AddRegion("SR1l4jMET350to450highMT2WHighPU","SR1l4jMET350to450highMT2WHighPU", "SR1l4jMET350to450highMT2W",highPU);
    AddRegion("SR1l4jMET450to550highMT2WHighPU","SR1l4jMET450to550highMT2WHighPU", "SR1l4jMET450to550highMT2W",highPU);
    AddRegion("SR1l4jMET550to650highMT2WHighPU","SR1l4jMET550to650highMT2WHighPU", "SR1l4jMET550to650highMT2W",highPU);
    AddRegion("SR1l4jMET650toInfhighMT2WHighPU","SR1l4jMET650toInfhighMT2WHighPU", "SR1l4jMET650toInfhighMT2W",highPU);

    AddRegion("Baseline","Baseline", &Baseline);
    AddRegion("BaselineHighPU","BaselineHighPU","Baseline",highPU);
    AddRegion("BaselineLowPU","BaselineLowPU","Baseline",lowPU);
    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(12900.);

    Create1DHistos();
    //Add2DHisto("LeptonPT","MET");
    //Add2DHisto("nJets","MET");

    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);
    bool useTriggerInfo = currentProcessType == "data" ? true: false;
    myEvent.crossSection  = GetDatasetCrossSection(currentDataset);
    //cout << "cross section: " << myEvent.crossSection << endl;

/*
    //--- Info about the type of dataset is transmitted to the WeightfFactory
    currentProcessType == "data" ? WFact.SetIsData(true): WFact.SetIsData(false);
    currentProcessType == "signal" ? WFact.SetIsFastSim(true): WFact.SetIsFastSim(false);
    //---------------------
*/

    TFile *fle = NULL;
    float weightSignal = -13;
    if (currentProcessType == "signal" && (myEvent.gen_stop_m.at(0) == 300 || myEvent.gen_stop_m.at(0) == 325 || myEvent.gen_stop_m.at(0) == 350 || myEvent.gen_stop_m.at(0) == 375 || myEvent.gen_stop_m.at(0) == 400 || myEvent.gen_stop_m.at(0) == 425 || myEvent.gen_stop_m.at(0) == 450 || myEvent.gen_stop_m.at(0) == 475 || myEvent.gen_stop_m.at(0) == 500 ))
    {
        if(currentDataset != storedDataset && h2 == NULL) //@MJ@ TODO this can work only with one signal dataset!!!
        {
            storedDataset = currentDataset;
            TString fName =  babyTuplePath+currentDataset+".root";
            fle = new TFile(fName);
            h2 = (TH2D*)fle->Get("hStopNeutralino")->Clone();
            xaxis = h2->GetXaxis();
            yaxis = h2->GetYaxis();
        }
        if (h2 == NULL) throw std::runtime_error("The histogram used for CS was not filled!");
        float neutralinoMass = myEvent.gen_neutralino_m.at(0);
        float stopMass = myEvent.gen_stop_m.at(0);
        float sigCrossSection = returnSigCS(stopMass);
        if(myEvent.gen_stop_m.at(0) == 600 || myEvent.gen_stop_m.at(0) == 800)
        {
            pair<uint32_t, uint32_t> yek = make_pair( stopMass,neutralinoMass );
            auto it = scanMap.find( yek );
        
            if ( it != scanMap.end() )
            {
                currentProcessClass = it->second;
            }
            else
            {
                cout << "process class not found ;stop" << stopMass <<" ,neutralino: " <<  neutralinoMass << endl;
            }
       
        }
        else if(stopMass-neutralinoMass == 175)
        {
            currentProcessClass = "grouped";
        }
        else
        {
        }
        Int_t binx = xaxis->FindBin(stopMass);
        Int_t biny = yaxis->FindBin(neutralinoMass);
        uint32_t totalNrOfEvents = h2->GetBinContent(binx, biny);
        myEvent.totalNumberOfInitialEvent = h2->GetEntries();
        weightSignal = sigCrossSection * GetLumi() * myEvent.mc_weight / totalNrOfEvents;
        //cout << "signal CS " << sigCrossSection << " mc weight " << myEvent.mc_weight << " nr of events " << totalNrOfEvents << endl;
     }

    //@MJ@ TODO do a method from this
    if (currentProcessClass == "test" && (myEvent.numberOfGeneratedLeptons >= 2))
    {
        currentProcessClass = "lostLepton";
    }
    else if (currentProcessClass == "test" && (myEvent.numberOfGeneratedLeptons < 2) && currentDataset.find("TTJets")!=std::string::npos)
    {
        currentProcessClass = "singleLeptonFromT";
    }
    else if (currentProcessClass == "test" && (myEvent.numberOfGeneratedLeptons < 2))
    {
        currentProcessClass = "singleLepton";
    }
    else if(currentProcessClass == "test") 
    {
        cout << "nr of leptons " << myEvent.numberOfGeneratedLeptons <<endl;
        throw std::runtime_error("This should not happen");
    }
    else
    {}

    recompute(useTriggerInfo, currentDataset);

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent; //@MJ@ TODO cross section form file?!

    //--- Compute Weights  ---//
    //we should use hadronFlavour and not partonFlavour but it is not available yet
    //WFact.BtagWeighComputor (myEvent.ak4pfjets_pt, myEvent.ak4pfjets_eta, myEvent.ak4pfjets_hadronFlavour, myEvent.ak4pfjets_CSV);// should be called once per event
   /* WFact.BtagWeighComputor (myEvent.ak4pfjets_pt, myEvent.ak4pfjets_eta, myEvent.ak4pfjets_hadronFlavour, myEvent.ak4pfjets_CSV);// should be called once per event
    double btagWeight = WFact.GetBtagW();
*/
    //cout<<"btagWeight = "<<btagWeight<<endl;

    //if(currentDataset.find("WJets")!=std::string::npos || currentDataset.find("W1Jets")!=std::string::npos || currentDataset.find("W2Jets")!=std::string::npos || currentDataset.find("W3Jets")!=std::string::npos || currentDataset.find("W4Jets")!=std::string::npos)
    //{
    //    weightLumi = weightLumi/2;
    //}
    //cout << "CS " << myEvent.crossSection << " lumi " << GetLumi() << " mc weight " << myEvent.mc_weight << " nr of events " << myEvent.totalNumberOfInitialEvent << endl;

    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    if (currentProcessType == "signal") weight = weightSignal;
    //if (currentProcessType == "signal") weight = 1;
    AutoFillProcessClass(currentProcessClass, weight);

    //cout << "weight for process " << currentProcessType << " is " << weight << endl;

    *process = currentProcessClass;
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

    vector<string> processClassLabelList;
    GetProcessClassLabelList(&process_table_name); //@MJ@ TODO why no single lepton, why table name with signal

    for(std::vector<string>::iterator it = processClassLabelList.begin(); it != processClassLabelList.end(); ++it)
    {
        if(*it=="test" || *it=="data" || *it=="throw")
            continue;
        else
        {
           if(*it == "lostLepton")
           {
              process_name.push_back("2l");
              process_syst.push_back(1.3);
              process_stat.push_back(1);
              process_table_name.push_back(*it);
           }
           else if(*it == "singleLepton")
           {
              process_name.push_back("1l");
              process_syst.push_back(1.3);
              process_stat.push_back(1);
              process_table_name.push_back(*it);
           }
           else if(*it == "singleLeptonFromT")
           {
              process_name.push_back("tto1l");
              process_syst.push_back(1.3);
              process_stat.push_back(1);
              process_table_name.push_back(*it);
           }
           else if(*it == "rare")
           {
              process_name.push_back("rare");
              process_syst.push_back(1.3);
              process_stat.push_back(1);
              process_table_name.push_back(*it);
           }
           else
           {
              continue;
           }
        }
    }

    // Schedule plots
    //

    SchedulePlots("1DSuperimposed");
    SchedulePlots("1DSuperimposedNoNorm");
    SchedulePlots("1DStack");
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

    vector<string> processClassTags;
    GetProcessClassTagList(&processClassTags);
    for(uint32_t p =0; p< processClassTags.size(); p++)
    {
        if(GetProcessClassType(processClassTags.at(p)) == "signal")
        {
            //cout << "filling signal: " << processClassTags.at(p) << " ,size: " << processClassTags.size() << endl;
            process_signal.push_back(processClassTags.at(p));
            process_signal_syst.push_back(1.3);
            process_signal_stat.push_back(1);
        }
    }

    string outFold = "datacards";
    system(string("mkdir -p "+outFold).c_str());

    cout << " sieze of signal points: " << process_signal.size() << endl;

    for(uint32_t i = 0; i<process_signal.size(); i++)
    {
          //cout << "in here:" << i << endl;
          ofstream myfile(outFold+"/"+process_signal.at(i)+".txt");
          if (myfile.is_open())
          {
              myfile << "sig " << process_signal.at(i) << " " << process_signal_syst.at(i) << " " << process_signal_stat.at(i) << endl ;
              //cout << "sig" << process_signal.at(i) << " " << process_signal_syst.at(i) << " " << process_signal_stat.at(i) << endl ;
              for(uint32_t j = 0; j<process_name.size(); j++)
              {
                  myfile << process_name.at(j) << " " << process_table_name.at(j) << " " << process_syst.at(j) << " " << process_stat.at(j) << endl ;
                  //cout << process_name.at(j) << " " << process_table_name.at(j) << " " << process_syst.at(j) << " " << process_stat.at(j) << endl ;
              }
              myfile.close();
          }
    }

    vector<string> signalReg = {"SR1l2jMET250to350","SR1l2jMET350to450","SR1l2jMET450toInf","SR1l3jMET250to350","SR1l3jMET350to450","SR1l3jMET450to550","SR1l3jMET550toInf","SR1l4jMET250to350lowMT2W","SR1l4jMET350to450lowMT2W","SR1l4jMET450toInflowMT2W","SR1l4jMET250to350highMT2W","SR1l4jMET350to450highMT2W","SR1l4jMET450to550highMT2W","SR1l4jMET550to650highMT2W","SR1l4jMET650toInfhighMT2W" };
    TableDataMC(this, signalReg,"lepChannel", "includeSignal" ).Print("signalReg.tab", 4);
    TableDataMC(this, signalReg,"lepChannel", "includeSignal" ).PrintLatex("signalReg.tex", 4);
    ofstream sigfile(outFold+"/signalReg.txt");
    if (sigfile.is_open())
    {
        for(uint32_t k = 0; k<signalReg.size(); k++)
        {
	    sigfile << signalReg.at(k)  << endl ;
        }
        sigfile.close();
    }


    vector<string> signalRegPU = {"SR1l2jMET250to350LowPU","SR1l2jMET350to450LowPU","SR1l2jMET450toInfLowPU","SR1l3jMET250to350LowPU","SR1l3jMET350to450LowPU","SR1l3jMET450to550LowPU","SR1l3jMET550toInfLowPU","SR1l4jMET250to350lowMT2WLowPU","SR1l4jMET350to450lowMT2WLowPU","SR1l4jMET450toInflowMT2WLowPU","SR1l4jMET250to350highMT2WLowPU","SR1l4jMET350to450highMT2WLowPU","SR1l4jMET450to550highMT2WLowPU","SR1l4jMET550to650highMT2WLowPU","SR1l4jMET650toInfhighMT2WLowPU", "SR1l2jMET250to350HighPU","SR1l2jMET350to450HighPU","SR1l2jMET450toInfHighPU","SR1l3jMET250to350HighPU","SR1l3jMET350to450HighPU","SR1l3jMET450to550HighPU","SR1l3jMET550toInfHighPU","SR1l4jMET250to350lowMT2WHighPU","SR1l4jMET350to450lowMT2WHighPU","SR1l4jMET450toInflowMT2WHighPU","SR1l4jMET250to350highMT2WHighPU","SR1l4jMET350to450highMT2WHighPU","SR1l4jMET450to550highMT2WHighPU","SR1l4jMET550to650highMT2WHighPU","SR1l4jMET650toInfhighMT2WHighPU", "SR1l2jMET250to350","SR1l2jMET350to450","SR1l2jMET450toInf","SR1l3jMET250to350","SR1l3jMET350to450","SR1l3jMET450to550","SR1l3jMET550toInf","SR1l4jMET250to350lowMT2W","SR1l4jMET350to450lowMT2W","SR1l4jMET450toInflowMT2W","SR1l4jMET250to350highMT2W","SR1l4jMET350to450highMT2W","SR1l4jMET450to550highMT2W","SR1l4jMET550to650highMT2W","SR1l4jMET650toInfhighMT2W" };
    TableDataMC(this, signalRegPU,"lepChannel", "includeSignal" ).Print("signalRegPU.tab", 4);
    TableDataMC(this, signalRegPU,"lepChannel", "includeSignal" ).PrintLatex("signalRegPU.tex", 4);
    
    vector<string> signalRegBlPU = { "BaselineHighPU" ,"BaselineLowPU" };
    TableDataMC(this, signalRegBlPU,"lepChannel", "includeSignal" ).Print("signalRegBlPU.tab", 4);
    TableDataMC(this, signalRegBlPU,"lepChannel", "includeSignal" ).PrintLatex("signalRegBlPU.tex", 4);


    vector<string> signalRegOneBinPU = { "allRegions" ,"allRegionsLowPU", "allRegionsHighPU" };
    TableDataMC(this, signalRegOneBinPU,"lepChannel", "includeSignal" ).Print("signalRegOneBinPU.tab", 4);
    TableDataMC(this, signalRegOneBinPU,"lepChannel", "includeSignal" ).PrintLatex("signalRegOneBinPU.tex", 4);
    
    vector<string> totYield = {"SR1l", "CR1l", "CR2l", "SR1l2jMET250to350","SR1l2jMET350to450","SR1l2jMET450toInf","SR1l3jMET250to350","SR1l3jMET350to450","SR1l3jMET450to550","SR1l3jMET550toInf","SR1l4jMET250to350lowMT2W","SR1l4jMET350to450lowMT2W","SR1l4jMET450toInflowMT2W","SR1l4jMET250to350highMT2W","SR1l4jMET350to450highMT2W","SR1l4jMET450to550highMT2W","SR1l4jMET550to650highMT2W","SR1l4jMET650toInfhighMT2W" , "CR2l2jMET250to350","CR2l2jMET350to450","CR2l2jMET450toInf","CR2l3jMET250to350","CR2l3jMET350to450","CR2l3jMET450to550","CR2l3jMET550toInf","CR2l4jMET250to350lowMT2W","CR2l4jMET350to450lowMT2W","CR2l4jMET450toInflowMT2W","CR2l4jMET250to350highMT2W","CR2l4jMET350to450highMT2W","CR2l4jMET450to550highMT2W","CR2l4jMET550to650highMT2W","CR2l4jMET650toInfhighMT2W" , "CR1l2jMET250to350","CR1l2jMET350to450","CR1l2jMET450toInf","CR1l3jMET250to350","CR1l3jMET350to450","CR1l3jMET450to550","CR1l3jMET550toInf","CR1l4jMET250to350lowMT2W","CR1l4jMET350to450lowMT2W","CR1l4jMET450toInflowMT2W","CR1l4jMET250to350highMT2W","CR1l4jMET350to450highMT2W","CR1l4jMET450to550highMT2W","CR1l4jMET550to650highMT2W","CR1l4jMET650toInfhighMT2W"};
    TableDataMC(this, totYield,"lepChannel" ).Print("yield.tab", 4);
    TableDataMC(this, totYield,"lepChannel" ).PrintLatex("yield.tex", 4);


    vector<string> yieldCR1l = {"CR1l", "CR1l2jMET250to350","CR1l2jMET350to450","CR1l2jMET450toInf","CR1l3jMET250to350","CR1l3jMET350to450","CR1l3jMET450to550","CR1l3jMET550toInf","CR1l4jMET250to350lowMT2W","CR1l4jMET350to450lowMT2W","CR1l4jMET450toInflowMT2W","CR1l4jMET250to350highMT2W","CR1l4jMET350to450highMT2W","CR1l4jMET450to550highMT2W","CR1l4jMET550to650highMT2W","CR1l4jMET650toInfhighMT2W"};
    TableDataMC(this, yieldCR1l,"lepChannel", "includeSignal" ).Print("yieldCR1l.tab", 4);
    TableDataMC(this, yieldCR1l,"lepChannel", "includeSignal" ).PrintLatex("yieldCR1l.tex", 4);
    
    vector<string> yieldCR2l = {"CR2l", "CR2l2jMET250to350","CR2l2jMET350to450","CR2l2jMET450toInf","CR2l3jMET250to350","CR2l3jMET350to450","CR2l3jMET450to550","CR2l3jMET550toInf","CR2l4jMET250to350lowMT2W","CR2l4jMET350to450lowMT2W","CR2l4jMET450toInflowMT2W","CR2l4jMET250to350highMT2W","CR2l4jMET350to450highMT2W","CR2l4jMET450to550highMT2W","CR2l4jMET550to650highMT2W","CR2l4jMET650toInfhighMT2W"};
    TableDataMC(this, yieldCR2l,"lepChannel", "includeSignal" ).Print("yieldCR2l.tab", 4);
    TableDataMC(this, yieldCR2l,"lepChannel", "includeSignal" ).PrintLatex("yieldCR2l.tex", 4);
    
    vector<string> tfreg = {"2jMET250to350","2jMET350to450","2jMET450toInf","3jMET250to350","3jMET350to450","3jMET450to550","3jMET550toInf","4jMET250to350lowMT2W","4jMET350to450lowMT2W","4jMET450toInflowMT2W","4jMET250to350highMT2W","4jMET350to450highMT2W","4jMET450to550highMT2W","4jMET550to650highMT2W","4jMET650toInfhighMT2W"};

    TFProducer prod(tfreg, "lostLepton", "CR2l");
    prod.produceTFTable("yield.tab", "lostLeptonTF");
    TFProducer prod2(tfreg, "singleLepton", "CR1l");
    prod2.produceTFTable("yield.tab", "singleLeptonTF");
    TFProducer prod3(tfreg, "singleLeptonFromT", "CR1l");
    prod3.produceTFTable("yield.tab", "singleLeptonFromTTF");

    cout << "end of processing" << endl;
 }
