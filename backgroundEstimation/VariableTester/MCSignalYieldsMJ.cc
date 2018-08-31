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
#define USE_NEW_VAR
#define USE_GEN_LOSTLEPTON
#include "../../common/Reader_CommonFormat.h"
babyEvent myEvent;

//#include "../../common/common.h"
#include "../../common/TFFactory.h"
#include "../../Selection/testHJFixed.h"
//#include "../../common/selection2016.h"
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



//----------------------
// temporary selection 
//----------------------
bool SR1l_2j() {return SR1l() && myEvent.ngoodjets == 2;}
bool SR1l_3j() {return SR1l() && myEvent.ngoodjets == 3;}
bool SR1l_4j() {return SR1l() && myEvent.ngoodjets >=4;}


bool Pass2LooseB(){
   if(myEvent.ngoodjets >=2){
     int nbtag = 0;
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
     	if(myEvent.ak4pfjets_CSV[i]>=0.46) nbtag++;
     }
     if (nbtag>=2) return true;
   }
   return false;
}

///*
int nMatched2LooseB;

bool SortByCSV(pair<float,int> i, pair<float,int> j) {return (i.first>j.first);}

void MatchBQuarks(){
	vector<pair<float,int> > jets;
	if(myEvent.ak4pfjets_CSV.size() != myEvent.ak4pfjets_hadronFlavour.size()){
		return ;
	}
	for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
		float CSV = myEvent.ak4pfjets_CSV[i];
		int flavour = myEvent.ak4pfjets_hadronFlavour[i];
		pair<float,int> pair;
		pair.first = CSV;
		pair.second = flavour;
		//jets.push_back(make_pair<float,int>(myEvent.ak4pfjets_CSV[i],myEvent.ak4pfjets_hadronFlavour[i]));
		//jets.push_back(make_pair<float,int>(CSV,flavour));
		jets.push_back(pair);
	}
	//sort the collection
	std::sort(jets.begin(), jets.end(), SortByCSV);
	//Check the pdgid of the 2 first jets
	nMatched2LooseB = 0;
	if(jets.size()>=2){
		cout<<jets[0].first<<" "<<jets[1].first<<endl;
		if(fabs(jets[0].second)==5) nMatched2LooseB++;
		if(fabs(jets[1].second)==5) nMatched2LooseB++;
	}
}

/*
bool SR1l_2lb() {return SR1l() && Pass2LooseB();}
bool SR1l_2j_2lb() {return SR1l_2j() && Pass2LooseB();}
bool SR1l_3j_2lb() {return SR1l_3j() && Pass2LooseB();}
bool SR1l_4j_2lb() {return SR1l_4j() && Pass2LooseB();}
bool CR1l_2lb() {return CR1l() && Pass2LooseB();}
*/
//----------------------

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


//new variable

float mStop = 0;
float mNeutralino = 0;
float METSig= 0;
float DeltaPhibb = 0;

bool lepChannel() 
{ 
    return true; 
}
    
//Add this as a global variable 
WeightFactory WFact;

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    //babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/0909/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/BabyTuples/14_09/";
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/BabyTuples/16_09/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/Synchro/CMSSW_8_0_5/src/store/babyTuples/TriggerStudy/";
    totalNumberOfWorkers = 6;


    //@EC@: to be done: put default files in the code itself or load info from a external config file
    string WPath = "../../Tools/Weighting/files/";
    //string WPath = "files/";
    //string CSVfileFullSim = WPath+"CSVv2_ichep_slimmed.csv";
    string CSVfileFullSim = WPath+"CSVv2_ichep.csv";
    string CSVfileFastSim = WPath+"CSV_13TEV_Combined_14_7_2016.csv";
    string BtagEffFullSimFilename = WPath+"btageff__ttbar_powheg_pythia8_25ns.root";
    string BtagEffFastSimFilename = WPath+"btageff__SMS-T1bbbb-T1qqqq_fastsim.root";
    WFact.InitializeBtagSFTool( CSVfileFullSim,  CSVfileFastSim,  BtagEffFullSimFilename,  BtagEffFastSimFilename);
    WFact.InitializeLeptonSFTool();
    PrintBoxedMessage("Weighting:Initialization ended");

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

    AddVariable("StopMass","stop mass", "GeV",400,200,1200, &mStop);
    AddVariable("NeutralinoMass","lsp mass", "GeV",400,0,1000, &mNeutralino);
    AddVariable("MET", "MET",  "MET", 20 ,200,1000,  &(myEvent.pfmet), "");

    
    //AddVariable("nMatched2LooseB","nMatched2LooseB","nMatched2LooseB",3,0,2,&nMatched2LooseB,"");


    // ------------------
    // Datasets
    // ------------------
    ///*
    AddProcessClass("rare", "rare", "background", kBlue);//@MJ@ TODO K-factor?
    	//AddDataset("ttZ","rare",0,0.7826);
    	AddDataset("tZq","rare",0,0.0758);
    	AddDataset("ZZ","rare",0,0.564);
    	//AddDataset("WZ","rare",0,3.06);

   //*/

    ///*
    AddProcessClass("throw", "throw", "signal", kBlue);
     	AddDataset("T2tt_400to1200", "throw", 0, 0 );
    //
    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    TString fNameTmp =  babyTuplePath+"T2tt_400to1200.root";
    ftmp = new TFile(fNameTmp);
    htmp = (TH2D*)ftmp->Get("hStopNeutralino")->Clone();
   
    
    int citer = 0;
    
    for(uint32_t bx = 0; bx < htmp->GetNbinsX(); bx++)
    {
        for(uint32_t by = 0; by < htmp->GetNbinsY(); by++)
        {
            if(htmp->GetBinContent(bx+1,by+1))
            {
                //if( htmp->GetXaxis()->GetBinLowEdge(bx+1) == 500 || htmp->GetXaxis()->GetBinLowEdge(bx+1) == 1000) //@MJ@ TODO avoid too many regions
                if( (htmp->GetXaxis()->GetBinLowEdge(bx+1) == 425 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 300) ||
                    (htmp->GetXaxis()->GetBinLowEdge(bx+1) == 650 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 350)  ||
		    (htmp->GetXaxis()->GetBinLowEdge(bx+1) == 900 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 50))
		{
                    std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinLowEdge(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinLowEdge(by+1) << std::endl;
                    pair<uint32_t, uint32_t> key = make_pair( htmp->GetXaxis()->GetBinLowEdge(bx+1), htmp->GetYaxis()->GetBinLowEdge(by+1));
                    string stops = to_string(htmp->GetXaxis()->GetBinLowEdge(bx+1));
                    string neutrs = to_string( htmp->GetYaxis()->GetBinLowEdge(by+1));
                    scanMap[key] = stops+"_"+neutrs;
                    int i_stop, i_neutrs;
		    stringstream ss; ss<<stops ; ss>>i_stop;
		    ss<<neutrs ; ss>>i_neutrs;
		    stringstream ss2;
		    cout<<ss.str()<<endl;
		    ss2<<"("<<i_stop<<","<<i_neutrs<<")"<<endl;
		    cout<<"name = "<<i_stop<<" "<<i_neutrs<<endl;
		    cout<<"ss = "<<ss2.str()<<endl;
		    //AddProcessClass( stops+"_"+neutrs, stops+"_"+neutrs, "signal", kViolet+citer++);
                    AddProcessClass( stops+"_"+neutrs, ss.str(), "signal", kViolet+citer++);
		}
            }
        }

    }

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

   
   /*

    AddProcessClass("data", "data", "data", COLORPLOT_BLACK); //kViolet);
    	AddDataset("SE_0", "data", 0, 0 );
    	AddDataset("SE_1", "data", 0, 0 );
    	AddDataset("SE_C", "data", 0, 0 );
    	AddDataset("SE_D", "data", 0, 0 );
        AddDataset("SM_0", "data", 0, 0 );
        AddDataset("SM_1", "data", 0, 0 );
        AddDataset("SM_C", "data", 0, 0 );
        AddDataset("SM_D", "data", 0, 0 );
        //AddDataset("MET_0", "data", 0, 0 );
        //AddDataset("MET_1", "data", 0, 0 );
    */ 
    //AddProcessClass("test", "other", "background",kBlack);
    ///*
    AddProcessClass("test", "1l", "background",kGreen);
    	AddDataset("ST_s","test",0,10.11*0.364176);
    	AddDataset("ST_tW_top","test",0,38.09*0.5135);
    	AddDataset("ST_tW_atop","test",0,38.09*0.5135);
    	AddDataset("ST_t","test",0,80.95*0.324);
    	AddDataset("W1JetsToLNuTune","test",0, 9493*1.238);
    	AddDataset("W2JetsToLNuTune","test",0, 3120*1.231);
    	AddDataset("W3JetsToLNuTune","test",0, 942.3*1.231);
    	AddDataset("W4JetsToLNuTune","test",0, 524.2*1.114);
    	AddDataset("TTWtoQQ","test",0,0.4062);
    	AddDataset("TTWtoLNu","test",0,0.2043);
    	AddDataset("TTT","test",0,1.0);
    	AddDataset("VV","test",0,12.05*0.9917);

    AddProcessClass("singleLeptonFromT", "1l from top", "background", kOrange);
	AddDataset("TTJetsSLtop", "singleLeptonFromT", 0, 114.6*1.594 );
    	AddDataset("TTJetsSLatopv1","singleLeptonFromT",0,114.6*1.594);

    AddProcessClass("lostLepton", "lost Lepton", "background", kPink);
    	AddDataset("TTJetsDLv0v4","lostLepton",0, 57.35*1.5225);
    	//AddDataset("WJetsToLNuTune","test",0,60781.5*1.01);

   ///* 
    //AddProcessClass("singleLepton", "1l", "background", kGreen);
   //*/

    // ------------------
    // Regions
    // ------------------
/*   
AddRegion("SR1l23jLowMlb_MET250to350","SR1l23jLowMlb_MET250to350",&SR1l23jLowMlb_MET250to350);
AddRegion("SR1l23jLowMlb_MET350to450","SR1l23jLowMlb_MET350to450",&SR1l23jLowMlb_MET350to450);
AddRegion("SR1l23jLowMlb_MET450to550","SR1l23jLowMlb_MET450to550",&SR1l23jLowMlb_MET450to550);
AddRegion("SR1l23jLowMlb_MET550toInf","SR1l23jLowMlb_MET550toInf",&SR1l23jLowMlb_MET550toInf);
AddRegion("SR1l23jHighMlb_MET250to350","SR1l23jHighMlb_MET250to350",&SR1l23jHighMlb_MET250to350);
AddRegion("SR1l23jHighMlb_MET350to450","SR1l23jHighMlb_MET350to450",&SR1l23jHighMlb_MET350to450);
AddRegion("SR1l23jHighMlb_MET450toInf","SR1l23jHighMlb_MET450toInf",&SR1l23jHighMlb_MET450toInf);
AddRegion("SR1l4jLowMT2WLowMlb_MET250to350","SR1l4jLowMT2WLowMlb_MET250to350",&SR1l4jLowMT2WLowMlb_MET250to350);
AddRegion("SR1l4jLowMT2WLowMlb_MET350to450","SR1l4jLowMT2WLowMlb_MET350to450",&SR1l4jLowMT2WLowMlb_MET350to450);
AddRegion("SR1l4jLowMT2WLowMlb_MET450to550","SR1l4jLowMT2WLowMlb_MET450to550",&SR1l4jLowMT2WLowMlb_MET450to550);
AddRegion("SR1l4jLowMT2WLowMlb_MET550to650","SR1l4jLowMT2WLowMlb_MET550to650",&SR1l4jLowMT2WLowMlb_MET550to650);
AddRegion("SR1l4jLowMT2WLowMlb_MET650toInf","SR1l4jLowMT2WLowMlb_MET650toInf",&SR1l4jLowMT2WLowMlb_MET650toInf);
AddRegion("SR1l4jLowMT2WHighMlb_MET250to350","SR1l4jLowMT2WHighMlb_MET250to350",&SR1l4jLowMT2WHighMlb_MET250to350);
AddRegion("SR1l4jLowMT2WHighMlb_MET350to450","SR1l4jLowMT2WHighMlb_MET350to450",&SR1l4jLowMT2WHighMlb_MET350to450);
AddRegion("SR1l4jLowMT2WHighMlb_MET450to550","SR1l4jLowMT2WHighMlb_MET450to550",&SR1l4jLowMT2WHighMlb_MET450to550);
AddRegion("SR1l4jLowMT2WHighMlb_MET550toInf","SR1l4jLowMT2WHighMlb_MET550toInf",&SR1l4jLowMT2WHighMlb_MET550toInf);
AddRegion("SR1l4jMidMT2WLowMlb_MET250to350","SR1l4jMidMT2WLowMlb_MET250to350",&SR1l4jMidMT2WLowMlb_MET250to350);
AddRegion("SR1l4jMidMT2WLowMlb_MET350toInf","SR1l4jMidMT2WLowMlb_MET350toInf",&SR1l4jMidMT2WLowMlb_MET350toInf);
AddRegion("SR1l4jMidMT2WHighMlb_MET250toInf","SR1l4jMidMT2WHighMlb_MET250toInf",&SR1l4jMidMT2WHighMlb_MET250toInf);
AddRegion("SR1l4jHighMT2WLowMlb_MET250to350","SR1l4jHighMT2WLowMlb_MET250to350",&SR1l4jHighMT2WLowMlb_MET250to350);
AddRegion("SR1l4jHighMT2WLowMlb_MET350to450","SR1l4jHighMT2WLowMlb_MET350to450",&SR1l4jHighMT2WLowMlb_MET350to450);
AddRegion("SR1l4jHighMT2WLowMlb_MET450to600","SR1l4jHighMT2WLowMlb_MET450to600",&SR1l4jHighMT2WLowMlb_MET450to600);
AddRegion("SR1l4jHighMT2WLowMlb_MET600toInf","SR1l4jHighMT2WLowMlb_MET600toInf",&SR1l4jHighMT2WLowMlb_MET600toInf);
AddRegion("SR1l4jHighMT2WHighMlb_MET250to400","SR1l4jHighMT2WHighMlb_MET250to400",&SR1l4jHighMT2WHighMlb_MET250to400);
AddRegion("SR1l4jHighMT2WHighMlb_MET400toInf","SR1l4jHighMT2WHighMlb_MET400toInf",&SR1l4jHighMT2WHighMlb_MET400toInf);
*/


AddRegion("SR1l23jLowMlb_MET250to350","SR1l23jLowMlb_MET250to350",&SR1l23jLowMlb_MET250to350);
AddRegion("SR1l23jLowMlb_MET350to450","SR1l23jLowMlb_MET350to450",&SR1l23jLowMlb_MET350to450);
AddRegion("SR1l23jLowMlb_MET450to550","SR1l23jLowMlb_MET450to550",&SR1l23jLowMlb_MET450to550);
AddRegion("SR1l23jLowMlb_MET550toInf","SR1l23jLowMlb_MET550toInf",&SR1l23jLowMlb_MET550toInf);
AddRegion("SR1l23jHighMlb_MET250to350","SR1l23jHighMlb_MET250to350",&SR1l23jHighMlb_MET250to350);
AddRegion("SR1l23jHighMlb_MET350to450","SR1l23jHighMlb_MET350to450",&SR1l23jHighMlb_MET350to450);
AddRegion("SR1l23jHighMlb_MET450to550","SR1l23jHighMlb_MET450to550",&SR1l23jHighMlb_MET450to550);
AddRegion("SR1l23jHighMlb_MET550toInf","SR1l23jHighMlb_MET550toInf",&SR1l23jHighMlb_MET550toInf);
AddRegion("SR1l4jLowMT2WLowMlb_MET250to350","SR1l4jLowMT2WLowMlb_MET250to350",&SR1l4jLowMT2WLowMlb_MET250to350);
AddRegion("SR1l4jLowMT2WLowMlb_MET350to450","SR1l4jLowMT2WLowMlb_MET350to450",&SR1l4jLowMT2WLowMlb_MET350to450);
AddRegion("SR1l4jLowMT2WLowMlb_MET450to550","SR1l4jLowMT2WLowMlb_MET450to550",&SR1l4jLowMT2WLowMlb_MET450to550);
AddRegion("SR1l4jLowMT2WLowMlb_MET550to650","SR1l4jLowMT2WLowMlb_MET550to650",&SR1l4jLowMT2WLowMlb_MET550to650);
AddRegion("SR1l4jLowMT2WLowMlb_MET650toInf","SR1l4jLowMT2WLowMlb_MET650toInf",&SR1l4jLowMT2WLowMlb_MET650toInf);
AddRegion("SR1l4jLowMT2WHighMlb_MET250to350","SR1l4jLowMT2WHighMlb_MET250to350",&SR1l4jLowMT2WHighMlb_MET250to350);
AddRegion("SR1l4jLowMT2WHighMlb_MET350to450","SR1l4jLowMT2WHighMlb_MET350to450",&SR1l4jLowMT2WHighMlb_MET350to450);
AddRegion("SR1l4jLowMT2WHighMlb_MET450to550","SR1l4jLowMT2WHighMlb_MET450to550",&SR1l4jLowMT2WHighMlb_MET450to550);
AddRegion("SR1l4jLowMT2WHighMlb_MET550toInf","SR1l4jLowMT2WHighMlb_MET550toInf",&SR1l4jLowMT2WHighMlb_MET550toInf);
AddRegion("SR1l4jMidMT2WLowMlb_MET250to350","SR1l4jMidMT2WLowMlb_MET250to350",&SR1l4jMidMT2WLowMlb_MET250to350);
AddRegion("SR1l4jMidMT2WLowMlb_MET350to450","SR1l4jMidMT2WLowMlb_MET350to450",&SR1l4jMidMT2WLowMlb_MET350to450);
AddRegion("SR1l4jMidMT2WLowMlb_MET450toInf","SR1l4jMidMT2WLowMlb_MET450toInf",&SR1l4jMidMT2WLowMlb_MET450toInf);
AddRegion("SR1l4jMidMT2WHighMlb_MET250to400","SR1l4jMidMT2WHighMlb_MET250to400",&SR1l4jMidMT2WHighMlb_MET250to400);
AddRegion("SR1l4jMidMT2WHighMlb_MET400toInf","SR1l4jMidMT2WHighMlb_MET400toInf",&SR1l4jMidMT2WHighMlb_MET400toInf);
AddRegion("SR1l4jHighMT2WLowMlb_MET250to350","SR1l4jHighMT2WLowMlb_MET250to350",&SR1l4jHighMT2WLowMlb_MET250to350);
AddRegion("SR1l4jHighMT2WLowMlb_MET350to450","SR1l4jHighMT2WLowMlb_MET350to450",&SR1l4jHighMT2WLowMlb_MET350to450);
AddRegion("SR1l4jHighMT2WLowMlb_MET450to600","SR1l4jHighMT2WLowMlb_MET450to600",&SR1l4jHighMT2WLowMlb_MET450to600);
AddRegion("SR1l4jHighMT2WLowMlb_MET600toInf","SR1l4jHighMT2WLowMlb_MET600toInf",&SR1l4jHighMT2WLowMlb_MET600toInf);
AddRegion("SR1l4jHighMT2WHighMlb_MET250to400","SR1l4jHighMT2WHighMlb_MET250to400",&SR1l4jHighMT2WHighMlb_MET250to400);
AddRegion("SR1l4jHighMT2WHighMlb_MET400to650","SR1l4jHighMT2WHighMlb_MET400to650",&SR1l4jHighMT2WHighMlb_MET400to650);
AddRegion("SR1l4jHighMT2WHighMlb_MET650toInf","SR1l4jHighMT2WHighMlb_MET650toInf",&SR1l4jHighMT2WHighMlb_MET650toInf);

    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(12900.);
    //SetLumi(3900.);

    Create1DHistos();
    //Add2DHisto("LeptonPT","MET");
    //Add2DHisto("nJets","MET");
    Add2DHisto("StopMass","NeutralinoMass");

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


    //--- Info about the type of dataset is transmitted to the WeightfFactory
    currentProcessType == "data" ? WFact.SetIsData(true): WFact.SetIsData(false);
    currentProcessType == "signal" ? WFact.SetIsFastSim(true): WFact.SetIsFastSim(false);
    //---------------------

    // --- Recompute variables event per event
    METSig = myEvent.pfmet/myEvent.ak4_HT;
    DeltaPhibb = fabs(myEvent.DeltaPhibb);
    MatchBQuarks(); // call this function to compute nMatched2LooseB 

    mStop = -999;
    mNeutralino = -99;
    if (currentProcessType == "signal"){
      mStop = myEvent.gen_stop_m.at(0);
      mNeutralino = myEvent.gen_neutralino_m.at(0); 
      //cout<<"Masses = "<<mStop<<" "<<mNeutralino<<endl;
    }

    TFile *file = NULL;
    float weightSignal = -13;
    if (currentProcessType == "signal")//  && (myEvent.gen_stop_m.at(0) == 500 || myEvent.gen_stop_m.at(0) == 425 || myEvent.gen_stop_m.at(0) == 900))
    {
        if(currentDataset != storedDataset && h2 == NULL) //@MJ@ TODO this can work only with one signal dataset!!!
        {
            storedDataset = currentDataset;
            TString fName =  babyTuplePath+currentDataset+".root";
            file = new TFile(fName);
            h2 = (TH2D*)file->Get("hStopNeutralino")->Clone();
            xaxis = h2->GetXaxis();
            yaxis = h2->GetYaxis();
        }
        if (h2 == NULL) throw std::runtime_error("The histogram used for CS was not filled!");
        float neutralinoMass = myEvent.gen_neutralino_m.at(0);
        float stopMass = myEvent.gen_stop_m.at(0);
        float sigCrossSection = returnSigCS(stopMass);
        pair<uint32_t, uint32_t> yek = make_pair( stopMass,neutralinoMass );
        auto it = scanMap.find( yek );
        
        if ( it != scanMap.end() )
        {
            currentProcessClass = it->second;
        }
        else
        {
        	//    cout << "process class not found ;stop" << stopMass <<" ,neutralino: " <<  neutralinoMass << endl;
        	//return;
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
    else if (currentProcessClass == "test" && (myEvent.numberOfGeneratedLeptons < 2))
    {
        //currentProcessClass = "singleLepton";
        currentProcessClass = "test";
    }
    else
    {}

    //recompute(useTriggerInfo, currentDataset);

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent; //@MJ@ TODO cross section form file?!

    //--- Compute Weights  ---//
    //we should use hadronFlavour and not partonFlavour but it is not available yet
    //cout<<myEvent.ak4pfjets_pt.size()<<" "<<myEvent.ak4pfjets_eta.size()<<" "<<myEvent.ak4pfjets_hadronFlavour.size()<<" "<<myEvent.ak4pfjets_CSV.size()<<endl;
    WFact.BtagWeighComputor (myEvent.ak4pfjets_pt, myEvent.ak4pfjets_eta, myEvent.ak4pfjets_hadronFlavour, myEvent.ak4pfjets_CSV);// should be called once per event
    //WFact.BtagWeighComputor (myEvent.ak4pfjets_pt, myEvent.ak4pfjets_eta, myEvent.ak4pfjets_partonFlavour, myEvent.ak4pfjets_CSV);// should be called once per event
    double btagWeight = WFact.GetBtagW();
    if(btagWeight==0) cout<<"btagWeight = "<<btagWeight<<endl;

    double lepWeight = 0;
    if(currentProcessType != "data"){
    	WFact.LeptonWeightComputor(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_pdgid, myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_pdgid, myEvent.nvetoleps, myEvent.numberOfSelectedLeptons, myEvent.genLostLeptons_pt.size(), myEvent.genLostLeptons_pt, myEvent.genLostLeptons_eta, myEvent.genLostLeptons_pdgid );
	lepWeight = WFact.GetLepW();
    	//cout<<"lepw = "<<lepWeight<<endl;
    }



    //if(currentDataset.find("WJets")!=std::string::npos || currentDataset.find("W1Jets")!=std::string::npos || currentDataset.find("W2Jets")!=std::string::npos || currentDataset.find("W3Jets")!=std::string::npos || currentDataset.find("W4Jets")!=std::string::npos)
    //{
    //    weightLumi = weightLumi/2;
    //}
    //cout << "CS " << myEvent.crossSection << " lumi " << GetLumi() << " mc weight " << myEvent.mc_weight << " nr of events " << myEvent.totalNumberOfInitialEvent << endl;

    float weight     = weightLumi*lepWeight*btagWeight;
    if (currentProcessType == "signal") weight = weightSignal*lepWeight*btagWeight;
    if (currentProcessType == "data") weight = 1.0;
    //else weight*=btagWeight*lepWeight;
    //else weight*=1;
    
    AutoFillProcessClass(currentProcessClass, weight);

    //cout << "weight for process " << currentProcessType << " is " << weight << endl;

    *process = currentProcessClass;
    if(counter % 10000 == 0)
    {
        cout << counter << endl;
    }

    //file->Delete();
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

    //SchedulePlots("1DSuperimposed");
    //SchedulePlots("1DStack");
    //SchedulePlots("1DDataMCComparison");
    //SchedulePlots("1DFigureOfMerit","var=DeltaPhibb,cutType=keepHighValues");
    //SchedulePlots("1DFigureOfMerit","var=DeltaPhibb,cutType=keepLowValues,type=sOverSqrtB,minBackground=0");
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

   /*
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
   */
vector<string> signalReg = { "SR1l23jLowMlb_MET250to350" , "SR1l23jLowMlb_MET350to450" , "SR1l23jLowMlb_MET450to550" , "SR1l23jLowMlb_MET550toInf" , "SR1l23jHighMlb_MET250to350" , "SR1l23jHighMlb_MET350to450" , "SR1l23jHighMlb_MET450to550" , "SR1l23jHighMlb_MET550toInf" , "SR1l4jLowMT2WLowMlb_MET250to350" , "SR1l4jLowMT2WLowMlb_MET350to450" , "SR1l4jLowMT2WLowMlb_MET450to550" , "SR1l4jLowMT2WLowMlb_MET550to650" , "SR1l4jLowMT2WLowMlb_MET650toInf" , "SR1l4jLowMT2WHighMlb_MET250to350" , "SR1l4jLowMT2WHighMlb_MET350to450" , "SR1l4jLowMT2WHighMlb_MET450to550" , "SR1l4jLowMT2WHighMlb_MET550toInf" , "SR1l4jMidMT2WLowMlb_MET250to350" , "SR1l4jMidMT2WLowMlb_MET350to450" , "SR1l4jMidMT2WLowMlb_MET450toInf" , "SR1l4jMidMT2WHighMlb_MET250to400" , "SR1l4jMidMT2WHighMlb_MET400toInf" , "SR1l4jHighMT2WLowMlb_MET250to350" , "SR1l4jHighMT2WLowMlb_MET350to450" , "SR1l4jHighMT2WLowMlb_MET450to600" , "SR1l4jHighMT2WLowMlb_MET600toInf" , "SR1l4jHighMT2WHighMlb_MET250to400" , "SR1l4jHighMT2WHighMlb_MET400to650" , "SR1l4jHighMT2WHighMlb_MET650toInf"};
//        vector<string> signalReg = { "SR1l23jLowMlb_MET250to350" , "SR1l23jLowMlb_MET350to450" , "SR1l23jLowMlb_MET450to550" , "SR1l23jLowMlb_MET550toInf" , "SR1l23jHighMlb_MET250to350" , "SR1l23jHighMlb_MET350to450" , "SR1l23jHighMlb_MET450toInf" , "SR1l4jLowMT2WLowMlb_MET250to350" , "SR1l4jLowMT2WLowMlb_MET350to450" , "SR1l4jLowMT2WLowMlb_MET450to550" , "SR1l4jLowMT2WLowMlb_MET550to650" , "SR1l4jLowMT2WLowMlb_MET650toInf" , "SR1l4jLowMT2WHighMlb_MET250to350" , "SR1l4jLowMT2WHighMlb_MET350to450" , "SR1l4jLowMT2WHighMlb_MET450to550" , "SR1l4jLowMT2WHighMlb_MET550toInf" , "SR1l4jMidMT2WLowMlb_MET250to350" , "SR1l4jMidMT2WLowMlb_MET350toInf" , "SR1l4jMidMT2WHighMlb_MET250toInf" , "SR1l4jHighMT2WLowMlb_MET250to350" , "SR1l4jHighMT2WLowMlb_MET350to450" , "SR1l4jHighMT2WLowMlb_MET450to600" , "SR1l4jHighMT2WLowMlb_MET600toInf" , "SR1l4jHighMT2WHighMlb_MET250to400" , "SR1l4jHighMT2WHighMlb_MET400toInf" }; 
    TableDataMC(this, signalReg,"lepChannel", "includeSignal" ).Print("yield.tab", 4);
    TableZbi(this, signalReg,"lepChannel").Print("Zbi.tab", 4);
    CombineCardMaker card(this, signalReg, "lepChannel", "throw", "StopMass", "NeutralinoMass");
    card.Print("card.tab", 4);
    card.ProduceCard("cards", 1.2, 1.5);
    
    cout << "end of processing" << endl;
 }
