#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <boost/math/constants/constants.hpp>

#define USE_VAR_BASELINE
#define USE_LEP1
#define USE_LEP2
#define USE_GEN_LOSTLEPTON
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

#include "../../common/Reader_CommonFormat.h"
babyEvent myEvent;

#include "../../common/TFFactory.h"
#include "../../Selection/test.h"
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

const double pi = boost::math::constants::pi<double>();
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

uint32_t nrLooseB(){
     uint32_t nbtag = 0;
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
     	if(myEvent.ak4pfjets_CSV[i]>=0.46) nbtag++;
     }
   return nbtag;
}

float differenceBpt(){
     uint32_t nbtag = 0;
     vector<float> bpt;
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) nbtag++;
     }
    if(nbtag >=2)
    {
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){ //@MJ@ TODO do this better and with the jets with highest CSVV2
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) bpt.push_back(myEvent.ak4pfjets_pt[i]);
     }
    
   return abs(bpt.at(0) - bpt.at(1));
    }
   return -999;
}

float dphiB(){
     uint32_t nbtag = 0;
     vector<float> bphi;
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) nbtag++;
     }
    if(nbtag >=2)
    {
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){ //@MJ@ TODO do this better and with the jets with highest CSVV2
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) bphi.push_back(myEvent.ak4pfjets_phi[i]);
     }

    
   return abs(bphi.at(0) - bphi.at(1)) > pi ? abs(bphi.at(0) - bphi.at(1)) - pi: abs(bphi.at(0) - bphi.at(1)) ;
    }
    return -999;
}


float dphiBl(){
     uint32_t nbtag = 0;
     vector<float> bphi;
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) nbtag++;
     }
    if(nbtag >=2)
    {
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){ //@MJ@ TODO do this better and with the jets with highest CSVV2
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) bphi.push_back(myEvent.ak4pfjets_phi[i]);
     }
    
   return abs(bphi.at(0) - myEvent.lep1_phi) > pi ? abs(bphi.at(0) - myEvent.lep1_phi) - pi:abs(bphi.at(0) - myEvent.lep1_phi) ;
    }
    return -999;
}


float dphiBMET(){
     uint32_t nbtag = 0;
     vector<float> bphi;
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) nbtag++;
     }
    if(nbtag >=2)
    {
     for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){ //@MJ@ TODO do this better and with the jets with highest CSVV2
     	if(myEvent.ak4pfjets_CSV[i]>=0.8) bphi.push_back(myEvent.ak4pfjets_phi[i]);
     }
    
   return abs(bphi.at(0) - myEvent.pfmet_phi) > pi ? abs(bphi.at(0) - myEvent.pfmet_phi) -pi :abs(bphi.at(0) - myEvent.pfmet_phi) ;
    }
    return -999;
}


int nMatched2LooseB = 0;
bool SortByCSV(pair<float,float> i, pair<float,float> j) {return (i.first>j.first);}
bool SortByPT(pair<float,float> i, pair<float,float> j) {return (i.second>j.second);}

float jetPtAssymetry()
{
	vector<pair<float,int> > jets;
	for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
		float CSV = myEvent.ak4pfjets_CSV[i];
		float pt = myEvent.ak4pfjets_pt[i];
		pair<float,float> pair;
		pair.first = CSV;
		pair.second = pt;
		jets.push_back(pair);
        }

	std::sort(jets.begin(), jets.end(), SortByCSV);

	if(jets.size()>=4){
            float num = abs( (jets[0].second + jets[1].second) - (jets[2].second + jets[3].second));
            float den = jets[0].second + jets[1].second;
            return num/den;
        }

       return -999;

}

bool descendingSort (int i,int j) { return (i>j); }

float jetPtAssymetryPT()
{

	std::sort(myEvent.ak4pfjets_pt.begin(), myEvent.ak4pfjets_pt.end(), descendingSort);

	if( myEvent.ak4pfjets_pt.size()>3){
            float num = abs( (myEvent.ak4pfjets_pt[0] + myEvent.ak4pfjets_pt[1]) - (myEvent.ak4pfjets_pt[2] + myEvent.ak4pfjets_pt[3]));
            float den = myEvent.ak4pfjets_pt[0] + myEvent.ak4pfjets_pt[1];
            return num/den;
        }

       return -999;

}

float jetPtAssymetryBPT()
{
	vector<pair<float,int> > jets;
	for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
		float CSV = myEvent.ak4pfjets_CSV[i];
		float pt = myEvent.ak4pfjets_pt[i];
		pair<float,float> pair;
		pair.first = CSV;
		pair.second = pt;
		jets.push_back(pair);
        }

	std::sort(jets.begin(), jets.end(), SortByCSV);

	if(jets.size()>=4){
            float nump1 = (jets[0].second + jets[1].second);
            float den = jets[0].second + jets[1].second;
            jets.erase(jets.begin(), jets.begin()+2);
	    std::sort(jets.begin(), jets.end(), SortByPT);
            float nump2 =  (jets[0].second + jets[1].second);
            return abs(nump1-nump2)/den;
        }

       return -999;

}
float ptRMS()
{
     
    float ptsum =  0;
    if (myEvent.ak4pfjets_pt.size() == 0) return -999;
    for( uint32_t p =0; p< myEvent.ak4pfjets_pt.size(); p++)
    {
        ptsum += (myEvent.ak4pfjets_pt.at(p) * myEvent.ak4pfjets_pt.at(p));
    }
    float rms = sqrt(ptsum/ myEvent.ak4pfjets_pt.size());
    return rms;
}

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
		jets.push_back(pair);
        }
	std::sort(jets.begin(), jets.end(), SortByCSV);
	nMatched2LooseB = 0;
	if(jets.size()>=2){
		cout<<jets[0].first<<" "<<jets[1].first<<endl;
		if(fabs(jets[0].second)==5) nMatched2LooseB++;
		if(fabs(jets[1].second)==5) nMatched2LooseB++;
	}
}

float jetEtaAssymetry()
{
	vector<pair<float,int> > jets;
	for(unsigned int i=0;i<myEvent.ak4pfjets_CSV.size();i++){
		float CSV = myEvent.ak4pfjets_CSV[i];
		float eta = myEvent.ak4pfjets_eta[i];
		pair<float,float> pair;
		pair.first = CSV;
		pair.second = eta;
		jets.push_back(pair);
        }

	std::sort(jets.begin(), jets.end(), SortByCSV);

	if(jets.size()>=4){
            float num = abs( (jets[0].second + jets[1].second) - (jets[2].second + jets[3].second));
            float den = jets[0].second + jets[1].second;
            return num/den;
        }

       return -999;

}

float jetEtaAssymetryPT()
{
	vector<pair<float,int> > jets;
	for(unsigned int i=0;i<myEvent.ak4pfjets_pt.size();i++){
		float pt = myEvent.ak4pfjets_pt[i];
		float eta = myEvent.ak4pfjets_eta[i];
		pair<float,float> pair;
		pair.first = pt;
		pair.second = eta;
		jets.push_back(pair);
        }

	std::sort(jets.begin(), jets.end(), SortByCSV);


	if( myEvent.ak4pfjets_pt.size()>3){
            float num = abs( (jets[0].second + jets[1].second) - (jets[2].second + jets[3].second));
            float den = jets[0].second + jets[1].second;
            return num/den;
        }

       return -999;

}


float etaOfSystem()
{
    TLorentzVector lv;
    lv.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
    //cout << "eta1:  " << lv.Eta() << " lepeta: " <<  myEvent.lep1_eta << endl;
    for(uint32_t j = 0; j< myEvent.ak4pfjets_pt.size(); j++)
    {
        //cout << myEvent.ak4pfjets_pt.size() << " " << myEvent.ak4pfjets_eta.size() << " " << myEvent.ak4pfjets_phi.size() << " " <<  myEvent.ak4pfjets_mass.size() << endl;
        TLorentzVector jv;
        jv.SetPtEtaPhiM(myEvent.ak4pfjets_pt.at(j), myEvent.ak4pfjets_eta.at(j), myEvent.ak4pfjets_phi.at(j), myEvent.ak4pfjets_mass.at(j));
        lv += jv;
    }
   //cout << "eta2:  " << lv.Eta() << endl;
   return lv.Eta();
//return 0;
}

float dphiMETl()
{
    float dphi = abs(myEvent.pfmet_phi - myEvent.lep1_phi) > pi ? abs(myEvent.pfmet_phi - myEvent.lep1_phi) - pi : abs(myEvent.pfmet_phi - myEvent.lep1_phi);
    float dphimodpi2 = abs(dphi > (pi/2) ? dphi - (pi/2) : dphi);
    return dphimodpi2;

}
/*
float getJetPt(uint32_t jetNr) //already defined, use for recomputation of ngood jets
{
    if(ak4pfjets_pt.size() >= jetNr)
        return ak4pfjets_pt.at(jetNr - 1);
    return -999;
}
*/
int looseB = -1;
float ptDiffB = -999;
float phiDiffB = -999;
float phiDiffBl = -999;
float phiDiffBMET = -999;
float lepPtToMET = -999;
float relLepPtToMET = -999;
float ptAssB = -999;
float ptAssPT = -999;
float ptAssBPT = -999;
float ptrms = -999;
float phiDiffMETlep = -999;    
float etaAssB = -999;
float etaAssPT = -999;
float etaSyst = -999;
float jet1Pt = -999;
float jet2Pt = -999;
float jet3Pt = -999;
float jet4Pt = -999;
float jet1PtToHT = -999;
float jet2PtToHT = -999;
float jet3PtToHT = -999;
float jet4PtToHT = -999;
float jet1PtToMET = -999;
float jet2PtToMET = -999;
float jet3PtToMET = -999;
float jet4PtToMET = -999;
float pt1PlusPt2 = -999;
float pt1MinPt2 = -999;
float pt12Ass = -999;
float METAndHT = -999;
float METAndHTAndLPT = -999;
float mlb = -999;
TH2D *correl = new TH2D("MTcor", "MTCor", 100, 20, 1000, 100, 0,7);

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
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/BabyTuples/16_09/";
    totalNumberOfWorkers = 6;


    //@EC@: to be done: put default files in the code itself or load info from a external config file
    string WPath = "../../Tools/Weighting/files/";
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
    // fixed binning
    AddVariable("MET", "MET",  "MET", 20 ,200,1000,  &(myEvent.pfmet), "");
    AddVariable("MT2W", "MT2W",  "MT2W", 20 ,0,500,  &(myEvent.MT2W), "");
    AddVariable("Mlb", "Mlb",  "Mlb", 20 ,0,500,  &(mlb ), "");
    AddVariable("MT", "MT",  "MT", 20 ,100,1000,  &(myEvent.mt_met_lep), "");
    AddVariable("topness","topness","topness",20,-20,20,&(myEvent.topness),"");
    AddVariable("ngoodbtags","ngoodbtags","ngoodbtags",10,0,10,&(myEvent.ngoodbtags),"");
    AddVariable("looseB","looseB","looseB",10,0,10,&(looseB),"noUnderflowInFirstBin");
    AddVariable("ptDiffB","ptDiffB","ptDiffB", 20, 0, 700,&(ptDiffB),"noUnderflowInFirstBin");
    AddVariable("phiDiffB","phiDiffB","phiDiffB", 20, 0, 4,&(phiDiffB),"noUnderflowInFirstBin");
    AddVariable("phiDiffMETlep","phiDiffMETlep","phiDiffMETlep", 20, 0, 4,&(phiDiffMETlep),"noUnderflowInFirstBin");
    AddVariable("lepPtToMET","lepPtToMET","lepPtToMET",20, 0,1.4,&(lepPtToMET),"noUnderflowInFirstBin");
    AddVariable("relLepPtToMET","relLepPtToMET","relLepPtToMET",20, 0,1,&(relLepPtToMET),"noUnderflowInFirstBin");
    AddVariable("ptrms","ptrms","ptrms",20, 0,500,&(ptrms),"noUnderflowInFirstBin");
    AddVariable("HT","HT","HT", 20, 0, 2000,&(myEvent.HT),"noUnderflowInFirstBin"); //@MJ@ TODO correct HT?!
    AddVariable("ptAssB","ptAssB","ptAssB",20, 0,6,&(ptAssB),"noUnderflowInFirstBin");
    AddVariable("ptAssPT","ptAssPT","ptAssPT",20, 0,1.5,&(ptAssPT),"noUnderflowInFirstBin");
    AddVariable("ptAssBPT","ptAssBPT","ptAssBPT",20, 0,10,&(ptAssBPT),"noUnderflowInFirstBin");
    AddVariable("etaAssB","etaAssB","etaAssB",20, 0,100,&(etaAssB),"noUnderflowInFirstBin");
    AddVariable("etaAssPT","etaAssPT","etaAssPT",20, 0,100,&(etaAssPT),"noUnderflowInFirstBin");
    AddVariable("etaSyst","etaSyst","etaSyst",20, 0,4,&(etaSyst),"noUnderflowInFirstBin");
    AddVariable("phiDiffBl","phiDiffBl","phiDiffBl", 20, 0, 4,&(phiDiffBl),"noUnderflowInFirstBin");
    AddVariable("phiDiffBMET","phiDiffBMET","phiDiffBMET", 20, 0, 4,&(phiDiffBMET),"noUnderflowInFirstBin");
    AddVariable("jet1Pt","jet1Pt","", 50, 0, 1000,&(jet1Pt),"noUnderflowInFirstBin");
    AddVariable("jet2Pt","jet2Pt","", 50, 0, 1000,&(jet2Pt),"noUnderflowInFirstBin");
    AddVariable("jet3Pt","jet3Pt","", 50, 0, 1000,&(jet3Pt),"noUnderflowInFirstBin");
    AddVariable("jet4Pt","jet4Pt","", 50, 0, 1000,&(jet4Pt),"noUnderflowInFirstBin");
    AddVariable("jet1PtToHT","jet1PtToHT","", 50, 0, 1,&(jet1PtToHT),"noUnderflowInFirstBin");
    AddVariable("jet2PtToHT","jet2PtToHT","", 50, 0, 1,&(jet2PtToHT),"noUnderflowInFirstBin");
    AddVariable("jet3PtToHT","jet3PtToHT","", 50, 0, 1,&(jet3PtToHT),"noUnderflowInFirstBin");
    AddVariable("jet4PtToHT","jet4PtToHT","", 50, 0, 1,&(jet4PtToHT),"noUnderflowInFirstBin");
    AddVariable("jet1PtToMET","jet1PtToMET","", 50, 0, 3,&(jet1PtToMET),"noUnderflowInFirstBin");
    AddVariable("jet2PtToMET","jet2PtToMET","", 50, 0, 3,&(jet2PtToMET),"noUnderflowInFirstBin");
    AddVariable("jet3PtToMET","jet3PtToMET","", 50, 0, 3,&(jet3PtToMET),"noUnderflowInFirstBin");
    AddVariable("jet4PtToMET","jet4PtToMET","", 50, 0, 3,&(jet4PtToMET),"noUnderflowInFirstBin");
    AddVariable("pt1PlusPt2","pt1PlusPt2","", 50, 0, 2000,&(pt1PlusPt2),"noUnderflowInFirstBin");
    AddVariable("pt1MinPt2","pt1MinPt2","", 50, 0, 1000,&(pt1MinPt2),"");
    AddVariable("pt12Ass","pt12Ass","", 20, 0,1,&(pt12Ass),"noUnderflowInFirstBin");
    AddVariable("METAndHT","METAndHT","", 50, 0,3000,&(METAndHT),"noUnderflowInFirstBin"); //MET+HT
    AddVariable("METAndHTAndLPT","METAndHTAndLPT","", 50, 0,3000,&(METAndHTAndLPT),"noUnderflowInFirstBin"); //MET+HT+PT of lep
    AddVariable("M3b", "M3b",  "M3b", 50 ,0,1000,  &(myEvent.M3b ), "noUnderflowInFirstBin");
    //AddVariable("","","",100,-20,20,&(myEvent.),"");

    // ------------------
    // Datasets
    // ------------------
    AddProcessClass("rare", "rare", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("ttZ","rare",0,0.7826);
    	AddDataset("tZq","rare",0,0.0758);
    	AddDataset("ZZ","rare",0,0.564);
   	AddDataset("WZ","rare",0,3.06);

    AddProcessClass("throw", "throw", "signal", kBlue);
     	//AddDataset("T2tt_400to1200", "throw", 0, 0 );
     	AddDataset("T2tt_400to1200", "throw", 0, 0 );
    //
    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    //TString fNameTmp =  babyTuplePath+"T2tt_400to1200.root";
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
                if( (htmp->GetXaxis()->GetBinLowEdge(bx+1) == 400 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 300 ) ||  ( htmp->GetXaxis()->GetBinLowEdge(bx+1) == 650 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 350 ) ||( htmp->GetXaxis()->GetBinLowEdge(bx+1) == 900 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 50 )) //@MJ@ TODO avoid too many regions
                {
                    std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinLowEdge(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinLowEdge(by+1) << " content: " << htmp->GetBinContent(bx+1,by+1) << std::endl;
                    pair<uint32_t, uint32_t> key = make_pair( htmp->GetXaxis()->GetBinLowEdge(bx+1), htmp->GetYaxis()->GetBinLowEdge(by+1));
                    string stops = to_string(htmp->GetXaxis()->GetBinLowEdge(bx+1));
                    string neutrs = to_string( htmp->GetYaxis()->GetBinLowEdge(by+1));
                    scanMap[key] = stops+"_"+neutrs;
                    AddProcessClass( stops+"_"+neutrs, stops+"_"+neutrs, "signal", kViolet+citer++ );
                }
            }
        }

    }

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

    //AddProcessClass("data", "data", "data", kViolet);
      //  AddDataset("SE_0", "data", 0, 0 );
      //  AddDataset("SE_1", "data", 0, 0 );
      //  AddDataset("SE_C", "data", 0, 0 );
      //  AddDataset("SE_D", "data", 0, 0 );
      //  AddDataset("SM_0", "data", 0, 0 );
      //  AddDataset("SM_1", "data", 0, 0 );
      //  AddDataset("SM_C", "data", 0, 0 );
      //  AddDataset("SM_D", "data", 0, 0 );
      //AddDataset("MET_0", "data", 0, 0 );
      //AddDataset("MET_1", "data", 0, 0 );
  
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
AddRegion("SR1lbaselineOnly_MET250toInf","SR1lbaselineOnly_MET250toInf",&SR1lbaselineOnly_MET250toInf);
AddRegion("SR1lbaselineOnlyLeadJetsRedifined_MET250toInf","SR1lbaselineOnlyLeadJetsRedifined_MET250toInf",&SR1lbaselineOnlyLeadJetsRedifined_MET250toInf);
AddRegion("SR1lbaselineOnlyGoodJetsRedifined4j_MET250toInf","SR1lbaselineOnlyGoodJetsRedifined4j_MET250toInf",&SR1lbaselineOnlyGoodJetsRedifined4j_MET250toInf);
AddRegion("SR1l23jLowMlb_MET250toInf","SR1l23jLowMlb_MET250toInf",&SR1l23jLowMlb_MET250toInf);
AddRegion("SR1l23jHighMlb_MET250toInf","SR1l23jHighMlb_MET250toInf",&SR1l23jHighMlb_MET250toInf);
AddRegion("SR1l4jLowMT2WLowMlb_MET250toInf","SR1l4jLowMT2WLowMlb_MET250toInf",&SR1l4jLowMT2WLowMlb_MET250toInf);
AddRegion("SR1l4jLowMT2WHighMlb_MET250toInf","SR1l4jLowMT2WHighMlb_MET250toInf",&SR1l4jLowMT2WHighMlb_MET250toInf);
AddRegion("SR1l4jMidMT2WLowMlb_MET250toInf","SR1l4jMidMT2WLowMlb_MET250toInf",&SR1l4jMidMT2WLowMlb_MET250toInf);
AddRegion("SR1l4jMidMT2WHighMlb_MET250toInf","SR1l4jMidMT2WHighMlb_MET250toInf",&SR1l4jMidMT2WHighMlb_MET250toInf);
AddRegion("SR1l4jHighMT2WLowMlb_MET250toInf","SR1l4jHighMT2WLowMlb_MET250toInf",&SR1l4jHighMT2WLowMlb_MET250toInf);
AddRegion("SR1l4jHighMT2WHighMlb_MET250toInf","SR1l4jHighMT2WHighMlb_MET250toInf",&SR1l4jHighMT2WHighMlb_MET250toInf);

 
    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(12900.);

    Create1DHistos();
    Add2DHisto("MT", "phiDiffMETlep");
    Add2DHisto("MET", "HT");
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
    if (currentProcessType == "signal" && ( (myEvent.gen_stop_m.at(0) == 400 && myEvent.gen_neutralino_m.at(0) ==300) || (myEvent.gen_stop_m.at(0) == 650 && myEvent.gen_neutralino_m.at(0) ==350) ||( myEvent.gen_stop_m.at(0) == 900 && myEvent.gen_neutralino_m.at(0) == 50) ))
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
        if(true)
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
        //else if(stopMass-neutralinoMass == 175)
        //{
        //    currentProcessClass = "grouped";
        //}
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
    else if (currentProcessClass == "test" && (myEvent.numberOfGeneratedLeptons < 2))
    {
        currentProcessClass = "test";
    }
    else
    {}

    //
    //recompute(useTriggerInfo, currentDataset);

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent; //@MJ@ TODO cross section form file?!

    //--- Compute Weights  ---//
    WFact.BtagWeighComputor (myEvent.ak4pfjets_pt, myEvent.ak4pfjets_eta, myEvent.ak4pfjets_hadronFlavour, myEvent.ak4pfjets_CSV);// should be called once per event
    double btagWeight = WFact.GetBtagW();
    if(btagWeight==0) cout<<"btagWeight = "<<btagWeight<<endl;

    double lepWeight = 0;
    if(currentProcessType != "data"){
        WFact.LeptonWeightComputor(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_pdgid, myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_pdgid, myEvent.nvetoleps, myEvent.numberOfSelectedLeptons, myEvent.genLostLeptons_pt.size(), myEvent.genLostLeptons_pt, myEvent.genLostLeptons_eta, myEvent.genLostLeptons_pdgid );
        lepWeight = WFact.GetLepW();
    }

    ///new vars
    looseB = nrLooseB();
    ptDiffB = differenceBpt();
    phiDiffB = dphiB();
    phiDiffBl = dphiBl();
    phiDiffBMET = dphiBMET();
    lepPtToMET = myEvent.lep1_pt / myEvent.pfmet;
    relLepPtToMET =(abs (myEvent.lep1_pt - myEvent.pfmet)) / myEvent.pfmet;
    ptAssB = jetPtAssymetry();
    ptAssPT = jetPtAssymetryPT();
    ptAssBPT = jetPtAssymetryBPT();
    etaAssB = jetEtaAssymetry();
    etaAssPT = jetEtaAssymetryPT();
    etaSyst = etaOfSystem();
    ptrms = ptRMS();
    phiDiffMETlep = dphiMETl();     
    jet1Pt = getJetPt(1);
    jet2Pt = getJetPt(2);
    jet3Pt = getJetPt(3);
    jet4Pt = getJetPt(4);
    jet1PtToHT = getJetPt(1)/myEvent.HT;
    jet2PtToHT = getJetPt(2)/myEvent.HT;
    jet3PtToHT = getJetPt(3)/myEvent.HT;
    jet4PtToHT = getJetPt(4)/myEvent.HT;
    jet1PtToMET = getJetPt(1)/myEvent.pfmet;
    jet2PtToMET = getJetPt(2)/myEvent.pfmet;
    jet3PtToMET = getJetPt(3)/myEvent.pfmet;
    jet4PtToMET = getJetPt(4)/myEvent.pfmet;
    pt1PlusPt2 = getJetPt(1)+ getJetPt(2);
    pt1MinPt2 = getJetPt(1)- getJetPt(2);
    pt12Ass = (getJetPt(1)- getJetPt(2)) / (getJetPt(1)+ getJetPt(2));
    METAndHT = myEvent.pfmet + myEvent.HT;
    METAndHTAndLPT = myEvent.pfmet + myEvent.HT + myEvent.lep1_pt ;
    mlb = minMlb();
    float weight     = weightLumi*btagWeight*lepWeight;
    if (currentProcessType == "data") weight = 1.0;
    if (currentProcessType == "signal") weight = weightSignal*btagWeight*lepWeight;
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
    SchedulePlots("1DFigureOfMerit","var=phiDiffB,cutType=keepLowValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=HT,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=phiDiffMETlep,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=jet1Pt,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=jet2Pt,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=jet3Pt,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=jet4Pt,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=looseB,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=METAndHT,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=METAndHTAndLPT,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("1DFigureOfMerit","var=M3b,cutType=keepHighValues,type=Zbi,minBackground=0");
    SchedulePlots("2D");

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
    WritePlots("./plotsTestMy/");

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

vector<string> signalReg = { "SR1lbaselineOnly_MET250toInf" , "SR1l23jLowMlb_MET250toInf" , "SR1l23jHighMlb_MET250toInf" , "SR1l4jLowMT2WLowMlb_MET250toInf" , "SR1l4jLowMT2WHighMlb_MET250toInf" , "SR1l4jMidMT2WLowMlb_MET250toInf" , "SR1l4jMidMT2WHighMlb_MET250toInf" , "SR1l4jHighMT2WLowMlb_MET250toInf" , "SR1l4jHighMT2WHighMlb_MET250toInf" ,  };
    TableDataMC(this, signalReg,"lepChannel", "includeSignal" ).Print("signalReg2.tab", 4);
    TableDataMC(this, signalReg,"lepChannel", "includeSignal" ).PrintLatex("signalReg2.tex", 4);
    ofstream sigfile(outFold+"/signalReg2.txt");
    if (sigfile.is_open())
    {
        for(uint32_t k = 0; k<signalReg.size(); k++)
        {
	    sigfile << signalReg.at(k)  << endl ;
        }
        sigfile.close();
    }
    
/*
    TableDataMC(this, yield,"lepChannel" ).Print("yield.tab", 4);
    TableDataMC(this, yield,"lepChannel" ).PrintLatex("yield.tex", 4);
*/
   
    vector<string> tfreg = {};
    TFProducer prod(tfreg, "lostLepton", "CR2l");
    prod.produceTFTable("yield.tab", "lostLeptonTF");
    TFProducer prod2(tfreg, "singleLepton", "CR1l");
    prod2.produceTFTable("yield.tab", "singleLeptonTF");
    TFProducer prod3(tfreg, "singleLeptonFromT", "CR1l");
    prod3.produceTFTable("yield.tab", "singleLeptonFromTTF");


    cout << "one lepton correlation: " << correl->GetCorrelationFactor() << endl;
    cout << "end of processing" << endl;
 }
