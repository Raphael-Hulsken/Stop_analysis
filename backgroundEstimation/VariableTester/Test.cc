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


bool SR1l_2lb() {return SR1l() && Pass2LooseB();}
bool SR1l_2j_2lb() {return SR1l_2j() && Pass2LooseB();}
bool SR1l_3j_2lb() {return SR1l_3j() && Pass2LooseB();}
bool SR1l_4j_2lb() {return SR1l_4j() && Pass2LooseB();}
bool CR1l_2lb() {return CR1l() && Pass2LooseB();}


//new variable
float mStop = 0;
float mNeutralino = 0;
float leadjet_pt = 0;
float Meff = 0;
float deta_lep = 0;
float dphi_lep = 0;
int dphi_deta = 0;

float super_dphi_lep = 0;
float super_deta_lep = 0;
float super_dphi_deta = 0;
int category = 0;




int ComputeDphiDeta(float dphi, float deta){
  int indice_phi;
  int indice_eta;
  int nbins = 5;
  for(int i=0;i<nbins;i++){
    if(dphi>=TMath::Pi()/nbins*i && dphi<TMath::Pi()/nbins*(i+1)){
      indice_phi = i;
      break;
     }
  }
  for(int i=0;i<nbins;i++){
    if(deta>=TMath::Pi()/nbins*i && deta<TMath::Pi()/nbins*(i+1)){
      indice_eta = i;
      break;
     }
  }
  return nbins*indice_phi+indice_eta;

}

int ComputeDphiDetaV(float dphi, float deta){
  int indice_phi;
  int indice_eta;
  vector<double> bins_phi = {0,TMath::Pi()*0.2,TMath::Pi()*0.4,TMath::Pi()*0.6,TMath::Pi()*0.8,TMath::Pi()};
  vector<double> bins_eta = {0,0.4,0.8,1.2,2.0,100};
  for(int i=0;i<bins_phi.size()-1;i++){
    if(dphi>=bins_phi[i]&& dphi<bins_phi[i+1]){
      indice_phi = i;
      break;
     }
  }
  for(int i=0;i<bins_eta.size();i++){
    if(deta>=bins_eta[i]&& deta<bins_eta[i+1]){
      indice_eta = i;
      break;
     }
  }
  return bins_eta.size()*indice_phi+indice_eta;

}


bool SR_2l(){
  
  // compute mll
  TLorentzVector l1;
  l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
  TLorentzVector l2;
  l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
  TLorentzVector ll = l1+l2;
  double mll = ll.M();
  //cout<<"mll = "<<mll<<endl;

  dphi_lep = fabs(l1.DeltaPhi(l2));
   
  if( myEvent.ngoodleps == 2 &&  ( (abs(myEvent.lep1_pdgid) == 11 && abs(myEvent.lep2_pdgid) == 13) || (abs(myEvent.lep1_pdgid) == 13 && abs(myEvent.lep2_pdgid) == 11) )   
     && mll > 20 
     && myEvent.ngoodjets >=2 && myEvent.ngoodbtags>=2 && myEvent.pfmet>0){
  //if (myEvent.ngoodleps == 2 && myEvent.ngoodjets >=2 && myEvent.ngoodbtags>=2 && myEvent.pfmet>30){
  //if (myEvent.ngoodleps == 2 && myEvent.ngoodjets >=2 && myEvent.ngoodbtags>=2 && myEvent.pfmet>0){
  	return true;
  }
  return false;
}


bool SR_2l_2j() { return SR_2l() && myEvent.ngoodjets == 2;}
bool SR_2l_3j() { return SR_2l() && myEvent.ngoodjets == 3 && myEvent.ak4pfjets_pt[0]<300 ;}
bool SR_2l_4j() { return SR_2l() && myEvent.ngoodjets >= 4 && myEvent.ak4pfjets_pt[0]<300 ;}
bool SR_2l_ISR() { return SR_2l() && myEvent.ngoodjets >= 3 && myEvent.ak4pfjets_pt[0]>=300;}


bool SR_2l_2j_lowMass(){
	return SR_2l_2j() && Meff<400;
}
bool SR_2l_lowMass(){
	return SR_2l() && Meff<400;
}

bool SR_2l_3j_incl() { return SR_2l() && myEvent.ngoodjets >= 3;}



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


bool lepChannel() 
{ 
    return true; 
}


bool SuggestedSelection(){
	if (myEvent.pfmet < 250) return false;
	if (myEvent.ak4_HT < 250) return false;
	if (myEvent.ngoodjets < 2) return false;
	if (myEvent.ngoodleps != 1) return false; // >=1 lepton
	if (myEvent.lep1_pt < 25) return false;
	return true;
	//if (myEvent.
	// lepton requirement
	// trigger requirement 
}

bool Sel0b(){
	return (SuggestedSelection() && myEvent.ngoodbtags==0);
}

bool Sel1b(){
	return (SuggestedSelection() && myEvent.ngoodbtags>=1);
}

bool Sel2b(){
	return (SuggestedSelection() && myEvent.ngoodbtags>=2);
}


//Add this as a global variable 
WeightFactory WFact;

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    //babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/0909/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/BabyTuples/14_09/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/BabyTuples/16_09/";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/BabyTuples/25_11/";
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11";
    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2016/Synchro/CMSSW_8_0_5/src/store/babyTuples/TriggerStudy/";
    totalNumberOfWorkers = 5;


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
    vector<float> METBins = {0,50,100,150,200,250,500};
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
    AddVariable("MET", "MET",  "GeV", 25 ,0,500,  &(myEvent.pfmet), "logY");
    AddVariable("MET2", "MET",  "GeV", (int) (METBins.size()-1), METBins.data(),  &(myEvent.pfmet), "logY");
    AddVariable("topness","topness","",20,-20,20,&(myEvent.topness),"logY");
    AddVariable("nJets","nJets","",8,1,8,&(myEvent.ngoodjets),"logY");
    AddVariable("nBJets","nBJets","", 4, 1, 4, &(myEvent.ngoodbtags),"logY");
    AddVariable("ak4_HT", "ak4_HT", "GeV",25,0,500,&(myEvent.ak4_HT),"logY");
    AddVariable("lep_pt","p_{T}(lep)", "GeV",50, 0,500,&(myEvent.lep1_pt),"logY");
    AddVariable("lep_eta","eta (lep)", "",50, -2.5,2.5,&(myEvent.lep1_eta),"logY");
    AddVariable("deta_lep","#Delta#eta(l,l)","",25,0,TMath::Pi(),&deta_lep,"");
    AddVariable("dphi_lep","#Delta#phi(l,l)","",25,0,TMath::Pi(),&dphi_lep,"");
    //AddVariable("deta_lep","#Delta#eta(l,l)","",10,0,TMath::Pi(),&deta_lep,"");
    //AddVariable("dphi_lep","#Delta#phi(l,l)","",10,0,TMath::Pi(),&dphi_lep,"");
    AddVariable("super_dphi_lep","#Delta#phi(l,l)","",100,0,TMath::Pi()*4,&super_dphi_lep,"");
    AddVariable("super_deta_lep","#Delta#eta(l,l)","",100,0,TMath::Pi()*4,&super_deta_lep,"");
    AddVariable("category", "category","",4,0,3,&category,""); 
    AddVariable("dphi_deta","","",25,0,25,&dphi_deta,"");
    AddVariable("super_dphi_deta","","",100,0,100,&super_dphi_deta,"");
    AddVariable("leadjet_pt","leading jet p_{T}", "GeV",10,0,500,&leadjet_pt,"logY");
    AddVariable("Meff","Mass effective", "GeV",50,0,1000,&Meff,"logY");

    AddVariable("StopMass","stop mass", "GeV",24,140,260, &mStop);
    AddVariable("NeutralinoMass","lsp mass", "GeV",20,0,100, &mNeutralino);


    // ------------------
    // Datasets
    // ------------------
    
    /*
    AddProcessClass("throw", "throw", "signal", kBlue);
     	AddDataset("T2tt_400to1200", "throw", 0, 0 );
    */
    //
    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    TString fNameTmp =  babyTuplePath+"T2tt.root";
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
                if( (htmp->GetXaxis()->GetBinLowEdge(bx+1) == 425 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 250) ||
                    (htmp->GetXaxis()->GetBinLowEdge(bx+1) == 500 && htmp->GetYaxis()->GetBinLowEdge(by+1) == 1)  ||
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

   //*/
   
  /*

    AddProcessClass("data", "data", "data", COLORPLOT_BLACK); //kViolet);
    */
    /*
	AddDataset("SE_0", "data", 0, 0 );
    	AddDataset("SE_1", "data", 0, 0 );
    	AddDataset("SE_C", "data", 0, 0 );
    	AddDataset("SE_D", "data", 0, 0 );
        AddDataset("SM_0", "data", 0, 0 );
        AddDataset("SM_1", "data", 0, 0 );
        AddDataset("SM_C", "data", 0, 0 );
        AddDataset("SM_D", "data", 0, 0 );
        */
/*	AddDataset("MET_0", "data", 0, 0 );
        AddDataset("MET_1", "data", 0, 0 );
        AddDataset("MET_C", "data", 0, 0 );
        AddDataset("MET_D", "data", 0, 0 );
 */   
    //AddProcessClass("test", "other", "background",kBlack);
    
    AddProcessClass("ttbar", "t#bar{t}", "background",kBlue);
	
    	//AddDataset("TTJetsDLv0v4","ttbar",0, 57.35*1.5225);
    	AddDataset("ttbar_diLept_madgraph_pythia8_25ns","ttbar",0, 57.35*1.5225);
	/*
	AddDataset("TTJetsSLtop", "ttbar", 0, 114.6*1.594 );
    	AddDataset("TTJetsSLatopv1","ttbar",0,114.6*1.594);
    	//AddDataset("TTJetsSLatopext","ttbar",0,114.6*1.594);
    	AddDataset("TTJetsDLv0v4","ttbar",0, 57.35*1.5225);
    	*/
	//AddDataset("TTJetsDLext","ttbar",0, 57.35*1.5225);
	
	/*
	//AddDataset("WJetsToLNuTune","test",0,60781.5*1.01);
    	AddDataset("W1JetsToLNuTune","test",0, 9493*1.238);
    	AddDataset("W2JetsToLNuTune","test",0, 3120*1.231);
    	AddDataset("W3JetsToLNuTune","test",0, 942.3*1.231);
    	AddDataset("W4JetsToLNuTune","test",0, 524.2*1.114);
        */

	/*
    AddProcessClass("Wjets", "W+jets", "background",kGreen+1);
	AddDataset("WJetsToLNuHT100To200ext","Wjets",0,1345*1.21);
	AddDataset("WJetsToLNuHT200To400v3v2","Wjets",0,359.7*1.21);
	//AddDataset("WJetsToLNuHT200To400ext","Wjets",0,359.7*1.21);
	AddDataset("WJetsToLNuHT400To600ext","Wjets",0,48.91*1.21);
	AddDataset("WJetsToLNuHT600To800v3v2","Wjets",0,12.1*1.21);
	AddDataset("WJetsToLNuHT800To1200v3v1","Wjets",0,5.50*1.21);
	AddDataset("WJetsToLNuHT1200To2500v3v2","Wjets",0,1.33*1.21);
	//AddDataset("WJetsToLNuHT1200To2500ext","Wjets",0,1.33*1.21);
	AddDataset("WJetsToLNuHT2500ToInfv3v1","Wjets",0,0.032*1.21);
	*/
	/*
    AddProcessClass("singleTop", "single top", "background",kRed);
	AddDataset("ST_s","singleTop",0,10.11*0.364176);
    	AddDataset("ST_tW_top","singleTop",0,38.09*0.5135);
    	AddDataset("ST_tW_atop","singleTop",0,38.09*0.5135);
    	//AddDataset("ST_t","singleTop",0,80.95*0.324);
    	AddDataset("ST_t","singleTop",0,70.71);
    
*/
/*	
    AddProcessClass("other", "other", "background",kBlack);
	AddDataset("TTWtoQQ","other",0,0.4062);
    	AddDataset("TTWtoLNu","other",0,0.2043);
    	AddDataset("TTT","other",0,1.0);
    	AddDataset("VV","other",0,12.05*0.9917);
  */
/*
    AddProcessClass("rare", "rare", "background", kBlue);//@MJ@ TODO K-factor?
    	AddDataset("ttZ","other",0,0.7826);
    	//AddDataset("tZq","rare",0,0.0758);
    	//AddDataset("ZZ","rare",0,0.564);
    	//AddDataset("WZ","rare",0,3.06);
	AddDataset("TTWtoQQ","rare",0,0.4062);
    	AddDataset("TTWtoLNu","rare",0,0.2043);

  */  
    //AddProcessClass("lostLepton", "lost Lepton", "background", kPink);
    //AddProcessClass("singleLepton", "1l", "background", kGreen);
    //AddProcessClass("singleLeptonFromT", "1l from top", "background", kOrange);

    AddProcessClass("T2tt", "T2tt", "signal", kRed);//@MJ@ TODO K-factor?
	//AddDataset("T2tt","T2tt", 1.);
	//AddDataset("T2tt_250to350","T2tt", 1.);
	AddDataset("T2tt_mStop_425_mLSP_325_25ns","T2tt", 1.);

    AddProcessClass("175-0", "175-0","signal", kRed-4);
    AddProcessClass("200-25", "200-25","signal", kRed-2);
    AddProcessClass("225-50", "225-50","signal", kRed);
    AddProcessClass("250-75", "250-75","signal", kRed+2);

    AddProcessClass("stop","stop","background",kGreen);

    // ------------------
    // Regions
    // ------------------
   

    
    AddRegion("SR_2l","",&SR_2l);
    AddRegion("SR_2l_3j_incl","",&SR_2l_3j_incl);
    AddRegion("SR_2l_lowMass", "", &SR_2l_lowMass); 
    AddRegion("SR_2l_2j_lowMass", "", &SR_2l_2j_lowMass); 

    AddRegion("SR_2l_2j","",&SR_2l_2j);
    AddRegion("SR_2l_3j","",&SR_2l_3j);
    AddRegion("SR_2l_4j","",&SR_2l_4j);
    AddRegion("SR_2l_ISR","",&SR_2l_ISR);

   
    /*
    AddRegion("SR1l","SR1l",&SR1l);
    AddRegion("SR1l_2j","SR1l_2j",&SR1l_2j);
    AddRegion("SR1l_3j","SR1l_3j",&SR1l_3j);
    AddRegion("SR1l_4j","SR1l_4j",&SR1l_4j);
    AddRegion("CR1l","CR1l",&CR1l);
    AddRegion("CR2l","CR2l",&CR2l);

    AddRegion("SR1l_2lb","SR1l_2lb",&SR1l_2lb);
    AddRegion("CR1l_2lb","CR1l_2lb",&CR1l_2lb);
    AddRegion("SR1l_2j_2lb","SR1l_2j_2lb", &SR1l_2j_2lb);
    AddRegion("SR1l_3j_2lb","SR1l_3j_2lb", &SR1l_3j_2lb);
    AddRegion("SR1l_4j_2lb","SR1l_4j_2lb", &SR1l_4j_2lb);
    */


    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","", &lepChannel);

    //SetLumi(12900.);
    SetLumi(37800.);
    //SetLumi(42800.);
    //SetLumi(3900.);

    Create1DHistos();
    Add2DHisto("StopMass","NeutralinoMass");
    Add2DHisto("dphi_lep","deta_lep");

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
    
    leadjet_pt = 0;
    if(myEvent.ak4pfjets_pt.size()>0) leadjet_pt = myEvent.ak4pfjets_pt[0];
    

    if(SR_2l() && myEvent.ak4pfjets_pt.size()>=2){
    	TLorentzVector l1;
    	TLorentzVector l2;
    	TLorentzVector j1;
    	TLorentzVector j2;
  	l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
  	l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
	///cout<<myEvent.jet_mass[0]<<" "<<myEvent.jet_mass[1]<<endl;
	///*
	j1.SetPtEtaPhiM(myEvent.ak4pfjets_pt[0], myEvent.ak4pfjets_eta[0], myEvent.ak4pfjets_phi[0], myEvent.jet_mass[0]); 
	j2.SetPtEtaPhiM(myEvent.ak4pfjets_pt[1], myEvent.ak4pfjets_eta[1], myEvent.ak4pfjets_phi[1], myEvent.jet_mass[0]); 

        Meff = (l1+l2+j1+j2).M(); 
    	//*/
    }



    //--- Info about the type of dataset is transmitted to the WeightfFactory
    currentProcessType == "data" ? WFact.SetIsData(true): WFact.SetIsData(false);
    currentProcessType == "signal" ? WFact.SetIsFastSim(true): WFact.SetIsFastSim(false);
    //---------------------

    // --- Recompute variables event per event


    mStop = -999;
    mNeutralino = -99;
    if (currentProcessType == "signal"){
      mStop = myEvent.gen_stop_m.at(0);
      mNeutralino = myEvent.gen_neutralino_m.at(0); 
      //cout<<"Masses = "<<mStop<<" "<<mNeutralino<<endl;
    }

    TFile *file = NULL;
    float weightSignal = -13;
    if (currentProcessType == "signal") //  && (myEvent.gen_stop_m.at(0) == 500 || myEvent.gen_stop_m.at(0) == 425 || myEvent.gen_stop_m.at(0) == 900))
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
        
	Int_t binx = xaxis->FindBin(stopMass);
        Int_t biny = yaxis->FindBin(neutralinoMass);
        uint32_t totalNrOfEvents = h2->GetBinContent(binx, biny);
        myEvent.totalNumberOfInitialEvent = h2->GetEntries();
        weightSignal = sigCrossSection * GetLumi() * myEvent.mc_weight / totalNrOfEvents;
        //cout <<weightSignal <<  "signal CS " << sigCrossSection << " mc weight " << myEvent.mc_weight << " nr of events " << totalNrOfEvents << endl;
        
        if ( it != scanMap.end() )
        {
            currentProcessClass = it->second;
        }
        else
        {
        	 //   cout << "process class not found ;stop" << stopMass <<" ,neutralino: " <<  neutralinoMass << endl;
        	//return;
	}
     }

    //@MJ@ TODO do a method from this
    /*
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
        //currentProcessClass = "singleLepton";
        currentProcessClass = "test";
    }
    else
    {}
    */

    if( currentProcessClass == "T2tt"){
        currentProcessClass = "T2tt";
      if(mStop == 275 && mNeutralino == 100){
	currentProcessClass = "stop";
      }
      /*
      if(mStop == 175 && mNeutralino == 0)
        currentProcessClass = "175-0";
      if(mStop == 200 && mNeutralino == 25)
        currentProcessClass = "200-25";
      if(mStop == 225 && mNeutralino == 50)
        currentProcessClass = "225-50";
      if(mStop == 250 && mNeutralino == 75)
        currentProcessClass = "250-75";
      */
    }

    //cout<<"n-tops: "<<myEvent.genTops_pt.size()<<endl;
    double top_pt = 0;
    if(myEvent.genTops_pt.size()!=0){
    	//cout<<"I found tops: "<<myEvent.genTops_pt[0]<<endl;
	for(unsigned int i=0;i<myEvent.genTops_pt.size();i++){
		if(myEvent.genTops_pdgid[i] == 6)
			top_pt = myEvent.genTops_pt[i];
	}
    }

    recompute(useTriggerInfo, currentDataset);

    //cout<<GetLumi()<<endl;
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

    float weight     = weightLumi;
    if (currentProcessType == "signal") weight = weightSignal;
    if (currentProcessType == "data") weight = 1.0;
    else weight*=btagWeight*lepWeight;
    //else weight*=1;
    //if(top_pt!=0) cout<<"top pt weight = "<< WFact.TopPTWeightComputor(top_pt)<<endl;
    //if(top_pt!=0) weight*= WFact.TopPTWeightComputor(top_pt);
   



    deta_lep = fabs(myEvent.lep1_eta - myEvent.lep2_eta);
    dphi_deta = ComputeDphiDetaV(dphi_lep, deta_lep);

    int categ = 0;
    if(SR_2l_2j()) categ = 0;
    if(SR_2l_3j()) categ = 1;
    if(SR_2l_4j()) categ = 2;
    if(SR_2l_ISR()) categ = 3;

    category = categ;
	
    /*
    if(SR_2l_ISR()) categ = 3;
    else{
    	if(myEvent.pfmet<150) categ = 0;
	if(myEvent.pfmet>=150 && myEvent.pfmet<250) categ = 1;
	else categ = 2;
    }
    */


    super_dphi_lep = categ*TMath::Pi()+dphi_lep; 
    super_deta_lep = categ*TMath::Pi()+deta_lep;
    super_dphi_deta = categ*25+dphi_deta;


    /*
    float phi_1 = myEvent.lep1_phi;
    float phi_2 = myEvent.lep2_phi;
    //go from 0 to 2 Pi
    while( phi_1 >= 2*TMath::Pi()) phi_1-=TMath::Pi();
    while( phi_1 <= 0) phi_1+=TMath::Pi();
    while( phi_2 >= 2*TMath::Pi()) phi_2-=TMath::Pi();
    while( phi_2 <= 0) phi_2+=TMath::Pi();


    //dphi_lep = fabs(myEvent.lep2_phi - myEvent.lep1_phi);
    dphi_lep = fabs(phi_1 - phi_2);
    while( dphi_lep >= TMath::Pi()) dphi_lep-=TMath::Pi();
    while( dphi_lep <= 0) dphi_lep+=TMath::Pi();
   */
    //cout<<myEvent.lep1_phi<<" "<<myEvent.lep2_phi<<" "<<dphi_lep<<endl;

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

    vector<string> inputRegionTags = {"SR_2l_2j","SR_2l_3j","SR_2l_4j","SR_2l_ISR"};
    string channel = "lepChannel";
    TableDataMC table(this, inputRegionTags, channel);//, string options = "")
    table.Print("yield.tab", 4);
    cout<<"CR0b: SF: "<<endl;
    Figure fig = table.Get("CR0b", "data")/table.Get("CR0b", "totalSM");
    cout<<fig.Print()<<endl;
    cout<<"CR1b: SF: "<<endl;
    fig = table.Get("CR1b", "data")/table.Get("CR1b","totalSM");
    cout<<fig.Print()<<endl;
    cout<<"CR2b: SF: "<<endl;
    fig = table.Get("CR2b", "data")/table.Get("CR2b","totalSM");
    cout<<fig.Print()<<endl;


    // Schedule plots
    //

    SchedulePlots("1DSuperimposed");
    SchedulePlots("1DStack");
    SchedulePlots("1DDataMCComparison");
    SchedulePlots("2D");
    //SchedulePlots("1DFigureOfMerit","var=DeltaPhibb,cutType=keepHighValues");
    //SchedulePlots("1DFigureOfMerit","var=DeltaPhibb,cutType=keepLowValues,type=sOverSqrtB,minBackground=0");
    //SchedulePlots("2D");

    // Config plots
    SetGlobalStringOption("Plot", "infoTopRight", "CMS Preliminary");
    SetGlobalStringOption("Plot", "infoTopLeft",  "12.9 fb^{-1} (#sqrt{s} = 13 TeV)");

    SetGlobalBoolOption("Plot", "exportPdf", false);
    SetGlobalBoolOption("Plot", "exportEps", true);
    SetGlobalBoolOption("Plot", "exportPng", false);
    SetGlobalFloatOption("DataMCRatio","min",0);
    SetGlobalFloatOption("DataMCRatio","max",2);

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plotsTest/");

    //CombineCardMaker card(this, signalReg, "lepChannel", "throw", "StopMass", "NeutralinoMass");
    CombineCardMaker card(this, inputRegionTags, "lepChannel", "T2tt", "StopMass", "NeutralinoMass");
    card.Print("card_yield.tab", 4);
    card.ProduceCard("cards");
    CombineCardMaker card_Zbi(this, inputRegionTags, "lepChannel", "T2tt", "StopMass", "NeutralinoMass", true);
    card_Zbi.Print("card_zbi.tab", 4);
    
    // ######################
    //  Tables and other stuff
    // ######################

    cout << "end of processing" << endl;
 }
