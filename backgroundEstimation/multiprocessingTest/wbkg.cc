#include <vector>
#include <iostream>
#include <map>
#include <string>
using namespace std;

#include "TTree.h"
//babyEvent myEvent;




//-----------------------------------
//-- Load only necessary branches
//-----------------------------------

#define USE_GEN_INFO
#define USE_GEN_INFO_EXT
//#define USE_SKIMMING_VAR
#define USE_JETS
#define USE_JETS_EXT
#define USE_GEN_INFO
#define USE_GEN_INFO_EXT
#define USE_LEP1
//#define USE_LEP1_EXT
#define USE_LEP2
//#define USE_LEP2_EXT
#define USE_JETS
//#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_VAR_BASELINE
#define USE_GLOBAL_VAR
//#define USE_OLD_VAR
#define USE_AK8_JETS
#define USE_AK10_JETS


//-----------------------------------
// -- load Babytuple Format --
// -- Need first to activate the branches
//-----------------------------------
//#include "../common/BabyTupleFormat.h"
#include "../../common/Reader_CommonFormat.h"

#include "../../common/selectionDefinitions13TeV.h"


//---------------------------------------
// -- Add variables computed on the flyt
// -- Need first to load the format
//--------------------------------------
#include "OnTheFlyVariables.h"
#include "chi2/chi2.h"

// Should be called after selectionDefinition 'cause it need myEvent varible ...
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"



//bool muonChannelSelector() { return true; }
//bool electronChannelSelector() { return true; }
bool combinedChannelSelector() { return true; }
uint64_t counting;
int nrOfW;
TString storedDataset = "";
TH2D *h2 = NULL;
float nSjak8 = -13;
float nSjak10 = -13;
float number = 7;
float number2[2];
double chi2;
double chi2_v0;
double chi2_v1;
double chi2_v2;
double chi2_v3;
double chi2_v4;
double chi2_v5;
map< pair<uint32_t,uint32_t>, string > scanMap;
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/2902/";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Oct21_postpatch8";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 10;


    //AddVariable("bla", "bla",  "", 10,   0, 10,  &(number), "noUnderflowInFirstBin");
    //AddVariable("bla2", "bla2",  "", 10,   0, 10,  number2, "noUnderflowInFirstBin");
    //AddVariable("chi2", "chi2",  "", 20,   0, 20,  &(chi2), "noUnderflowInFirstBin");

    AddVariable("chi2_v0", "chi2",  "", 20,   0, 100,  &(chi2_v0), "noUnderflowInFirstBin");
    AddVariable("chi2_v1", "chi2",  "", 20,   0, 100,  &(chi2_v1), "noUnderflowInFirstBin");
    AddVariable("chi2_v2", "chi2",  "", 20,   0, 100,  &(chi2_v2), "noUnderflowInFirstBin");
    AddVariable("chi2_v3", "chi2",  "", 20,   0, 100,  &(chi2_v3), "noUnderflowInFirstBin");
    AddVariable("chi2_v4", "chi2",  "", 20,   0, 100,  &(chi2_v4), "noUnderflowInFirstBin");
    AddVariable("chi2_v5", "chi2",  "", 20,   0, 100,  &(chi2_v5), "noUnderflowInFirstBin");
    /*
    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 300,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 800,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    */
// ...

    // background
    //AddProcessClass("SingleTop", "Single top", "background", kGreen);
    //	AddDataset("ST_tW-atop", "SingleTop", 0, 35.85 );
//    	AddDataset("ST_tW-top", "SingleTop", 0, 35.85 );
  //  	AddDataset("ST_s", "SingleTop", 0, 3.38 );
    
  /*
  AddProcessClass("T2tt_600-950_1to450", "T2tt_600-950_1to450", "signal", kBlue);
    	AddDataset("600-950_1to450", "T2tt_600-950_1to450", 0, 0 );
 */
    //Process classes for all scans in one sample
    //fead all scans and create map of names
    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    TString fNameTmp =  babyTuplePath+"600-950_1to450.root";
    ftmp = new TFile(fNameTmp);
    htmp = (TH2D*)ftmp->Get("hStopNeutralino")->Clone();
   
    AddProcessClass("T2tt_600-950_1to450", "T2tt_600-950_1to450", "signal", kBlue);
    	AddDataset("600-950_1to450", "T2tt_600-950_1to450", 0, 0 );
    
    for(uint32_t bx = 0; bx < htmp->GetNbinsX(); bx++)
    {
        for(uint32_t by = 0; by < htmp->GetNbinsY(); by++)
        {
            if(htmp->GetBinContent(bx+1,by+1))
            {
                //std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinLowEdge(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinLowEdge(by+1) << std::endl;
                pair<uint32_t, uint32_t> key = make_pair( htmp->GetXaxis()->GetBinLowEdge(bx+1), htmp->GetYaxis()->GetBinLowEdge(by+1));
                string stops = to_string(htmp->GetXaxis()->GetBinLowEdge(bx+1));
                string neutrs = to_string( htmp->GetYaxis()->GetBinLowEdge(by+1));
                scanMap[key] = stops+"_"+neutrs;
                cout<<stops<<" "<<neutrs<<endl;
		if(int(htmp->GetXaxis()->GetBinLowEdge(bx+1)) == 600 && int(htmp->GetYaxis()->GetBinLowEdge(by+1)) == 100) 
		AddProcessClass( stops+"_"+neutrs, stops+"_"+neutrs, "signal", kBlue);
            }
        }

    }

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

 
    /*
    // Signal(s)
  
  //  AddProcessClass("SingleTop", "Single top", "background", kGreen);
  //  	AddDataset("ST_s", "SingleTop", 0, 3.38 );
  //  	AddDataset("ST_t-atop", "SingleTop", 0, 26.49 );
  //  	AddDataset("ST_t-top", "SingleTop", 0, 44.51 );
  //  	AddDataset("ST_tW-atop", "SingleTop", 0, 35.85 );
   // 	AddDataset("ST_tW-top", "SingleTop", 0 , 35.85 );

   // AddProcessClass("ttV", "tt+boson", "background", kMagenta);
    //	AddDataset("TTW_ln", "ttV", 0,0.70*0.32);
   // 	AddDataset("TTW_qq", "ttV", 0, 0.70*0.675);
   // 	AddDataset("TTZ_ll", "ttV", 0, 0.62*0.2);
   // 	AddDataset("TTZ_qq", "ttV", 0, 0.62*0.7);

  //  AddProcessClass("VV", "di-bosons", "background", kMagenta+1);
  //  	AddDataset("WW_aMC", "VV", 0 , 48.4);
  //  	AddDataset("ZZ", "VV", 0, 16.523);
   // 	AddDataset("WZ", "VV", 0 , 47.13);
 */ 
   AddProcessClass("TTJets", "TT+jets - 1l", "background", kRed);
        AddDataset("TTjets_M5_1", "TTJets", 0, 831.76);
   AddDataset("TTjets_M5_2", "TTJets", 0, 831.76);
      AddProcessClass("TT_2l", "TT+jets - 2l", "background", kAzure+1);
     // AddProcessClass("TT_0l", "TT+jets - 0l", "background", kViolet+6);
//*/ 
   ///*
   AddProcessClass("WJets", "W+jetsNew", "background", kBlack);
    	AddDataset("WJetsToLNu_HT-100To200_1", "WJets", 0, 1345*1.21); 
    	AddDataset("WJetsToLNu_HT-100To200_2", "WJets", 0, 1345*1.21); 
    	AddDataset("WJetsToLNu_HT-200To400", "WJets", 0, 359.7*1.21 ); 
    	AddDataset("WJetsToLNu_HT-400To600", "WJets", 0, 48.91*1.21 ); 
    	AddDataset("WJetsToLNu_HT-600To800", "WJets", 0, 12.05*1.21); 
    	AddDataset("WJetsToLNu_HT-800To1200", "WJets", 0, 5.501*1.21 ); 
    	AddDataset("WJetsToLNu_HT-1200To2500", "WJets", 0, 1.329*1.21); 
    	AddDataset("WJetsToLNu_HT-2500ToInf", "WJets", 0, 0.03216*1.21); 
//*/
  
    AddRegion("defaultSearchAllCuts", "defaultSearchAllCuts", &defaultSearchAllCuts);
    ///*
    AddRegion("DefaultBin1", "DefaultBin1", &DefaultBin1);
    AddRegion("DefaultBin2", "DefaultBin2", &DefaultBin2);
    AddRegion("DefaultBin3", "DefaultBin3", &DefaultBin3);
    AddRegion("DefaultBin4", "DefaultBin4", &DefaultBin4);
    AddRegion("DefaultBin5", "DefaultBin5", &DefaultBin5);
    //AddRegion("DefaultBin6", "DefaultBin6", &DefaultBin6); this one is incorrect
    AddRegion("DefaultBin7", "DefaultBin7", &DefaultBin7);
    //*/

    
     // ...
    //AddChannel("muon", "#mu channel", &muonChannelSelector);
    //AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    // ...

    SetLumi(2440.);

    Create1DHistos();
    //Add2DHisto("MET","nrOfW");
    /*Add2DHisto("genWPt","genWdR");
    Add2DHisto("genTopPt","genTopdR");
    Add2DHisto("genBarTopPt","genBarTopdR");
*/
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    //store previous dataset and if it does not equal to other clone histogram from file

    counting++;
    //number = 2;
    //number2[0] = 3;
    //number2[1] = 5;
    
    //SetSizeOfVarArray("bla2", sizeof(number2)/sizeof(float));
    //cout << "processing event" << endl;
    
    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);
    TFile *fle = NULL;
    float weightSignal = -13;
    if (currentProcessType == "signal")
    {
        if(currentDataset != storedDataset)
        {
            storedDataset = currentDataset;
            TString fName =  babyTuplePath+currentDataset+".root";
            fle = new TFile(fName);
            h2 = (TH2D*)fle->Get("hStopNeutralino")->Clone();
        }
        //std::cout << "reading stop mass from file: "<<babyTuplePath+currentDataset <<".root"<< std::endl;
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
        TAxis *xaxis = h2->GetXaxis();
        TAxis *yaxis = h2->GetYaxis();
        Int_t binx = xaxis->FindBin(stopMass);
        Int_t biny = yaxis->FindBin(neutralinoMass);
        uint32_t totalNrOfEvents = h2->GetBinContent(binx, biny);
        myEvent.totalNumberOfInitialEvent = h2->GetEntries();
        weightSignal = sigCrossSection * GetLumi() * myEvent.mc_weight / totalNrOfEvents;
        //cout << "signal CS " << sigCrossSection << " mc weight " << myEvent.mc_weight << " nr of events " << totalNrOfEvents << endl;
    }
    
    bool new_formula = true; 
    int njets_max = 6;
    //cout<<"@@@ chi2_v0"<<endl;
    chi2_v0 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass, 1, njets_max, false);
    //cout<<"@@@ chi2_v1"<<endl;
    chi2_v1 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass, 1, njets_max, new_formula);
    ///*
    chi2_v2 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass, 2, njets_max, new_formula);
    chi2_v3 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass, 3, njets_max, new_formula);
    chi2_v4 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass, 4, njets_max, new_formula);
    //*/
    ////cout<<"@@@ chi2_v5"<<endl;
    chi2_v5 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass, 5, njets_max, new_formula);

    //chi2 = Chi2(myEvent.jet_pt, myEvent.jet_eta, myEvent.jet_phi, myEvent.jet_mass);
    //cout << "chi2: " << chi2 << endl;

    /*if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons == 2))
    {
        currentProcessClass = "TT_2l";
    }

    if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons == 0))
    {
        currentProcessClass = "TT_0l";
    }*/

    //cout << "after signal things" << endl;

    //@MJ@ TODO
    //some fake info due to missing info in babytuple
   //myEvent.ngoodbtags = 1;
    
/*
    cout << "size " << myEvent.ak8pfjets_trimmed_mass.size() << endl;
    cout << "size pt " << myEvent.ak8pfjets_pt.size() << endl;

i
    for(uint32_t v = 0; v <myEvent.ak8pfjets_trimmed_mass.size(); v++)
    {    
    if(counting % 2 == 0)
    {
        myEvent.ak8pfjets_trimmed_mass.push_back(80);
        myEvent.ak8pfjets_pt.push_back(350);
    }
    } 
 */    
   // cout << "size 2 " << myEvent.ak8pfjets_trimmed_mass.size() << endl;
    //cout << "size pt 2 " << myEvent.ak8pfjets_pt.size() << endl;
    //cout<<myEvent.mt_met_lep<<endl;
    //cout << " sel: "<< goesInPreselection() << endl;
    myEvent.crossSection  = GetDatasetCrossSection(currentDataset);
    //cout<<"xsextion : "<<myEvent.crossSection<<endl;
    //cout<<"tot weight: "<<myEvent.totalNumberOfInitialEvent<<endl;
    //cout<<"weight: "<<myEvent.mc_weight<<endl;
    // Compute on the fly variables if needed
/*
     
    ComputeOnTheFlyVariables();
    
    //find generated particles and their properties
    //Initial clean up
    onTheFlyVariables.m_genWPt = -1;
    onTheFlyVariables.m_genWdR = -1;
    onTheFlyVariables.m_genTopPt = -1;
    onTheFlyVariables.m_genTopdR = -1;
    onTheFlyVariables.m_genBarTopPt = -1;
    onTheFlyVariables.m_genBarTopdR = -1;
   
    myEvent.dphi_ak4pfjets_met = calculateDPhi();
 
    //find W decaying to two quarks
    findGenParticleProps(24, &(onTheFlyVariables.m_genWPt), &(onTheFlyVariables.m_genWdR)); //W
    //cout << "1) back in main; W: gen pt is: " << onTheFlyVariables.m_genWPt << " gen dR is: " <<  onTheFlyVariables.m_genWdR << " end of one event" << endl;
    
    //find t decaying to W+ and quark
    findGenParticleProps(6, &(onTheFlyVariables.m_genTopPt), &(onTheFlyVariables.m_genTopdR)); //top
    //cout << "2) back in main; top:  gen pt is: " << onTheFlyVariables.m_genTopPt << " gen dR is: " <<  onTheFlyVariables.m_genTopdR << " end of one event" << endl;

    //find bar t decaying to W- and quark
    findGenParticleProps(-6, &(onTheFlyVariables.m_genBarTopPt), &(onTheFlyVariables.m_genBarTopdR)); //top
    //cout << "3) back in main; bar top:  gen pt is: " << onTheFlyVariables.m_genBarTopPt << " gen dR is: " <<  onTheFlyVariables.m_genBarTopdR << " end of one event" << endl;

    minJetdRToGenParticle(24, &(onTheFlyVariables.m_mindR));
    
    //count efficiency of W tagging
    countEfficiency(24);
    //cout << "W ak8 pt" << onTheFlyVariables.m_ak8recoWPt << "ake ak8 W pt" << onTheFlyVariables.m_ak8recoWFakePt << endl;
*/

    //if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons == 2))
    //               currentProcessClass = "TT_2l";
   
     // Compute weight for current event

    //@MJ@ TODO if signal do differently
    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;

    //cout << "bkg  CS " << myEvent.crossSection << " mc weight " << myEvent.mc_weight << " nr of events " << myEvent.totalNumberOfInitialEvent << endl;
    float weight     = weightLumi;
    //float weight     = 1;
    //float weight     = 1;
    if (currentProcessType == "data") weight = 1.0;
    if (currentProcessType == "signal") weight = weightSignal;


    //cout << "weight " << weight << endl;
    // Fill this event in the histo collections

    //if(chi2 < 2)
        AutoFillProcessClass(currentProcessClass, weight);
    //cout << "process " << currentProcessClass << " weight " << weight << endl;
    //AutoFillProcessClass(currentProcessClass, 0.001);

    // @MJ@ TODO push back sth
    //myEvent.ak8pfjets_trimmed_mass.clear();
    //myEvent.ak8pfjets_pt.clear();
    
    if(counting % 5000 == 0)
    {
        cout << counting << endl;

    }

    //cout << "up to here" << std::endl;
   
    //delete fle; //@MJ@ TODO be aware of deleting in here, may not be good with null
    
}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    // ######################
    //  Plot configuration and production
    // ######################

    // Schedule plots

    SchedulePlots("1DSuperimposed");
    //SchedulePlots("1DStack");
    //SchedulePlots("2D");

    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", true);
    SetGlobalBoolOption("Plot", "exportEps", false);
    SetGlobalBoolOption("Plot", "exportPng", false);

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plots2/");

    // ######################
    //  Tables and other stuff
    // ######################

   vector<string> regionsDef = {"DefaultBin1", "DefaultBin2", "DefaultBin3", "DefaultBin4", "DefaultBin5", "DefaultBin7"};
   TableBackgroundSignal(this, regionsDef,"combinedChannel" ).Print("wbkg.tab", 4);
   TableBackgroundSignal(this, regionsDef,"combinedChannel" ).PrintLatex("wbkg.tex", 4);
   TableToBackgroundRatio(this, regionsDef,"combinedChannel" ).Print("wbkgToB.tab", 4);
   TableToBackgroundRatio(this, regionsDef,"combinedChannel" ).PrintLatex("wbkgToB.tex", 4);
   TableZbi(this, regionsDef,"combinedChannel" ).Print("wbkgZbi.tab", 4);
   TableZbi(this, regionsDef,"combinedChannel" ).PrintLatex("wbkgZbi.tex", 4);

   
}

