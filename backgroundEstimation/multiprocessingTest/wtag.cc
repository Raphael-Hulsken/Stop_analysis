#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <cstring>
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
map< pair<uint32_t,uint32_t>, string > scanMap;
int lepJetOverlap;
float ak10pt[100]; //= new float();
float ak10prunedMass[100]; //= new float();
float ak10nSub[100]; //= new float();

// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/mjansova/store/tmp/2902/";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Oct21_postpatch8";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 8;


    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 300,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 800,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("genW4Ak10", "genW4Ak10",  "",1 , 0  , 2,  &(onTheFlyVariables.m_genW4Ak10), "noUnderflowInFirstBin");
    AddVariable("recW4Ak10", "recW4Ak10",  "",1  ,0 ,2,  &(onTheFlyVariables.m_recW4Ak10), "noUnderflowInFirstBin");
    AddVariable("recMatchedW4Ak10", "recMatchedW4Ak10",  "",1 ,  0 ,2 ,  &(onTheFlyVariables.m_recMatchedW4Ak10), "noUnderflowInFirstBin");
    AddVariable("recMatchedFakeW4Ak10", "recMatchedFakeW4Ak10",  "",1 , 0  ,2 ,  &(onTheFlyVariables.m_recMatchedFakeW4Ak10), "noUnderflowInFirstBin");
    AddVariable("lepJetOverlap", "lepJetOverlap",  "", 5,  0, 5,  &(lepJetOverlap), "noUnderflowInFirstBin");
    AddVariable("ak10pt", "ak10pt",  "", 50,  0, 700 ,  ak10pt , "noUnderflowInFirstBin");
    AddVariable("ak10prunedMass", "ak10prunedMass",  "", 50,  50, 130,  ak10prunedMass, "noUnderflowInFirstBin");
    AddVariable("ak10nSub", "ak10nSub",  "", 20,  0, 1 ,  ak10nSub , "noUnderflowInFirstBin");
// ...

    // background
    AddProcessClass("ttV", "tt+boson", "background", kMagenta);
    //	AddDataset("TTW_ln", "ttV", 0,0.70*0.32);
   // 	AddDataset("TTW_qq", "ttV", 0, 0.70*0.675);
    	AddDataset("TTZ_ll", "ttV", 0, 0.62*0.2);
   // 	AddDataset("TTZ_qq", "ttV", 0, 0.62*0.7);

    
  AddProcessClass("T2tt_600-950_1to450", "T2tt_600-950_1to450", "signal", kBlue);
    	AddDataset("600-950_1to450", "T2tt_600-950_1to450", 0, 0 );
 
    //Process classes for all scans in one sample
    //fad all scans and create map of names
    TFile *ftmp = NULL;
    TH2D *htmp = NULL;
    TString fNameTmp =  babyTuplePath+"600-950_1to450.root";
    ftmp = new TFile(fNameTmp);
    htmp = (TH2D*)ftmp->Get("hStopNeutralino")->Clone();
   
    
    for(uint32_t bx = 0; bx < htmp->GetNbinsX(); bx++)
    {
        for(uint32_t by = 0; by < htmp->GetNbinsY(); by++)
        {
            if(htmp->GetBinContent(bx+1,by+1))
            {
                std::cout << "bin Xedge: " << htmp->GetXaxis()->GetBinLowEdge(bx+1) << " bin Y edge " << htmp->GetYaxis()->GetBinLowEdge(by+1) << std::endl;
                pair<uint32_t, uint32_t> key = make_pair( htmp->GetXaxis()->GetBinLowEdge(bx+1), htmp->GetYaxis()->GetBinLowEdge(by+1));
                string stops = to_string(htmp->GetXaxis()->GetBinLowEdge(bx+1));
                string neutrs = to_string( htmp->GetYaxis()->GetBinLowEdge(by+1));
                scanMap[key] = stops+"_"+neutrs;
                AddProcessClass( stops+"_"+neutrs, stops+"_"+neutrs, "signal", kBlue);
            }
        }

    }

    delete htmp;
    delete ftmp;
    htmp =NULL;
    ftmp =NULL;

  

    AddProcessClass("SingleTop", "Single top", "background", kGreen);
    	AddDataset("ST_tW-atop", "SingleTop", 0, 35.6 );
    	AddDataset("ST_tW-top", "SingleTop", 0, 35.6 );
    	AddDataset("ST_t_1", "SingleTop", 0, 70.69 );
    	AddDataset("ST_t_2", "SingleTop", 0, 70.69 );
    	AddDataset("ST_s", "SingleTop", 0, 3.682 );

  
   AddProcessClass("TTJets", "TT+jets - 1l", "background", kRed);
        AddDataset("TTjets_M5_1", "TTJets", 0, 831.76);
        AddDataset("TTjets_M5_2", "TTJets", 0, 831.76);
     // AddProcessClass("TT_2l", "TT+jets - 2l", "background", kAzure+1);
     // AddProcessClass("TT_0l", "TT+jets - 0l", "background", kViolet+6);
 
   AddProcessClass("WJets", "W+jetsNew", "background", kBlack);
    	AddDataset("WJetsToLNu_HT-100To200_1", "WJets", 0, 1345*1.21); 
    	AddDataset("WJetsToLNu_HT-100To200_2", "WJets", 0, 1345*1.21); 
    	AddDataset("WJetsToLNu_HT-200To400", "WJets", 0, 359.7*1.21 ); 
    	AddDataset("WJetsToLNu_HT-400To600", "WJets", 0, 48.91*1.21 ); 
    	AddDataset("WJetsToLNu_HT-600To800", "WJets", 0, 12.1*1.21); 
    	AddDataset("WJetsToLNu_HT-800To1200", "WJets", 0, 5.50*1.21 ); 
    	AddDataset("WJetsToLNu_HT-1200To2500", "WJets", 0, 1.33*1.21); 
    	AddDataset("WJetsToLNu_HT-2500ToInf", "WJets", 0, 0.032*1.21); 

    AddRegion("defaultSearchSpec3", "defaultSearchSpec3", &defaultSearchSpec3);
    AddRegion("defaultSearchSpec4", "defaultSearchSpec4", &defaultSearchSpec4);
    AddRegion("defaultSearchSpec4250", "defaultSearchSpec4250", &defaultSearchSpec4250);
  
    AddRegion("DefaultBin1", "DefaultBin1", &DefaultBin1);
    AddRegion("DefaultBin2", "DefaultBin2", &DefaultBin2);
    AddRegion("DefaultBin3", "DefaultBin3", &DefaultBin3);
    AddRegion("DefaultBin4", "DefaultBin4", &DefaultBin4);
    AddRegion("DefaultBin5", "DefaultBin5", &DefaultBin5);
    AddRegion("DefaultBin6", "DefaultBin6", &DefaultBin6);
    AddRegion("DefaultBin7", "DefaultBin7", &DefaultBin7);

    AddRegion("Default2Bin3", "Default2Bin3", &Default2Bin3);
    AddRegion("Default2Bin4", "Default2Bin4", &Default2Bin4);
    AddRegion("Default2Bin5", "Default2Bin5", &Default2Bin5);

    AddRegion("Default3Bin3", "Default3Bin3", &Default3Bin3);
    AddRegion("Default3Bin4", "Default3Bin4", &Default3Bin4);
    AddRegion("Default3Bin5", "Default3Bin5", &Default3Bin5);
    
    //AddRegion("lowDmNoAk8JetsBin", "lowDmNoAk8JetsBin", &lowDmNoAk8JetsBin);
    //AddRegion("lowDmOnePlusAk8JetBin", "lowDmOnePlusAk8JetBin", &lowDmOnePlusAk8JetBin);
    
    AddRegion("NoAk8JetsBin", "NoAk8JetsBin", &NoAk8JetsBin);
    AddRegion("OnePlusAk8JetBin", "OnePlusAk8JetBin", &OnePlusAk8JetBin);

    AddRegion("threeJetsNoAk8JetsBin", "threeJetsNoAk8JetsBin", &threeJetsNoAk8JetsBin);
    AddRegion("threeJetsOnePlusAk8JetBin", "threeJetsOnePlusAk8JetBin", &threeJetsOnePlusAk8JetBin);

    //AddRegion("lowDmNoAk10JetsBin", "lowDmNoAk10JetsBin", &lowDmNoAk10JetsBin);
    //AddRegion("lowDmOnePlusAk10JetBin", "lowDmOnePlusAk10JetBin", &lowDmOnePlusAk10JetBin);
    
    AddRegion("NoAk10JetsBin", "NoAk10JetsBin", &NoAk10JetsBin);
    AddRegion("OnePlusAk10JetBin", "OnePlusAk10JetBin", &OnePlusAk10JetBin);

    AddRegion("threeJetsNoAk10JetsBin", "threeJetsNoAk10JetsBin", &threeJetsNoAk10JetsBin);
    AddRegion("threeJetsOnePlusAk10JetBin", "threeJetsOnePlusAk10JetBin", &threeJetsOnePlusAk10JetBin);
    
     // ...
    //AddChannel("muon", "#mu channel", &muonChannelSelector);
    //AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    // ...

    SetLumi(2260.);

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
    lepJetOverlap = 4;

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

    /*if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons == 2))
    {
        currentProcessClass = "TT_2l";
    }

    if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons == 0))
    {
        currentProcessClass = "TT_0l";
    }*/

    //cout << "after signal things" << endl;
    //at m_genW4Ak10;
    //  float m_recW4Ak10;
    //    float m_recMatchedW4Ak10;
    //      float m_recMatchedFakeW4Ak10;

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
    //
    if(myEvent.ak10pfjets_tau2.size() != myEvent.ak10pfjets_tau1.size())
        throw std::runtime_error("tau1 and tau2 vectors used to comute nsub have different leghts");
    vector<float> subjetiness;

    for(uint32_t v=0; v< myEvent.ak10pfjets_tau2.size(); v++ )
    {
        subjetiness.push_back(myEvent.ak10pfjets_tau2.at(v)/myEvent.ak10pfjets_tau1.at(v));
    }

    if(myEvent.ak10pfjets_pt.size() != 0)
    {
        memcpy ( ak10pt, &myEvent.ak10pfjets_pt[0], myEvent.ak10pfjets_pt.size()*sizeof(float) );
        memcpy ( ak10prunedMass, &myEvent.ak10pfjets_pruned_mass[0], myEvent.ak10pfjets_pruned_mass.size()*sizeof(float) );
        memcpy ( ak10nSub, &subjetiness[0], subjetiness.size()*sizeof(float) );
        SetSizeOfVarArray("ak10pt", myEvent.ak10pfjets_pt.size());
        SetSizeOfVarArray("ak10prunedMass", myEvent.ak10pfjets_pruned_mass.size());
        SetSizeOfVarArray("ak10nSub", subjetiness.size());
    }
    else
    {
        SetSizeOfVarArray("ak10pt", 0);
        SetSizeOfVarArray("ak10prunedMass", 0);
        SetSizeOfVarArray("ak10nSub", 0);
    }


     countEffAndFRAK10(24);
     lepJetOverlap = static_cast<int>(evaluateLeptonAndAk10Overlap());
    //nrOfW = nrOfWTaggs();
    

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


    //onTheFlyVariables.m_genW4Ak10 /= weight;
    //onTheFlyVariables.m_recW4Ak10 /= weight;
    //onTheFlyVariables.m_recMatchedW4Ak10 /= weight;
    //onTheFlyVariables.m_recMatchedFakeW4Ak10 /= weight;
    //cout << "weight " << weight << endl;
    // Fill this event in the histo collections

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
    SchedulePlots("1DStack");
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
    WritePlots("./plots5/"); //was 4

    // ######################
    //  Tables and other stuff
    // ######################


   vector<string> totYield = {"DefaultBin1", "DefaultBin2", "DefaultBin3", "DefaultBin4", "DefaultBin5", "DefaultBin6", "DefaultBin7", "NoAk8JetsBin", "OnePlusAk8JetBin", "threeJetsNoAk8JetsBin", "threeJetsOnePlusAk8JetBin", "NoAk10JetsBin", "OnePlusAk10JetBin", "threeJetsNoAk10JetsBin", "threeJetsOnePlusAk10JetBin" , "Default2Bin3", "Default2Bin4", "Default2Bin5", "Default3Bin3", "Default3Bin4", "Default3Bin5" };
   TableBackgroundSignal(this, totYield,"combinedChannel" ).Print("yield.tab", 4);
   TableBackgroundSignal(this, totYield,"combinedChannel" ).PrintLatex("yield.tex", 4);
  
   vector<string> regionsDef = {"DefaultBin1", "DefaultBin2", "DefaultBin3", "DefaultBin4", "DefaultBin5", "DefaultBin6", "DefaultBin7"};
   TableBackgroundSignal(this, regionsDef,"combinedChannel" ).Print("defReg.tab", 4);
   TableBackgroundSignal(this, regionsDef,"combinedChannel" ).PrintLatex("defReg.tex", 4);
   TableToBackgroundRatio(this, regionsDef,"combinedChannel" ).Print("defRegToB.tab", 4);
   TableToBackgroundRatio(this, regionsDef,"combinedChannel" ).PrintLatex("defRegBToB.tex", 4);
   TableZbi(this, regionsDef,"combinedChannel" ).Print("defRegZbi.tab", 4);
   TableZbi(this, regionsDef,"combinedChannel" ).PrintLatex("defRegZbi.tex", 4);


   vector<string> regionsDef2 = {"DefaultBin1", "DefaultBin2", "Default2Bin3", "Default2Bin4", "Default2Bin5","DefaultBin6", "DefaultBin7"};
   TableBackgroundSignal(this, regionsDef2,"combinedChannel" ).Print("defReg2.tab", 4);
   TableBackgroundSignal(this, regionsDef2,"combinedChannel" ).PrintLatex("defReg2.tex", 4);
   TableToBackgroundRatio(this, regionsDef2,"combinedChannel" ).Print("defRegToB2.tab", 4);
   TableToBackgroundRatio(this, regionsDef2,"combinedChannel" ).PrintLatex("defRegBToB2.tex", 4);
   TableZbi(this, regionsDef2,"combinedChannel" ).Print("defRegZbi2.tab", 4);
   TableZbi(this, regionsDef2,"combinedChannel" ).PrintLatex("defRegZbi2.tex", 4);
  

   vector<string> regionsDef3 = {"DefaultBin1", "DefaultBin2", "Default3Bin3", "Default3Bin4", "Default3Bin5", "DefaultBin6","DefaultBin7"};
   TableBackgroundSignal(this, regionsDef3,"combinedChannel" ).Print("defReg3.tab", 4);
   TableBackgroundSignal(this, regionsDef3,"combinedChannel" ).PrintLatex("defReg3.tex", 4);
   TableToBackgroundRatio(this, regionsDef3,"combinedChannel" ).Print("defRegToB3.tab", 4);
   TableToBackgroundRatio(this, regionsDef3,"combinedChannel" ).PrintLatex("defRegBToB3.tex", 4);
   TableZbi(this, regionsDef3,"combinedChannel" ).Print("defRegZbi3.tab", 4);
   TableZbi(this, regionsDef3,"combinedChannel" ).PrintLatex("defRegZbi3.tex", 4);
 
   vector<string> regions4Ak8 = {"DefaultBin3", "NoAk8JetsBin", "OnePlusAk8JetBin"};
   TableBackgroundSignal(this, regions4Ak8,"combinedChannel" ).Print("4ak8.tab", 4);
   TableBackgroundSignal(this, regions4Ak8,"combinedChannel" ).PrintLatex("4ak8.tex", 4);
   TableToBackgroundRatio(this, regions4Ak8,"combinedChannel" ).Print("4ak8ToB.tab", 4);
   TableToBackgroundRatio(this, regions4Ak8,"combinedChannel" ).PrintLatex("4ak8ToB.tex", 4);
   TableZbi(this, regions4Ak8,"combinedChannel" ).Print("4ak8Zbi.tab", 4);
   TableZbi(this, regions4Ak8,"combinedChannel" ).PrintLatex("4ak8Zbi.tex", 4);

   vector<string> regions3Ak8 = {"threeJetsNoAk8JetsBin", "threeJetsOnePlusAk8JetBin"};
   TableBackgroundSignal(this, regions3Ak8,"combinedChannel" ).Print("3ak8.tab", 4);
   TableBackgroundSignal(this, regions3Ak8,"combinedChannel" ).PrintLatex("3ak8.tex", 4);
   TableToBackgroundRatio(this, regions3Ak8,"combinedChannel" ).Print("3ak8ToB.tab", 4);
   TableToBackgroundRatio(this, regions3Ak8,"combinedChannel" ).PrintLatex("3ak8ToB.tex", 4);
   TableZbi(this, regions3Ak8,"combinedChannel" ).Print("3ak8Zbi.tab", 4);
   TableZbi(this, regions3Ak8,"combinedChannel" ).PrintLatex("3ak8Zbi.tex", 4);
   
   //vector<string> regionslowDmAk8 = {"lowDmNoAk8JetsBin", "lowDmOnePlusAk8JetBin"};
   //TableBackgroundSignal(this, regionslowDmAk8,"combinedChannel" ).Print("lowDmak8.tab", 4);
   //TableBackgroundSignal(this, regionslowDmAk8,"combinedChannel" ).PrintLatex("lowDmak8.tex", 4);
   //TableToBackgroundRatio(this, regionslowDmAk8,"combinedChannel" ).Print("lowDmak8ToB.tab", 4);
   //TableToBackgroundRatio(this, regionslowDmAk8,"combinedChannel" ).PrintLatex("lowDmak8ToB.tex", 4);
   //TableZbi(this, regionslowDmAk8,"combinedChannel" ).Print("lowDmak8Zbi.tab", 4);
   //TableZbi(this, regionslowDmAk8,"combinedChannel" ).PrintLatex("lowDmak8Zbi.tex", 4);

   vector<string> regions4Ak10 = {"DefaultBin3", "NoAk10JetsBin", "OnePlusAk10JetBin"};
   TableBackgroundSignal(this, regions4Ak10,"combinedChannel" ).Print("4ak10.tab", 4);
   TableBackgroundSignal(this, regions4Ak10,"combinedChannel" ).PrintLatex("4ak10.tex", 4);
   TableToBackgroundRatio(this, regions4Ak10,"combinedChannel" ).Print("4ak10ToB.tab", 4);
   TableToBackgroundRatio(this, regions4Ak10,"combinedChannel" ).PrintLatex("4ak10ToB.tex", 4);
   TableZbi(this, regions4Ak10,"combinedChannel" ).Print("4ak10Zbi.tab", 4);
   TableZbi(this, regions4Ak10,"combinedChannel" ).PrintLatex("4ak10Zbi.tex", 4);

   vector<string> regions3Ak10 = {"threeJetsNoAk10JetsBin", "threeJetsOnePlusAk10JetBin"};
   TableBackgroundSignal(this, regions3Ak10,"combinedChannel" ).Print("3ak10.tab", 4);
   TableBackgroundSignal(this, regions3Ak10,"combinedChannel" ).PrintLatex("3ak10.tex", 4);
   TableToBackgroundRatio(this, regions3Ak10,"combinedChannel" ).Print("3ak10ToB.tab", 4);
   TableToBackgroundRatio(this, regions3Ak10,"combinedChannel" ).Print("3ak10ToB.tab", 4);
   TableZbi(this, regions3Ak10,"combinedChannel" ).PrintLatex("3ak10Zbi.tex", 4);
   TableZbi(this, regions3Ak10,"combinedChannel" ).PrintLatex("3ak10Zbi.tex", 4);

   //vector<string> regionslowDmAk10 = {"lowDmNoAk10JetsBin", "lowDmOnePlusAk10JetBin"};
   //TableBackgroundSignal(this, regionslowDmAk10,"combinedChannel" ).Print("lowDmak10.tab", 4);
   //TableBackgroundSignal(this, regionslowDmAk10,"combinedChannel" ).PrintLatex("lowDmak10.tex", 4);
   //TableBackgroundSignal(this, regionslowDmAk10,"combinedChannel" ).Print("lowDmak10ToB.tab", 4);
   //TableBackgroundSignal(this, regionslowDmAk10,"combinedChannel" ).PrintLatex("lowDmak10ToB.tex", 4);
   //TableZbi(this, regionslowDmAk10,"combinedChannel" ).Print("lowDmak10Zbi.tab", 4);
   //TableZbi(this, regionslowDmAk10,"combinedChannel" ).PrintLatex("lowDmak10Zbi.tex", 4);
  
}

