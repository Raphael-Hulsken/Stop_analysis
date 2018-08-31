#include <vector>
#include <iostream>
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
//#define USE_GEN_INFO
//#define USE_GEN_INFO_EXT
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



bool muonChannelSelector() { return true; }
bool electronChannelSelector() { return true; }
uint64_t counting;
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 1;



    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 300,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 300,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets), "noUnderflowInFirstBin");
    AddVariable("Nbjets", "N_{b-jets}",  "", 6,   0, 5,  &(myEvent.ngoodbtags), "noUnderflowInFirstBin");
    AddVariable("NSubJets", "N_{sub-jets}",  "", 6,   0, 5,  &(onTheFlyVariables.m_ak8WSubJets), "noUnderflowInFirstBin");
    AddVariable("NFakeSubJets", "N_{fakesub-jets}",  "", 6,   0, 5,  &(onTheFlyVariables.m_ak8WFakeSubJets), "noUnderflowInFirstBin");
    AddVariable("SRBins", "Signal Region box",  "", 11,   -5, 5,  &(onTheFlyVariables.SRbins), "noUnderflowInFirstBin");
    AddVariable("genWPt", "Pt of gen W",  "", 50,   0, 1000,  &(onTheFlyVariables.m_genWPt), "noUnderflowInFirstBin");
    AddVariable("genTopPt", "Pt of gen top quark",  "", 50,   0, 1000,  &(onTheFlyVariables.m_genTopPt), "noUnderflowInFirstBin");
    AddVariable("genBarTopPt", "Pt of gen bar top quark",  "", 50,   0, 1000,  &(onTheFlyVariables.m_genBarTopPt), "noUnderflowInFirstBin");
    AddVariable("genWdR", "#delta R related to gen W",  "", 50,   0, 4,  &(onTheFlyVariables.m_genWdR), "noUnderflowInFirstBin");
    AddVariable("genTopdR", "#deltar R related to gen top quark",  "", 50,   0, 4,  &(onTheFlyVariables.m_genTopdR), "noUnderflowInFirstBin");
    AddVariable("genBarTopdR", "#deltar R related to gen bar top quark",  "", 50,   0, 4,  &(onTheFlyVariables.m_genBarTopdR), "noUnderflowInFirstBin");
    AddVariable("recWPt", "Pt of rec W",  "", 50,   0, 1000,  &(onTheFlyVariables.m_recWPt), "noUnderflowInFirstBin");
    AddVariable("ak8WMass", "mass of W ak8 jet",  "",50,   0, 400,  &(onTheFlyVariables.m_ak8WMass), "noUnderflowInFirstBin");
    AddVariable("ak8WPrunedMass", "pruned mass of W ak8 jet",  "", 50,   0, 400,  &(onTheFlyVariables.m_ak8WPrunedMass), "noUnderflowInFirstBin");
    AddVariable("ak8WTrimmedMass", "trimmed mass of W ak8 jet",  "", 50,   0, 400,  &(onTheFlyVariables.m_ak8WTrimmedMass), "noUnderflowInFirstBin");
    AddVariable("ak8WSoftDropMass", "softdrop mass of W ak8 jet",  "",50,   0, 400,  &(onTheFlyVariables.m_ak8WSoftDropMass), "noUnderflowInFirstBin");
    AddVariable("ak8recoWPt", "pt of W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8recoWPt), "noUnderflowInFirstBin");
    AddVariable("ak8MET", "MET of W ak8 jet",  "", 50,   0, 300,  &(onTheFlyVariables.m_ak8MET), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeMass", "mass of fake W ak8 jet",  "",50,   0, 400,  &(onTheFlyVariables.m_ak8WFakeMass), "noUnderflowInFirstBin");
    AddVariable("ak8WFakePrunedMass", "pruned mass of fake W ak8 jet",  "",50,   0, 400,  &(onTheFlyVariables.m_ak8WFakePrunedMass), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeTrimmedMass", "trimmed mass of fake W ak8 jet",  "", 50,   0, 400,  &(onTheFlyVariables.m_ak8WFakeTrimmedMass), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeSoftDropMass", "softdrop mass of fake W ak8 jet",  "", 50,   0, 400,  &(onTheFlyVariables.m_ak8WFakeSoftDropMass), "noUnderflowInFirstBin");
    AddVariable("ak8recoWFakePt", "pt of fake W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8recoWFakePt), "noUnderflowInFirstBin");
    AddVariable("ak8FakeMET", "MET of fake W ak8 jet",  "", 50,   0, 300,  &(onTheFlyVariables.m_ak8FakeMET), "noUnderflowInFirstBin");
    // ...

    AddProcessClass("T2tt_425_325", "T2tt_425_325", "background", kPink);
        AddDataset("T2tt_425_325", "T2tt_425_325", 0, 0);
    
   //AddProcessClass("T2tt_500_325", "T2tt_500_325", "background", kMagenta+3);
   //     AddDataset("T2tt_500_325", "T2tt_500_325", 0, 0);
   
    //AddProcessClass("T2tt_650_325", "T2tt_650_325", "background", kViolet+6);
    //    AddDataset("T2tt_650_325", "T2tt_650_325", 0, 0);

    //AddProcessClass("T2tt_850_325", "T2tt_850_325", "background", kOrange+7);
    //    AddDataset("T2tt_850_325", "T2tt_850_325", 0, 0);

    // ...

    AddRegion("preselectionTT", "PreselectionTT", &goesInPreselectionNoVetoNoMetCut);
    //AddRegion("BaselineSelection", "BaselineSelection", &goesInBaselineSearchSR);
    // ...

    AddChannel("muon", "#mu channel", &muonChannelSelector);
    AddChannel("electron", "e channel", &electronChannelSelector);
    // ...

    SetLumi(2169005.);  //@MJ@ TODO correct lumi

    Create1DHistos();
    Add2DHisto("genWPt","genWdR");
    Add2DHisto("genTopPt","genTopdR");
    Add2DHisto("genBarTopPt","genBarTopdR");
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
     
    //cout<<myEvent.mt_met_lep<<endl;
    //cout << " sel: "<< goesInPreselection() << endl;
    myEvent.crossSection  = CrossSection(currentDataset);
    cout<<"xsextion : "<<myEvent.crossSection<<endl;
    cout<<"tot weight: "<<myEvent.totalNumberOfInitialEvent<<endl;
    //cout<<"weight: "<<myEvent.mc_weight<<endl;
    // Compute on the fly variables if needed

   // myEvent.totalNumberOfInitialEvent = 50000;
     
    ComputeOnTheFlyVariables();
    
    //find generated particles and their properties
    //Initial clean up
    onTheFlyVariables.m_genWPt = -1;
    onTheFlyVariables.m_genWdR = -1;
    onTheFlyVariables.m_genTopPt = -1;
    onTheFlyVariables.m_genTopdR = -1;
    onTheFlyVariables.m_genBarTopPt = -1;
    onTheFlyVariables.m_genBarTopdR = -1;
    
    //find W decaying to two quarks
    findGenParticleProps(24, &(onTheFlyVariables.m_genWPt), &(onTheFlyVariables.m_genWdR)); //W
    //cout << "1) back in main; W: gen pt is: " << onTheFlyVariables.m_genWPt << " gen dR is: " <<  onTheFlyVariables.m_genWdR << " end of one event" << endl;
    
    //find t decaying to W+ and quark
    findGenParticleProps(6, &(onTheFlyVariables.m_genTopPt), &(onTheFlyVariables.m_genTopdR)); //top
    //cout << "2) back in main; top:  gen pt is: " << onTheFlyVariables.m_genTopPt << " gen dR is: " <<  onTheFlyVariables.m_genTopdR << " end of one event" << endl;

    //find bar t decaying to W- and quark
    findGenParticleProps(-6, &(onTheFlyVariables.m_genBarTopPt), &(onTheFlyVariables.m_genBarTopdR)); //top
    //cout << "3) back in main; bar top:  gen pt is: " << onTheFlyVariables.m_genBarTopPt << " gen dR is: " <<  onTheFlyVariables.m_genBarTopdR << " end of one event" << endl;

    //count efficiency of W tagging
    countEfficiency(24);
    //cout << "W ak8 pt" << onTheFlyVariables.m_ak8recoWPt << "ake ak8 W pt" << onTheFlyVariables.m_ak8recoWFakePt << endl;
    
    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    // Compute weight for current event

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;

    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;

    // Fill this event in the histo collections

    AutoFillProcessClass(currentProcessClass, weight);
    //AutoFillProcessClass(currentProcessClass, 0.001);

    counting++;
    if(counting % 1000 == 0)
    {
        cout << counting << endl;
    }
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
    SchedulePlots("2D");

    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", true);
    SetGlobalBoolOption("Plot", "exportEps", false);
    SetGlobalBoolOption("Plot", "exportPng", false);

    //SetGlobalBoolOption("1DSuperimposed", "includeSignal", true);
    //SetGlobalStringOption("1DStack", "includeSignal", "stack");

	    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plotsSig/");

    // ######################
    //  Tables and other stuff
    // ######################

    TableDataMC(this,{"preselection"},"muon", "").Print("table.tab",4);

}



