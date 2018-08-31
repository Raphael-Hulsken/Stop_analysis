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
bool combinedChannelSelector() { return true; }
uint64_t counting;
int nrOfW;
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_all_3j_lowpt_dphi_v2";
    //babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Oct21_postpatch8";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 1;


    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 300,  &(myEvent.mt_met_lep), "noUnderflowInFirstBin");
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 800,  &(myEvent.pfmet), "noUnderflowInFirstBin");
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W), "noUnderflowInFirstBin");
    /*AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets), "noUnderflowInFirstBin");
    AddVariable("Nbjets", "N_{b-jets}",  "", 5,   0, 5,  &(myEvent.ngoodbtags), "noUnderflowInFirstBin");
    AddVariable("NSubJets", "N_{sub-jets}",  "", 5,   0, 5,  &(onTheFlyVariables.m_ak8WSubJets), "noUnderflowInFirstBin");
    AddVariable("NFakeSubJets", "N_{fakesub-jets}",  "", 5,   0, 5,  &(onTheFlyVariables.m_ak8WFakeSubJets), "noUnderflowInFirstBin");
    AddVariable("SRBins", "Signal Region box",  "", 11,   -5, 5,  &(onTheFlyVariables.SRbins), "noUnderflowInFirstBin");
    AddVariable("effN1", "effN1",  "", 1, 0, 2,  &(onTheFlyVariables.m_recEffDef1), "noUnderflowInFirstBin");
    AddVariable("effD1", "effD1",  "", 1, 0, 2,  &(onTheFlyVariables.m_genEffDef1), "noUnderflowInFirstBin");
    AddVariable("effN2", "effN2",  "", 1, 0, 2,  &(onTheFlyVariables.m_recEffDef2), "noUnderflowInFirstBin");
    AddVariable("effD2", "effD2",  "", 1, 0, 2,  &(onTheFlyVariables.m_genEffDef2), "noUnderflowInFirstBin");
    AddVariable("effFakeN", "effFakeN",  "", 1, 0, 2,  &(onTheFlyVariables.m_matchedAndTaggedFakes), "noUnderflowInFirstBin");
    AddVariable("effFakeD", "efffakeD",  "", 1, 0, 2,  &(onTheFlyVariables.m_matchedFakes), "noUnderflowInFirstBin");
    AddVariable("NSJreals", "NSJreals",  "", 30, 0, 1.5,  &(onTheFlyVariables.m_subjetinessReals), "noUnderflowInFirstBin");
    AddVariable("NSJfakes", "NSJfakes",  "", 30, 0, 1.5,  &(onTheFlyVariables.m_subjetinessFakes), "noUnderflowInFirstBin");
    */
    AddVariable("nrOfW", "nrOfW",  "", 10, 0, 10,  &(nrOfW), "noUnderflowInFirstBin");
    /*AddVariable("genWPt", "Pt of gen W",  "", 50,   0, 1000,  &(onTheFlyVariables.m_genWPt), "noUnderflowInFirstBin");
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
    AddVariable("ak8WGenPt", "gen pt of W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8WGenPt), "noUnderflowInFirstBin");
    AddVariable("ak8MET", "MET of W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8MET), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeMass", "mass of fake W ak8 jet",  "",50,   0, 400,  &(onTheFlyVariables.m_ak8WFakeMass), "noUnderflowInFirstBin");
    AddVariable("ak8WFakePrunedMass", "pruned mass of fake W ak8 jet",  "",50,   0, 400,  &(onTheFlyVariables.m_ak8WFakePrunedMass), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeTrimmedMass", "trimmed mass of fake W ak8 jet",  "", 50,   0, 400,  &(onTheFlyVariables.m_ak8WFakeTrimmedMass), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeSoftDropMass", "softdrop mass of fake W ak8 jet",  "", 50,   0, 400,  &(onTheFlyVariables.m_ak8WFakeSoftDropMass), "noUnderflowInFirstBin");
    AddVariable("ak8recoWFakePt", "pt of fake W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8recoWFakePt), "noUnderflowInFirstBin");
    AddVariable("ak8WFakeGenPt", "gen pt of fake W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8WFakeGenPt), "noUnderflowInFirstBin");
    AddVariable("ak8FakeMET", "MET of fake W ak8 jet",  "", 50,   0, 1000,  &(onTheFlyVariables.m_ak8FakeMET), "noUnderflowInFirstBin");
    AddVariable("dphimet", "dphi met",  "", 50,   0, 7,  &(myEvent.dphi_ak4pfjets_met), "noUnderflowInFirstBin");
    AddVariable("mindR", "minimal dR between W and ak8 jet",  "", 50,   0, 2,  &(onTheFlyVariables.m_mindR), "noUnderflowInFirstBin");
 */
// ...

    // background
    
    // Signal(s)
  
    AddProcessClass("SingleTop", "Single top", "background", kGreen);
    	AddDataset("ST_s", "SingleTop", 0, 3.38 );
    	AddDataset("ST_t-atop", "SingleTop", 0, 26.49 );
    	AddDataset("ST_t-top", "SingleTop", 0, 44.51 );
    	AddDataset("ST_tW-atop", "SingleTop", 0, 35.85 );
    	AddDataset("ST_tW-top", "SingleTop", 0 , 35.85 );

    AddProcessClass("ttV", "tt+boson", "background", kMagenta);
    	AddDataset("TTW_ln", "ttV", 0,0.70*0.32);
    	AddDataset("TTW_qq", "ttV", 0, 0.70*0.675);
    	AddDataset("TTZ_ll", "ttV", 0, 0.62*0.2);
    	AddDataset("TTZ_qq", "ttV", 0, 0.62*0.7);

    AddProcessClass("VV", "di-bosons", "background", kMagenta+1);
    	AddDataset("WW_aMC", "VV", 0 , 48.4);
    	AddDataset("ZZ", "VV", 0, 16.523);
    	AddDataset("WZ", "VV", 0 , 47.13);
    
    AddProcessClass("T2tt_850_325", "T2tt_850_325", "background",  kMagenta+3);
        AddDataset("T2tt_850_325","T2tt_850_325", 0, 0.01896);
    
    AddProcessClass("T2tt_650_325", "T2tt_650_325", "signal",  kViolet+6);
        AddDataset("T2tt_650_325","T2tt_650_325", 0, 0.107045);
    
    AddProcessClass("T2tt_500_325", "T2tt_500_325", "signal",  kOrange+7);
        AddDataset("T2tt_500_325","T2tt_500_325", 0, 0.51848);
  
   AddProcessClass("TTJets", "TT+jets", "background", kRed);
        AddDataset("TTjets_M5", "TTJets", 0, 831.76);
    //AddProcessClass("TT_2l", "TT+jets - 2l", "background", kCyan);

 
//    AddProcessClass("WJets14", "W+jets14", "background", kBlue);
//    	AddDataset("WJetsToLNu_HT-100to200_P14", "WJets14", 0, 1345*1.21 ); 
//    	AddDataset("WJetsToLNu_HT-200to400_14", "WJets14", 0, 359.7*1.21 ); 
//    	AddDataset("WJetsToLNu_HT-400to600_14", "WJets14", 0, 48.91*1.21 ); 
//    	AddDataset("WJetsToLNu_HT-600toInf_P14", "WJets14", 0, 18.77*1.21 ); 
 
   AddProcessClass("WJetsNew", "W+jetsNew", "background", kBlack);
    	AddDataset("WJetsToLNu_HT-100To200", "WJetsNew", 0, 1345*1.21); 
    	AddDataset("WJetsToLNu_HT-200To400", "WJetsNew", 0, 359.7*1.21 ); 
    	AddDataset("WJetsToLNu_HT-400To600", "WJetsNew", 0, 48.91*1.21 ); 
    	AddDataset("WJetsToLNu_HT-600To800", "WJetsNew", 0, 12.05*1.21); 
    	AddDataset("WJetsToLNu_HT-600ToInf", "WJetsNew", 0, 18.77*1.21 ); 
    	AddDataset("WJetsToLNu_HT-800To1200", "WJetsNew", 0, 5.501*1.21 ); 
    	AddDataset("WJetsToLNu_HT-1200To2500", "WJetsNew", 0, 1.329*1.21); 
    	AddDataset("WJetsToLNu_HT-2500ToInf", "WJetsNew", 0, 0.03216*1.21); 

    AddRegion("defaultSearch", "defaultSearch", &defaultSearch);
    AddRegion("boostedSearch", "boostedSearch", &boostedSearch);
/*   
    AddRegion("DefaultBin1", "DefaultBin1", &DefaultBin1);
    AddRegion("DefaultBin2", "DefaultBin2", &DefaultBin2);
    AddRegion("DefaultBin3", "DefaultBin3", &DefaultBin3);
    AddRegion("DefaultBin4", "DefaultBin4", &DefaultBin4);
    AddRegion("DefaultBin5", "DefaultBin5", &DefaultBin5);
    
    AddRegion("NoAk8JetsBin1", "NoAk8JetsBin1", &NoAk8JetsBin1);
    AddRegion("NoAk8JetsBin2", "NoAk8JetsBin2", &NoAk8JetsBin2);
    AddRegion("NoAk8JetsBin3", "NoAk8JetsBin3", &NoAk8JetsBin3);
    AddRegion("NoAk8JetsBin4", "NoAk8JetsBin4", &NoAk8JetsBin4);
    AddRegion("NoAk8JetsBin5", "NoAk8JetsBin5", &NoAk8JetsBin5);

    AddRegion("AtLeastOneAk8JetBin1", "AtLeastOneAk8JetBin1", &AtLeastOneAk8JetBin1);
    AddRegion("AtLeastOneAk8JetBin2", "AtLeastOneAk8JetBin2", &AtLeastOneAk8JetBin2);
    AddRegion("AtLeastOneAk8JetBin3", "AtLeastOneAk8JetBin3", &AtLeastOneAk8JetBin3);
    AddRegion("AtLeastOneAk8JetBin4", "AtLeastOneAk8JetBin4", &AtLeastOneAk8JetBin4);
    AddRegion("AtLeastOneAk8JetBin5", "AtLeastOneAk8JetBin5", &AtLeastOneAk8JetBin5);
    
    AddRegion("AtLeastOneAk8JetNoCutBin1", "AtLeastOneAk8JetNoCutBin1", &AtLeastOneAk8JetNoCutBin1);
    AddRegion("AtLeastOneAk8JetNoCutBin2", "AtLeastOneAk8JetNoCutBin2", &AtLeastOneAk8JetNoCutBin2);
    AddRegion("AtLeastOneAk8JetNoCutBin3", "AtLeastOneAk8JetNoCutBin3", &AtLeastOneAk8JetNoCutBin3);
    AddRegion("AtLeastOneAk8JetNoCutBin4", "AtLeastOneAk8JetNoCutBin4", &AtLeastOneAk8JetNoCutBin4);
    AddRegion("AtLeastOneAk8JetNoCutBin5", "AtLeastOneAk8JetNoCutBin5", &AtLeastOneAk8JetNoCutBin5);

    AddRegion("ThreeJetsDefaultBin1", "ThreeJetsDefaultBin1", &ThreeJetsDefaultBin1);
    AddRegion("ThreeJetsDefaultBin2", "ThreeJetsDefaultBin2", &ThreeJetsDefaultBin2);
    AddRegion("ThreeJetsDefaultBin3", "ThreeJetsDefaultBin3", &ThreeJetsDefaultBin3);
    AddRegion("ThreeJetsDefaultBin4", "ThreeJetsDefaultBin4", &ThreeJetsDefaultBin4);
    AddRegion("ThreeJetsDefaultBin5", "ThreeJetsDefaultBin5", &ThreeJetsDefaultBin5);
    
    AddRegion("ThreeJetsDefaultMiasBin", "ThreeJetsDefaultMiasBin", &ThreeJetsDefaultMiasBin);
    AddRegion("Default2", "Default2", &Default2);
   
    AddRegion("ThreeJetsInEventNoAk8Bin1", "ThreeJetsInEventNoAk8Bin1", &ThreeJetsInEventNoAk8Bin1);
    AddRegion("ThreeJetsInEventNoAk8Bin2", "ThreeJetsInEventNoAk8Bin2", &ThreeJetsInEventNoAk8Bin2);
    AddRegion("ThreeJetsInEventNoAk8Bin3", "ThreeJetsInEventNoAk8Bin3", &ThreeJetsInEventNoAk8Bin3);
    AddRegion("ThreeJetsInEventNoAk8Bin4", "ThreeJetsInEventNoAk8Bin4", &ThreeJetsInEventNoAk8Bin4);
    AddRegion("ThreeJetsInEventNoAk8Bin5", "ThreeJetsInEventNoAk8Bin5", &ThreeJetsInEventNoAk8Bin5);

    AddRegion("ThreeJetsInEventNoAk8MiasBin", "ThreeJetsInEventNoAk8MiasBin", &ThreeJetsInEventNoAk8MiasBin);
    
    AddRegion("ThreeJetsInEventOneOfThemAk8Bin1", "ThreeJetsInEventOneOfThemAk8Bin1", &ThreeJetsInEventOneOfThemAk8Bin1);
    AddRegion("ThreeJetsInEventOneOfThemAk8Bin2", "ThreeJetsInEventOneOfThemAk8Bin2", &ThreeJetsInEventOneOfThemAk8Bin2);
    AddRegion("ThreeJetsInEventOneOfThemAk8Bin3", "ThreeJetsInEventOneOfThemAk8Bin3", &ThreeJetsInEventOneOfThemAk8Bin3);
    AddRegion("ThreeJetsInEventOneOfThemAk8Bin4", "ThreeJetsInEventOneOfThemAk8Bin4", &ThreeJetsInEventOneOfThemAk8Bin4);
    AddRegion("ThreeJetsInEventOneOfThemAk8Bin5", "ThreeJetsInEventOneOfThemAk8Bin5", &ThreeJetsInEventOneOfThemAk8Bin5);
    
    AddRegion("ThreeJetsInEventOneOfThemAk8MiasBin", "ThreeJetsInEventOneOfThemAk8MiasBin", &ThreeJetsInEventOneOfThemAk8MiasBin);

    AddRegion("lowMT2WDefaultBin1", "lowMT2WDefaultBin1", &lowMT2WDefaultBin1);
    AddRegion("lowMT2WDefaultBin2", "lowMT2WDefaultBin2", &lowMT2WDefaultBin2);
    AddRegion("lowMT2WDefaultBin3", "lowMT2WDefaultBin3", &lowMT2WDefaultBin3);
    AddRegion("lowMT2WDefaultBin4", "lowMT2WDefaultBin4", &lowMT2WDefaultBin4);
    AddRegion("lowMT2WDefaultBin5", "lowMT2WDefaultBin5", &lowMT2WDefaultBin5);
    
    AddRegion("lowMT2WThreeJetsDefaultMiasBin", "lowMT2WThreeJetsDefaultMiasBin", &lowMT2WThreeJetsDefaultMiasBin);


    AddRegion("NewSRDefaultForLowMETBin1", "NewSRDefaultForLowMETBin1", &NewSRDefaultForLowMETBin1);
    AddRegion("NewSR0ak8Bin2", "NewSR0ak8Bin2", &NewSR0ak8Bin2);
    AddRegion("NewSR1andMoreak8Bin3", "NewSR1andMoreak8Bin3", &NewSR1andMoreak8Bin3);

    AddRegion("NewSR0ak8Bin2WOCuts", "NewSR0ak8Bin2WOCuts", &NewSR0ak8Bin2WOCuts);
    AddRegion("NewSR1andMoreak8Bin3WOCuts", "NewSR1andMoreak8Bin3WOCuts", &NewSR1andMoreak8Bin3WOCuts);
    
    AddRegion("Default2Bin1", "Default2Bin1", &Default2Bin1);
    AddRegion("Default2Bin2", "Default2Bin2", &Default2Bin2);
    AddRegion("Default2Bin3", "Default2Bin3", &Default2Bin3);
    
    AddRegion("Default2OneAk8Bin1", "Default2OneAk8Bin1", &Default2OneAk8Bin1);
    AddRegion("Default2OneAk8Bin2", "Default2OneAk8Bin2", &Default2OneAk8Bin2);
    AddRegion("Default2OneAk8Bin3", "Default2OneAk8Bin3", &Default2OneAk8Bin3);
    
    AddRegion("Default2NoAk8Bin1", "Default2NoAk8Bin1", &Default2NoAk8Bin1);
    AddRegion("Default2NoAk8Bin2", "Default2NoAk8Bin2", &Default2NoAk8Bin2);
    AddRegion("Default2NoAk8Bin3", "Default2NoAk8Bin3", &Default2NoAk8Bin3);
    */
    //AddRegion("NewThreeJetsInEventOneOfThemAk8MiasBin", "NewThreeJetsInEventOneOfThemAk8MiasBin", &NewThreeJetsInEventOneOfThemAk8MiasBin);
    //AddRegion("NewThreeJetsInEventNoAk8MiasBin", "NewThreeJetsInEventNoAk8MiasBin", &NewThreeJetsInEventNoAk8MiasBin);
    
     //AddRegion("preselectionTT", "preselectionTT", &goesInPreselectionNoVetoNoMetCut);
    // ...

    //AddChannel("muon", "#mu channel", &muonChannelSelector);
    //AddChannel("electron", "e channel", &electronChannelSelector);
    AddChannel("combinedChannel","e/#mu-channel",&combinedChannelSelector);
    // ...

    SetLumi(2000.);  //@MJ@ TODO correct lumi

    Create1DHistos();
    Add2DHisto("MET","nrOfW");
    /*Add2DHisto("genWPt","genWdR");
    Add2DHisto("genTopPt","genTopdR");
    Add2DHisto("genBarTopPt","genBarTopdR");
*/
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{

    counting++;
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

    nrOfW = nrOfWTaggs();

    
    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    //if (currentProcessClass == "TT_1l" && (myEvent.numberOfGeneratedLeptons == 2))
    //               currentProcessClass = "TT_2l";
   
     // Compute weight for current event

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;

    float weight     = weightLumi;
    //float weight     = 1;
    //float weight     = 1;
    if (currentProcessType == "data") weight = 1.0;

    //cout << "weight " << weight << endl;
    // Fill this event in the histo collections

    AutoFillProcessClass(currentProcessClass, weight);
    //cout << "process " << currentProcessClass << " weight " << weight << endl;
    //AutoFillProcessClass(currentProcessClass, 0.001);

    // @MJ@ TODO push back sth
    //myEvent.ak8pfjets_trimmed_mass.clear();
    //myEvent.ak8pfjets_pt.clear();
    
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

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plots3/");

    // ######################
    //  Tables and other stuff
    // ######################
    /*
   vector<string> regions = {"boostedSearch","DefaultBin1", "DefaultBin2", "DefaultBin3", "DefaultBin4", "DefaultBin5", "NoAk8JetsBin1", "NoAk8JetsBin2", "NoAk8JetsBin3", "NoAk8JetsBin4", "NoAk8JetsBin5", "AtLeastOneAk8JetBin1", "AtLeastOneAk8JetBin2", "AtLeastOneAk8JetBin3", "AtLeastOneAk8JetBin4", "AtLeastOneAk8JetBin5", "AtLeastOneAk8JetNoCutBin1", "AtLeastOneAk8JetNoCutBin2", "AtLeastOneAk8JetNoCutBin3", "AtLeastOneAk8JetNoCutBin4", "AtLeastOneAk8JetNoCutBin5", "ThreeJetsDefaultBin1", "ThreeJetsDefaultBin2", "ThreeJetsDefaultBin3", "ThreeJetsDefaultBin4", "ThreeJetsDefaultBin5", "ThreeJetsDefaultMiasBin",  "ThreeJetsInEventNoAk8Bin1", "ThreeJetsInEventNoAk8Bin2", "ThreeJetsInEventNoAk8Bin3", "ThreeJetsInEventNoAk8Bin4", "ThreeJetsInEventNoAk8Bin5", "ThreeJetsInEventNoAk8MiasBin", "ThreeJetsInEventOneOfThemAk8Bin1", "ThreeJetsInEventOneOfThemAk8Bin2", "ThreeJetsInEventOneOfThemAk8Bin3", "ThreeJetsInEventOneOfThemAk8Bin4", "ThreeJetsInEventOneOfThemAk8Bin5", "ThreeJetsInEventOneOfThemAk8MiasBin", "lowMT2WDefaultBin1",  "lowMT2WDefaultBin2", "lowMT2WDefaultBin3", "lowMT2WDefaultBin4", "lowMT2WDefaultBin5", "lowMT2WThreeJetsDefaultMiasBin", "NewSRDefaultForLowMETBin1", "NewSR0ak8Bin2", "NewSR1andMoreak8Bin3", "Default2Bin1", "Default2Bin2", "Default2Bin3", "NewThreeJetsInEventOneOfThemAk8MiasBin", "NewThreeJetsInEventNoAk8MiasBin","NewSR0ak8Bin2WOCuts", "NewSR1andMoreak8Bin3WOCuts"};
    
   TableBackgroundSignal(this, regions,"combinedChannel" ).Print("yieldTable3fb.tab", 4);
   
   vector<string> regionsDef = {"DefaultBin1", "DefaultBin2", "DefaultBin3", "DefaultBin4", "DefaultBin5"};
   TableBackgroundSignal(this, regionsDef,"combinedChannel" ).Print("def3fb.tab", 4);
   TableBackgroundSignal(this, regionsDef,"combinedChannel" ).PrintLatex("def3fb.tex", 4);
   
   vector<string> regionsNoAk8 = {"NoAk8JetsBin1", "NoAk8JetsBin2", "NoAk8JetsBin3", "NoAk8JetsBin4", "NoAk8JetsBin5"};
   TableBackgroundSignal(this, regionsNoAk8,"combinedChannel" ).Print("noak83fb.tab", 4);
   TableBackgroundSignal(this, regionsNoAk8,"combinedChannel" ).PrintLatex("noak83fb.tex", 4);

   vector<string> regionsAk8 = {"AtLeastOneAk8JetBin1", "AtLeastOneAk8JetBin2", "AtLeastOneAk8JetBin3", "AtLeastOneAk8JetBin4", "AtLeastOneAk8JetBin5"};
   TableBackgroundSignal(this, regionsAk8,"combinedChannel" ).Print("ak83fb.tab", 4);
   TableBackgroundSignal(this, regionsAk8,"combinedChannel" ).PrintLatex("ak83fb.tex", 4);

   vector<string> regions3Jets = {"ThreeJetsDefaultMiasBin", "ThreeJetsInEventNoAk8MiasBin", "ThreeJetsInEventOneOfThemAk8MiasBin", "NewThreeJetsInEventOneOfThemAk8MiasBin", "NewThreeJetsInEventNoAk8MiasBin"};
   TableBackgroundSignal(this, regions3Jets,"combinedChannel" ).Print("3jets3fb.tab", 4);
   TableBackgroundSignal(this, regions3Jets,"combinedChannel" ).PrintLatex("3jets3fb.tex", 4);

   vector<string> regionsLowMT2W = {"lowMT2WDefaultBin1",  "lowMT2WDefaultBin2", "lowMT2WDefaultBin3", "lowMT2WDefaultBin4", "lowMT2WDefaultBin5", "lowMT2WThreeJetsDefaultMiasBin"};
   TableBackgroundSignal(this, regionsLowMT2W,"combinedChannel" ).Print("lowMT2W3fb.tab", 4);
   TableBackgroundSignal(this, regionsLowMT2W,"combinedChannel" ).PrintLatex("lowMT2W23fb.tex", 4);
   
   vector<string> regions4JetsNew = {"NewSRDefaultForLowMETBin1", "NewSR0ak8Bin2", "NewSR1andMoreak8Bin3", "NewSR0ak8Bin2WOCuts", "NewSR1andMoreak8Bin3WOCuts"};
   TableBackgroundSignal(this, regions4JetsNew,"combinedChannel" ).Print("newSR4Jets3fb.tab", 4);
   TableBackgroundSignal(this, regions4JetsNew,"combinedChannel" ).PrintLatex("newSR4Jets3fb.tex", 4);
   
   vector<string> default2 = {"Default2Bin1", "Default2Bin2", "Default2Bin3"};
   TableBackgroundSignal(this, default2, "combinedChannel" ).Print("default23fb.tab", 4);
   TableBackgroundSignal(this, default2, "combinedChannel" ).PrintLatex("default23fb.tex", 4);


   vector<string> newBLregionsNoAk8 = {"Default2NoAk8Bin1", "Default2NoAk8Bin2", "Default2NoAk8Bin3"};
   TableBackgroundSignal(this, newBLregionsNoAk8,"combinedChannel" ).Print("newBLnoak82fb.tab", 4);
   TableBackgroundSignal(this, newBLregionsNoAk8,"combinedChannel" ).PrintLatex("newBLnoak82fb.tex", 4);

   vector<string> newBLregionsAk8 = {"Default2OneAk8Bin1", "Default2OneAk8Bin2", "Default2OneAk8Bin3"};
   TableBackgroundSignal(this, newBLregionsAk8,"combinedChannel" ).Print("newBLak82fb.tab", 4);
   TableBackgroundSignal(this, newBLregionsAk8,"combinedChannel" ).PrintLatex("newBLak82fb.tex", 4);
*/
}

