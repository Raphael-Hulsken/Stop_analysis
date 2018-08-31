
#include <vector>
using namespace std;

#include "../common/common.h"
/*
#include "TTree.h"
#include "BabyTupleFormat.h"
babyEvent myEvent;
*/


#include "OnTheFlyVariables.h"

//#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"
//#include "interface/BabyScrewdriver.h"

bool muonChannelSelector() { return true; }
bool goesInAnyChannel()                             { return (goesInSingleLeptonChannel() || goesInDoubleLeptonChannel());                  }



/*
double getSumOfWeight(const string& path, const string& filename){

    // Open file and get mcweight
    TFile*  theFile = TFile::Open((path+filename).c_str(),"READ");
    TTree*  theTree = (TTree*) theFile->Get("babyTuple");
    double  sumOfWeight = 0;	       
    TH1D* h = theTree->Draw("Sum$(mc_weight)");
    sumOfWeight = h->Integral();			           
    theFile->Close();
    delete h;
    return sumOfWeight;
}
*/

int getNofEvents(const string& path, const string& filename){
    // Open file and get mcweight
    TFile*  theFile = TFile::Open((path+filename+".root").c_str(),"READ");
    TTree*  theTree = (TTree*) theFile->Get("babyTuple");
    return theTree->GetEntries();
}

// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop2015/BabyTuples/v1_30_06/";
    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/MantaRay-patch7-pfcand/";
    
    totalNumberOfWorkers = 7;

    AddVariable("MT", "whatever",  "", 10,   0, 100,  &(myEvent.MT));
    //AddVariable("x", "whatever2", "", 20, -40,  40,  &(onTheFlyVariables.someCoolVariable));
    // ...

    AddProcessClass("W+jets",   "W+jets",                          "background",kYellow);
        AddDataset("WJets-MC@NLO",    "W+jets", 100, 1);
	//double sumW = getSumOfWeight(babyTuplePath,"WJets-MC@NLO");
	//cout<<"sum weights ="<<sumW<<endl;
	cout<<getNofEvents(babyTuplePath,"WJets-MC@NLO")<<endl;

    //AddProcessClass("ttbar",   "ttbar",                          "background",kBlue);
     //   AddDataset("TTJets-FXFX",    "ttbar", 0, 0);

    /*
    AddProcessClass("singleTop",   "singleTop",                          "background",kGreen);
        AddDataset("ST-t",    "singleTop", 0, 0);
    
    AddProcessClass("WW",   "WW",                          "background",kOrange);
        AddDataset("WW",    "WW", 0, 0);

    AddProcessClass("TTW",   "TTW",                          "background",kOrange+1);
        AddDataset("TTW",    "TTW", 0, 0);
    
    AddProcessClass("TTZ",   "TTZ",                          "background",kOrange+1);
        AddDataset("TTZ-lep",    "TTZ", 0, 0);
        AddDataset("TTZ-had",    "TTZ", 0, 0);
   */
    // ...

    AddRegion("preselection", "Preselection", { Cut("MT", '>', 50) });
    AddRegion("presel_MTpeak",          "Preselection (MT peak)",      &goesInMTpeak);
    AddRegion("presel_MTtail",          "Preselection (MT peak)",      &goesInMTtail);
    // ...

    //AddChannel("muon", "#mu channel", &muonChannelSelector);
    AddChannel("singleElec",   "e-channel",       &goesInSingleElecChannel  );
    AddChannel("singleMuon",   "#mu-channel",     &goesInSingleMuonChannel  );
    AddChannel("allChannels",  "",                &goesInAnyChannel         );
    // ...

    SetLumi(3.14);

    Create1DHistos();
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    // Compute on the fly variables if needed

    ComputeOnTheFlyVariables();

    cout<<"print"<<endl;
    cout<<"MT: "<<myEvent.mt_met_lep<<endl;
    cout<<"goesInAnyChannel: "<<goesInAnyChannel()<<endl;
    cout<<"goesInMTpeak: "<<goesInMTpeak()<<endl;
    cout<<"goesInMTtail: "<<goesInMTtail()<<endl;

    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    // Compute weight for current event

    float weightLumi = myEvent.crossSection * GetLumi() / myEvent.totalNumberOfInitialEvent;

    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;

    // Fill this event in the histo collections

    //AutoFillProcessClass(currentProcessClass, weight);
    AutoFillProcessClass(currentProcessClass, 0.001);
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

    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 8 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", true);
    SetGlobalBoolOption("Plot", "exportEps", false);
    SetGlobalBoolOption("Plot", "exportPng", false);

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plots/");

    // ######################
    //  Tables and other stuff
    // ######################

    TableDataMC(this,{"preselection","presel_MTpeak","presel_MTtail"},"muon", "").Print("table.tab",4);

}



