//
// ########################################################################
// #  Thanks to Wouter Verkerke, awesome stat teacher, for original code  #
// #                                                                      #
// # 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #501      #
// # Using simultaneous p.d.f.s to describe simultaneous fits to multiple #
// # datasets                                                             #
// #                                                                      #
// # 07/2008                                                              #
// ########################################################################

// ROOT headers

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TRandom.h"
#include "TMath.h"

using namespace RooFit;

// STL headers

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// Sonic screwdriver headers

#include "interface/Table.h"
#include "interface/Figure.h"

using namespace theDoctor;

// ####################################
// # MT tail correction configuration #
// ####################################

#define INPUT_FOLDER  "./plots/"
#define OUTPUT_FOLDER "./results/"

#define OBSERVABLE_FOR_FIT "Mlb_leadb_bin2"
//#define OBSERVABLE_FOR_FIT "DeltaRlj"


#define PROCESS_NAME_TT_1L "TT_1l"
#define PROCESS_NAME_TT_2L "TT_2l"
#define PROCESS_NAME_WJETS "WJets"
#define PROCESS_NAME_RARE "OTHER"

//#define CHANNEL_NAME singleLepton
#define CHANNEL_NAME "muon"


TRandom* randomnessGenerator;

// Uncertainty of the template fit method itself,
// coming from test with MCstat, JES, algorithms, ...

#define TEMPLATE_FIT_METHOD_UNCERTAINTY 0.2



// ########################################
// # Regions for cut-and-count correction #
// ########################################

vector<string> listIndividualCuts =
{
    "no MT cut",
    "CR_MET250_lowMT",
    "CR_MET250_lowMT_3j",
    "CR0b_MET250_lowMT",
    "CR0b_MET250_lowMT_3j",
    "50$\\lt$MET$\\lt$100",
    "100$\\lt$MET$\\lt$150",
    "150$\\lt$MET$\\lt$200",
    "200$\\lt$MET$\\lt$250",
    "250$\\lt$MET$\\lt$300",
    "$MET$\\gt$300",
    "MT2W>200",
    "200<MT2W>250",
    "MT2W>250"
    /*
    "$M_{T}$$\\ge$80",
    "$M_{T}$$\\ge$90",
    "$M_{T}$$\\ge$100",
    "$M_{T}$$\\ge$110",
    "$M_{T}$$\\ge$120",
    "$M_{T}$$\\ge$130",
    "$M_{T}$$\\ge$140",
    "80$\\le$$M_{T}$$\\le$90",
    "90$\\le$$M_{T}$$\\le$100",
    "100$\\le$$M_{T}$$\\le$110",
    "110$\\le$$M_{T}$$\\le$120",
    "120$\\le$$M_{T}$$\\le$130",
    "130$\\le$$M_{T}$$\\le$140",
    */
    //"MT_150"
    //"2j",
    //"3j"
    //"4j"
};

vector<string> listIndividualCuts_MTtail =
{
    "CR0b_presel",
    "CR_MET250_lowMT",
    "CR_MET250_lowMT_3j",
    "CR0b_MET250_lowMT",
    "CR0b_MET250_lowMT_3j",
    "CR0b_presel_MET50",
    "CR0b_presel_MET50",
    "CR0b_presel_MET100",
    "CR0b_presel_MET150",
    "CR0b_presel_MET200",
    "CR0b_presel_MET250",
    "CR0b_presel_MET300",
    "CR0b_presel_MT2Wtail",
    "CR0b_presel_MT2W200",
    "CR0b_presel_MT2W250"
    /*
    "CR0b_presel_MTtail_80",
    "CR0b_presel_MTtail_90",
    "CR0b_presel_MTtail_100",
    "CR0b_presel_MTtail_110",
    "CR0b_presel_MTtail_120",
    "CR0b_presel_MTtail_130",
    "CR0b_presel_MTtail_140",
    "CR0b_presel_MTtail_80_ex",
    "CR0b_presel_MTtail_90_ex",
    "CR0b_presel_MTtail_100_ex",
    "CR0b_presel_MTtail_110_ex",
    "CR0b_presel_MTtail_120_ex",
    "CR0b_presel_MTtail_130_ex",
    */
    //"CR0b_presel_MTtail_150"
    //"CR0b_presel_2j_MTtail",
    //"CR0b_presel_3j_MTtail"
    //"CR0b_presel_4j_MTtail"
};

vector<string>  listIndividualCuts_MTpeak =
{
    "CR0b_presel",
    "CR_MET250_lowMT",
    "CR_MET250_lowMT_3j",
    "CR0b_MET250_lowMT",
    "CR0b_MET250_lowMT_3j",
    "CR0b_presel_MET50",
    "CR0b_presel_MET100",
    "CR0b_presel_MET150",
    "CR0b_presel_MET200",
    "CR0b_presel_MET250",
    "CR0b_presel_MET300",
    "CR0b_presel",
    "CR0b_presel",
    "CR0b_presel"
    /*
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    */
    //"CR0b_presel_MTpeak"
    //"CR0b_presel_2j_MTpeak",
    //"CR0b_presel_3j_MTpeak"
    //"CR0b_presel_4j_MTpeak"
};

vector<string> listCutAndCounts =
{
    //"preselection"
};

map<string,vector<string> > listCutAndCounts_cuts;

void initCutAndCountCuts()
{
    //listCutAndCounts_cuts["preselection"]                     = { "MT_100"                     };
}
