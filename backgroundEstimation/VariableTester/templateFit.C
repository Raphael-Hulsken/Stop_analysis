
// ########################################################################
// # Credits goes to Wouter Verkerke for original code                    #
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

#define INPUT_FOLDER  "./plotsTest/"
#define OUTPUT_FOLDER "./results/"

//#define OBSERVABLE_FOR_FIT "dphi_lep"
//#define OBSERVABLE_FOR_FIT "deta_lep"
//#define OBSERVABLE_FOR_FIT "dphi_deta"
//#define OBSERVABLE_FOR_FIT "super_dphi_lep"
#define OBSERVABLE_FOR_FIT "super_deta_lep"
//#define OBSERVABLE_FOR_FIT "super_dphi_deta"



//#define OBSERVABLE_FOR_FIT "lep_pt"

#define REGION_FOR_FIT "lepChannel"


TRandom* randomnessGenerator;

// Uncertainty of the template fit method itself,
// coming from test with MCstat, JES, algorithms, ...

#define TEMPLATE_FIT_METHOD_UNCERTAINTY 0.2


// ##############################
// # Regions for BDT correction #
// ##############################
/*
vector<string> listBDTSignalRegions =
{
    "BDT_T2tt_1",
    "BDT_T2tt_2",
    "BDT_T2tt_5",
    "BDT_T2bw075_1",
    "BDT_T2bw075_2",
    "BDT_T2bw075_3",
    "BDT_T2bw075_5",
    "BDT_T2bw050_1",
    "BDT_T2bw050_3",
    "BDT_T2bw050_4",
    "BDT_T2bw050_5",
    "BDT_T2bw050_6",
    "BDT_T2bw025_1",
    "BDT_T2bw025_3",
    "BDT_T2bw025_4",
    "BDT_T2bw025_6"
};

vector<string> listBDTSignalRegions_MTtail =
{
    "BDT_MTtail_T2tt_1",
    "BDT_MTtail_T2tt_2",
    "BDT_MTtail_T2tt_5",
    "BDT_MTtail_T2bw075_1",
    "BDT_MTtail_T2bw075_2",
    "BDT_MTtail_T2bw075_3",
    "BDT_MTtail_T2bw075_5",
    "BDT_MTtail_T2bw050_1",
    "BDT_MTtail_T2bw050_3",
    "BDT_MTtail_T2bw050_4",
    "BDT_MTtail_T2bw050_5",
    "BDT_MTtail_T2bw050_6",
    "BDT_MTtail_T2bw025_1",
    "BDT_MTtail_T2bw025_3",
    "BDT_MTtail_T2bw025_4",
    "BDT_MTtail_T2bw025_6"
};

vector<string> listBDTSignalRegions_MTpeak =
{
    "BDT_MTpeak_T2tt_1",
    "BDT_MTpeak_T2tt_2",
    "BDT_MTpeak_T2tt_5",
    "BDT_MTpeak_T2bw075_1",
    "BDT_MTpeak_T2bw075_2",
    "BDT_MTpeak_T2bw075_3",
    "BDT_MTpeak_T2bw075_5",
    "BDT_MTpeak_T2bw050_1",
    "BDT_MTpeak_T2bw050_3",
    "BDT_MTpeak_T2bw050_4",
    "BDT_MTpeak_T2bw050_5",
    "BDT_MTpeak_T2bw050_6",
    "BDT_MTpeak_T2bw025_1",
    "BDT_MTpeak_T2bw025_3",
    "BDT_MTpeak_T2bw025_4",
    "BDT_MTpeak_T2bw025_6"
};

vector<string> listBDTSignalRegions_MTpeak_NoBtag =
{
    "BDT_MTPeakNoBtag_T2tt_1",
    "BDT_MTPeakNoBtag_T2tt_2",
    "BDT_MTPeakNoBtag_T2tt_5",
    "BDT_MTPeakNoBtag_T2bw075_1",
    "BDT_MTPeakNoBtag_T2bw075_2",
    "BDT_MTPeakNoBtag_T2bw075_3",
    "BDT_MTPeakNoBtag_T2bw075_5",
    "BDT_MTPeakNoBtag_T2bw050_1",
    "BDT_MTPeakNoBtag_T2bw050_3",
    "BDT_MTPeakNoBtag_T2bw050_4",
    "BDT_MTPeakNoBtag_T2bw050_5",
    "BDT_MTPeakNoBtag_T2bw050_6",
    "BDT_MTPeakNoBtag_T2bw025_1",
    "BDT_MTPeakNoBtag_T2bw025_3",
    "BDT_MTPeakNoBtag_T2bw025_4",
    "BDT_MTPeakNoBtag_T2bw025_6"
};

vector<string> listBDTSignalRegions_MTpeak_OneBtag =
{
    "BDT_MTPeakOneBtag_T2tt_1",
    "BDT_MTPeakOneBtag_T2tt_2",
    "BDT_MTPeakOneBtag_T2tt_5",
    "BDT_MTPeakOneBtag_T2bw075_1",
    "BDT_MTPeakOneBtag_T2bw075_2",
    "BDT_MTPeakOneBtag_T2bw075_3",
    "BDT_MTPeakOneBtag_T2bw075_5",
    "BDT_MTPeakOneBtag_T2bw050_1",
    "BDT_MTPeakOneBtag_T2bw050_3",
    "BDT_MTPeakOneBtag_T2bw050_4",
    "BDT_MTPeakOneBtag_T2bw050_5",
    "BDT_MTPeakOneBtag_T2bw050_6",
    "BDT_MTPeakOneBtag_T2bw025_1",
    "BDT_MTPeakOneBtag_T2bw025_3",
    "BDT_MTPeakOneBtag_T2bw025_4",
    "BDT_MTPeakOneBtag_T2bw025_6"
};

// ########################################
// # Regions for cut-and-count correction #
// ########################################

vector<string> listIndividualCuts =
{
    "MT_100",
    "MT_120",
    "MT_125",
    "MT_130",
    "MT_135",
    "MT_140",
    "MET_200",
    "MET_250",
    "MET_300",
    "MET_350",
    "METoverSqrtHT_6",
    "METoverSqrtHT_7",
    "METoverSqrtHT_8",
    "METoverSqrtHT_9",
    "METoverSqrtHT_10",
    "METoverSqrtHT_12",
    "BPt_100",
    "BPt_180",
    "DPhi_02",
    "DPhi_08",
    "ISRJet",
    "MT2W_180",
    "MT2W_190",
    "MT2W_200"
};

vector<string> listIndividualCuts_MTtail =
{
    "CR0btag_MTtail_MT_100",
    "CR0btag_MTtail_MT_120",
    "CR0btag_MTtail_MT_125",
    "CR0btag_MTtail_MT_130",
    "CR0btag_MTtail_MT_135",
    "CR0btag_MTtail_MT_140",

    "CR0btag_MTtail_MET_200",
    "CR0btag_MTtail_MET_250",
    "CR0btag_MTtail_MET_300",
    "CR0btag_MTtail_MET_350",

    "CR0btag_MTtail_METoverSqrtHT_6",
    "CR0btag_MTtail_METoverSqrtHT_7",
    "CR0btag_MTtail_METoverSqrtHT_8",
    "CR0btag_MTtail_METoverSqrtHT_9",
    "CR0btag_MTtail_METoverSqrtHT_10",
    "CR0btag_MTtail_METoverSqrtHT_12",

    "CR0btag_MTtail_BPt_100",
    "CR0btag_MTtail_BPt_180",

    "CR0btag_MTtail_DPhi_02",
    "CR0btag_MTtail_DPhi_08",

    "CR0btag_MTtail_ISRJet",

    "CR0btag_MTtail_MT2W_180",
    "CR0btag_MTtail_MT2W_190",
    "CR0btag_MTtail_MT2W_200"
};

vector<string>  listIndividualCuts_MTpeak =
{
    "0btag_MTpeak",
    "0btag_MTpeak",
    "0btag_MTpeak",
    "0btag_MTpeak",
    "0btag_MTpeak",
    "0btag_MTpeak",

    "CR0btag_MTpeak_MET_200",
    "CR0btag_MTpeak_MET_250",
    "CR0btag_MTpeak_MET_300",
    "CR0btag_MTpeak_MET_350",

    "CR0btag_MTpeak_METoverSqrtHT_6",
    "CR0btag_MTpeak_METoverSqrtHT_7",
    "CR0btag_MTpeak_METoverSqrtHT_8",
    "CR0btag_MTpeak_METoverSqrtHT_9",
    "CR0btag_MTpeak_METoverSqrtHT_10",
    "CR0btag_MTpeak_METoverSqrtHT_12",

    "CR0btag_MTpeak_BPt_100",
    "CR0btag_MTpeak_BPt_180",

    "CR0btag_MTpeak_DPhi_02",
    "CR0btag_MTpeak_DPhi_08",

    "CR0btag_MTpeak_ISRJet",

    "CR0btag_MTpeak_MT2W_180",
    "CR0btag_MTpeak_MT2W_190",
    "CR0btag_MTpeak_MT2W_200"
};

vector<string> listCutAndCounts =
{
    "preselection",
    "cutAndCount_T2tt_offShellLoose",
    "cutAndCount_T2tt_offShellTight",
    "cutAndCount_T2tt_lowDeltaM",
    "cutAndCount_T2tt_mediumDeltaM",
    "cutAndCount_T2tt_highDeltaM",

    "cutAndCount_T2bw025_offShell",
    "cutAndCount_T2bw025_lowMasses",
    "cutAndCount_T2bw025_highMasses",

    "cutAndCount_T2bw050_offShell",
    "cutAndCount_T2bw050_lowMasses",
    "cutAndCount_T2bw050_mediumDeltaM",
    "cutAndCount_T2bw050_highDeltaM",

    "cutAndCount_T2bw075_lowDeltaM",
    "cutAndCount_T2bw075_mediumDeltaM",
    "cutAndCount_T2bw075_highDeltaM"
};

map<string,vector<string> > listCutAndCounts_cuts;

void initCutAndCountCuts()
{
    listCutAndCounts_cuts["preselection"]                     = { "MT_100"                     };
    listCutAndCounts_cuts["cutAndCount_T2tt_offShellLoose"]   = { "MT_125"                     };
    listCutAndCounts_cuts["cutAndCount_T2tt_offShellTight"]   = { "MT_130", "MET_300"          };
    listCutAndCounts_cuts["cutAndCount_T2tt_lowDeltaM"]       = { "MT_140", "METoverSqrtHT_8"  };
    listCutAndCounts_cuts["cutAndCount_T2tt_mediumDeltaM"]    = { "MT_140", "MET_200"          };
    listCutAndCounts_cuts["cutAndCount_T2tt_highDeltaM"]      = { "MT_130", "MET_300"          }; // (FIXME) NB : real MET cut here is 350 ...

    listCutAndCounts_cuts["cutAndCount_T2bw025_offShell"]     = { "MT_120", "METoverSqrtHT_9"  };
    listCutAndCounts_cuts["cutAndCount_T2bw025_lowMasses"]    = { "MT_120", "METoverSqrtHT_6"  };
    listCutAndCounts_cuts["cutAndCount_T2bw025_highMasses"]   = { "MT_120", "MET_300"          };

    listCutAndCounts_cuts["cutAndCount_T2bw050_offShell"]     = { "MT_120", "METoverSqrtHT_9"  };
    listCutAndCounts_cuts["cutAndCount_T2bw050_lowMasses"]    = { "MT_135", "METoverSqrtHT_6"  };
    listCutAndCounts_cuts["cutAndCount_T2bw050_mediumDeltaM"] = { "MT_140", "METoverSqrtHT_7"  };
    listCutAndCounts_cuts["cutAndCount_T2bw050_highDeltaM"]   = { "MT_120", "MET_300"          };

    listCutAndCounts_cuts["cutAndCount_T2bw075_lowDeltaM"]    = { "MT_120", "METoverSqrtHT_12" };
    listCutAndCounts_cuts["cutAndCount_T2bw075_mediumDeltaM"] = { "MT_130", "METoverSqrtHT_10" };
    listCutAndCounts_cuts["cutAndCount_T2bw075_highDeltaM"]   = { "MT_140", "MET_300"          };
}
*/

std::pair<double,double>  GetSF(RooFitResult* res, string param)
{

    RooRealVar* par_init = (RooRealVar*) res->floatParsInit().find(param.c_str());
    RooRealVar* par_final = (RooRealVar*) res->floatParsFinal().find(param.c_str());
    double vinit = par_init->getVal();
    double vfinal = par_final->getVal();
    double SF = 1.;
    if(par_init!=0) SF = vfinal/vinit;

    double SFerror = 0.;
    if(vfinal!=0) SFerror = SF*par_final->getError()/vfinal;
    cout << "vfinal, parfinalGetErro, SF = " << vfinal << " ; " << par_final->getError() << " ; " << SF << endl;
    cout << "SF  = " << SF << endl;
    cout << "SF error = " << SFerror << endl;
    std::pair<double,double> pSF(SF,SFerror);
    return pSF;
}

//----------  Retrieve pdf (norm distrib) -----------------//

TH1F* GetHisto(TFile* fin, string region, string process, string varname, float& norm)
{
    string cname = REGION_FOR_FIT+string("/")+region+"/"+varname;
    TCanvas* c = (TCanvas*) fin->Get(cname.c_str());
    string hname = "v:"+varname+"|p:"+process+"|r:"+region+"|c:"+REGION_FOR_FIT+"|t:1DEntries";
    TH1F* h = 0;
    TList* l = c->GetListOfPrimitives();
    TPad* pad = (TPad*) l->At(0);
    THStack* stack = (THStack*) pad->GetPrimitive("");
    h = (TH1F*) stack->GetHists()->FindObject(hname.c_str());
    cout<<hname<<endl;
    cout<<"h"<<endl;
    norm = h->Integral();
    // NEW temporary
    ///*
    TH1F* htmp = (TH1F*) h->Clone();
    htmp->SetName(varname.c_str());
    TFile* fout = new TFile((process+string(".root")).c_str(),"RECREATE");
    fout->cd();
    htmp->Write();
    fout->Write();
    fout->Close();
    //*/
    //----
    return (TH1F*) h->Clone();
}

RooHistPdf* GetRooHistPdf(TFile* fin, string region, string process, string varname, RooRealVar* var, float& norm, bool do_mcstat)
{
    TH1F* h  =  GetHisto(fin,region,process,varname, norm);
    if(do_mcstat)
    {
        //randomisation of the histo ...
        for(int i=1;i<=h->GetNbinsX();i++)
        {
            h->SetBinContent(i,randomnessGenerator->Gaus(h->GetBinContent(i),h->GetBinError(i)));
        }
    }
    string rdhname = "rdh_"+region+"_"+process;
    RooDataHist *rdh  = new RooDataHist(rdhname.c_str(),rdhname.c_str(),RooArgList(*var),Import(*h));
    string pdfname = "pdf_"+region+"_"+process;
    RooHistPdf *pdf  = new  RooHistPdf(pdfname.c_str(),pdfname.c_str(),RooArgSet(*var),*rdh);
    return pdf;
}
//---------------------------------------------------------//

//----------  Retrieve data histo  -----------------
TH1F* GetData(TFile* fin, string region, string varname)
{
    string cname = REGION_FOR_FIT+string("/")+region+"/"+varname;
    TCanvas* c = (TCanvas*) fin->Get(cname.c_str());
    TList* l = c->GetListOfPrimitives();
    TPad* pad = (TPad*) l->At(0);
    string hname = "v:"+varname+"|r:"+region+"|c:"+REGION_FOR_FIT+"|t:1DSumData";
    TH1F* h = (TH1F*) pad->GetPrimitive(hname.c_str());
    return (TH1F*) h->Clone();
}

RooDataHist* GetRooData(TFile* fin, string region, string varname, RooRealVar* var)
{
    TH1F* h = GetData(fin, region, varname);
    string dname = "data_"+region;
    RooDataHist *datah = new RooDataHist(dname.c_str(),dname.c_str(), RooArgList(*var), Import(*h));
    return datah;
}

//processes and norm should have the same size
RooDataHist* GetRooToyData(TFile* fin, string region, string varname, RooRealVar* var, vector<string> processes, vector<float> norm)
{
    TH1F* h = 0;
    for(unsigned int i=0;i<processes.size();i++){
    	TH1* htmp  =  GetHisto(fin,region,processes[i],varname, norm[i]);
    	if(i==0) h = (TH1F*) htmp->Clone();
	else h->Add(htmp);
    }
    string dname = "data_"+region;
    RooDataHist *datah = new RooDataHist(dname.c_str(),dname.c_str(), RooArgList(*var), Import(*h));
    return datah;
}

struct FitSetup
{
    //--name of the root file --//
    string filename;
    //--- variable used to the fit --//
    string varname;
    float  varMin;
    float  varMax;
    //--- region used to the fit ---//
    string region;
    //--- Xsection uncertainties --//
    float  xs_sysfactor;
    //bool   do_xs_tt2l_sys;
    bool   do_xs_DY_sys;
    //--- MC stat uncertainties --//
    bool do_mcstat;
    //-- algorithm used to the fit --//
    string type;
    string algo;
    //-- Init value uncert --//
    bool do_init_uncert;
    float init_ttbar;
    float init_stop;

    void Reset()
    {
        filename=string(INPUT_FOLDER)+"/1DDataMCComparison.root";
        varname=OBSERVABLE_FOR_FIT;
        varMin = 0;
        varMax = 1000;
        region = "0btag_MTtail";

        //--- Xsection uncertainties --// (or at least, the one who store them later since these variables don't say anything at all anyway...)
        xs_sysfactor=1.;
        //do_xs_tt2l_sys = false;
        do_xs_DY_sys = false;

        //--- MC stat uncertainties --//
        do_mcstat = false;

        //-- algorithm used to the fit --//
        type = "Minuit2";
        algo = "MIGRAD";

        //-- Init value uncert --//
        do_init_uncert = false;
        init_ttbar=1.;
        init_stop=1.;
    }
};

struct FitResult
{
    string conditions;
    pair<float,float> SF_ttbar;
    pair<float,float> SF_stop;
    float norm_ttbar;
    float correlation;
    float edm;

    void Reset()
    {
        norm_ttbar = 0;
        edm = 0;
        SF_ttbar = pair<float,float>(0,0);
        SF_stop = pair<float,float>(0,0);
        correlation = 0;
    }
    void Print()
    {
        cout<<"#################"<<endl;
        cout<<"# "<<conditions<<"\t SF_ttbar = "<<SF_ttbar.first<<"+/-"<<SF_ttbar.second<<"\t SF_stop = "<<SF_stop.first<<"+/-"<<SF_stop.second<<endl;
        cout<<"#\t edm = "<<edm<<" correlation = "<<correlation<<endl;
        cout<<"#\t init ttbar: "<<norm_ttbar<<" \t fitted ttbar: "<<norm_ttbar*SF_ttbar.first<<endl;
        cout<<"#################"<<endl;
    }
};


FitResult  doFit(const FitSetup& setup, string conditions, string fname=string(""))
{

    string varname = setup.varname;
    RooRealVar var(varname.c_str(),varname.c_str(),setup.varMin, setup.varMax);

    //string region="0btag_MTtail";
    string region= setup.region;

    //should it be an argument ?
    TFile* fin = 0;
    if(fname=="") fin = TFile::Open(setup.filename.c_str());
    //else fin = TFile::Open(fname.c_str());
    else fin = TFile::Open(fname.c_str());

    //-- normalisation in the MC --//
    float mc_norm_ttbar = 0;
    //float mc_norm_tt2l = 0;
    float mc_norm_stop = 0;
    float mc_norm_DY = 0;

    // C r e a t e   m o d e l   f o r  CR1_peak_lowM3b
    // -------------------------------------------------------------
    // Construct pdfs for ttbar, tt2l, stop and DY
    RooHistPdf *pdf_ttbar  = GetRooHistPdf(fin,region,"ttbar",varname,&var,mc_norm_ttbar, setup.do_mcstat);
    //RooHistPdf *pdf_tt2l   = GetRooHistPdf(fin,region,"ttbar_2l",varname,&var,mc_norm_tt2l, setup.do_mcstat);
    RooHistPdf *pdf_stop  = GetRooHistPdf(fin,region,"stop",varname,&var,mc_norm_stop, setup.do_mcstat);
    //RooHistPdf *pdf_DY   = GetRooHistPdf(fin,region,"DY",varname,&var,mc_norm_DY, setup.do_mcstat);

    // normalization factors (RooRealVar)
    float val_ttbar = mc_norm_ttbar;
    float val_stop = mc_norm_stop;
    if(setup.do_init_uncert)
    {
        val_ttbar = setup.init_ttbar*mc_norm_ttbar;
        val_stop = setup.init_stop*mc_norm_stop;
    }
    cout<<"normalization: ttbar = "<<mc_norm_ttbar<<" stop = "<<mc_norm_stop<<endl;
    RooRealVar norm_ttbar("norm_ttbar","norm_ttbar",val_ttbar,0.25*mc_norm_ttbar,10.*mc_norm_ttbar);
    RooRealVar norm_stop("norm_stop","norm_stop",val_stop,0.25*mc_norm_stop,10.*mc_norm_stop);
    //RooRealVar norm_tt2l("norm_tt2l","norm_tt2l",mc_norm_tt2l,0.25*mc_norm_tt2l,2*mc_norm_tt2l);
    //RooRealVar norm_DY("norm_DY","norm_DY",mc_norm_DY,0.25*mc_norm_DY,2*mc_norm_DY);
    // possibility to study a systematic on it
    //if(setup.do_xs_tt2l_sys) mc_norm_tt2l*=setup.xs_sysfactor;
    //if(setup.do_xs_DY_sys) mc_norm_DY*=setup.xs_sysfactor;
    //RooConstVar norm_tt2l("norm_tt2l","norm_tt2l",mc_norm_tt2l);
    //RooConstVar norm_DY("norm_DY","norm_DY",mc_norm_DY);

    /*
    RooAddPdf model("model","model",
            RooArgList(*pdf_ttbar,*pdf_tt2l,*pdf_stop,*pdf_DY),
            RooArgList(norm_ttbar,norm_tt2l,norm_stop,norm_DY)) ;
    */
    /*
    RooAddPdf model("model","model",
            RooArgList(*pdf_ttbar,*pdf_stop,*pdf_DY),
            RooArgList(norm_ttbar,norm_stop,norm_DY)) ;
    */
    RooAddPdf model("model","model",
            RooArgList(*pdf_ttbar,*pdf_stop),
            RooArgList(norm_ttbar,norm_stop)) ;


    //RooDataHist *data_CR1_peak_lowM3b = GetRooData(fin,region,varname,&var);
    
    vector<string> processes = {"ttbar","stop"};
    vector<float> norms = {mc_norm_ttbar,mc_norm_stop};
    RooDataHist* data_CR1_peak_lowM3b =  GetRooToyData(fin,region,varname,&var, processes, norms);

    fin->Close();


    //--  Constraints on single top and DY --//
    float RelUncert = 0.2;
    // Construct another Gaussian constraint p.d.f on "DY" bkg
    //RooGaussian constr_DY("constr_DY","constr_DY",norm_DY,RooConst(mc_norm_DY),RooConst(RelUncert*mc_norm_DY)) ;
    // Construct another Gaussian constraint p.d.f on "tt2l" bkg
    //RooGaussian constr_tt2l("constr_tt2l","constr_tt2l",norm_tt2l,RooConst(mc_norm_tt2l),RooConst(RelUncert*mc_norm_tt2l)) ;

    // P e r f o r m   t em p l a t e   f i t
    // ---------------------------------------------------

    //Minimizer(type,algo) -- Choose minimization package and algorithm to use. Default is MINUIT/MIGRAD through the RooMinimizer
    //                       interface, but DY can be specified (through RooMinimizer interface). Select OldMinuit to use
    //                       MINUIT through the old RooMinuit interface
    //
    //     Type         Algorithm
    //     ------       ---------
    //     OldMinuit    migrad, simplex, minimize (=migrad+simplex), migradimproved (=migrad+improve)
    //     Minuit       migrad, simplex, minimize (=migrad+simplex), migradimproved (=migrad+improve)
    //     Minuit2      migrad, simplex, minimize, scan
    //     GSLMultiMin  conjugatefr, conjugatepr, bfgs, bfgs2, steepestdescent
    //     GSLSimAn     -


    // ---  Perform simultaneous fit of model to data and model_ctl to data_ctl --//
    //RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,Save());
    //RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,ExternalConstraints(constr_DY),ExternalConstraints(constr_tt2l),PrintLevel(-1),Save(),
    //RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,ExternalConstraints(constr_DY),PrintLevel(-1),Save(),
    RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,PrintLevel(-1),Save(),
            Minimizer(setup.type.c_str(),setup.algo.c_str()));

    //--- Writing the results ---///
    FitResult fitRes;
    fitRes.Reset();
    fitRes.norm_ttbar = mc_norm_ttbar;
    fitRes.SF_ttbar = GetSF(res,"norm_ttbar");
    fitRes.SF_stop = GetSF(res,"norm_stop");
    fitRes.edm = res->edm();
    fitRes.correlation = res->correlationMatrix()[0][1];
    fitRes.conditions = conditions;

    return fitRes;

}

int main()
{
    randomnessGenerator = new TRandom();

    system((string("mkdir -p ")+OUTPUT_FOLDER).c_str());

    // ###########################
    // # Prepare final SFR table #
    // ###########################

    /*
    vector<string> columns = { "SFR_ttbar", "SFR_stop" };

    
    vector<string> listAllSignalRegion = listBDTSignalRegions;
                   listAllSignalRegion.insert(listAllSignalRegion.end(), listCutAndCounts.begin(), listCutAndCounts.end());

    vector<string> listRawRegion = listBDTSignalRegions;
                   listRawRegion.push_back("BDT_average");
                   listRawRegion.insert(listRawRegion.end(), listIndividualCuts.begin(), listIndividualCuts.end());

    Table tableSFRToBeUsed(columns,listAllSignalRegion);
    Table tableRawSFR(columns,listRawRegion);
   */
    // Create observables
    //RooRealVar var(OBSERVABLE_FOR_FIT,OBSERVABLE_FOR_FIT,0,TMath::Pi()) ;
    RooRealVar var(OBSERVABLE_FOR_FIT,OBSERVABLE_FOR_FIT,0,25) ;
    string varname(OBSERVABLE_FOR_FIT);

    FitSetup setup;
    FitResult res;
    string conditions;

    // ########################
    // #  ____________ _____  #
    // #  | ___ \  _  \_   _| #
    // #  | |_/ / | | | | |   #
    // #  | ___ \ | | | | |   #
    // #  | |_/ / |/ /  | |   #
    // #  \____/|___/   \_/   #
    // #                      #
    // ########################

    // Create histos
    /*
    TH1F h_SF_MTpeak_BDT_ttbar ("h_SF_MTpeak_BDT_ttbar", "", listBDTSignalRegions_MTtail.size(),0,listBDTSignalRegions_MTtail.size());
    TH1F h_SF_MTtail_BDT_ttbar ("h_SF_MTtail_BDT_ttbar", "", listBDTSignalRegions_MTtail.size(),0,listBDTSignalRegions_MTtail.size());
    TH1F h_SFR_BDT_ttbar       ("h_SFR_BDT_ttbar",       "", listBDTSignalRegions_MTtail.size(),0,listBDTSignalRegions_MTtail.size());
    TH1F h_SF_MTpeak_BDT_stop ("h_SF_MTpeak_BDT_stop", "", listBDTSignalRegions_MTtail.size(),0,listBDTSignalRegions_MTtail.size());
    TH1F h_SF_MTtail_BDT_stop ("h_SF_MTtail_BDT_stop", "", listBDTSignalRegions_MTtail.size(),0,listBDTSignalRegions_MTtail.size());
    TH1F h_SFR_BDT_stop       ("h_SFR_BDT_stop",       "", listBDTSignalRegions_MTtail.size(),0,listBDTSignalRegions_MTtail.size());
    
    float mean_SFttbar_value = 0;
    float mean_SFttbar_error = 0;
    float rms_SFttbar_value = 0;

    float mean_SFstop_value = 0;
    float mean_SFstop_error = 0;
    float rms_SFstop_value = 0;
    */

    /*
    for(unsigned int i=0;i<listBDTSignalRegions_MTtail.size();i++)
    {
        cout<<"%%%%%%%%%%%%%%%%%% "<<listBDTSignalRegions_MTtail[i]<<endl;

        string label = listBDTSignalRegions[i];
        */
	Figure SFR_ttbar;
        Figure SFR_stop;

        //MT tail
        setup.Reset(); conditions="sigRegions_tail";
        setup.region = "SR_2l";
	//setup.region=listBDTSignalRegions_MTtail[i];
        setup.varname=varname;
        setup.varMin=0;
        //setup.varMax=TMath::Pi();
        setup.varMax=25;

        res = doFit(setup,conditions);

        Figure SF_ttbar_tail = Figure(res.SF_ttbar.first,res.SF_ttbar.second);
        Figure SF_stop_tail = Figure(res.SF_stop.first,res.SF_stop.second);

	//Loop over the regions
	/*
	vector<string> regions = {"SR_2l","SR_2l_2j","SR_2l_3j","SR_2l_4j","SR_2l_ISR"};
	for(int r = 0; r<regions.size();r++){
		setup.region = regions[r];
        	res = doFit(setup,conditions);

	}
	*/


        /*
	h_SF_MTtail_BDT_ttbar.SetBinContent(i+1,res.SF_ttbar.first);
        h_SF_MTtail_BDT_ttbar.SetBinError(i+1,res.SF_ttbar.second);
        h_SF_MTtail_BDT_ttbar.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SF_MTtail_BDT_stop.SetBinContent(i+1,res.SF_stop.first);
        h_SF_MTtail_BDT_stop.SetBinError(i+1,res.SF_stop.second);
        h_SF_MTtail_BDT_stop.GetXaxis()->SetBinLabel(i+1,label.c_str());
	*/
	
        //MT peak
        /*
	setup.Reset(); conditions="sigRegions_peak"; setup.region=listBDTSignalRegions_MTpeak[i];
        setup.varname=varname;
        setup.varMin=0;
        setup.varMax=600;

        res = doFit(setup,conditions);

        Figure SF_ttbar_peak = Figure(res.SF_ttbar.first,res.SF_ttbar.second);
        Figure SF_stop_peak = Figure(res.SF_stop.first,res.SF_stop.second);

        h_SF_MTpeak_BDT_ttbar.SetBinContent(i+1,res.SF_ttbar.first);
        h_SF_MTpeak_BDT_ttbar.SetBinError(i+1,res.SF_ttbar.second);
        h_SF_MTpeak_BDT_ttbar.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SF_MTpeak_BDT_stop.SetBinContent(i+1,res.SF_stop.first);
        h_SF_MTpeak_BDT_stop.SetBinError(i+1,res.SF_stop.second);
        h_SF_MTpeak_BDT_stop.GetXaxis()->SetBinLabel(i+1,label.c_str());

        //Now compute the ration : SF_tail/SF_peak
        SFR_ttbar = SF_ttbar_tail / SF_ttbar_peak;
        SFR_stop = SF_stop_tail / SF_stop_peak;

        //-- do some additionnal test as function of the b-tag multiplicity

        //MT peak (no btag req)

        setup.Reset(); conditions="sigRegions_peak_NoBtag";
        setup.region=listBDTSignalRegions_MTpeak_NoBtag[i];
        setup.varMin=0;
        setup.varMax=600;

        res = doFit(setup,conditions);

        //MT peak (one btag req)

        setup.Reset(); conditions="sigRegions_peak_OneBtag";
        setup.region=listBDTSignalRegions_MTpeak_OneBtag[i];
        setup.varname=varname;
        setup.varMin=0;
        setup.varMax=600;

        res = doFit(setup,conditions);

        //Computation of mean/rms/ ..
        //It is based on  the ratio SF_tail/SF_peak

        //-- ttbar
        mean_SFttbar_value += SFR_ttbar.value();
        mean_SFttbar_error += SFR_ttbar.error();
        rms_SFttbar_value  += (SFR_ttbar.value()*SFR_ttbar.value());

        //-- stop
        mean_SFstop_value += SFR_stop.value();
        mean_SFstop_error += SFR_stop.error();
        rms_SFstop_value  += (SFR_stop.value()*SFR_stop.value());

        //SFR
        h_SFR_BDT_ttbar.SetBinContent(i+1,SFR_ttbar.value());
        h_SFR_BDT_ttbar.SetBinError(i+1,SFR_ttbar.error());
        h_SFR_BDT_ttbar.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SFR_BDT_stop.SetBinContent(i+1,SFR_stop.value());
        h_SFR_BDT_stop.SetBinError(i+1,SFR_stop.error());
        h_SFR_BDT_stop.GetXaxis()->SetBinLabel(i+1,label.c_str());

        cout << listBDTSignalRegions[i] << " : SF_ttbar_peak: " << SF_ttbar_peak.Print() <<" SF_stop_peak: "<<SF_stop_tail.Print()<<endl;
        cout << listBDTSignalRegions[i] << " : SF_ttbar_tail: " << SF_ttbar_tail.Print() <<" SF_stop_tail: "<<SF_stop_peak.Print()<<endl;
        cout << listBDTSignalRegions[i] << " : SFR_ttbar: "     << SFR_ttbar.Print()     <<" SFR_stop: "<<SFR_stop.Print()<<endl;

        tableRawSFR.Set("SFR_ttbar",listBDTSignalRegions[i],SFR_ttbar);
        tableRawSFR.Set("SFR_stop",listBDTSignalRegions[i],SFR_stop);
    }

    //Save plots in roofile
    TFile fCR1_BDT((string(OUTPUT_FOLDER)+"/results_BDT.root").c_str(),"RECREATE");
    h_SF_MTpeak_BDT_ttbar.Write();
    h_SF_MTpeak_BDT_stop.Write();
    h_SF_MTtail_BDT_ttbar.Write();
    h_SF_MTtail_BDT_stop.Write();
    h_SFR_BDT_ttbar.Write();
    h_SFR_BDT_stop.Write();


    //---------------------------------------//

    mean_SFttbar_value /= listBDTSignalRegions_MTtail.size();
    mean_SFttbar_error /= listBDTSignalRegions_MTtail.size();

    rms_SFttbar_value  /= listBDTSignalRegions_MTtail.size();
    rms_SFttbar_value  -= (mean_SFttbar_value*mean_SFttbar_value);

    Figure BDT_SFttbar(mean_SFttbar_value,sqrt(rms_SFttbar_value
                                              +mean_SFttbar_error*mean_SFttbar_error
                                              +pow(TEMPLATE_FIT_METHOD_UNCERTAINTY*mean_SFttbar_value,2)));
    //---------------------------------------//

    mean_SFstop_value /= listBDTSignalRegions_MTtail.size();
    mean_SFstop_error /= listBDTSignalRegions_MTtail.size();

    rms_SFstop_value  /= listBDTSignalRegions_MTtail.size();
    rms_SFstop_value  -= (mean_SFstop_value*mean_SFstop_value);

    Figure BDT_SFstop(mean_SFstop_value,sqrt(rms_SFstop_value
                                              +mean_SFstop_error*mean_SFstop_error
                                              +pow(TEMPLATE_FIT_METHOD_UNCERTAINTY*mean_SFstop_value,2)));
    //---------------------------------------//

    tableRawSFR.Set("SFR_ttbar","BDT_average",BDT_SFttbar);
    tableRawSFR.Set("SFR_stop","BDT_average",BDT_SFstop);
    for(unsigned int i=0;i<listBDTSignalRegions.size();i++)
    {
        tableSFRToBeUsed.Set("SFR_ttbar",listBDTSignalRegions[i],BDT_SFttbar);
        tableSFRToBeUsed.Set("SFR_stop",listBDTSignalRegions[i],BDT_SFstop);
    }

    // ########################
    // #  _____        _____  #
    // # /  __ \ ___  /  __ \ #
    // # | /  \/( _ ) | /  \/ #
    // # | |    / _ \/\ |     #
    // # | \__/\ (_>  < \__/\ #
    // #  \____/\___/\/\____/ #
    // #                      #
    // ########################

    std::map<string,Figure> SFR_CC_ttbar_map;
    std::map<string,Figure> SFR_CC_stop_map;

    //Create histos
    TH1F h_SF_MTpeak_CC_ttbar ("h_SF_MTpeak_CC_ttbar",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SF_MTtail_CC_ttbar ("h_SF_MTtail_CC_ttbar",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SFR_CC_ttbar       ("h_SFR_CC_ttbar",        "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SF_MTpeak_CC_stop ("h_SF_MTpeak_CC_stop",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SF_MTtail_CC_stop ("h_SF_MTtail_CC_stop",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SFR_CC_stop       ("h_SFR_CC_stop",        "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());

    for(unsigned int i=0;i<listIndividualCuts.size();i++)
    {
        cout<<"%%%%%%%%%%%%%%%%%% "<<listIndividualCuts_MTtail[i]<<endl;

        string label = listIndividualCuts[i];
        Figure SFR_ttbar;
        Figure SFR_stop;

        //MT tail
        setup.Reset();
        conditions="sigRegions_CC_tail";
        setup.region=listIndividualCuts_MTtail[i];
        setup.varname=varname;
        setup.varMin=0;
        setup.varMax=600;

        res = doFit(setup,conditions);
        SFR_ttbar=Figure(res.SF_ttbar.first,res.SF_ttbar.second);
        SFR_stop=Figure(res.SF_stop.first,res.SF_stop.second);
        h_SF_MTtail_CC_ttbar.SetBinContent(i+1,res.SF_ttbar.first);
        h_SF_MTtail_CC_ttbar.SetBinError(i+1,res.SF_ttbar.second);
        h_SF_MTtail_CC_ttbar.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SF_MTtail_CC_stop.SetBinContent(i+1,res.SF_stop.first);
        h_SF_MTtail_CC_stop.SetBinError(i+1,res.SF_stop.second);
        h_SF_MTtail_CC_stop.GetXaxis()->SetBinLabel(i+1,label.c_str());

        //MT peak
        setup.Reset();
        conditions="sigRegions_CC_peak";
        setup.region=listIndividualCuts_MTpeak[i];
        setup.varname=varname;
        setup.varMin=0;
        setup.varMax=600;

        res = doFit(setup,conditions);
        h_SF_MTpeak_CC_ttbar.SetBinContent(i+1,res.SF_ttbar.first);
        h_SF_MTpeak_CC_ttbar.SetBinError(i+1,res.SF_ttbar.second);
        h_SF_MTpeak_CC_ttbar.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SF_MTpeak_CC_stop.SetBinContent(i+1,res.SF_stop.first);
        h_SF_MTpeak_CC_stop.SetBinError(i+1,res.SF_stop.second);
        h_SF_MTpeak_CC_stop.GetXaxis()->SetBinLabel(i+1,label.c_str());

        // Now compute the ration : SF_tail/SF_peak
        SFR_ttbar /= Figure(res.SF_ttbar.first, res.SF_ttbar.second);
        SFR_stop /= Figure(res.SF_stop.first, res.SF_stop.second);

        // It is based on the ratio SF_tail/SF_peak
        SFR_CC_ttbar_map[label] = SFR_ttbar;
        SFR_CC_stop_map[label] = SFR_stop;

        // SFR
        h_SFR_CC_ttbar.SetBinContent(i+1,SFR_ttbar.value());
        h_SFR_CC_ttbar.SetBinError(i+1,SFR_ttbar.error());
        h_SFR_CC_ttbar.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SFR_CC_stop.SetBinContent(i+1,SFR_stop.value());
        h_SFR_CC_stop.SetBinError(i+1,SFR_stop.error());
        h_SFR_CC_stop.GetXaxis()->SetBinLabel(i+1,label.c_str());

        cout<<"individual cut: " << listIndividualCuts[i] << " ; SFR_ttbar: "<<SFR_ttbar.Print()<<" SFR_stop: "<<SFR_stop.Print()<<endl;

        tableRawSFR.Set("SFR_ttbar",listIndividualCuts[i],SFR_ttbar);
        tableRawSFR.Set("SFR_stop",listIndividualCuts[i],SFR_stop);
    }

    //Save plots in roofile
    TFile fCR1_CC((string(OUTPUT_FOLDER)+"/results_CC.root").c_str(),"RECREATE");
    h_SF_MTpeak_CC_ttbar.Write();
    h_SF_MTpeak_CC_stop.Write();
    h_SF_MTtail_CC_ttbar.Write();
    h_SF_MTtail_CC_stop.Write();
    h_SFR_CC_ttbar.Write();
    h_SFR_CC_stop.Write();

    //---------------------------------------------
    // Results for C&C
    //---------------------------------------------
    initCutAndCountCuts();
    vector<string> cuts;

    for(unsigned int r=0;r<listCutAndCounts.size();r++)
    {
        cuts = listCutAndCounts_cuts[listCutAndCounts[r]];
        Figure SFR_CC_ttbar;
        Figure SFR_CC_stop;
        for(unsigned i=0;i<cuts.size();i++)
        {
            cout<<cuts[i]<<" "<<SFR_CC_ttbar_map[cuts[i]].Print()<<endl;
            if(i == 0) SFR_CC_ttbar = SFR_CC_ttbar_map[cuts[i]];
            else SFR_CC_ttbar = Figure(SFR_CC_ttbar.value(),
                                       sqrt(pow(SFR_CC_ttbar.error(),2)
                                           +pow(SFR_CC_ttbar.value()
                                                -SFR_CC_ttbar_map[cuts[i].c_str()].value(),2)
                                           ));

            if(i == 0) SFR_CC_stop = SFR_CC_stop_map[cuts[i]];
            else SFR_CC_stop = Figure(SFR_CC_stop.value(),
                                       sqrt(pow(SFR_CC_stop.error(),2)
                                           +pow(SFR_CC_stop.value()
                                               -SFR_CC_stop_map[cuts[i].c_str()].value(),2)
                                           ));
        }

        // Add 20% of uncertainty for fit itself (JES, MC stat. ...)
        SFR_CC_ttbar *= Figure(1.0,TEMPLATE_FIT_METHOD_UNCERTAINTY);
        SFR_CC_stop *= Figure(1.0,TEMPLATE_FIT_METHOD_UNCERTAINTY);

        tableSFRToBeUsed.Set("SFR_ttbar",listCutAndCounts[r],SFR_CC_ttbar);
        tableSFRToBeUsed.Set("SFR_stop",listCutAndCounts[r],SFR_CC_stop);
    }

    tableRawSFR     .Print(string(OUTPUT_FOLDER)+"/rawSFR.tab"   ,4);
    tableSFRToBeUsed.Print(string(OUTPUT_FOLDER)+"/SF_MTtail.tab",4);

   */
}
