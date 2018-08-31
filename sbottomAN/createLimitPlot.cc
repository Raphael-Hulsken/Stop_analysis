#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
#include <cmath>
#include <exception>
#include <ctime>

#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THashList.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"

using namespace std;

int main(int argc, char *argv[])
{

    if(argc != 2)
        throw std::runtime_error("Bad number of arguments!");

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kMint);
    gROOT->ForceStyle();

    TString file1 = argv[1];

    TFile* f1 = NULL;
    TTree* t1 = NULL;
    f1 = TFile::Open(file1);
    if(f1==NULL)
        throw std::runtime_error("File 1 address not set");
    t1 = dynamic_cast<TTree*>(f1->Get("limit"));
    if(t1==NULL)
        throw std::runtime_error("Tree 1 address not set");

    double limit = -13;
    float mSbot = -13;
    float mNeutr = -13;
    t1->SetBranchAddress("limit", &limit );
    t1->SetBranchAddress("sbot_mass_point", &mSbot );
    t1->SetBranchAddress("neutr_mass_point", &mNeutr );

    TH2F* limitPlot = new TH2F("limitPlot", "limitPlot", 150, 0 , 1500, 130, 0, 1300); //@MJ@ TODO do smething more general
   
   Int_t nentries = (Int_t)t1->GetEntries();
   for (Int_t e=0; e<nentries; e++)
   {
       t1->GetEntry(e);

       if(e==2 || (e-2)%6==0)
       {
           cout << "entry " << e << " limit " <<limit << " for point: sbottom " << mSbot << " neutralino " << mNeutr << endl;
           //if(limit < 1)
           limitPlot->Fill(mSbot,mNeutr,limit);
 
       }

   }

   TCanvas can("limits","limits");
   limitPlot->Draw("colz");
   can.SaveAs("limitPlot.root");

} 
