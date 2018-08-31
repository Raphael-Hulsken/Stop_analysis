#include <iostream>
#include <exception>
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TRint.h"
using namespace std;

int main(int argc, char *argv[]){

    if(argc != 4)
        throw std::runtime_error("Bad number of arguments!");



    TString file1 = argv[1];
    TString file2 = argv[2];
    TString fileo = argv[3];


    TFile *f1 = new TFile(file1);
    f1->ls();
    TH2D * h1 = (TH2D*)f1->Get("sig");

    TFile *f2 = new TFile(file2);
    f2->ls();
    TH2D * h2 = (TH2D*)f2->Get("sig")->Clone();

    TCanvas* c = new TCanvas("c", "c");
    TH2D* output = (TH2D*)h1->Clone();
    output->Add(h2, -1);

/*   Int_t ncol = 100;
   Int_t colors[ncol];
   TColor *col;
   Double_t dg=1/(Double_t)ncol;
   Double_t grey=0;
   for (Int_t i=0; i<ncol; i++) {
      colors[i]= i+100;
      col = gROOT->GetColor(colors[i]);
      col->SetRGB(grey, grey, grey);
      grey = grey+dg;
   }
   output->SetContour(ncol);
   gStyle->SetPalette(100,colors);
  */ 


    output->Draw("colz");
    TFile fi2(fileo,"RECREATE");
    c->Write();
    output->Write();
    fi2.Close();

}
