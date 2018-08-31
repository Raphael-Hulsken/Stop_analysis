//use
// root -l -q AddDivideMultiply.C++ xx.steer nameOfProcess

#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
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

#include <fstream>
#include <iostream>
#include <vector>

void AddDivideMultiply()
{

    //steering file with information
    TString input = gApplication->Argv(5);
    //information if addition, division or multiplication required
    TString process = gApplication->Argv(6);
    
    vector<TString> in_stream;

    ifstream is;
    is.open(input);
    if( is.is_open() )
    {
	char * sline;
	TString tsline("");
	Int_t line = 0;
	string strline;
	while(getline(is, strline)) {
		tsline = strline;
		if(tsline.BeginsWith(" ")) continue;
		if(tsline.Length()==0)     continue;
		in_stream.push_back(tsline);
		line++;
	}
	is.close();
    }

    //number of vars should be same as nr of files!
    uint32_t nrOfVars = in_stream.at(0).Atof();
    TString dir = in_stream.at(1+3*nrOfVars);    


    //this will be read from steering
    vector<TString> inputRootFile;
    //this will be read from steering
    vector<TString> variable;
    //this will be read from steering
    vector<TString> inputCanvasName;

    for(uint32_t f = 0; f < nrOfVars; f++)
    {
        variable.push_back(in_stream.at(1+f));
        inputRootFile.push_back(in_stream.at(1+nrOfVars+f));
        inputCanvasName.push_back(in_stream.at(1+2*nrOfVars+f));
    }

    vector<TFile*> inputTFiles;
    vector<TH1D*> inputHist;
    vector<TCanvas*> inputCanvas;

    for(uint32_t files = 0; files < inputRootFile.size(); files++)
    {
        cout << "input file: " <<  inputRootFile.at(files) << endl;
 
        inputTFiles.push_back(new TFile(inputRootFile.at(files)));
        
        cout << "variable: " <<  variable.at(files) << endl;
        
        //inputTFiles.at(files)->cd("muon/preselection/");
       
        cout << "input canvas: " <<  inputCanvasName.at(files) << endl;
        cout << "dir " << dir << endl;
 
        inputCanvas.push_back(dynamic_cast<TCanvas*>(inputTFiles.at(files)->Get(dir + inputCanvasName.at(files))));

        inputHist.push_back(dynamic_cast<TH1D*>(inputCanvas.at(files)->GetPrimitive(variable.at(files))->Clone()));
    
    }

    cout << "nr of input files; " << inputRootFile.size() << endl;
    
    TCanvas* can = new TCanvas("can", "can");
    can->Divide(1,1);
    can->cd(1);
    
    //normalize and add
    if(process == "naddition" && inputRootFile.size() > 1)
    {   
        //normalize the first one
        Double_t norm1 = inputHist.at(0)->Integral(); 
        inputHist.at(0)->Scale(1/norm1, "width");
        for(uint32_t h = 0; h < inputHist.size(); h++)
        {   
            //normalize
            Double_t norm = inputHist.at(h+1)->Integral(); 
            inputHist.at(h+1)->Scale(1/norm, "width");
            //add
            inputHist.at(0)->Add(inputHist.at(1+h), 1);
        }

        //inputHist.at(0)->GetYaxis()->SetTitle("whatever");
        inputHist.at(0)->Draw();
    }

    //merge two histograms
    if(process == "merge")
    {  
        cout << "process merge" << endl;
        
        Double_t max = inputHist.at(0)->GetMaximum();
        inputHist.at(0)->SetMaximum(3*max);
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {   
            inputHist.at(h)->SetLineColor(h+1);
            if(h == 0)
                inputHist.at(h)->Draw("hist");
            else
                inputHist.at(h)->Draw("same hist");
        }
    }

    //efficency plus distributions
    if(process == "distrandeff" && inputRootFile.size() == 2)
    { 
        cout << "process distrandeff" << endl; 
        TH1D* division = dynamic_cast<TH1D*>(inputHist.at(0)->Clone());
        division->Divide(dynamic_cast<TH1D*>(inputHist.at(1)->Clone()));
        division->SetMaximum(2);
        division->Draw("hist");

        //inputHist.at(0)->GetYaxis()->SetTitle("whatever");
        Double_t n2 = inputHist.at(1)->Integral(); 
        cout << "gen events: "  << n2 << endl;
        inputHist.at(1)->SetLineColor(kGreen);
        inputHist.at(1)->Scale(1/n2, "width");
        inputHist.at(1)->Draw("same hist");

        Double_t n1 = inputHist.at(0)->Integral(); 
        cout << "reco events: "  << n1 << endl;
        inputHist.at(0)->SetLineColor(kBlue);
        inputHist.at(0)->Scale(1/n1, "width");
        inputHist.at(0)->Draw("hist same");
    }
  
    //ROC
    //warning -  histograms must have same amount of bins!
    if(process == "bsefficiency")
    {
        if(inputRootFile.size() == 2)
        {
            cout << "process bsefficiency" << endl;
           
            Int_t nrOfBins = inputHist.at(0)->GetNbinsX();
            TGraph* output = new TGraph(nrOfBins);
            
            //initial number of events
            Double_t norm1 = inputHist.at(0)->Integral();
            Double_t norm2 = inputHist.at(1)->Integral();

            for(uint32_t bin = 1; bin<(nrOfBins+1); bin++)
            {
                output->SetPoint(bin, (inputHist.at(0)->Integral(bin, nrOfBins))/norm1, (inputHist.at(1)->Integral(bin, nrOfBins))/norm2);
            }
            
            output->GetXaxis()->SetTitle(inputCanvasName.at(0));
            output->GetYaxis()->SetTitle(inputCanvasName.at(1));
            output->Draw("A*"); 
        }
    }
    
    //count mean and rms
    if(process == "countmean")
    {
        ofstream myfile ("TableOfMeanValues.txt");
        if (myfile.is_open())
        { 
            myfile << "process countmean" << endl;
        
            myfile << "-------mean and RMS table--------" << endl;
            myfile << "variable" << endl;      
            myfile << "        mean          RMS" << endl << endl;
            for(uint32_t h = 0; h <inputRootFile.size(); h++ )
            {   
                Double_t mean = inputHist.at(h)->GetMean(1);
                Double_t RMS = inputHist.at(h)->GetRMS();
                myfile << variable.at(h) << endl;
                myfile << "      " << mean << "        " << RMS << endl << endl;;
            }
      }
      else
          cout << "unable to open a file!" << endl;
    }
    
    can->SaveAs("./output/" + process + inputCanvasName.at(0) + ".eps");
    cout << "drawing" << endl;
    
    
    TString rootOut = "./output/" + process + inputCanvasName.at(0) + ".root";
    TFile fi(rootOut,"RECREATE");
    inputHist.at(0)->Write();
    fi.Close();
}
