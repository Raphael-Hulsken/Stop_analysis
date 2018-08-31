//use
// root -l -q AddDivideMultiply.C++ xx.steer nameOfProcess

#include <fstream>
#include <iostream>
#include <vector>

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
#include "TMath.h"
#include "TLegend.h"

#include "FakeDataMakerForNJets.h"

#include <../../sonicScrewdriver/interface/Table.h>
#include <../../sonicScrewdriver/interface/Figure.h>

using namespace std;
//using namespace theDoctor;

int main(int argc, char *argv[])
{
  
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();

    if(argc != 2)
    {
        cout << "incorrect number of coommand line parametres!" << endl;
        return -1;
    }

    //steering file with information
    TString input = argv[1];
    //information if addition, division or multiplication required
    //TString process = argv[2];
    
    //@MJ@ TODO DataReader::readData(input)
    //@MJ@ TODO DataReader::fillInfo()


    vector<TFile*> inputTFiles;
    vector<TH1D*> inputHist;
    vector<TCanvas*> inputCanvas;

    //@MJ@ TODO use object before all members of Data Reader

    for(uint32_t files = 0; files < inputRootFile.size(); files++)
    {
        cout << "input file: " <<  inputRootFile.at(files) << endl;
 
        inputTFiles.push_back(new TFile(inputRootFile.at(files)));
        
        cout << "variable: " <<  variable.at(files) << endl;
        cout << "input canvas: " <<  inputCanvasName.at(files) << endl;
        cout << "dir " << inputDir.at(files) << endl;
 
        inputCanvas.push_back(dynamic_cast<TCanvas*>(inputTFiles.at(files)->Get(inputDir.at(files) + inputCanvasName.at(files))));
        inputHist.push_back(dynamic_cast<TH1D*>(inputCanvas.at(files)->GetPrimitive(variable.at(files))->Clone()));
    }
    
    TString rootOut = "./outputiBkgEst/" + process + inputCanvasName.at(0) + ".root";
    TFile fi(rootOut,"RECREATE");
   
    TLegend *leg = new TLegend(0.7,0.5,0.9,0.9);
    TCanvas* can = new TCanvas("can", "can");
    can->Divide(1,1);
    can->cd(1);

    //first rescale the NJet (fake data) -> use method...
    //here the MET should be binned
    //steer fake data, MC
    //N Jet make factors 3/2 a 4/2 + n bTAGGS -> use method.... for fake data and mc
    //compute scale factors mc/data
    //rescale mc with this hist->scale
    // use weighted mc and data
    // compute the prediction from reweighted mc and data for each bin -> use method....
    //create table with met sr from fake data and prediction
   
	    //vector<>
 
    can->SaveAs("./outputBkgEst/" + process + inputCanvasName.at(0) + ".eps");
    
    fi.Close();

    //normalize and add
    if(process == "naddition" && inputRootFile.size() > 1)
    {   
        //normalize the first one
        Double_t norm1 = inputHist.at(0)->Integral(); 
        inputHist.at(0)->Scale(1/norm1);
        for(uint32_t h = 0; h < inputHist.size(); h++)
        {   
            //normalize
            Double_t norm = inputHist.at(h+1)->Integral(); 
            inputHist.at(h+1)->Scale(1/norm);
            //add
            inputHist.at(0)->Add(inputHist.at(1+h), 1);
        }

        //inputHist.at(0)->GetYaxis()->SetTitle("whatever");
        inputHist.at(0)->Draw();
        inputHist.at(0)->Write();
    }

    //merge two histograms
    if(process == "merge")
    {  
        cout << "process merge" << endl;
        
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {   
            Double_t norm = inputHist.at(h)->Integral(); 
            inputHist.at(h)->Scale(1/norm);
            inputHist.at(h)->SetLineColor(h+1);
            inputHist.at(h)->SetLineStyle(1);
            if(inputCanvasName.at(0) == "ak8recoWPt")
                inputHist.at(h)->GetXaxis()->SetTitle("Pt of ak8 jet");
            else if(inputCanvasName.at(0) == "NSJreals")
                inputHist.at(h)->GetXaxis()->SetTitle("#tau_{2} / #tau_{1}");
            else
                inputHist.at(h)->GetXaxis()->SetTitle("Raw mass");

            inputHist.at(h)->GetYaxis()->SetTitle("Entries");
            inputHist.at(h)->SetTitle("");
            if(h == 0)
            {
                Double_t max = inputHist.at(h)->GetMaximum();
                inputHist.at(h)->SetMaximum(2*max);
                inputHist.at(h)->Draw("hist");
            }
            else
            {
                inputHist.at(h)->Draw("same hist");
            }

           
           inputHist.at(h)->Write();
        }
    }
    
    //get nr of entries from 1-bin distribution
    if(process == "entries")
    { 
        cout << "entries" << endl;
        ofstream myfile;
        myfile.open (input +".txt");
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {
            myfile << variable.at(h) << " nr of entries: " << inputHist.at(h)->GetBinContent(1) << endl;
        } 
        myfile.close(); 
    }

    //get integral
    //@MJ@ TODO not general anymore!
    if(process == "integral")
    { 
        cout << "integral" << endl;
        ofstream intFile;
        intFile.open (input +"Integral.txt");
        vector<float> data;
        vector<string> names;
        
        vector<string> colI = {"baseline", "all cuts"};
        vector<string> rowI = {"s-channel", "t-channel", "tW", "0 gen lep", "1 gen lep", "2 gen lep", "tW 0 gen lep", "tW 1 gen lep", "tW 2 gen lep", "tt 1l", "tt 2l"};
        uint32_t rowIdI = 0;
        vector<Double_t> error;
        Double_t result;
 
        Table tI(colI, rowI);
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {
            error.push_back(0);
            if(h == 0 || h%2 == 0) //@MJ@ TODO be aware of this, not really safe!
            { 
                rowIdI = h/2;
                cout << "just to be sure, the var should be in baseline region: "<< variable.at(h) << endl;
                result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
                cout << "error " << error.at(h);
                tI.Set(colI.at(0), rowI.at(rowIdI), Figure(result, error.at(h)));
            }
            else if(h == 1 || h%2 == 1)
            {
                result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
                tI.Set(colI.at(1), rowI.at(rowIdI), Figure(result, error.at(h)));
            }
            intFile << variable.at(h) << " nr of entries: " << inputHist.at(h)->Integral() << endl;
            data.push_back(inputHist.at(h)->Integral());
            names.push_back(static_cast<string>(inputLegend.at(h)));
        } 
        intFile.close();
        tab.printTable(data, names, static_cast<string>(input +"Integral.tex"));
        tI.PrintLatex(static_cast<string>(input +"table.tex"));
    }
    //get integral
    //@MJ@ TODO not general anymore!
    if(process == "integralv2")
    { 
        cout << "integral" << endl;
        ofstream intFile;
        intFile.open (input +"Integral.txt");
        vector<float> data;
        vector<string> names;
        
        vector<string> colI = {"1 jets", "2 jets", "1 or 2 jets"};
        vector<string> rowI = {"tW 2 gen lep", "tt 2l"};
        uint32_t rowIdI = 0;
        vector<Double_t> error;
        Double_t result;
 
        Table tI(colI, rowI);
        for(uint32_t h = 0; h <inputRootFile.size(); h++ ) //@MJ@ TODO
        {
            error.push_back(0);
            if(h == 0 || h%3 == 0) //@MJ@ TODO be aware of this, not really safe!
            { 
                rowIdI = h/3;
                cout << "just to be sure, the var should be in baseline region: "<< variable.at(h) << endl;
                result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
                cout << "error " << error.at(h);
                tI.Set(colI.at(0), rowI.at(rowIdI), Figure(result, error.at(h)));
            }
            else if(h == 1 || h%3 == 1)
            {
                result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
                tI.Set(colI.at(1), rowI.at(rowIdI), Figure(result, error.at(h)));
            }
            else if(h == 2 || h%3 == 2)
            {
                result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
                tI.Set(colI.at(2), rowI.at(rowIdI), Figure(result, error.at(h)));
            }
            intFile << variable.at(h) << " nr of entries: " << inputHist.at(h)->Integral() << endl;
            data.push_back(inputHist.at(h)->Integral());
            names.push_back(static_cast<string>(inputLegend.at(h)));
        } 
        intFile.close();
        tab.printTable(data, names, static_cast<string>(input +"Integralv2.tex"));
        tI.PrintLatex(static_cast<string>(input +"tablev2.tex"));
    }

    //get integral and efficiency
    //@MJ@ TODO not general anymore!
    if(process == "intandeff")
    { 
        cout << "intandeff" << endl;
        ofstream intFile;
        intFile.open (input +"eff.txt");
        vector<float> data;
        vector<string> names;
        vector<string> col = {"MT cut", "MET cut", "MT2W cut", "all cuts"};
        vector<string> row = {"s-channel", "t-channel", "tW", "0 gen lep", "1 gen lep", "2 gen lep", "tW 0 gen lep", "tW 1 gen lep", "tW 2 gen lep", "tt 1l", "tt 2l"};
        
        Table t(col, row);
        Double_t def = 0;
        uint32_t rowId = 0;
        vector<Double_t> error;
        Double_t dError;
        Double_t result;
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {
            error.push_back(0);
            if(h == 0 || h%5 == 0) //@MJ@ TODO be aware of this, not really safe!
            { 
                def = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(),dError);
                rowId = h/5;
                cout << "just to be sure, the var should be in baseline region: "<< variable.at(h) << endl;
            }
            else if(h == 1 || h%5 == 1)
            {
               result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
               error.at(h) = sqrt(((error.at(h)/result)*(error.at(h)/result)) + ((dError/def)*(dError/def)))*(result/def);
               t.Set(col.at(0), row.at(rowId), Figure((result/def)*100, error.at(h)*100));
            }
            else if(h == 2 || h%5 == 2)
            {
               result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
               error.at(h) = sqrt(((error.at(h)/result)*(error.at(h)/result)) + ((dError/def)*(dError/def)))*(result/def);
               t.Set(col.at(1), row.at(rowId), Figure((result/def)*100, error.at(h)*100));
            }
            else if(h == 3 || h%5 == 3)
            {
               result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
               error.at(h) = sqrt(((error.at(h)/result)*(error.at(h)/result)) + ((dError/def)*(dError/def)))*(result/def);
               t.Set(col.at(2), row.at(rowId), Figure((result/def)*100, error.at(h)*100));
            }
            else if(h == 4 || h%5 == 4)
            {
               result = inputHist.at(h)->IntegralAndError(1,inputHist.at(h)->GetNbinsX(), error.at(h) );
               error.at(h) = sqrt(((error.at(h)/result)*(error.at(h)/result)) + ((dError/def)*(dError/def)))*(result/def);
               t.Set(col.at(3), row.at(rowId), Figure((result/def)*100, error.at(h)*100));
            }
            intFile << variable.at(h) << " nr of entries: " << (inputHist.at(h)->Integral())/def << endl;
            data.push_back((inputHist.at(h)->Integral())/def);
            names.push_back(static_cast<string>(inputLegend.at(h)));
        }
        intFile.close();
        tab.printTable(data, names, static_cast<string>(input +"Integral.tex"));
        t.PrintLatex(static_cast<string>(input +"table.tex"), 6, "");
    }

    //normalize and plot distribution
    if(process == "normalize")
    { 
        cout << "normalize" << endl;
 
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {
            Double_t norma =  inputHist.at(h)->Integral();
            inputHist.at(h)->Scale(1/norma);
            inputHist.at(h)->SetLineColor(1+h);
            inputHist.at(h)->GetXaxis()->SetTitle(inputCanvasName.at(h));
            //inputHist.at(h)->SetMaximum(0.05);
            inputHist.at(h)->SetLineWidth(3);
            leg->AddEntry(inputHist.at(h), inputLegend.at(h));
            if(h == 0)
            {
                Double_t max = inputHist.at(h)->GetMaximum();
                inputHist.at(h)->SetMaximum((1.5)*(max));
                inputHist.at(h)->Draw("hist e");
            }
            else
                inputHist.at(h)->Draw("hist same e");
        
            inputHist.at(h)->Write();
        }
        leg->Draw();
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
        inputHist.at(1)->Scale(1/n2);
        inputHist.at(1)->Draw("same hist");

        Double_t n1 = inputHist.at(0)->Integral(); 
        cout << "gen reals: "  << n1 << endl;
        inputHist.at(0)->SetLineColor(kBlue);
        inputHist.at(0)->Scale(1/n1);
        inputHist.at(0)->Draw("hist same");
        
        division->Write();
        inputHist.at(1)->Write();
        inputHist.at(0)->Write();

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
            
            output->SetTitle("");
            output->GetXaxis()->SetTitle("efficiency");
            output->GetYaxis()->SetTitle("efficiency on fake");
            output->Draw("A*"); 
            output->Write();
        }
    }
    
    //efficiency for ROC curves
    if(process == "eff4ROC")
    {
            cout << "process bsefficiency" << endl;
           
            Int_t nrOfBins = inputHist.at(0)->GetNbinsX();
            TGraph* output = new TGraph(nrOfBins);
            Double_t norm = inputHist.at(0)->Integral();
            
            //initial number of events
            Double_t norm1 = inputHist.at(0)->Integral();

            for(uint32_t bin = 1; bin<(nrOfBins+1); bin++)
            {
                output->SetPoint(bin, inputHist.at(0)->GetXaxis()->GetBinLowEdge(bin), (inputHist.at(0)->Integral(bin, nrOfBins))/norm);
            }
            
            output->GetXaxis()->SetTitle("cut efficiency");
            output->GetYaxis()->SetTitle(inputCanvasName.at(0));
            output->Draw("A*"); 
    }
    
    //signal/background ratio
    if(process == "sbratio")
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
                output->SetPoint(bin, inputHist.at(0)->GetXaxis()->GetBinLowEdge(bin), ((inputHist.at(0)->Integral(bin, nrOfBins))/norm1)/(TMath::Sqrt((inputHist.at(1)->Integral(bin, nrOfBins))/norm2)));
            }
            
            output->GetXaxis()->SetTitle(inputCanvasName.at(0));
            output->GetYaxis()->SetTitle(inputCanvasName.at(1));
            output->Draw("A*"); 
        }
    }
    
    //count mean and rms
    if(process == "countmean")
    {
        cout << "process countmean" << endl;
        ofstream myfile ("./output/TableOfMeanValues.txt");
        if (myfile.is_open())
        { 
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

    return 0;
}
