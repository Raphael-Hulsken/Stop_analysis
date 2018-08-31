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

#include "../../sonicScrewdriver/interface/Table.h"
#include "../../sonicScrewdriver/interface/Figure.h"
#include "latexTable.h"

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[])
{
  
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();

    if(argc != 3)
    {
        cout << "incorrect number of coommand line parametres!" << endl;
        return -1;
    }

    //steering file with information
    //TString input = gApplication->Argv(5);
    TString input = argv[1];
    //information if addition, division or multiplication required
    //TString process = gApplication->Argv(6);
    TString process = argv[2];
    
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
    //TString dir = in_stream.at(1+3*nrOfVars);    


    //this will be read from steering
    vector<TString> inputRootFile;
    //this will be read from steering
    vector<TString> variable;
    //this will be read from steering
    vector<TString> inputCanvasName;
    //this will be read from steering
    vector<TString> inputDir;
    //this will be read from steering
    vector<TString> inputLegend;


    latexTable tab;

    for(uint32_t f = 0; f < nrOfVars; f++)
    {
        variable.push_back(in_stream.at(1+f));
        inputRootFile.push_back(in_stream.at(1+nrOfVars+f));
        inputCanvasName.push_back(in_stream.at(1+2*nrOfVars+f));
        inputDir.push_back(in_stream.at(1+3*nrOfVars+f)); //@MJ@ TODO all steerings have to be broadened!
        inputLegend.push_back(in_stream.at(1+4*nrOfVars+f)); //@MJ@ TODO all steerings have to be broadened!
    }

    vector<TFile*> inputTFiles;
    vector<TH1D*> inputHist;
    vector<TCanvas*> inputCanvas;
    TH1::SetDefaultSumw2();

    for(uint32_t files = 0; files < inputRootFile.size(); files++)
    {
        cout << "input file: " <<  inputRootFile.at(files) << endl;
 
        inputTFiles.push_back(new TFile(inputRootFile.at(files)));
        
        cout << "variable: " <<  variable.at(files) << endl;
        
        //inputTFiles.at(files)->cd("muon/preselection/");
       
        cout << "input canvas: " <<  inputCanvasName.at(files) << endl;
        cout << "dir " << inputDir.at(files) << endl;
 
        inputCanvas.push_back(dynamic_cast<TCanvas*>(inputTFiles.at(files)->Get(inputDir.at(files) + inputCanvasName.at(files))));

        inputHist.push_back(dynamic_cast<TH1D*>(inputCanvas.at(files)->GetPrimitive(variable.at(files))->Clone()));
    
    }

    cout << "nr of input files; " << inputRootFile.size() << endl;
    
    TString rootOut = "./output/" + process + inputCanvasName.at(0) +  inputLegend.at(0)+ ".root";
    TFile fi(rootOut,"RECREATE");
   
    TLegend *leg = new TLegend(0.7,0.5,0.9,0.9);
 
    TCanvas* can = new TCanvas("can", "can");
    can->Divide(1,1);
    can->cd(1);
    
    //normalize and add
    if(process == "naddition" && inputRootFile.size() > 1)
    {   
        //normalize the first one
        cout << "in n addition" << endl;
        Double_t norm1 = inputHist.at(0)->Integral(); 
        inputHist.at(0)->Scale(1/norm1);
        for(uint32_t h = 1; h < inputHist.size()+1; h++)
        {   
            //normalize
            Double_t norm = inputHist.at(h)->Integral(); 
            inputHist.at(h)->Scale(1/norm);
            //add
            inputHist.at(0)->Add(inputHist.at(h), 1);
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
                //Double_t max = inputHist.at(h)->GetMean(2);
                inputHist.at(h)->SetMaximum((1.6)*(max));
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
    
    //division of two histos
    if(process == "division")
    { 
        cout << "process division" << endl; 
        TH1D* division = dynamic_cast<TH1D*>(inputHist.at(0)->Clone());
        division->Divide(dynamic_cast<TH1D*>(inputHist.at(1)->Clone()));
        division->SetLineColor(kPink);
        division->SetMaximum(0.2);
        //division->SetMaximum(0.2);
        //leg->AddEntry(division, inputLegend.at(0));
      
        //TH1D* division1 = dynamic_cast<TH1D*>(inputHist.at(2)->Clone());
        //division1->Divide(dynamic_cast<TH1D*>(inputHist.at(3)->Clone()));
        //division1->SetLineColor(kPink);
        //leg->AddEntry(division1, inputLegend.at(2));

        //TH1D* division2 = dynamic_cast<TH1D*>(inputHist.at(4)->Clone());
        //division2->Divide(dynamic_cast<TH1D*>(inputHist.at(5)->Clone()));
        //division2->SetLineColor(kViolet);
        //leg->AddEntry(division2, inputLegend.at(4));
        
        Int_t nrOfBins = division->GetNbinsX();
        Double_t ex = 0;
        Double_t exx = 0;
        Double_t cont;
        for(uint32_t bin = 3; bin<(nrOfBins-4); bin++)
        {
            cont = division->GetBinContent(bin);
            ex += cont;
            exx += (cont*cont);
            cout << "ex " << ex << endl;
            cout << "exx " << exx << endl;
        }
        cout << "variance " << (exx/(8)) - ((ex/(8))*(ex/(8))) << endl;
        cout << "RMS: " << sqrt( (exx/(8)) - ((ex/(8))*(ex/(8)))) << endl;

        division->SetTitle("");
        division->GetXaxis()->SetTitle("nvertex");
        division->GetYaxis()->SetTitle("efsionMTmcprofile.epsficiency");
        //division->SetMaximum(2);
        
        division->Draw("hist e");
        //division1->Draw("same hist");
        //division2->Draw("same hist");
        division->Write();
        //division1->Write();
        //division2->Write();
        //leg->Draw();
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
            output->GetXaxis()->SetTitle(inputLegend.at(0));
            output->GetYaxis()->SetTitle(inputLegend.at(1));
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
    
    //signal/background ratio efficiency
    if(process == "sbratio")
    {
        if(inputRootFile.size() == 2)
        {
            cout << "process sbratio" << endl;
           
            Int_t nrOfBins = inputHist.at(0)->GetNbinsX();
            TGraph* output = new TGraph(nrOfBins);
            
            //initial number of events
            Double_t norm1 = 1;//inputHist.at(0)->Integral();
            Double_t norm2 = 1;//inputHist.at(1)->Integral();

            for(uint32_t bin = 1; bin<(nrOfBins+1); bin++)
            {
                output->SetPoint(bin, inputHist.at(0)->GetXaxis()->GetBinCenter(bin), ((inputHist.at(0)->Integral(bin, nrOfBins))/norm1)/(((inputHist.at(1)->Integral(bin, nrOfBins))/norm2))); //square root
            }
            
            output->GetXaxis()->SetTitle(inputCanvasName.at(0));
            output->GetYaxis()->SetTitle("signal/background efficiency");
            output->Draw("A*"); 
        }
    }
    
    //signal efficiency
    if(process == "sratio")
    {
        if(true)
        {
            cout << "process sratio" << endl;
           
            Int_t nrOfBins = inputHist.at(0)->GetNbinsX();
            TGraph* output = new TGraph(nrOfBins);
            
            //initial number of events
            Double_t norm1 = 1;//inputHist.at(0)->Integral();

            for(uint32_t bin = 1; bin<(nrOfBins+1); bin++)
            {
                output->SetPoint(bin, inputHist.at(0)->GetXaxis()->GetBinCenter(bin), ((inputHist.at(0)->Integral(bin, nrOfBins))/norm1));
            }
            
            output->GetXaxis()->SetTitle(inputCanvasName.at(0));
            output->GetYaxis()->SetTitle("signal cut efficiency");
            output->Draw("A*"); 
        }
    }
    
    //signal/background double ratio  efficiency
    if(process == "sbdoubleratio")
    {
        if(inputRootFile.size() == 4)
        {
            cout << "process bsefficiency" << endl;
           
            Int_t nrOfBins = inputHist.at(0)->GetNbinsX();
            TGraph* output = new TGraph(nrOfBins);
            
            //initial number of events
            Double_t norm1 = inputHist.at(0)->Integral();
            Double_t norm2 = inputHist.at(1)->Integral();
            Double_t norm3 = inputHist.at(2)->Integral();
            Double_t norm4 = inputHist.at(3)->Integral();

            for(uint32_t bin = 1; bin<(nrOfBins+1); bin++)
            {
                Double_t r1 = ((inputHist.at(0)->Integral(bin, nrOfBins))/norm1)/(((inputHist.at(1)->Integral(bin, nrOfBins))/norm2));
                Double_t r2 = ((inputHist.at(2)->Integral(bin, nrOfBins))/norm3)/(((inputHist.at(3)->Integral(bin, nrOfBins))/norm4));
                output->SetPoint(bin, bin, r1/r2); //square root
            }
            
            output->GetXaxis()->SetTitle("bin");
            output->GetYaxis()->SetTitle("signal/background efficiency double ratio: " + inputCanvasName.at(0) + "/" + inputCanvasName.at(3) );
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
                Double_t mean =inputHist.at(h)->GetMean(1);
                Double_t RMS = inputHist.at(h)->GetRMS();
                myfile << variable.at(h) << endl;
                myfile << "      " << mean << "        " << RMS << endl << endl;;
            }
      }
      else
          cout << "unable to open a file!" << endl;
    }
    
    //hist mean
    if(process == "meanhist")
    {
     
        cout << "mean hist" << endl;
        TH1F *output = new TH1F("output", "", inputHist.size(), 0, 5*inputHist.size());
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
            {   
                Double_t mean = inputHist.at(h)->GetMean(1);
                cout << "mean: " << mean << endl;
                Double_t err = inputHist.at(h)->GetMeanError(1);
                output->SetBinContent(h+1, mean);
                output->SetBinError(h+1, err);
                //Double_t RMS = inputHist.at(h)->GetRMS();
            }
            output->SetTitle("");
            output->SetLineColor(kPink);
            output->SetMinimum(50);
            output->GetXaxis()->SetTitle("nvertex");
            output->GetYaxis()->SetTitle("mean " + inputCanvasName.at(0));
            output->Draw("hist e"); 
            output->Write();
    }
    
    //hist RMS
    if(process == "rmshist")
    {
     
        cout << "rms hist" << endl;
        TH1F *output = new TH1F("output", "", inputHist.size(), 0, 5*inputHist.size());
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
            {   
                Double_t RMS = inputHist.at(h)->GetRMS();
                Double_t err = inputHist.at(h)->GetRMSError(1);
                output->SetBinContent(h+1, RMS);
                output->SetBinError(h+1, err);
            }
            output->SetTitle("");
            output->SetLineColor(kPink);
            output->SetMinimum(30);
            output->GetXaxis()->SetTitle("nvertex");
            output->GetYaxis()->SetTitle("RMS " + inputCanvasName.at(0));
            output->Draw("hist e"); 
            output->Write();
    }
    
    //hist tail
    if(process == "tailhist") //@MJ@ TODO make the binning same, in order to use one value for MT. MET and MT2W cut!!!!!!!!!!
    {
     
        cout << "tail hist" << endl;
        TH1F *output = new TH1F("output", "", inputHist.size(), 0, 5*inputHist.size());
        uint32_t n = 0;
        if(inputCanvasName.at(0) == "MET") n = 10;
        else if(inputCanvasName.at(0) == "MT2W") n = 8;
        else if(inputCanvasName.at(0) == "MT") n = 6;
        Double_t err1;
        Double_t err2;
        Double_t err;
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
            {   
                Int_t bins = inputHist.at(h)->GetNbinsX(); 
                Double_t den = inputHist.at(h)->IntegralAndError(1,bins, err1);
                Double_t num = inputHist.at(h)->IntegralAndError(n+1, bins, err2);
                Double_t eff = num/den;
                //err = sqrt((eff*(1-eff))/(den));
                err = sqrt((eff*(1-eff))/(den/0.181212)); //4 tt1l
                output->SetBinContent(h+1,eff);
                //output->SetBinError(h+1, sqrt((err1*err1) + (err2*err2)));
                output->SetBinError(h+1, err);
                
                //Double_t tail = inputHist.at(h)->IntegralAndError(n+1, bins, err1);
                //output->SetBinContent(h+1,tail);
                //output->SetBinError(h+1, err1);
            }
            output->SetTitle("");
            output->SetLineColor(kPink);
            output->SetMinimum(0);
            //output->SetMaximum(0.02);
            output->GetXaxis()->SetTitle("nvertex");
            output->GetYaxis()->SetTitle("fraction in tail " + inputCanvasName.at(0));
            output->Draw("hist e"); 
            output->Write();
    }
    
    //tail number
    if(process == "tailnr")
    {
     
        cout << "tail nr" << endl;
        TH1F *output = new TH1F("output", "", inputHist.size(), 0, 5*inputHist.size());
        uint32_t n = 0;
        if(inputCanvasName.at(0) == "MET") n = 10;
        else if(inputCanvasName.at(0) == "MT2W") n = 8;
        else if(inputCanvasName.at(0) == "MT") n = 6;
        Double_t err1;
        Double_t err2;
        Double_t err3;
        Double_t err;
        for(uint32_t h = 0; h <inputRootFile.size()-1; h++ )
            {   
                Int_t bins = inputHist.at(h)->GetNbinsX(); 
                Double_t den = inputHist.at(h)->IntegralAndError(1,bins, err1);
                if(den == 0) continue;
                cout << "den " << den << endl;
                Double_t num = inputHist.at(h)->IntegralAndError(n+1, bins, err2);
                cout << "num " << num << endl;
                Double_t norm = inputHist.at(7)->IntegralAndError(1, bins, err3);
                cout << "norm " << norm << endl;
                Double_t eff = (num/den)*norm;
                cout << "efficiency " << eff << endl;
                cout << "error " << (err2/num)*eff << endl;
                //err = sqrt((eff*(1-eff))/(den));
                //err = sqrt((eff*(1-eff))/(den/0.181212)); //4 tt1l
                output->SetBinContent(h+1,eff);
                output->SetBinError(h+1, (err2/num)*eff);
                //output->SetBinError(h+1, err);
                
                //Double_t tail = inputHist.at(h)->IntegralAndError(n+1, bins, err1);
                //output->SetBinContent(h+1,tail);
                //output->SetBinError(h+1, err1);
            }
            output->SetTitle("");
            output->SetLineColor(kPink);
            output->SetMinimum(0);
            Double_t max = output->GetBinContent(4);
            output->SetMaximum(2*max);
            cout << "uncertainty " << ((output->GetBinContent(6) - output->GetBinContent(2))/2) / output->GetBinContent(4) << endl;
            cout << "uncertainty 2 " << (( ((output->GetBinContent(6) + output->GetBinContent(5))/2) -((output->GetBinContent(2) + output->GetBinContent(3))/2) ) / (2*output->GetBinContent(4)));
            //output->SetMaximum(0.02);
            output->GetXaxis()->SetTitle("nvertex");
            output->GetYaxis()->SetTitle(" yield tail * N/N_{i} " + inputCanvasName.at(0));
            output->Draw("hist e"); 
            output->Write();
    }
    
    if(process == "mcpuhist")
    {
        TH1F *h1 = new TH1F("mcpu", "mcpu", 53, 0, 53);
        h1->SetBinContent(1,4.8551E-07 );
        h1->SetBinContent(2,1.74806E-06 );
        h1->SetBinContent(3,3.30868E-06 );
        h1->SetBinContent(4,1.62972E-05 );
        h1->SetBinContent(5,4.95667E-05 );
        h1->SetBinContent(6,0.000606966 );
        h1->SetBinContent(7,0.003307249 );
        h1->SetBinContent(8,0.010340741 );
        h1->SetBinContent(9,0.022852296 );
        h1->SetBinContent(10,0.041948781 );
        h1->SetBinContent(11,0.058609363 );
        h1->SetBinContent(12, 0.067475755);
        h1->SetBinContent(13, 0.072817826);
        h1->SetBinContent(14, 0.075931405);
        h1->SetBinContent(15, 0.076782504);
        h1->SetBinContent(16, 0.076202319);
        h1->SetBinContent(17, 0.074502547);
        h1->SetBinContent(18, 0.072355135);
        h1->SetBinContent(19, 0.069642102);
        h1->SetBinContent(20, 0.064920999);
        h1->SetBinContent(21, 0.05725576);
        h1->SetBinContent(22, 0.047289348);
        h1->SetBinContent(23, 0.036528446);
        h1->SetBinContent(24, 0.026376131);
        h1->SetBinContent(25, 0.017806872);
        h1->SetBinContent(26, 0.011249422);
        h1->SetBinContent(27, 0.006643385);
        h1->SetBinContent(28, 0.003662904);
        h1->SetBinContent(29, 0.001899681);
        h1->SetBinContent(30, 0.00095614);
        h1->SetBinContent(31, 0.00050028);
        h1->SetBinContent(32, 0.000297353);
        h1->SetBinContent(33, 0.000208717);
        h1->SetBinContent(34, 0.000165856);
        h1->SetBinContent(35, 0.000139974);
        h1->SetBinContent(36, 0.000120481);
        h1->SetBinContent(37, 0.000103826);
        h1->SetBinContent(38, 8.88868E-05);
        h1->SetBinContent(39, 7.53323E-05);
        h1->SetBinContent(40, 6.30863E-05);
        h1->SetBinContent(41, 5.21356E-05);
        h1->SetBinContent(42, 4.24754E-05);
        h1->SetBinContent(43, 3.40876E-05);
        h1->SetBinContent(44, 2.69282E-05);
        h1->SetBinContent(45, 2.09267E-05);
        h1->SetBinContent(46, 1.5989E-05);
        h1->SetBinContent(47, 4.8551E-06);
        h1->SetBinContent(48, 2.42755E-06);
        h1->SetBinContent(49, 4.8551E-07);
        h1->SetBinContent(50, 2.42755E-07);
        h1->SetBinContent(51, 1.21378E-07);
        h1->SetBinContent(52, 4.8551E-08);
        
        h1->Draw("hist e"); 
        h1->Write();
    }


    //normalize and add
    if(process == "ndivision")
    {  

         TFile *f = new TFile("MyDataPileupHistogram.root");
         f->ls();
         TH1D * h1 = (TH1D*)f->Get("pileup");
         
         TFile *d = new TFile("mcpuhistMTmcprofile.root");
         d->ls();
         TH1D * h2 = (TH1D*)d->Get("mcpu");
  
        Int_t nrOfBins1 = h1->GetNbinsX();    //normalize the first one
        Int_t nrOfBins2 = h2->GetNbinsX();    //normalize the first one
        cout << "nbins 1 " << nrOfBins1 << " nbins 2 " << nrOfBins2 << endl;
        Double_t norm1 = h1->Integral(); 
        h1->Scale(1/norm1);
        Double_t norm2 = h2->Integral(); 
        h2->Scale(1/norm2);
        cout << "norm 1 " << norm1 << "norm 2" << norm2 << endl;
        TH1D* division = h1;
        division->Divide(h2);
        TH1D* division2 = (TH1D*)division->Clone();


        TFile fi2(rootOut,"RECREATE");
        
        //inputHist.at(0)->GetYaxis()->SetTitle("whatever");
        division2->Draw();
        division2->Write();
        fi2.Close();
    }

    can->SaveAs("./output/" + process + inputCanvasName.at(0)+ inputLegend.at(0) + ".eps");
    cout << "drawing" << endl;
    
    
    fi.Close();

    return 0;
}
