#include <exception>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//vim uncertainties.txt
//./ttZUnctRootFile tableUncZnunu.tab  realregions.txt uncertainties.txt

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 4)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTab = argv[1];
        string inputFile = argv[2];
        string uncNames = argv[3];

        vector<string> regions;
        //specify dataset to be read

        //read signal regions 
	string line;
        ifstream regfile(inputFile);
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                regions.push_back(line);
            }
            regfile.close();
        }
     
        vector<string> regLabels = regions;
        const std::string ss = "less";
        const std::string t = "<";
        
        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ss, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ss.size(), t );
                n += t.size();
            } 
        }              

       TH1::SetDefaultSumw2();
       Table tab(inputTab);


        //read names of sysematics 
        vector<string> colI;
	string line2;
        ifstream systfile(uncNames);
        if (systfile.is_open())
        {
            colI.push_back("yield");
            colI.push_back("mcStatsDN");
            colI.push_back("mcStatsUP");
            while ( getline (systfile,line2) )
            {
                if(line2 == "topPtmodeling2") 
                    continue;
                colI.push_back(line2);
            }
            systfile.close();
        }
       
        vector<TH1D*> histo;
        for(uint32_t h =0; h<colI.size(); h++)
        {
               histo.push_back(new TH1D(colI.at(h).c_str(), colI.at(h).c_str(), regions.size(), 0, regions.size()) );
        }

       
       //tables for all realtive and absolute uncertainties

       TFile fi2("bkgEst_Znunu.root","RECREATE");
       //loop over all regions
       for(uint32_t s=0; s<colI.size();s++)
       {
           for(uint32_t r=0; r< regions.size(); r++)
           {
               if(colI.at(s) == "mcStatsDN")
               {
                   theDoctor::Figure resall = tab.Get("yield", regions.at(r));
                   histo.at(s)->SetBinContent(r+1,resall.value()-resall.error() );      
                   
               }
               else if(colI.at(s) == "mcStatsUP")
               {
                   theDoctor::Figure resall = tab.Get("yield", regions.at(r));
                   histo.at(s)->SetBinContent(r+1,resall.value()+resall.error() );      
               }
               else
               {
                   theDoctor::Figure resall = tab.Get(colI.at(s), regions.at(r));
                   histo.at(s)->SetBinContent(r+1,resall.value() );      
                   histo.at(s)->SetBinError(r+1,resall.error() );      
               
               }
               histo.at(s)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
            }
            histo.at(s)->Write();
         }
           

       fi2.Close();
}
