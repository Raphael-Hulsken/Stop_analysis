#include <exception>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//vim uncertainties.txt
//./ZnunuUnctRootFileSummed groupsummedprocessesZnunu.tab groupsummedprocessesRelHigherZnunu.tab realregions.txt uncertainties.txt

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 5)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTab = argv[1];
        string inputTabRel = argv[2];
        string inputFile = argv[3];
        string uncNames = argv[4];

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
                cout << "row " << line << endl;
            }
            regfile.close();
        }
     
       cout << "in here 0.1 " << endl;
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
       cout << "in here 0.2 " << endl;

       TH1::SetDefaultSumw2();
       Table tab(inputTab);
       Table tabRel(inputTabRel);


       cout << "in here 0.3 " << endl;
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
                if(line2 == "topPtModeling") 
                    continue;
                colI.push_back(line2);
            }
            //colI.push_back("lumiDN");
            //colI.push_back("lumiUP");
            colI.push_back("total");
            systfile.close();
        }
       cout << "in here 1 " << endl;
        vector<TH1D*> histo;
        for(uint32_t h =0; h<colI.size(); h++)
        {
               histo.push_back(new TH1D(colI.at(h).c_str(), colI.at(h).c_str(), regions.size(), 0, regions.size()) );
        }

       cout << "in here 2 " << endl;
       
       //tables for all realtive and absolute uncertainties

       TFile fi2("ZNuNu_BkgEst.root","RECREATE");
       cout << "in here 3 " << endl;
       //loop over all regions
       for(uint32_t s=0; s<colI.size();s++)
       {
       cout << "in here 4 " << endl;
           for(uint32_t r=0; r< regions.size(); r++)
           {
               if(colI.at(s) == "mcStatsDN")
               {
       cout << "in here 5 " << endl;
                   theDoctor::Figure resall = tab.Get("yield", regions.at(r));
                   if(resall.value() == 0)
                       histo.at(s)->SetBinContent(r+1, 0 );
                   else
                       histo.at(s)->SetBinContent(r+1,resall.value()-resall.error() );     

                  cout << "MC stat unc " << (resall.error()/resall.value())*100 << endl; 
                   
               }
               else if(colI.at(s) == "mcStatsUP")
               {
       cout << "in here 6 " << endl;
                   theDoctor::Figure resall = tab.Get("yield", regions.at(r));
                   if(resall.value() == 0)
                       histo.at(s)->SetBinContent(r+1, 0 );
                   else
                       histo.at(s)->SetBinContent(r+1,resall.value()+resall.error() );      
               }
               /*else if(colI.at(s) == "lumiDN" || colI.at(s) == "lumiUP") //@MJ@ TODO redo
               {
       cout << "in here 7 " << endl;
                   theDoctor::Figure relUnc = tabRel.Get("luminosity", regions.at(r));
                   theDoctor::Figure resall = tab.Get("yield", regions.at(r));
                   if(resall.value() == 0)
                       histo.at(s)->SetBinContent(r+1, 0 );
                   else
                   {
                      Figure lumiUD(0,0);
                      if(colI.at(s) == "lumiDN")
                          lumiUD = resall  - (resall*(relUnc/100));
                      else
                          lumiUD = resall  + (resall*(relUnc/100));

                       histo.at(s)->SetBinContent(r+1, lumiUD.value() );
                       histo.at(s)->SetBinError(r+1, lumiUD.error() );
                   }      
                   
               }*/
               else if(colI.at(s) == "total")//@MJ@ TODO redo
               {
       cout << "in here 8 " << endl;
                   theDoctor::Figure relUnc = tabRel.Get("total", regions.at(r));
                   theDoctor::Figure resall = tab.Get("yield", regions.at(r));
                   if(resall.value() == 0)
                       histo.at(s)->SetBinContent(r+1, 0 );
                   else
                   {
                       histo.at(s)->SetBinContent(r+1, resall.value() );
                       histo.at(s)->SetBinError(r+1, resall.value()*(relUnc.value()/100) );
                   }      
                   
               }
               else
               {
       cout << "in here 9 " << endl;
                   theDoctor::Figure resall = tab.Get(colI.at(s), regions.at(r));
                   if(resall.value() == 0)
                   {
                       histo.at(s)->SetBinContent(r+1, 0 );
                       histo.at(s)->SetBinError(r+1, 0 );
                   }
                   else
                   {
                       histo.at(s)->SetBinContent(r+1,resall.value() );      
                       histo.at(s)->SetBinError(r+1,resall.error() );   
                   }   
               
               }
               histo.at(s)->GetXaxis()->SetBinLabel(r+1,regLabels.at(r).c_str());
            }
            histo.at(s)->Write();
         }
           

       fi2.Close();
}
