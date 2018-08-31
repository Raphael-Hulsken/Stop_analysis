#include <exception>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//vim uncertainties.txt
//./ZnunuSumGroupedTables ttZgrouptableUnc.tab WZgrouptableUnc.tab ZZgrouptableUnc.tab grouprealregions.txt groupuncertainties.txt groupttZWZZZ
using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 7)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTab1 = argv[1];
        string inputTab2 = argv[2];
        string inputTab3 = argv[3];
        string inputFile = argv[4];
        string uncNames = argv[5];
        string group = argv[6];

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
                cout << "regions " << line << endl;
            }
            regfile.close();
        }
     

       Table tab1(inputTab1);
       Table tab2(inputTab2);
       //Table tab3(inputTab3);

        //read names of sysematics 
        vector<string> colI;
	string line2;
        ifstream systfile(uncNames);
        if (systfile.is_open())
        {
            colI.push_back("yield");
       
            while ( getline (systfile,line2) )
            {
                colI.push_back(line2);
            }
            systfile.close();
        }
      
       vector<string> systOutNames;
       for(uint32_t s=0; s<colI.size(); s++)
       {
               if(s==0 || s%2 == 1)
               {
                   systOutNames.push_back(colI.at(s));
                   cout << "syst out anmes " << colI.at(s) << endl;
               }
       }
 

        cout << "in here 1" << endl; 
        Table tabSum( colI, regions );
        cout << "in here 2" << endl; 
        Table thigher(systOutNames, regions);
        cout << "in here 3" << endl; 

        for(uint32_t r = 0; r<regions.size(); r++)
        {
            float bef = 0;
            for(uint32_t c = 0; c< colI.size(); c++)
            {
                Figure f1 = tab1.Get(colI.at(c), regions.at(r));
                Figure f2 = tab2.Get(colI.at(c), regions.at(r));
                //Figure f3 = tab3.Get(colI.at(c), regions.at(r));
                Figure f3 = Figure(0,0);

                Figure f = f1+f2+f3;
                tabSum.Set(colI.at(c), regions.at(r),f );

                if(c==0)
                {
                    thigher.Set(colI.at(0), regions.at(r), Figure(f.error()/f.value(),0));

                }
                else if(c%2 != 0)
                {
                    bef = f.value();
                }
                else if(c%2 ==0)
                {
                       float rel= 0 ;
                       if(bef+f.value() != 0)
                           rel = (abs(bef-f.value())) /(bef+f.value()) ;
                       thigher.Set(systOutNames.at(c/2), regions.at(r), Figure(rel*100,0));


                }
            }
        }

 
        tabSum.Print(static_cast<string>(group+"summedprocessesZnunu.tab"),4);
        tabSum.PrintLatex(static_cast<string>(group+"summedprocessesZnunu.tex"),4); 
        thigher.Print(static_cast<string>(group+"summedprocessesRelHigherZnunu.tab"),4);
        thigher.PrintLatex(static_cast<string>(group+"summedprocessesRelHigherZnunu.tex"),4, "noError"); 


}
