#include <exception>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//vim uncertainties.txt
//./ZnunuUnctTablePerProcessTheoreticalUnc ttZtableUnc.tab ttZtableUncRelHigher.tab realregions.txt uncertainties.txt ttZgrouptableUncRelHigher.tab grouprealregions.txt groupuncertainties.txt ttZ

using namespace std;
using namespace theDoctor;


int main(int argc, char *argv[]){

        if(argc != 9)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTabAbs = argv[1];
        string inputTabRel = argv[2];
        string inputRegions = argv[3];
        string uncNames = argv[4];
        string theoTabRel = argv[5];
        string theoRegions = argv[6];
        string theouncNames = argv[7];
        string currentProcess = argv[8];
        //specify dataset to be read
        string group = "groupTheoretical";

        vector<string> regions;
        vector<string> regionsTheo;

        cout << "in here 1" << endl;
        //read signal regions 
	string line;
        ifstream regfile(inputRegions);
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                if (line.find("_I_") == std::string::npos)
                    regions.push_back(line);
          
            }
            regfile.close();
        }

        cout << "in here 2" << endl;
        vector<string> regionsTheoA;
        vector<string> regionsTheoC;
        //read the other regions
	string line2;
        ifstream regfile2(theoRegions);
        if (regfile2.is_open())
        {
            while ( getline (regfile2,line2) )
            {
                cout << " theo region " << line2 << endl;
                if (line2.find("_I_") == std::string::npos)
                {
                cout << " theo region not i " << line2 << endl;
                    regionsTheo.push_back(line2);
                    if (line2.find("_AB_") != std::string::npos)
                        regionsTheoA.push_back(line2);
                    if (line2.find("_CDEFGH_") != std::string::npos)
                        regionsTheoC.push_back(line2);
                }
            }
            regfile2.close();
        }
     

        cout << "in here 3" << endl;
        Table tabAbs(inputTabAbs);
        cout << "in here 3.1" << endl;
        Table tabRel(inputTabRel);
        cout << "in here 3.2" << endl;
        Table tabTheo(theoTabRel);

        cout << "in here 3.5" << endl;
        //read names of sysematics 
        vector<string> colI;
	string line3;
        ifstream systfile(theouncNames);
        if (systfile.is_open())
        {
            colI.push_back("yield");
       
            while ( getline (systfile,line3) )
            {
                colI.push_back(line3);
                cout << "line " << line3 << endl;
            }
            systfile.close();
        }
      
        cout << "in here 4" << endl;
       vector<string> systOutNames;
       for(uint32_t s=0; s<colI.size(); s++)
       {
               if(s==0 || s%2 == 1)
                   systOutNames.push_back(colI.at(s));
       }
 
 
        cout << "in here 5" << endl;
        for(uint32_t r = 0; r<regions.size(); r++)
        {

            if (regions.at(r).find("_A_") != std::string::npos || regions.at(r).find("_B_") != std::string::npos)
            {
        cout << "in here 6" << endl;
                string METBin = regions.at(r);
                METBin.erase(0, 7);
                cout << "MET bin AB " <<METBin << endl;
                for(uint32_t ab=0; ab<regionsTheoA.size(); ab++ )
                {
        cout << "in here 7" << endl;
                    if(regionsTheoA.at(ab).find(METBin) != std::string::npos )
                    {
                        cout << "MET bin AB " << METBin << " ,real region " << regions.at(r) << " ,teho region " << regionsTheoA.at(ab) << endl;
                        Figure yieldAB(0,0);
                        for(uint32_t c= 0; c<colI.size(); c++)
                        {
                            if(c==0)
                            {
                                yieldAB = tabAbs.Get(colI.at(c),regions.at(r) );
                            }
                            if(c%2 != 0)
                            {
                                Figure relTheoreticalInPer = tabTheo.Get(colI.at(c), regionsTheoA.at(ab) );
                                Figure relTheoretical = Figure( relTheoreticalInPer.value()/100, 0);
                                if(yieldAB.value()!= 0)
                                {
                                    tabAbs.Set(colI.at(c),regions.at(r), yieldAB-(yieldAB*relTheoretical) );
                                    tabAbs.Set(colI.at(c+1),regions.at(r), yieldAB+(yieldAB*relTheoretical) );
                                    tabRel.Set(colI.at(c),regions.at(r),relTheoreticalInPer);
                                }
                                else
                                {
                                    tabAbs.Set(colI.at(c),regions.at(r), Figure(0,0) );
                                    tabAbs.Set(colI.at(c+1),regions.at(r), Figure(0,0) );
                                    tabRel.Set(colI.at(c),regions.at(r),Figure(0,0));
 
                                }
                                 

                            }
                        }
                        
                    }

                }
            }
            else if (regions.at(r).find("_C_") != std::string::npos || regions.at(r).find("_D_") != std::string::npos || regions.at(r).find("_E_") != std::string::npos || regions.at(r).find("_F_") != std::string::npos ||regions.at(r).find("_G_") != std::string::npos || regions.at(r).find("_H_") != std::string::npos)
            {
        cout << "in here 8" << endl;
                string METBin = regions.at(r);
                METBin.erase(0, 7);
                cout << "MET bin CDEFGH " <<METBin << endl;
                for(uint32_t cd=0; cd<regionsTheoC.size(); cd++ )
                {
        cout << "in here 9" << endl;
                    if(regionsTheoC.at(cd).find(METBin) != std::string::npos )
                    {
                        cout << "MET bin CDEFGH " << METBin << " ,real region " << regions.at(r) << " ,teho region " << regionsTheoC.at(cd) << endl;
                        Figure yieldCD(0,0);
                        for(uint32_t c= 0; c<colI.size(); c++)
                        {
                            if(c==0)
                            {
                                yieldCD = tabAbs.Get(colI.at(c),regions.at(r) );
                            }
                            if(c%2 != 0)
                            {
                                Figure relTheoreticalInPer = tabTheo.Get(colI.at(c), regionsTheoC.at(cd) );
                                Figure relTheoretical = Figure( relTheoreticalInPer.value()/100, 0);
                                if(yieldCD.value() !=0)
                                {
                                    tabAbs.Set(colI.at(c),regions.at(r), yieldCD-(yieldCD*relTheoretical) );
                                    tabAbs.Set(colI.at(c+1),regions.at(r), yieldCD+(yieldCD*relTheoretical) );
                                    tabRel.Set(colI.at(c),regions.at(r),relTheoreticalInPer);
                                }
                                else
                                {
                                    tabAbs.Set(colI.at(c),regions.at(r), Figure(0,0) );
                                    tabAbs.Set(colI.at(c+1),regions.at(r), Figure(0,0) );
                                    tabRel.Set(colI.at(c),regions.at(r),Figure(0,0));
 
                                }
                                 

                            }
                        }
                    }

                }
            }
            else
            {
               cout << "The issue, the region was not found" << endl;
            }
         

        cout << "in here 10" << endl;
 
            /*for(uint32_t c = 0; c< colI.size(); c++)
            {
                Figure f1 = tab1.Get(colI.at(c), regions.at(r));
                Figure f2 = tab2.Get(colI.at(c), regions.at(r));
                Figure f3 = tab3.Get(colI.at(c), regions.at(r));

                Figure f = f1+f2+f3;
                tabSum.Set(colI.at(c), regions.at(r),f );

                if(c==0)
                {
                    y1 = f1;
                    y2 = f2;
                    y3 = f3;
                    unc1 = 0;
                    unc2 = 0;
                    unc3 = 0;
                    unc1Stat = 0;
                    unc2Stat = 0;
                    unc3Stat = 0;
                    unc1Syst = 0;
                    unc2Syst = 0;
                    unc3Syst = 0;

                    unc1 += y1.error()*y1.error();
                    unc2 += y2.error()*y2.error();
                    unc3 += y3.error()*y3.error();
                    unc1Stat = y1.error()*y1.error();
                    unc2Stat = y2.error()*y2.error();
                    unc3Stat = y3.error()*y3.error();

                }
                else if(c%2 != 0)
                {
                    bef1 = f1.value();
                    bef2 = f2.value();
                    bef3 = f3.value();
                }
                else if(c%2 ==0)
                {
                   uncHigher1 = abs(bef1-y1.value()) > abs(f1.value() - y1.value()) ? abs(bef1-y1.value()): abs(f1.value() - y1.value());
                   uncHigher2 = abs(bef2-y2.value()) > abs(f2.value() - y2.value()) ? abs(bef2-y2.value()): abs(f2.value() - y2.value());
                   uncHigher3 = abs(bef3-y3.value()) > abs(f3.value() - y3.value()) ? abs(bef3-y3.value()): abs(f3.value() - y3.value());

                   if(colI.at(c) == "lepSFDN" || colI.at(c) == "lepSFUP"|| colI.at(c) == "topPtModeling" || colI.at(c) == "topPtmodeling2")
                   {
                       cout << "in top pt 1" << endl;
                       unc1 += uncHigher1*uncHigher1;
                       unc2 += uncHigher2*uncHigher2;
                       unc3 += uncHigher3*uncHigher3;
                       unc1Syst += uncHigher1*uncHigher1;
                       unc2Syst += uncHigher2*uncHigher2;
                       unc3Syst += uncHigher3*uncHigher3;
                   }
                   else
                   {
                       unc1 += (abs(bef1-f1.value())/2) * (abs(bef1-f1.value())/2) ;
                       unc2 += (abs(bef2-f2.value())/2) * (abs(bef2-f2.value())/2) ;
                       unc3 += (abs(bef3-f3.value())/2) * (abs(bef3-f3.value())/2) ;
                       unc1Syst += (abs(bef1-f1.value())/2) * (abs(bef1-f1.value())/2) ;
                       unc2Syst += (abs(bef2-f2.value())/2) * (abs(bef2-f2.value())/2) ;
                       unc3Syst += (abs(bef3-f3.value())/2) * (abs(bef3-f3.value())/2) ;

                   }

                }
                else
                {
                    cout << "something went terribly wrong" << endl;
                }

                if(c == colI.size() -1)
                {
                    unc1 += (lumi*y1.value())*(lumi*y1.value());
                    unc2 += (lumi*y2.value())*(lumi*y2.value());
                    unc3 += (lumi*y3.value())*(lumi*y3.value());
                    unc1Syst += (lumi*y1.value())*(lumi*y1.value());
                    unc2Syst += (lumi*y2.value())*(lumi*y2.value());
                    unc3Syst += (lumi*y3.value())*(lumi*y3.value());

                    totYields.Set(processes.at(0), regions.at(r), Figure(y1.value(),sqrt(unc1)));
                    totYields.Set(processes.at(1), regions.at(r), Figure(y2.value(),sqrt(unc2)));
                    totYields.Set(processes.at(2), regions.at(r), Figure(y3.value(),sqrt(unc3)));
                    cout << "stat " << sqrt(unc1Stat) << " syst " << sqrt(unc1Syst) << endl;
                    totYieldsStatSyst.Set(processes.at(0), regions.at(r), Figure(y1.value(),sqrt(unc1Stat),sqrt(unc1Syst)));
                    totYieldsStatSyst.Set(processes.at(1), regions.at(r), Figure(y2.value(),sqrt(unc2Stat),sqrt(unc2Syst)));
                    totYieldsStatSyst.Set(processes.at(2), regions.at(r), Figure(y3.value(),sqrt(unc3Stat),sqrt(unc3Syst)));
                }
            }
        }

        float uncBef = 0;
        float uncHigher = 0;
        for(uint32_t r=0; r<regions.size();r++)
        {
            Figure yieldVal;
            float uncTot = 0;
            float uncTotStat = 0;
            float uncTotSyst = 0;
            for(uint32_t c = 0; c< colI.size(); c++)
            {
                Figure resall = tabSum.Get(colI.at(c), regions.at(r));
                
                if(c==0)
                {
                    yieldVal = resall;
                    uncTot += yieldVal.error()*yieldVal.error() ;
                    uncTotStat = yieldVal.error()*yieldVal.error() ;
                    if(resall.value() !=0)
                        thigher.Set(systOutNames.at(0), regions.at(r), Figure((resall.error()/resall.value())*100, 0));
                    else
                        thigher.Set(systOutNames.at(0), regions.at(r), Figure(-1, 0));

                }
                else if(c%2 != 0)
                {
                    uncBef = resall.value();
                }
                else if(c%2 == 0)
                {
                   uncHigher = abs(uncBef-yieldVal.value()) > abs(resall.value() - yieldVal.value()) ? abs(uncBef-yieldVal.value()): abs(resall.value() - yieldVal.value());


                   if(colI.at(c) == "lepSFDN" || colI.at(c) == "lepSFUP"|| colI.at(c) == "topPtModeling" || colI.at(c) == "topPtmodeling2")
                   {
                       cout << "in top pt 2" << colI.at(c) << endl;
                       uncTot += uncHigher*uncHigher ;
                       uncTotSyst += uncHigher*uncHigher ;
                       if(yieldVal.value() != 0)
                           thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((uncHigher/yieldVal.value())*100, 0));
                       else
                           thigher.Set(systOutNames.at(c/2), regions.at(r), Figure(0, 0));

                   }
                   else
                   {

                       float center;
                       if(uncBef+resall.value() != 0)
                           center = (abs(uncBef-resall.value()))/(uncBef+resall.value());
                       else
                           center = 0;
                       
                       uncTot += (abs(uncBef-resall.value())/2) * (abs(uncBef-resall.value())/2) ;
                       uncTotSyst += (abs(uncBef-resall.value())/2) * (abs(uncBef-resall.value())/2) ;

                       thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((center)*100, 0));

                   }

                   if(c==colI.size() -1)
                       cout << "yield" << yieldVal.value() << "uncertainty is" << sqrt(uncTot) << endl;

                }
                else
                {
                    cout << "there is a mistake somewhere " << endl;
                }
                if(c == colI.size() -1)
                {
                    uncTot += (lumi*yieldVal.value())*(lumi*yieldVal.value());
                    uncTotSyst += (lumi*yieldVal.value())*(lumi*yieldVal.value());
                    totYields.Set(processes.at(3), regions.at(r), Figure(yieldVal.value(),sqrt(uncTot)));
                    totYieldsStatSyst.Set(processes.at(3), regions.at(r), Figure(yieldVal.value(),sqrt(uncTotStat),sqrt(uncTotSyst)));
                }

            }*/

        }

       //summary plots
       /* for(uint32_t s=0; s<systOutNames.size(); s++)
        {
            vector<float> oneSysts;
            for(uint32_t x=0; x<regions.size(); x++)
            {
                Figure oneSyst = thigher.Get(systOutNames.at(s), regions.at(x));
                oneSysts.push_back(oneSyst.value());
            }
            sort(oneSysts.begin(), oneSysts.end());
            for(uint32_t z=0; z<oneSysts.size(); z++)
            {
                if(oneSysts.at(z) != 0 )
                {
                    tsummary.Set(uncRange.at(0), systOutNamesAll.at(s), oneSysts.at(z));
                    break;
                }
            }
            tsummary.Set(uncRange.at(1), systOutNamesAll.at(s), oneSysts.at(regions.size()-1));
        }
        tsummary.Set(uncRange.at(1), systOutNamesAll.at(systOutNames.size()), Figure(6.2,0));
*/
        //@MJ@ TODO add lumi, CS

        tabAbs.Print(static_cast<string>("recomputed"+currentProcess+"tableUnc.tab"),4);
        tabAbs.PrintLatex(static_cast<string>("recomputed"+currentProcess+"tableUnc.tex"),4);
        tabRel.Print(static_cast<string>("recomputed"+currentProcess+"tableUncRelHigher.tab"),4, "noError");
        tabRel.PrintLatex(static_cast<string>("recomputed"+currentProcess+"tableUncRelHigher.tex"),4, "noError");
 
/*        tabSum.Print(static_cast<string>(group+"summedprocessesZnunu.tab"),4);
        tabSum.PrintLatex(static_cast<string>(group+"summedprocessesZnunu.tex"),4); 
        thigher.Print(static_cast<string>(group+"summedprocessesRelHigherZnunu.tab"),4, "noError");
        thigher.PrintLatex(static_cast<string>(group+"summedprocessesRelHigherZnunu.tex"),4, "noError"); 
        tsummary.Print(static_cast<string>(group+"summedprocessesRelHigherZnunuSummary.tab"),4, "noError");
        tsummary.PrintLatex(static_cast<string>(group+"summedprocessesRelHigherZnunuSummary.tex"),4, "noError"); 
        totYields.Print(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrors.tab"),4);
        totYields.PrintLatex(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrors.tex"),4); 
        totYieldsStatSyst.Print(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrorsStatSystBreakdown.tab"),4, "systError");
        totYieldsStatSyst.PrintLatex(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrorsStatSystBreakdown.tex"),4, "systError"); 
*/

}
