#include <exception>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//vim uncertainties.txt
//./ZnunuUnctSummedTables recomputedttZtableUnc.tab recomputedWZtableUnc.tab recomputedZZtableUnc.tab realregions.txt uncertainties.txt group
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
        //float lumi = 0.062;
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
            //colI.push_back("mcStatsDN");
            //colI.push_back("mcStatsUP");
            //colI.push_back("lumiDN");
            //colI.push_back("lumiUP");
            //if(currentProcess == "ZZ")
            //{
            //    colI.push_back("crossSectionDN");
            //    colI.push_back("crossSectionUP");
            //}
       
            while ( getline (systfile,line2) )
            {
                if(line2 == "topPtModeling" ) 
                    continue;
                if(line2 == "topPtmodeling2") 
                    continue;
                colI.push_back(line2);
            }
            systfile.close();
        }
      
       vector<string> systOutNames;
       for(uint32_t s=0; s<colI.size(); s++)
       {
               if(s==0 || s%2 == 1)
                   systOutNames.push_back(colI.at(s));
       }
        //systOutNames.push_back("luminosity");
        systOutNames.push_back("total");
 
        vector<string> systOutNamesAll = systOutNames;
        //systOutNamesAll.push_back("luminosity");
        //systOutNamesAll.push_back("pile-up");
        //systOutNamesAll.push_back("ZZ normalization");
 
        Table tabSum( colI, regions );
        Table thigher(systOutNames, regions);
        vector<string> uncRange = {"syst low (%)", "syst high (%)"};
        Table tsummary( uncRange, systOutNames);
        vector<string> processes = {"ttZ", "WZ", "Total"};
        Table totYields(processes, regions );
        Table totYieldsStatSyst(processes, regions );


	//recompute PU and set the normalization unc
	Figure PUDNtot(0.0);
	Figure PUUPtot(0.0);
	/*Figure PUyieldtot(0.0);
	for(uint32_t p=0; p<regions.size(); p++)
	{
	    PUDNtot += ( tab1.Get("PUdown", regions.at(p)) + tab2.Get("PUdown", regions.at(p))); //+ tab3.Get("PUdown", regions.at(p)));
	    PUUPtot += ( tab1.Get("PUup", regions.at(p)) +  tab2.Get("PUup", regions.at(p))); //+ tab3.Get("PUup", regions.at(p)));
	    PUyieldtot += ( tab1.Get("yield", regions.at(p)) +  tab2.Get("yield", regions.at(p))); //+ tab3.Get("PUup", regions.at(p)));
	    //PUyieldtot += tI.Get("yield", regions.at(p));
	}
	Figure two(2,0);
	//Figure perDiff = (PUUPtot-PUDNtot)/(PUUPtot+PUDNtot);
	Figure perDiff1 = (PUUPtot-PUyieldtot)/(PUyieldtot);
	Figure perDiff2 = (PUDNtot-PUyieldtot)/(PUyieldtot);
	Figure fPerDiff;
        if(abs(perDiff1.value()) > abs(perDiff2.value()))
        { 
            fPerDiff= Figure(abs(perDiff1.value()), perDiff1.error());
        }
        else
            fPerDiff= Figure(abs(perDiff2.value()), perDiff2.error());
	Figure PUY1(0.0);
	Figure PUY2(0.0);
	Figure PUY3(0.0);
	for(uint32_t p=0; p<regions.size(); p++)
	{
	    
	    PUY1 =  tab1.Get("yield", regions.at(p))  ;
	    PUY2 =  tab2.Get("yield", regions.at(p))  ;
	    //PUY3 =  tab3.Get("yield", regions.at(p))  ;
	    //if(PUY.value() !=0)
	    //{
		tab1.Set("PUdown", regions.at(p), PUY1 - (PUY1*fPerDiff));
		tab1.Set("PUup", regions.at(p), PUY1 + (PUY1*fPerDiff));
		tab2.Set("PUdown", regions.at(p), PUY2 - (PUY2*fPerDiff));
		tab2.Set("PUup", regions.at(p), PUY2 + (PUY2*fPerDiff));
		//tab3.Set("PUdown", regions.at(p), PUY3 - (PUY3*fPerDiff));
		//tab3.Set("PUup", regions.at(p), PUY3 + (PUY3*fPerDiff));

	}*/
        //@MJ@ TODO update one day thigher with new PU

        for(uint32_t r = 0; r<regions.size(); r++)
        {
            Figure y1;
            Figure y2;
            Figure y;
            //Figure y3;
            float unc1;
            float unc2;
            //float unc3;
            float bef1;
            float bef2;
            //float bef3;
            float uncHigher1;
            float uncHigher2;
            //float uncHigher3;
            float unc1Stat;
            float unc2Stat;
            //float unc3Stat;
            float unc1Syst;
            float unc2Syst;
            //float unc3Syst;
            float relJES = 0.04; //4%
            float relPU = 0.01; //4%
            Figure frelJES(0.04,0); //4%
            Figure frelPU(0.01,0); //4%
            for(uint32_t c = 0; c< colI.size(); c++)
            {
                Figure f1 = tab1.Get(colI.at(c), regions.at(r));
                Figure f2 = tab2.Get(colI.at(c), regions.at(r));
                //Figure f3 = tab3.Get(colI.at(c), regions.at(r));

                Figure f = f1+f2;//+f3;
                tabSum.Set(colI.at(c), regions.at(r),f );

                if(c==0)
                {
                    y1 = f1;
                    y2 = f2;
                    y=f;
                    //y3 = f3;
                    unc1 = 0;
                    unc2 = 0;
                    //unc3 = 0;
                    unc1Stat = 0;
                    unc2Stat = 0;
                    //unc3Stat = 0;
                    unc1Syst = 0;
                    unc2Syst = 0;
                    //unc3Syst = 0;

                    unc1 += y1.error()*y1.error();
                    unc2 += y2.error()*y2.error();
                    //unc3 += y3.error()*y3.error();
                    unc1Stat = y1.error()*y1.error();
                    unc2Stat = y2.error()*y2.error();
                    //unc3Stat = y3.error()*y3.error();

                }
                //@MJ@ fixed JES
                /*else if(colI.at(c) == "jesDN")
                {
                    tabSum.Set(colI.at(c), regions.at(r),y-(y*frelJES) ); 
                       unc1 += (y1.value()*relJES)* (y1.value()*relJES) ;
                       unc2 += (y2.value()*relJES)* (y2.value()*relJES) ;
                     
                       unc1Syst += (y1.value()*relJES)* (y1.value()*relJES) ;
                       unc2Syst += (y2.value()*relJES)* (y2.value()*relJES) ;
                }
                else if(colI.at(c) == "jesUP")
                {
                    tabSum.Set(colI.at(c), regions.at(r),y+(y*frelJES) ); 
                }*/
                //@MJ@ PU=1%
                /*else if(colI.at(c) == "PUdown")
                {
                    tabSum.Set(colI.at(c), regions.at(r),y-(y*frelPU) ); 
                       unc1 += (y1.value()*relPU)* (y1.value()*relPU) ;
                       unc2 += (y2.value()*relPU)* (y2.value()*relPU) ;
                     
                       unc1Syst += (y1.value()*relPU)* (y1.value()*relPU) ;
                       unc2Syst += (y2.value()*relPU)* (y2.value()*relPU) ;
                }
                else if(colI.at(c) == "PUup")
                {
                    tabSum.Set(colI.at(c), regions.at(r),y+(y*frelPU) ); 
                }*/
                else if(c%2 != 0)
                {
                    bef1 = f1.value();
                    bef2 = f2.value();
                    //bef3 = f3.value();
                }
                else if(c%2 ==0)
                {
                   uncHigher1 = abs(bef1-y1.value()) > abs(f1.value() - y1.value()) ? abs(bef1-y1.value()): abs(f1.value() - y1.value());
                   uncHigher2 = abs(bef2-y2.value()) > abs(f2.value() - y2.value()) ? abs(bef2-y2.value()): abs(f2.value() - y2.value());
                   //uncHigher3 = abs(bef3-y3.value()) > abs(f3.value() - y3.value()) ? abs(bef3-y3.value()): abs(f3.value() - y3.value());

                   if(colI.at(c) == "lepSFDN" || colI.at(c) == "lepSFUP" ||  colI.at(c) == "ISRnjetsDown" || colI.at(c) == "ISRnjetsUp"||  colI.at(c) == "PUdown" || colI.at(c) == "PUup")
                   {
                       cout << "in top pt 1" << endl;
                       unc1 += uncHigher1*uncHigher1;
                       unc2 += uncHigher2*uncHigher2;
                     //  unc3 += uncHigher3*uncHigher3;
                       unc1Syst += uncHigher1*uncHigher1;
                       unc2Syst += uncHigher2*uncHigher2;
                       //unc3Syst += uncHigher3*uncHigher3;
                   }
                   else
                   {
                       float rel1 = 0 ;
                       float rel2 = 0 ;
                       //float rel3= 0 ;
                       if(bef1+f1.value() != 0)
                           rel1 = (abs(bef1-f1.value())) /(bef1+f1.value()) ;
                       if(bef2+f2.value() != 0)
                           rel2 = (abs(bef2-f2.value())) /(bef2+f2.value()) ;
                       //if(bef3+f3.value() != 0)
                       //    rel3= (abs(bef3-f3.value())) / (bef3+f3.value()) ;

                       unc1 += (y1.value()*rel1)* (y1.value()*rel1) ;
                       unc2 += (y2.value()*rel2)* (y2.value()*rel2) ;
                       //unc3 += (y3.value()*rel3)* (y3.value()*rel3) ;
                     
                       unc1Syst += (y1.value()*rel1)* (y1.value()*rel1) ;
                       unc2Syst += (y2.value()*rel2)* (y2.value()*rel2) ;
                       //unc3Syst += (y3.value()*rel3)* (y3.value()*rel3)  ;

                   }

                }
                else
                {
                    cout << "something went terribly wrong" << endl;
                }

                if(c == colI.size() -1)
                {
                    //unc1 += (lumi*y1.value())*(lumi*y1.value()); //luminosity
                    //unc2 += (lumi*y2.value())*(lumi*y2.value());
                    //unc3 += (lumi*y3.value())*(lumi*y3.value());
                    //unc1Syst += (lumi*y1.value())*(lumi*y1.value());
                    //unc2Syst += (lumi*y2.value())*(lumi*y2.value());
                    //unc3Syst += (lumi*y3.value())*(lumi*y3.value());

                    totYields.Set(processes.at(0), regions.at(r), Figure(y1.value(),sqrt(unc1)));
                    totYields.Set(processes.at(1), regions.at(r), Figure(y2.value(),sqrt(unc2)));
                    //totYields.Set(processes.at(2), regions.at(r), Figure(y3.value(),sqrt(unc3)));
                    cout << "stat " << sqrt(unc1Stat) << " syst " << sqrt(unc1Syst) << endl;
                    totYieldsStatSyst.Set(processes.at(0), regions.at(r), Figure(y1.value(),sqrt(unc1Stat),sqrt(unc1Syst)));
                    totYieldsStatSyst.Set(processes.at(1), regions.at(r), Figure(y2.value(),sqrt(unc2Stat),sqrt(unc2Syst)));
                    //totYieldsStatSyst.Set(processes.at(2), regions.at(r), Figure(y3.value(),sqrt(unc3Stat),sqrt(unc3Syst)));
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


                   if(colI.at(c) == "lepSFDN" || colI.at(c) == "lepSFUP" ||  colI.at(c) == "ISRnjetsDown" || colI.at(c) == "ISRnjetsUp")
                   {
                       cout << "in top pt 2" << colI.at(c) << endl;
                       uncTot += uncHigher*uncHigher ;
                       uncTotSyst += uncHigher*uncHigher ;
                       if(yieldVal.value() != 0)
                           thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((uncHigher/yieldVal.value())*100, 0));
                       else
                           thigher.Set(systOutNames.at(c/2), regions.at(r), Figure(0, 0));

                   }
                   //@MJ@ fixed jes
                   /*else if(colI.at(c) == "jesDN" || colI.at(c) == "jesUP")
                   {
                       float trelJES = 0.04; //4%
                       uncTot += (trelJES*yieldVal.value()) * (trelJES*yieldVal.value()) ;
                       uncTotSyst +=  (trelJES*yieldVal.value()) * (trelJES*yieldVal.value()) ;

                       thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((trelJES)*100, 0));
                   }*/
                   //@MJ@ PU=3%
                   else if(colI.at(c) == "PUdown" || colI.at(c) == "PUup")
                   {
                       cout << "in here PU2" << endl;
                       float trelPU = 0.03; //1%
                       uncTot += (trelPU*yieldVal.value()) * (trelPU*yieldVal.value()) ;
                       uncTotSyst +=  (trelPU*yieldVal.value()) * (trelPU*yieldVal.value()) ;

                       cout << "region " << regions.at(r) << "uncPU" << (trelPU*yieldVal.value()) * (trelPU*yieldVal.value()) << endl;
                       thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((trelPU)*100, 0));
                   }
                   else
                   {

                       float center;
                       if(uncBef+resall.value() != 0)
                           center = (abs(uncBef-resall.value()))/(uncBef+resall.value());
                       else
                           center = 0;
                       
                       uncTot += (center*yieldVal.value()) * (center*yieldVal.value()) ;
                       uncTotSyst +=  (center*yieldVal.value()) * (center*yieldVal.value()) ;

                       thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((center)*100, 0));

                   }

                   if(c==colI.size()-1)
                       cout << "yield" << yieldVal.value() << "uncertainty is" << sqrt(uncTot) << endl;

                }
                else
                {
                    cout << "there is a mistake somewhere " << endl;
                }
                if(c == colI.size() -1)
                {
                    uint32_t lastEl = systOutNames.size();
                    //uncTot += (lumi*yieldVal.value())*(lumi*yieldVal.value());
                    //uncTotSyst += (lumi*yieldVal.value())*(lumi*yieldVal.value());


                    //thigher.Set(systOutNames.at(lastEl-2), regions.at(r),Figure(lumi*100,0));
                    thigher.Set(systOutNames.at(lastEl-1), regions.at(r), Figure((sqrt(uncTot)/yieldVal.value())*100,0));

                    totYields.Set(processes.at(2), regions.at(r), Figure(yieldVal.value(),sqrt(uncTot)));
                    totYieldsStatSyst.Set(processes.at(2), regions.at(r), Figure(yieldVal.value(),sqrt(uncTotStat),sqrt(uncTotSyst)));
                }

            }

        }

       //summary plots
        for(uint32_t s=0; s<systOutNames.size(); s++)
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
                    tsummary.Set(uncRange.at(0), systOutNames.at(s), oneSysts.at(z));
                    break;
                }
            }
            tsummary.Set(uncRange.at(1), systOutNames.at(s), oneSysts.at(regions.size()-1));
        }
        //tsummary.Set(uncRange.at(1), systOutNamesAll.at(systOutNames.size()), Figure(6.2,0));

        //@MJ@ TODO add lumi, CS

 
        tabSum.Print(static_cast<string>(group+"summedprocessesZnunu.tab"),4);
        tabSum.PrintLatex(static_cast<string>(group+"summedprocessesZnunu.tex"),4); 
        thigher.Print(static_cast<string>(group+"summedprocessesRelHigherZnunu.tab"),4);
        thigher.PrintLatex(static_cast<string>(group+"summedprocessesRelHigherZnunu.tex"),4, "noError"); 
        tsummary.Print(static_cast<string>(group+"summedprocessesRelHigherZnunuSummary.tab"),4, "noError");
        tsummary.PrintLatex(static_cast<string>(group+"summedprocessesRelHigherZnunuSummary.tex"),4, "noError"); 
        totYields.Print(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrors.tab"),4);
        totYields.PrintLatex(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrors.tex"),4); 
        totYieldsStatSyst.Print(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrorsStatSystBreakdown.tab"),4, "systError");
        totYieldsStatSyst.PrintLatex(static_cast<string>(group+"summedprocessesZnunuYieldsAndErrorsStatSystBreakdown.tex"),4, "systError"); 


}
