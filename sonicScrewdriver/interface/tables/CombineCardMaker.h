#ifndef CombineCardMaker_h
#define CombineCardMaker_h

#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "interface/Figure.h"
#include "interface/Table.h"
#include "interface/SonicScrewdriver.h"

using namespace std;

namespace theDoctor
{

    class CombineCardMaker : public Table
    {

     public:

            CombineCardMaker(){};

            CombineCardMaker(SonicScrewdriver* screwdriver, vector<string> inputRegionTags, string channel, string SigProcessClass, string sigVarX, string sigVarY)
	    {
                vector<string> rawProcessesTags;
                vector<string> rawProcessesLabels;
                vector<string> processesTags;
                vector<string> processesLabels;
                vector<string> regionsTags;
                vector<string> regionsLabels;

                screwdriver->GetProcessClassTagList  (&rawProcessesTags  );
                screwdriver->GetProcessClassLabelList(&rawProcessesLabels);

                screwdriver->GetRegionTagList  (&regionsTags  );
                screwdriver->GetRegionLabelList(&regionsLabels);
      

                // Sort backgrounds, signal, and add "total background" line
                for (unsigned int i = 0 ; i < rawProcessesTags.size() ; i++)
                {
                    string type = screwdriver->GetProcessClassType(rawProcessesTags[i]);
                    if (type == "background")
                    {
                        processesTags.push_back(rawProcessesTags[i]);
                        processesLabels.push_back(rawProcessesTags[i]);
                    }
                }

                processesTags.push_back("totalSM");
                processesLabels.push_back("total SM");

                bool dataFilled = false;
                for (unsigned int i = 0 ; i < rawProcessesTags.size() ; i++)
                {
                    string type2 = screwdriver->GetProcessClassType(rawProcessesTags[i]);
                    if (type2 == "data")
                    {
                        cout << "filling data in card maker " << endl;
                        if(dataFilled == true)
                            throw std::runtime_error("data are already present!");
                        processesTags.push_back(rawProcessesTags[i]);
                        processesLabels.push_back(rawProcessesTags[i]);
                        dataFilled = true;
                    }
                }

                /*
		for (unsigned int i = 0 ; i < rawProcessesTags.size() ; i++)
                {
                    string type = screwdriver->GetProcessClassType(rawProcessesTags[i]);
                    if (type == "signal")
                    {
                        processesTags.push_back(rawProcessesTags[i]);
                        processesLabels.push_back(rawProcessesTags[i]);
                    }
                }
		*/
		//Add signal labels by looking to non empty bins (for the moment in the first region only)
		 if(regionsTags.size()==0){
		 	cerr<<"CombinedCardMaker::no regions have been defined !"<<endl;
		 	return;
		 }
 
                 if(SigProcessClass != "none")
                 {
		 TH2D* histo2Dtmp = screwdriver-> Get2DHistoClone(sigVarX, sigVarY, SigProcessClass, regionsTags[0], channel);
		 for(int binX = 1 ; binX<histo2Dtmp->GetNbinsX(); binX++){

			for(int binY = 1 ; binY<histo2Dtmp->GetNbinsX(); binY++){
				if( histo2Dtmp->GetBinContent(binX, binY) !=0 ){
                        		stringstream ss;
					//rouding the masses by step of 5 GeV

					ss<<"("<<(int)histo2Dtmp->GetXaxis()->GetBinCenter(binX)<<","<<(int)histo2Dtmp->GetYaxis()->GetBinCenter(binY)<<")";
					processesTags.push_back(ss.str());
                        		processesLabels.push_back(ss.str());
				}
                        	

			}

		}
                }


                // Get labels for input regions
                if (inputRegionTags.size() != 0)
                {
                    vector<string> inputRegionLabels;

                    for (unsigned int i = 0 ; i < inputRegionTags.size() ; i++)
                    {
		      bool found = false;
                        for (unsigned int j = 0 ; j < regionsTags.size()  ; j++)
                        {
                            if (inputRegionTags[i] == regionsTags[j]) { inputRegionLabels.push_back(regionsLabels[j]); found = true; break; }
                        }
                        if (!found) WARNING_MSG << "Region " << inputRegionTags[i] << " was not found." << endl;
                    }

                    regionsTags   = inputRegionTags;
                    regionsLabels = inputRegionLabels;

                }

                Init(regionsTags,processesTags);
                SetLabels(regionsLabels,processesLabels);

                // Actually fill the table signal/bkg and bkg/bkg
                for (unsigned int r = 0 ; r < regionsTags.size()          ; r++)
                {
      		
		    //retrive signal yield from 2D plot 
		    TH2D* histo2DSigYields = screwdriver-> Get2DHistoClone(sigVarX, sigVarY, SigProcessClass, regionsTags[r], channel);
                    
		    Figure tmpTotal(0.0,0.0);
                    Figure tmpZero(0.0,0.0);
                    for (unsigned int p = 0 ; p < processesTags.size() ; p++)
                    {
                        if (processesTags[p] == "totalSM") continue;
                        string type = screwdriver->GetProcessClassType(processesTags[p]);
                        Figure currentYield = screwdriver->GetYieldAndError(processesTags[p],
                                                                            regionsTags[r],
                                                                            channel);
                       
			Set(regionsTags[r], processesTags[p], currentYield);
                        if (type != "background") continue;
                        tmpTotal += currentYield;

                    }

                    if(SigProcessClass != "none")
                    {
		    //Scan the Signal 2D plan
		    for(int binX = 1 ; binX<histo2DSigYields->GetNbinsX(); binX++){
			for(int binY = 1 ; binY<histo2DSigYields->GetNbinsX(); binY++){
				Figure sigYield (histo2DSigYields->GetBinContent(binX, binY), histo2DSigYields->GetBinError(binX, binY));
				Figure zbi = Zbi(sigYield.value(), tmpTotal.value());
                        	
				stringstream ss;
				ss<<"("<<(int)histo2DSigYields->GetXaxis()->GetBinCenter(binX)<<","<<(int)histo2DSigYields->GetYaxis()->GetBinCenter(binY)<<")";
				string label = ss.str();
                        	cout<<"yield: "<<ss.str()<<" "<<histo2DSigYields->GetBinContent(binX, binY);//<<endl;
				cout<<" "<<sigYield.value()<<" "<<tmpTotal.value();
				cout<<" "<<zbi.value()<<endl;

                            	//Set(regionsTags[r], label, zbi);
                            	Set(regionsTags[r], label, sigYield);
			}

		    }
                    }
		    ///*
                    for (unsigned int p = 0 ; p < processesTags.size() ; p++)
                    {
                        if (processesTags[p] == "totalSM"){
                        	/*
				Figure currentYield = screwdriver->GetYieldAndError(processesTags[p],
                                                                            regionsTags[r],
                                                                            channel);
                            	*/
				//put the sum of bkg
				cout<<"I found totalSM - "<<tmpTotal.value()<<endl;
				Set(regionsTags[r], processesTags[p], tmpTotal);
			
			}
		    	/*
                        if (processesTags[p] == "totalSM") continue;
                        string type = screwdriver->GetProcessClassType(processesTags[p]);
                        Figure currentYield = screwdriver->GetYieldAndError(processesTags[p],
                                                                            regionsTags[r],
                                                                            channel);
                        Figure zbi = Zbi(currentYield.value(), tmpTotal.value());

                        if (type == "background") 
                            Set(regionsTags[r], processesTags[p], tmpZero);
                        else 
                            Set(regionsTags[r], processesTags[p], zbi);
                    	*/
		    }
		    //*/
		    //Loop over signal points


                    //Set(regionsTags[r],"totalSM", tmpZero);
                }

            };

            double Zbi(double n_sig, double n_b, double rel_uncert_b = 0.3 )
            {
                double n_on = n_sig+n_b;
                double mu_b_hat = n_b;
                double sigma_b = rel_uncert_b*n_b;
                double tau = mu_b_hat/(sigma_b*sigma_b);
                double n_off = tau*mu_b_hat;
                double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
                double Z_Bi = sqrt(2)*TMath::ErfInverse(1 - 2*P_Bi);
                return Z_Bi;
            };
            ~CombineCardMaker()
            {
            };

	    void UpdateCardTable(string inputTab, string TFtab, string bkgType){

                Table inTable(inputTab);
                vector<string> processList = inTable.rowTags; 
                vector<string> regionsList = inTable.colTags;

                Table TFTable(TFtab);
                vector<string> TFRegionsList = TFTable.rowTags; 
                vector<string> TFList = TFTable.colTags;

                if(regionsList.size() != TFRegionsList.size())
                    throw std::runtime_error("Tables do not have same number of signal regions");
		
                uint16_t upColumnNr = 11111;
                for (unsigned int k = 0 ; k < TFList.size() ; k++)
                {
                   if(TFList.at(k).find(bkgType)!=std::string::npos)
                   {
                       upColumnNr = k;
                       break;
                   }
                }
                if( upColumnNr == 11111)
                    throw std::runtime_error("Background estimation column not found");

		for (unsigned int i = 0 ; i < processList.size() ; i++)
                {
			if( processList.at(i).find(bkgType)!=std::string::npos ) 
                        {
				for(unsigned int j = 0; j < regionsList.size() ; j++){ 
				    Figure y = TFTable.Get(upColumnNr, j);
                                    inTable.Set(regionsList.at(j), processList.at(i), y);
                                }
                        }
                }

		for (unsigned int i = 0 ; i < regionsList.size() ; i++)
                {
                    Figure totSMVal(0,0);
		    for (unsigned int j = 0 ; j < processList.size() ; j++)
                    {
	                if( !(processList.at(j).find("total SM")!=std::string::npos || processList.at(j).find("totalSM")!=std::string::npos) )
                        {
		            Figure z = inTable.Get(i,j);
                            totSMVal+=z;
                        }
                        else
                        {
                            inTable.Set(regionsList.at(i), processList.at(j), totSMVal); 
                            break;
                        }
                    }
                }
                inTable.Print(inputTab,4);
            };  

	    void ProduceCard(string inputTab, string outputdir,float sigRelUncert = 1.1, float bkgRelUncert = 1.2){

                Table inTable(inputTab);
                vector<string> processList = inTable.rowTags; 
                vector<string> regionsList = inTable.colTags; 
		vector<int> iBkgLine;
		for (unsigned int i = 0 ; i < processList.size() ; i++){
			cout<<processList.at(i)<<endl;
			if(processList.at(i).find("bkg")!=std::string::npos) iBkgLine.push_back(i);
		}
		if(iBkgLine.size()==0){
			cerr<<"CombinedCardMaker::backgrounds not found !"<<endl;
			return;
		}
		int iDataLine=-1;
		for (unsigned int i = 0 ; i < processList.size() ; i++){
			if(processList.at(i).find("data")!=std::string::npos) iDataLine =i;
		}
		if(iDataLine==-1){
			cout<<"CombinedCardMaker::no data"<<endl;
		}

                ofstream shfile;
                shfile.open (outputdir+"/combineCommands.sh");
                shfile << "rm allLimits.root" << endl;
		shfile << "rm allLimits2.root" << endl;
                shfile << "rm higgsCombine*.root" << endl;
		//looking for signal
		for (unsigned int i = 0 ; i < processList.size() ; i++){
			cout<<processList.at(i)<<endl;
			if(processList.at(i).find("(")!=std::string::npos && processList.at(i).find(")")!=std::string::npos){
				cout<<"selected !"<<endl;
				string filename;
				size_t found = processList.at(i).find(",");
				string val1 = processList.at(i).substr(processList.at(i).find("(")+1,found-1);
				string val2 = processList.at(i).substr(found+1,processList.at(i).find(")")-found-1);
				filename = outputdir+string("/")+val1+string("_")+val2+string(".tab");
                                shfile << "combine -M Asymptotic " << val1+string("_")+val2+string(".tab") << " -n " << val1+string("_")+val2 << endl;
				shfile <<"combine -M ProfileLikelihood --significance " << val1+string("_")+val2+string(".tab") <<" -t -1 --expectSignal=1 --toysFreq -n "<< val1+string("_")+val2 <<endl;
		



				ofstream ofile(filename.c_str());
				ofile<<"imax "<<regionsList.size()<<" number of bins"<<endl;
				ofile<<"jmax "<< iBkgLine.size() <<" number of processes"<<endl; // bkg + signal only
				ofile<<"kmax "<<regionsList.size()*2+3<<" number of nuisance parameters"<<endl; // only a systematic on the bkg
					//	ofile<<"kmax "<<5<<" number of nuisance parameters"<<endl; // only a systematic on the bkg
				ofile<<"----------------------------------------------------"<<endl;
				ofile<<"bin\t";
				for(unsigned int j = 0; j < regionsList.size() ; j++){ ofile<<"ch"<<j<<"\t"; }
				ofile<<endl;
				ofile<<"----------------------------------------------------"<<endl;
				ofile<<"observation\t";
				for(unsigned int j = 0; j < regionsList.size() ; j++){ ofile<<inTable.Get(j,iDataLine).value()<<"\t"; }
				ofile<<endl;
				ofile<<"----------------------------------------------------"<<endl;
				ofile<<"bin\t";
				for(unsigned int j = 0; j < regionsList.size() ; j++)
                                {
                                    ofile<<"ch"<<j<<"\t"; 
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
                                        ofile<<"ch"<<j<<"\t"; 
                                    }
                                };
				ofile<<endl;
				ofile<<"process\t";
				//for(unsigned int j = 0; j < nCol ; j++){ ofile<<"bkg\tsignal\t"; }
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
                                { 
                                    ofile<<"sig\t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
                                        ofile<<processList.at(l) <<"\t"; 
                                    }
                                }
				ofile<<endl;
				ofile<<"process\t";
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
                                { 
                                    ofile<<"0\t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
                                        ofile<<l+1<<"\t"; 
                                    }
                                }
				ofile<<endl;
				ofile<<"rate\t";
				for(unsigned int j = 0; j < regionsList.size() ; j++){ 
					//signal yield (current one)
			            ofile<<inTable.Get(j,i).value()<<"\t";
					//background yield
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
					ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<"\t";
                                    }
				}
				ofile<<endl;
				ofile<<"----------------------------------------------------"<<endl;


				if(regionsList.size()==4)

				  {
				    


				//------------------------For I with table incertenties -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------			

				    
				    ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	
				    
				    
				    vector<double> table_bkg_lep_top;
				    string ligne2;
				    ifstream fichier2("./incertitude_lep_from_top.tab");			
				    while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg1\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromTop")
					      {
						
						ofile<<1+table_bkg_lep_top[j+27]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				    fichier2.close();
				    
				    
				    vector<double> table_bkg_Z_nunu;
				    string ligne3;
				    ifstream fichier3("./incertitude_Z_nunu.tab");
				    
				    while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg2\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgZnunu")
					      {
						
						ofile<<1+table_bkg_Z_nunu[j+27]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				    fichier3.close();
				    
				    
				    
				    
				    vector<double> table_bkg_lep_W;
				    string ligne;
				    ifstream fichier("./incertitude_lep_from_W.tab");
				    
				    while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					
				      }
				    
				    
				    
				    
				    for(unsigned int k = 0; k  < regionsList.size() ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					double x;
					x=1/(table_bkg_lep_W[k+27]*table_bkg_lep_W[k+27]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier.close();
				    
				    
				    
				    vector<double> table_bkg_lost_lep;
				    string ligne1;
				    ifstream fichier1("./incertitude_lost_lepton.tab");
				    while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					
				      }
				    
				    
				    
				    for(unsigned int k = 0; k  < regionsList.size() ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					double x;
					x=1/(table_bkg_lost_lep[k+27]*table_bkg_lost_lep[k+27]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier1.close();
				    
				    
				    

				  

				//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




				//---------------------------------- For I 0.1 signal and 0.2 bkg --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				    //	    	ofile<<"uncert_bkg	 lnN 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2"<<endl;   




					//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


					//--------------- 0.2 decorelated 5 incertitude-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


				    /*		ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	"<<endl; 	
				
				ofile<<"uncert_bkg1\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgOneLepFromW")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg2\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgLostLepton")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg3\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgOneLepFromTop")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg4\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgZnunu")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;*/


				//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

					//------------------ same mais lost lep et from W decorelé -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


				    /*					
	ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	"<<endl; 	
				
			



				ofile<<"uncert_bkg1\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgOneLepFromTop")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg2\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgZnunu")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;

				
				for(unsigned int k = 0; k  < regionsList.size() ; k++)
				  { 
				    ofile<<"uncert_bkg"<<3+k<<"\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<" \t";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromW" && k==j)
					      {
						ofile<<bkgRelUncert<<" \t";	
					      }
					    else
					      { 
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				  }

				
				for(unsigned int k = 0; k  < regionsList.size() ; k++)
				  { 
				    ofile<<"uncert_bkg2"<<34+k<<"\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<" \t";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgLostLepton" && k==j)
					      {
						ofile<<bkgRelUncert<<" \t";	
					      }
					    else
					      { 
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				    }*/
					

				//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



				//------------------------------------------- vrai incertitude sur tout correlé --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			
				    /*		vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					    
					  }
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					

					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					

					ofile<<"uncert_bkg\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"1.1"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j+27]<<" \t";
						  }
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j+27]<<" \t";
						  }
						if (processList[l]=="bkgOneLepFromW" )//&& j==k)
						  {
						    
						    ofile<<1+table_bkg_lep_W[j+27]<<" \t";
						  }
						if (processList[l]=="bkgLostLepton")// && j==k)
						  {
						    
						    ofile<<1+table_bkg_lost_lep[j+27]<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					fichier3.close();
					fichier.close();
					fichier1.close();
				*/


				//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


				//-------------------------------- incertitude papier 5 incertitude correlé ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


				    /*
					ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	


					
					vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					
					{
					
					table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					
					}
					
					ofile<<"uncert_bkg1\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					{ 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j+27]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg2\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j+27]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier3.close();
					



					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg3\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" )//&& j==k)
						  {
						    
						    ofile<<1+table_bkg_lep_W[j+27]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier.close();


					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					
					
					ofile<<"uncert_bkg4\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton")// && j==k)
						  {
						    
						    ofile<<1+table_bkg_lost_lep[j+27]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier1.close();
				
*/


				//------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- 



				





				//-------------- -------------  vrai incertitude mais from W et lost lep décorelé ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- ------------- 

					
				    /*
				
					ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- "<<endl; 	


					vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg1\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg2\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier3.close();
					



					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					
					for(unsigned int k = 0; k  < regionsList.size() ; k++)
					  { 
					    ofile<<"uncert_bkg"<<3+k<<"\t lnN \t";
					    //for(unsigned int j = 0; j < nCol ; j++)
					    for(unsigned int j = 0; j  < regionsList.size() ; j++)
					      { 
						ofile<<"-"<<"\t ";
						for(uint32_t l=0; l<iBkgLine.size(); l++ )
						  {
						    if (processList[l]=="bkgOneLepFromW" && j==k)
						      {
							
							ofile<<1+table_bkg_lep_W[j]<<" \t";
						      }
						    else
						      {
							ofile<<"-"<<" \t";
						      }
						  }	
					      }
					    ofile<<endl;
					  }
					fichier.close();


					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					
					
					for(unsigned int k = 0; k  < regionsList.size() ; k++)
					  { 
					    ofile<<"uncert_bkg"<<34+k<<"\t lnN \t";
					    //for(unsigned int j = 0; j < nCol ; j++)
					    for(unsigned int j = 0; j  < regionsList.size() ; j++)
					      { 
						ofile<<"-"<<"\t ";
						for(uint32_t l=0; l<iBkgLine.size(); l++ )
						  {
						    if (processList[l]=="bkgLostLepton" && j==k)
						      {
							
							ofile<<1+table_bkg_lost_lep[j]<<" \t";
						      }
						    else
						      {
							ofile<<"-"<<" \t";
						      }
						  }	
					      }
					    ofile<<endl;
					  }
					fichier1.close();
					
				    */
				//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				  }
				
				else if(regionsList.size()==8)

				  {



				    ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-     1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	
				    
				    
				    vector<double> table_bkg_lep_top;
				    string ligne2;
				    ifstream fichier2("./incertitude_lep_from_top.tab");			
				    while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg1\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				
				    for(unsigned int j = 0; j  < 4 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromTop")
					      {
						
						ofile<<1+table_bkg_lep_top[j+27]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    for(unsigned int j = 4; j  < 8 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromTop")
					      {
						
						ofile<<1+table_bkg_lep_top[j+23]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				      
				    ofile<<endl;
				    fichier2.close();
				    
				    
				    vector<double> table_bkg_Z_nunu;
				    string ligne3;
				    ifstream fichier3("./incertitude_Z_nunu.tab");
				    
				    while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg2\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < 4 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgZnunu")
					      {
						
						ofile<<1+table_bkg_Z_nunu[j+27]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				  for(unsigned int j = 4; j  < 8 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgZnunu")
					      {
						
						ofile<<1+table_bkg_Z_nunu[j+23]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				  
				    ofile<<endl;
				    fichier3.close();
				    
				    
				    
				    
				    vector<double> table_bkg_lep_W;
				    string ligne;
				    ifstream fichier("./incertitude_lep_from_W.tab");
				    
				    while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					
				      }
				    
				    
				    
				    
				    for(unsigned int k = 0; k  < 4 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					double x;
					x=1/(table_bkg_lep_W[k+27]*table_bkg_lep_W[k+27]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }

   
				    for(unsigned int k = 4; k  < 8 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					double x;
					x=1/(table_bkg_lep_W[k+23]*table_bkg_lep_W[k+23]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier.close();
				    
				    
				    
				    vector<double> table_bkg_lost_lep;
				    string ligne1;
				    ifstream fichier1("./incertitude_lost_lepton.tab");
				    while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					
				      }
				    
				    
				    
				    for(unsigned int k = 0; k  < 4 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					double x;
					x=1/(table_bkg_lost_lep[k+27]*table_bkg_lost_lep[k+27]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    

				    for(unsigned int k = 4; k  < 8 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					double x;
					x=1/(table_bkg_lost_lep[k+23]*table_bkg_lost_lep[k+23]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier1.close();
				    

				  }


				//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

				else if(regionsList.size()==9)
				  
				  {
				    
				    
				    
				    ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-     1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-      1.1	-      - 	 -	-"<<endl; 	
				    
				    
				    vector<double> table_bkg_lep_top;
				    string ligne2;
				    ifstream fichier2("./incertitude_lep_from_top.tab");			
				    while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg1\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				
				    for(unsigned int j = 0; j  < 5 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromTop")
					      {
						
						ofile<<1+table_bkg_lep_top[j+7]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    for(unsigned int j = 5; j  < 9 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromTop")
					      {
						
						ofile<<1+table_bkg_lep_top[j+22]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				      
				    ofile<<endl;
				    fichier2.close();
				    
				    
				    vector<double> table_bkg_Z_nunu;
				    string ligne3;
				    ifstream fichier3("./incertitude_Z_nunu.tab");
				    
				    while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg2\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < 5 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgZnunu")
					      {
						
						ofile<<1+table_bkg_Z_nunu[j+7]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				  for(unsigned int j = 5; j  < 9 ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgZnunu")
					      {
						
						ofile<<1+table_bkg_Z_nunu[j+22]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				  
				    ofile<<endl;
				    fichier3.close();
				    
				    
				    
				    
				    vector<double> table_bkg_lep_W;
				    string ligne;
				    ifstream fichier("./incertitude_lep_from_W.tab");
				    
				    while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					
				      }
				    
				    
				    
				    
				    for(unsigned int k = 0; k  < 5 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					double x;
					x=1/(table_bkg_lep_W[k+7]*table_bkg_lep_W[k+7]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }

   
				    for(unsigned int k = 5; k  < 9 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					double x;
					x=1/(table_bkg_lep_W[k+22]*table_bkg_lep_W[k+22]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier.close();
				    
				    
				    
				    vector<double> table_bkg_lost_lep;
				    string ligne1;
				    ifstream fichier1("./incertitude_lost_lepton.tab");
				    while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					
				      }
				    
				    
				    
				    for(unsigned int k = 0; k  < 5 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					double x;
					x=1/(table_bkg_lost_lep[k+7]*table_bkg_lost_lep[k+7]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    

				    for(unsigned int k = 5; k  < 9 ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					double x;
					x=1/(table_bkg_lost_lep[k+22]*table_bkg_lost_lep[k+22]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier1.close();
				    

				  }


				//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



				
				else if(regionsList.size()==5)
				  
				  {
				    
				    ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-      1.1	 - 	- 	- 	-"<<endl; 	
				    
				    
				    vector<double> table_bkg_lep_top;
				    string ligne2;
				    ifstream fichier2("./incertitude_lep_from_top.tab");			
				    while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg1\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromTop")
					      {
						
						ofile<<1+table_bkg_lep_top[j+7]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				    fichier2.close();
				    
				    
				    vector<double> table_bkg_Z_nunu;
				    string ligne3;
				    ifstream fichier3("./incertitude_Z_nunu.tab");
				    
				    while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					
				      }
				    
				    ofile<<"uncert_bkg2\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<"\t ";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgZnunu")
					      {
						
						ofile<<1+table_bkg_Z_nunu[j+7]<<" \t";
					      }
					    else
					      {
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				    fichier3.close();
				    
				    
				    
				    
				    vector<double> table_bkg_lep_W;
				    string ligne;
				    ifstream fichier("./incertitude_lep_from_W.tab");
				    
				    while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					
				      }
				    
				    
				    
				    
				    for(unsigned int k = 0; k  < regionsList.size() ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					double x;
					x=1/(table_bkg_lep_W[k+7]*table_bkg_lep_W[k+7]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier.close();
				    
				    
				    
				    vector<double> table_bkg_lost_lep;
				    string ligne1;
				    ifstream fichier1("./incertitude_lost_lepton.tab");
				    while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
				      
				      {
					
					table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					
				      }
				    
				    
				    
				    for(unsigned int k = 0; k  < regionsList.size() ; k++)
				      { 
					cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					double x;
					x=1/(table_bkg_lost_lep[k+7]*table_bkg_lost_lep[k+7]);
					int y;
					y= (int) x;
					if (y==0)
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					  }
					else
					  {
					    ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					  }
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton" && j==k)
						  {
						    if (y==0)
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
						      }
						    else
						      {
							ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
						      }
						  }
						else
						  { 
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					
					
				      }
				    fichier1.close();
				  
			}

  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

				else
				  {
				    
				    //--------------------------  une incertitude tout correlé signal :0.1 bkg: 0.2----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				    /*	ofile<<"uncert_bkg	 lnN 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	1.2 	1.2	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2      1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2    1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2 	1.1	 1.2 	1.2 	1.2 	1.2"<<endl; 	*/
				//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


				//--------------------   5 incertitudes chacune correlée---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				/*	ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	
				
				ofile<<"uncert_bkg1\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgOneLepFromW")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg2\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgLostLepton")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg3\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgOneLepFromTop")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg4\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgZnunu")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;*/
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




//pareil mais lost lepton et from W sont décorelé ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

				/*	ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	
				
			



				ofile<<"uncert_bkg1\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgOneLepFromTop")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;


				ofile<<"uncert_bkg2\t lnN \t";
				//for(unsigned int j = 0; j < nCol ; j++)
				for(unsigned int j = 0; j  < regionsList.size() ; j++)
				  { 
				    ofile<<"-"<<" \t";
                                    for(uint32_t l=0; l<iBkgLine.size(); l++ )
                                    {
				      if (processList[l]=="bkgZnunu")
					{
					  ofile<<bkgRelUncert<<" \t";	
					}
				      else
					{ 
					  ofile<<"-"<<" \t";
		                  	}
                                }	
				}
				ofile<<endl;

				
				for(unsigned int k = 0; k  < regionsList.size() ; k++)
				  { 
				    ofile<<"uncert_bkg"<<3+k<<"\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<" \t";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgOneLepFromW" && k==j)
					      {
						ofile<<bkgRelUncert<<" \t";	
					      }
					    else
					      { 
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				  }

				
				for(unsigned int k = 0; k  < regionsList.size() ; k++)
				  { 
				    ofile<<"uncert_bkg2"<<34+k<<"\t lnN \t";
				    //for(unsigned int j = 0; j < nCol ; j++)
				    for(unsigned int j = 0; j  < regionsList.size() ; j++)
				      { 
					ofile<<"-"<<" \t";
					for(uint32_t l=0; l<iBkgLine.size(); l++ )
					  {
					    if (processList[l]=="bkgLostLepton" && k==j)
					      {
						ofile<<bkgRelUncert<<" \t";	
					      }
					    else
					      { 
						ofile<<"-"<<" \t";
					      }
					  }	
				      }
				    ofile<<endl;
				  }

				*/
				//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





					//------------------ 5 incertitudes avec incertitudes du papier ln N---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

				
				/*	ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	



					vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg1\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg2\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier3.close();
					



					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg3\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromW" )//&& j==k)
						  {
						    
						    ofile<<1+table_bkg_lep_W[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier.close();


					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					
					
					ofile<<"uncert_bkg4\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgLostLepton")// && j==k)
						  {
						    
						    ofile<<1+table_bkg_lost_lep[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier1.close();
				*/
				//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------








////-----------------------incertitude du papier mais from W et lost lepton decorelé-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


				/*
					ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	


					vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg1\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg2\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier3.close();
					



					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					
					for(unsigned int k = 0; k  < regionsList.size() ; k++)
					  { 
					    ofile<<"uncert_bkg"<<3+k<<"\t lnN \t";
					    //for(unsigned int j = 0; j < nCol ; j++)
					    for(unsigned int j = 0; j  < regionsList.size() ; j++)
					      { 
						ofile<<"-"<<"\t ";
						for(uint32_t l=0; l<iBkgLine.size(); l++ )
						  {
						    if (processList[l]=="bkgOneLepFromW" && j==k)
						      {
							
							ofile<<1+table_bkg_lep_W[j]<<" \t";
						      }
						    else
						      {
							ofile<<"-"<<" \t";
						      }
						  }	
					      }
					    ofile<<endl;
					  }
					fichier.close();


					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					
					
					for(unsigned int k = 0; k  < regionsList.size() ; k++)
					  { 
					    ofile<<"uncert_bkg"<<34+k<<"\t lnN \t";
					    //for(unsigned int j = 0; j < nCol ; j++)
					    for(unsigned int j = 0; j  < regionsList.size() ; j++)
					      { 
						ofile<<"-"<<"\t ";
						for(uint32_t l=0; l<iBkgLine.size(); l++ )
						  {
						    if (processList[l]=="bkgLostLepton" && j==k)
						      {
							
							ofile<<1+table_bkg_lost_lep[j]<<" \t";
						      }
						    else
						      {
							ofile<<"-"<<" \t";
						      }
						  }	
					      }
					    ofile<<endl;
					  }
					fichier1.close();
				
				
				*/

				//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




//-------------------------------------------------------background from paper all corelated ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
				/*		
					vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					    
					  }
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					

					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					

					ofile<<"uncert_bkg\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"1.1"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j]<<" \t";
						  }
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j]<<" \t";
						  }
						if (processList[l]=="bkgOneLepFromW" )//&& j==k)
						  {
						    
						    ofile<<1+table_bkg_lep_W[j]<<" \t";
						  }
						if (processList[l]=="bkgLostLepton")// && j==k)
						  {
						    
						    ofile<<1+table_bkg_lost_lep[j]<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					fichier3.close();
					fichier.close();
					fichier1.close();
				*/

					///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



				
					ofile<<"uncert_bkg	 lnN 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	- 	1.1	 - 	- 	- 	-"<<endl; 	


					vector<double> table_bkg_lep_top;
					string ligne2;
					ifstream fichier2("./incertitude_lep_from_top.tab");			
					while(getline(fichier2, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_top.push_back(stod(ligne2));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg1\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgOneLepFromTop")
						  {
						    
						    ofile<<1+table_bkg_lep_top[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier2.close();
					
					
					vector<double> table_bkg_Z_nunu;
					string ligne3;
					ifstream fichier3("./incertitude_Z_nunu.tab");
					
					while(getline(fichier3, ligne3))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_Z_nunu.push_back(stod(ligne3));  // on l'affiche
					    
					  }
					
					ofile<<"uncert_bkg2\t lnN \t";
					//for(unsigned int j = 0; j < nCol ; j++)
					for(unsigned int j = 0; j  < regionsList.size() ; j++)
					  { 
					    ofile<<"-"<<"\t ";
					    for(uint32_t l=0; l<iBkgLine.size(); l++ )
					      {
						if (processList[l]=="bkgZnunu")
						  {
						    
						    ofile<<1+table_bkg_Z_nunu[j]<<" \t";
						  }
						else
						  {
						    ofile<<"-"<<" \t";
						  }
					      }	
					  }
					ofile<<endl;
					fichier3.close();
					



					vector<double> table_bkg_lep_W;
					string ligne;
					ifstream fichier("./incertitude_lep_from_W.tab");
					
					while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lep_W.push_back(stod(ligne));  // on l'affiche
					    
					  }
					



					for(unsigned int k = 0; k  < regionsList.size() ; k++)
					  { 
					    cout<<regionsList.size()<<" "<<table_bkg_lep_W.size()<<endl;
					    double x;
					    x=1/(table_bkg_lep_W[k]*table_bkg_lep_W[k]);
					    int y;
					    y= (int) x;
					    if (y==0)
					      {
						ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<1<<" \t";
					      }
					    else
					      {
						ofile<<"uncert_bkg"<<k+3<<"\t gmN  "<<y<<" \t";
					      }
					    //for(unsigned int j = 0; j < nCol ; j++)
					    for(unsigned int j = 0; j  < regionsList.size() ; j++)
					      { 
						ofile<<"-"<<"\t ";
						for(uint32_t l=0; l<iBkgLine.size(); l++ )
						  {
						    if (processList[l]=="bkgOneLepFromW" && j==k)
						      {
							if (y==0)
							  {
							    ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
							  }
							else
							  {
							    ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
							  }
						      }
						    else
						      { 
							ofile<<"-"<<" \t";
						      }
						  }	
					      }
					    ofile<<endl;
					    

					  }
					fichier.close();


					
					vector<double> table_bkg_lost_lep;
					string ligne1;
					ifstream fichier1("./incertitude_lost_lepton.tab");
					while(getline(fichier1, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
					  
					  {
					    
					    table_bkg_lost_lep.push_back(stod(ligne1));  // on l'affiche
					    
					  }
					
				
				   
					for(unsigned int k = 0; k  < regionsList.size() ; k++)
					  { 
					    cout<<regionsList.size()<<" "<<table_bkg_lost_lep.size()<<endl;
					    double x;
					    x=1/(table_bkg_lost_lep[k]*table_bkg_lost_lep[k]);
					    int y;
					    y= (int) x;
					    if (y==0)
					      {
						ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<1<<" \t";
					      }
					    else
					      {
						ofile<<"uncert_bkg"<<k+34<<"\t gmN  "<<y<<" \t";
					      }
					    //for(unsigned int j = 0; j < nCol ; j++)
					    for(unsigned int j = 0; j  < regionsList.size() ; j++)
					      { 
						ofile<<"-"<<"\t ";
						for(uint32_t l=0; l<iBkgLine.size(); l++ )
						  {
						    if (processList[l]=="bkgLostLepton" && j==k)
						      {
							if (y==0)
							  {
							    ofile<<inTable.Get(j,iBkgLine.at(l)).value()<<" \t";
							  }
							else
							  {
							    ofile<<inTable.Get(j,iBkgLine.at(l)).value()/y<<" \t";
							  }
						      }
						    else
						      { 
							ofile<<"-"<<" \t";
						      }
						  }	
					      }
					    ofile<<endl;
					    

					  }
				      	fichier1.close();




				  }
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

				ofile.close();

			}
		}
                shfile << "hadd allLimits.root higgsCombine*.Asymptotic.mH120.root" << endl;
		shfile << "hadd allLimits2.root higgsCombine*.ProfileLikelihood.mH120.root" << endl;
		shfile.close();
	    	
	    }

    };

}

#endif
