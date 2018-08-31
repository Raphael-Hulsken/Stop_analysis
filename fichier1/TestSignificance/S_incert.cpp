#include <exception>
#include <iostream>
#include "stdlib.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TStyle.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"
#include "interface/Figure.h"

using namespace std;
using namespace theDoctor;
//------------------
//Function to compute the significance
//------------------
double Zbi(double n_sig, double n_b, double rel_uncert_b = 0.3 )
{
  //if (n_sig==0){return 0;}
  // else 
  //{
  double n_on = n_sig+n_b;
  //cout<<"n_on ="<<n_on<<endl;
  double mu_b_hat = n_b;
  //cout<<"mu_b_hat ="<<mu_b_hat<<endl;
  double sigma_b = rel_uncert_b*n_b;
  //cout<<"sigma_b ="<<sigma_b<<endl;
  double tau = mu_b_hat/(sigma_b*sigma_b);
  // cout<<"tau ="<<tau<<endl;
  double n_off = tau*mu_b_hat;
  //cout<<"n_off ="<<n_off<<endl;
  double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);//TMath::Beta(n_on,1+n_off);
  //  cout<<"P_Bi ="<<P_Bi<<endl;
  double Z_Bi = sqrt(2)*TMath::ErfInverse(1 - 2*P_Bi);
  //  cout<<"Z_Bi ="<<Z_Bi<<endl;
  if (Z_Bi<=0)
    {
      return 0;
    }
  else 
    { 
      return Z_Bi;
    }
  //}
}


//int main(int argc, char *argv[])
int main()
{

  vector<double> table;
  string ligne;
  string ligne1;
  string ligne2;

  ofstream fichier2("value2.dat");
  // fichier2<<"b,Zbi,combine"<<endl;
  


 for(double i=0.1; i<1;i+=0.1)
    {
      //cout<< Zbi(2,i,0.2)<<endl;
      for (double j =0;j<1;j+=0.2)
	{
      
	  ofstream fichier1("1900_50.tab");
	  fichier1<<"imax 1 number of bins"<<endl;
	  fichier1<<"jmax 1 number of processes"<<endl;
	  fichier1<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"observation	0		"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	ch0		"<<endl;
	  fichier1<<"process	sig	totbackg	"<<endl;
	  fichier1<<"process	0	1		"<<endl;
	  fichier1<<"rate	2	"<<i<<"	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"uncert_bkg	 lnN 	"<<1+j<<"	 1.2 	"<<endl;
	  fichier1.close();
	  
	  ofstream fichier3("1600_50.tab");
	  fichier3<<"imax 1 number of bins"<<endl;
	  fichier3<<"jmax 1 number of processes"<<endl;
	  fichier3<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"bin	ch0	"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"observation	0		"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"bin	ch0	ch0		"<<endl;
	  fichier3<<"process	sig	totbackg	"<<endl;
	  fichier3<<"process	0	1		"<<endl;
	  fichier3<<"rate	5	"<<i<<"	"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"uncert_bkg	 lnN 	"<<1+j<<"      	 1.2 	"<<endl;
	  fichier3.close();
	  
	  ofstream fichier4("1500_50.tab");
	  fichier4<<"imax 1 number of bins"<<endl;
	  fichier4<<"jmax 1 number of processes"<<endl;
	  fichier4<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"bin	ch0	"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"observation	0		"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"bin	ch0	ch0		"<<endl;
	  fichier4<<"process	sig	totbackg	"<<endl;
	  fichier4<<"process	0	1		"<<endl;
	  fichier4<<"rate	10	"<<i<<"	"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"uncert_bkg	 lnN 	"<<1+j<<"	 1.2 	"<<endl;
	  fichier4.close();
	  
	  system("sh combineCommands3.sh >> essaie3.txt");
	  system("sh combineCommands4.sh >> essaie4.txt");
	  system("sh combineCommands5.sh >> essaie5.txt");
	  ifstream fichier("essaie3.txt");
	  ifstream fichier5("essaie4.txt");
	  ifstream fichier6("essaie5.txt");
	  int x=0;
	  int y=0;
	  int z=0;
	  while( getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
	    
	    {
	      x++;
	      // cout<<"AAAA"<<endl;
	      if (x==3)
		{
		  while( getline(fichier5, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
		    
		    {
		      y++;
		      // cout<<"BBBB"<<endl;
		      if (y==3)
			{
			  while( getline(fichier6, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
			    
			    {
			      z++;
			      //cout<<"CCCC"<<endl;
			      if(x==3 && y==3 && z==3)
				{
				  size_t found = ligne.find(" ");
				  string val_neutr = ligne.substr(found+1);
				  size_t found1 = ligne1.find(" ");
				  string val_neutr1 = ligne1.substr(found1+1);
				  size_t found2 = ligne2.find(" ");
				  string val_neutr2 = ligne2.substr(found2+1);
				  //float spart_mass_point = stof(val_spart);
				  // float neutr_mass_point = stof(val_neutr);
				  fichier2<<i<<" "<<val_neutr<<" "<<val_neutr1<<" "<<val_neutr2<<" " ; // on l'affiche
				}
			    }
			}
		    }
		}
	    }
	  fichier.close();
	  fichier5.close();
	  fichier6.close();
	
	}
        fichier2<<endl;
    }



 for(double i=1; i<10;i++)
    {
      //cout<< Zbi(2,i,0.2)<<endl;
      for (double j =0;j<1;j+=0.2)
	{
      
	  ofstream fichier1("1900_50.tab");
	  fichier1<<"imax 1 number of bins"<<endl;
	  fichier1<<"jmax 1 number of processes"<<endl;
	  fichier1<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"observation	0		"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	ch0		"<<endl;
	  fichier1<<"process	sig	totbackg	"<<endl;
	  fichier1<<"process	0	1		"<<endl;
	  fichier1<<"rate	2	"<<i<<"	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"uncert_bkg	 lnN 	"<<1+j<<"	 1.2 	"<<endl;
	  fichier1.close();
	  
	  ofstream fichier3("1600_50.tab");
	  fichier3<<"imax 1 number of bins"<<endl;
	  fichier3<<"jmax 1 number of processes"<<endl;
	  fichier3<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"bin	ch0	"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"observation	0		"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"bin	ch0	ch0		"<<endl;
	  fichier3<<"process	sig	totbackg	"<<endl;
	  fichier3<<"process	0	1		"<<endl;
	  fichier3<<"rate	5	"<<i<<"	"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"uncert_bkg	 lnN 	"<<1+j<<"      	 1.2 	"<<endl;
	  fichier3.close();
	  
	  ofstream fichier4("1500_50.tab");
	  fichier4<<"imax 1 number of bins"<<endl;
	  fichier4<<"jmax 1 number of processes"<<endl;
	  fichier4<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"bin	ch0	"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"observation	0		"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"bin	ch0	ch0		"<<endl;
	  fichier4<<"process	sig	totbackg	"<<endl;
	  fichier4<<"process	0	1		"<<endl;
	  fichier4<<"rate	10	"<<i<<"	"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"uncert_bkg	 lnN 	"<<1+j<<"	 1.2 	"<<endl;
	  fichier4.close();
	  
	  system("sh combineCommands3.sh >> essaie3.txt");
	  system("sh combineCommands4.sh >> essaie4.txt");
	  system("sh combineCommands5.sh >> essaie5.txt");
	  ifstream fichier("essaie3.txt");
	  ifstream fichier5("essaie4.txt");
	  ifstream fichier6("essaie5.txt");
	  int x=0;
	  int y=0;
	  int z=0;
	  while( getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
	    
	    {
	      x++;
	      // cout<<"AAAA"<<endl;
	      if (x==3)
		{
		  while( getline(fichier5, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
		    
		    {
		      y++;
		      // cout<<"BBBB"<<endl;
		      if (y==3)
			{
			  while( getline(fichier6, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
			    
			    {
			      z++;
			      //cout<<"CCCC"<<endl;
			      if(x==3 && y==3 && z==3)
				{
				  size_t found = ligne.find(" ");
				  string val_neutr = ligne.substr(found+1);
				  size_t found1 = ligne1.find(" ");
				  string val_neutr1 = ligne1.substr(found1+1);
				  size_t found2 = ligne2.find(" ");
				  string val_neutr2 = ligne2.substr(found2+1);
				  //float spart_mass_point = stof(val_spart);
				  // float neutr_mass_point = stof(val_neutr);
				  fichier2<<i<<" "<<val_neutr<<" "<<val_neutr1<<" "<<val_neutr2<<" " ; // on l'affiche
				}
			    }
			}
		    }
		}
	    }
	  fichier.close();
	  fichier5.close();
	  fichier6.close();
	
	}
        fichier2<<endl;
    }




 for(double i=10; i<100;i+=10)
    {
      //cout<< Zbi(2,i,0.2)<<endl;
      for (double j =0;j<1;j+=0.2)
	{
      
	  ofstream fichier1("1900_50.tab");
	  fichier1<<"imax 1 number of bins"<<endl;
	  fichier1<<"jmax 1 number of processes"<<endl;
	  fichier1<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"observation	0		"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	ch0		"<<endl;
	  fichier1<<"process	sig	totbackg	"<<endl;
	  fichier1<<"process	0	1		"<<endl;
	  fichier1<<"rate	2	"<<i<<"	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"uncert_bkg	 lnN 	"<<1+j<<"	 1.2 	"<<endl;
	  fichier1.close();
	  
	  ofstream fichier3("1600_50.tab");
	  fichier3<<"imax 1 number of bins"<<endl;
	  fichier3<<"jmax 1 number of processes"<<endl;
	  fichier3<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"bin	ch0	"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"observation	0		"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"bin	ch0	ch0		"<<endl;
	  fichier3<<"process	sig	totbackg	"<<endl;
	  fichier3<<"process	0	1		"<<endl;
	  fichier3<<"rate	5	"<<i<<"	"<<endl;
	  fichier3<<"----------------------------------------------------"<<endl;
	  fichier3<<"uncert_bkg	 lnN 	"<<1+j<<"      	 1.2 	"<<endl;
	  fichier3.close();
	  
	  ofstream fichier4("1500_50.tab");
	  fichier4<<"imax 1 number of bins"<<endl;
	  fichier4<<"jmax 1 number of processes"<<endl;
	  fichier4<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"bin	ch0	"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"observation	0		"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"bin	ch0	ch0		"<<endl;
	  fichier4<<"process	sig	totbackg	"<<endl;
	  fichier4<<"process	0	1		"<<endl;
	  fichier4<<"rate	10	"<<i<<"	"<<endl;
	  fichier4<<"----------------------------------------------------"<<endl;
	  fichier4<<"uncert_bkg	 lnN 	"<<1+j<<"	 1.2 	"<<endl;
	  fichier4.close();
	  
	  system("sh combineCommands3.sh >> essaie3.txt");
	  system("sh combineCommands4.sh >> essaie4.txt");
	  system("sh combineCommands5.sh >> essaie5.txt");
	  ifstream fichier("essaie3.txt");
	  ifstream fichier5("essaie4.txt");
	  ifstream fichier6("essaie5.txt");
	  int x=0;
	  int y=0;
	  int z=0;
	  while( getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
	    
	    {
	      x++;
	      // cout<<"AAAA"<<endl;
	      if (x==3)
		{
		  while( getline(fichier5, ligne1))  // tant que l'on peut mettre la ligne dans "contenu"
		    
		    {
		      y++;
		      // cout<<"BBBB"<<endl;
		      if (y==3)
			{
			  while( getline(fichier6, ligne2))  // tant que l'on peut mettre la ligne dans "contenu"
			    
			    {
			      z++;
			      //cout<<"CCCC"<<endl;
			      if(x==3 && y==3 && z==3)
				{
				  size_t found = ligne.find(" ");
				  string val_neutr = ligne.substr(found+1);
				  size_t found1 = ligne1.find(" ");
				  string val_neutr1 = ligne1.substr(found1+1);
				  size_t found2 = ligne2.find(" ");
				  string val_neutr2 = ligne2.substr(found2+1);
				  //float spart_mass_point = stof(val_spart);
				  // float neutr_mass_point = stof(val_neutr);
				  fichier2<<i<<" "<<val_neutr<<" "<<val_neutr1<<" "<<val_neutr2<<" " ; // on l'affiche
				}
			    }
			}
		    }
		}
	    }
	  fichier.close();
	  fichier5.close();
	  fichier6.close();
	
	}
        fichier2<<endl;
    }
 
 //ifstream fichier1("Significance.dat");
 fichier2.close();
 
 
 
}



