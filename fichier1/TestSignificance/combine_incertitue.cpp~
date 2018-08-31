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
  string ligne;
  ofstream fichier2("value5.dat");
  double x;
  //std::string str;
  double k;
  //std::ostringstream strs;
  int nbr;
 
  for(double i=0.1; i<1;i+=0.1)
    {
      fichier2<<i;
      //cout<< Zbi(2,i,0.2)<<endl;
      for (double j =0;j<5;j++)
	{
	  cout<<j<<endl;
	  if (j==0)
	    {
	      x=1.0001;
	      k=1.0001;
	    }
	  
	  if (j==1)
	    {
	      x=1.1;
	      k=1.1;
	    }
	  
	  if (j==2)
	    {
	      x=0.9;
	      k=0.9;
	    }


	  if (j==3)
	    {
	      x=1.1;
	      k=0.9;
	    }
	  
	  if (j==4)
	    {
	      x=0.9;
	      k=1.1;
	    }
	  
	  
	  for (int pff=0 ;pff<3 ;pff++)
	    
	    {
	      if (pff==0)
		{
		  nbr=5;
		}
	      
	      if (pff==1)
		{
		  nbr=10;
		}
	      
	      if (pff==2)
		{
		  nbr=15;
		}
	      
	      
	      
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
		  fichier1<<"rate	"<<nbr<<"	"<<i<<"	"<<endl;
		  fichier1<<"----------------------------------------------------"<<endl;
		  fichier1<<"uncert_bkg   lnN   "<<x<<"    "<<k<<" 	"<<endl;
		  fichier1.close();
		
	      
	      

	      
	      system("sh combineCommands3.sh >> essaie3.txt");
	      
	      ifstream fichier("essaie3.txt");
	      int y =0;
	      while (getline(fichier, ligne))
		{
		  y++;
		  if (y==3)
		    {
		      size_t found = ligne.find(" ");
		      string val_neutr = ligne.substr(found+1);
		      fichier2<<" "<<val_neutr;
		    }
		}
	    
	      fichier.close();
	    }
	}
      fichier2<<endl;
    }



  for(double i=1; i<10;i++)
    {
      fichier2<<i;
      //cout<< Zbi(2,i,0.2)<<endl;
      for (double j =0;j<5;j++)
	{
	  cout<<j<<endl;
	  if (j==0)
	    {
	      x=1.0001;
	      k=1.0001;
	    }
	  
	  if (j==1)
	    {
	      x=1.1;
	      k=1.1;
	    }
	  
	  if (j==2)
	    {
	      x=0.9;
	      k=0.9;
	    }


	  if (j==3)
	    {
	      x=1.1;
	      k=0.9;
	    }
	  
	  if (j==4)
	    {
	      x=0.9;
	      k=1.1;
	    }
	  
	  
	  for (int pff=0 ;pff<3 ;pff++)
	    
	    {
	      if (pff==0)
		{
		  nbr=5;
		}
	      
	      if (pff==1)
		{
		  nbr=10;
		}
	      
	      if (pff==2)
		{
		  nbr=15;
		}
	      
	      
	      
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
		  fichier1<<"rate	"<<nbr<<"	"<<i<<"	"<<endl;
		  fichier1<<"----------------------------------------------------"<<endl;
		  fichier1<<"uncert_bkg   lnN   "<<x<<"    "<<k<<" 	"<<endl;
		  fichier1.close();
	       
	      
	      

	      
	      system("sh combineCommands3.sh >> essaie3.txt");
	      
	      ifstream fichier("essaie3.txt");
	      int y =0;
	      while (getline(fichier, ligne))
		{
		  y++;
		  if (y==3)
		    {
		      size_t found = ligne.find(" ");
		      string val_neutr = ligne.substr(found+1);
		      fichier2<<" "<<val_neutr;
		    }
		}
	    
	      fichier.close();
	    }
	}
      fichier2<<endl;
    }



  for(double i=10; i<100;i+=10)
    {
      fichier2<<i;
      //cout<< Zbi(2,i,0.2)<<endl;
      for (double j =0;j<5;j++)
	{
	  cout<<j<<endl;
	  if (j==0)
	    {
	      x=1.0001;
	      k=1.0001;
	    }
	  
	  if (j==1)
	    {
	      x=1.1;
	      k=1.1;
	    }
	  
	  if (j==2)
	    {
	      x=0.9;
	      k=0.9;
	    }


	  if (j==3)
	    {
	      x=1.1;
	      k=0.9;
	    }
	  
	  if (j==4)
	    {
	      x=0.9;
	      k=1.1;
	    }
	  
	  
	  for (int pff=0 ;pff<3 ;pff++)
	    
	    {
	      if (pff==0)
		{
		  nbr=5;
		}
	      
	      if (pff==1)
		{
		  nbr=10;
		}
	      
	      if (pff==2)
		{
		  nbr=15;
		}
	      
	      
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
		  fichier1<<"rate	"<<nbr<<"	"<<i<<"	"<<endl;
		  fichier1<<"----------------------------------------------------"<<endl;
		  fichier1<<"uncert_bkg   lnN     "<<x<<"    "<<k<<"	"<<endl;
		  fichier1.close();
		
	      
	      

	      
	      system("sh combineCommands3.sh >> essaie3.txt");
	      
	      ifstream fichier("essaie3.txt");
	      int y =0;
	      while (getline(fichier, ligne))
		{
		  y++;
		  if (y==3)
		    {
		      size_t found = ligne.find(" ");
		      string val_neutr = ligne.substr(found+1);
		      fichier2<<" "<<val_neutr;
		    }
		}
	    
	      fichier.close();
	    }
	}
      fichier2<<endl;
    }

  fichier2.close();
}
