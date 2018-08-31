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


//int main(int argc, char *argv[])
int main()
{

  vector<double> table;
  string ligne;
  string ligne1;
  string ligne2;

  ofstream fichier2("value6.dat");
  // fichier2<<"b,Zbi,combine"<<endl;
  

  double dat=0;
 for(double i=0.1; i<1;i+=0.1)
    {
      fichier2<<i;
      //cout<< Zbi(2,i,0.2)<<endl;
      
      for (int k=0; k<5;k++)

	{

	  if(k==0)
	    {
	      dat=1;
	    }

	  if(k==1)
	    {
	      dat=10;
	    }

	  if(k==2)
	    {
	      dat=30;
	    }

	  if(k==3)
	    {
	      dat=50;
	    }

	  if(k==4)
	    {
	      dat=100;
	    }

	  
	  ofstream fichier1("900_50.tab");
	  fichier1<<"imax 1 number of bins"<<endl;
	  fichier1<<"jmax 1 number of processes"<<endl;
	  fichier1<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"observation	"<<i+30<<"		"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	ch0		"<<endl;
	  fichier1<<"process	sig	totbackg	"<<endl;
	  fichier1<<"process	0	1		"<<endl;
	  fichier1<<"rate	"<<dat<<"	"<<i<<"	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"uncert_bkg	 lnN 	1	 1.2 	"<<endl;
	  fichier1.close();
	  
	  system("sh combineCommands.sh >> essaie.txt");
	  ifstream fichier("essaie.txt");
	  int x=0;
	  while( getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
	    
	    {
	      x++;
	      // cout<<"AAAA"<<endl;
	      if (x==3)
		{
		  size_t found = ligne.find(" ");
		  string val_neutr = ligne.substr(found+1);
		  fichier2<<" "<<val_neutr;
		} 
	    } 
	  fichier.close(); 
	}
      

      fichier2<<endl;  // on l'affiche
     
    }







 for(double i=1; i<10;i++)
    {
      fichier2<<i;
      //cout<< Zbi(2,i,0.2)<<endl;
      
      for (int k=0; k<5;k++)

	{

	  if(k==0)
	    {
	      dat=1;
	    }

	  if(k==1)
	    {
	      dat=10;
	    }

	  if(k==2)
	    {
	      dat=30;
	    }

	  if(k==3)
	    {
	      dat=50;
	    }

	  if(k==4)
	    {
	      dat=100;
	    }

	  
	  ofstream fichier1("900_50.tab");
	  fichier1<<"imax 1 number of bins"<<endl;
	  fichier1<<"jmax 1 number of processes"<<endl;
	  fichier1<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"observation	"<<i+30<<"		"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	ch0		"<<endl;
	  fichier1<<"process	sig	totbackg	"<<endl;
	  fichier1<<"process	0	1		"<<endl;
	  fichier1<<"rate	"<<dat<<"	"<<i<<"	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"uncert_bkg	 lnN 	1	 1.2 	"<<endl;
	  fichier1.close();
	  
	  system("sh combineCommands.sh >> essaie.txt");
	  ifstream fichier("essaie.txt");
	  int x=0;
	  while( getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
	    
	    {
	      x++;
	      // cout<<"AAAA"<<endl;
	      if (x==3)
		{
		  size_t found = ligne.find(" ");
		  string val_neutr = ligne.substr(found+1);
		  fichier2<<" "<<val_neutr;
		} 
	    } 
	  fichier.close(); 
	}
      

      fichier2<<endl;  // on l'affiche
     
    }
 











for(double i=10; i<100;i+=10)
    {
      fichier2<<i;
      //cout<< Zbi(2,i,0.2)<<endl;
      
      for (int k=0; k<5;k++)

	{

	  if(k==0)
	    {
	      dat=1;
	    }

	  if(k==1)
	    {
	      dat=10;
	    }

	  if(k==2)
	    {
	      dat=30;
	    }

	  if(k==3)
	    {
	      dat=50;
	    }

	  if(k==4)
	    {
	      dat=100;
	    }

	  
	  ofstream fichier1("900_50.tab");
	  fichier1<<"imax 1 number of bins"<<endl;
	  fichier1<<"jmax 1 number of processes"<<endl;
	  fichier1<<"kmax 1 number of nuisance parameters"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"observation	"<<i+30<<"		"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"bin	ch0	ch0		"<<endl;
	  fichier1<<"process	sig	totbackg	"<<endl;
	  fichier1<<"process	0	1		"<<endl;
	  fichier1<<"rate	"<<dat<<"	"<<i<<"	"<<endl;
	  fichier1<<"----------------------------------------------------"<<endl;
	  fichier1<<"uncert_bkg	 lnN 	1	 1.2 	"<<endl;
	  fichier1.close();
	  
	  system("sh combineCommands.sh >> essaie.txt");
	  ifstream fichier("essaie.txt");
	  int x=0;
	  while( getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
	    
	    {
	      x++;
	      // cout<<"AAAA"<<endl;
	      if (x==3)
		{
		  size_t found = ligne.find(" ");
		  string val_neutr = ligne.substr(found+1);
		  fichier2<<" "<<val_neutr;
		} 
	    } 
	  fichier.close(); 
	}
      

      fichier2<<endl;  // on l'affiche
     
    }





 //ifstream fichier1("Significance.dat");
 fichier2.close();
  
  
  
}



