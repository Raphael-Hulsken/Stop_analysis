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
#include "Math/TRandomEngine.h"
#include <TRandom3.h>
#include <TGraph.h>


using namespace std;
using namespace theDoctor;
using namespace TMath;

//int main(int argc, char *argv[])
int main()
{

  /* TF1 *mygaus_b_s = new TF1("mygaus_b_s","TMath::Gaus(x,8,.5)",-100,100);
  mygaus_b_s->Draw();

  TF1 *mygaus_b = new TF1("mygaus_b","TMath::Gaus(x,3,.5)",-100,100);
  mygaus_b->Draw();
  */
  double mu1=14;
  double mu2=10;
  double sig1=1.8;
  double sig2=1.8;
  double value_dat=13;
  double value_dat_CLS;
  int nb_event=1000000;
  int nb_bin=10000;

   TH1* histo_b = new TH1I("histo_b", "histo_b_s", nb_bin, 0, 300);

   TH1* histo_b_s = new TH1I("histo_b_s", "histo_b_s", nb_bin,0, 300);
  
   TH1* histo_CLS_b = new TH1I("histo_CLS_b", "histo_CLS_b", nb_bin, -20, 20);

   TH1* histo_CLS_b_s = new TH1I("histo_CLS_b_s", "histo_CLS_b_s", nb_bin,-20, 20);
  

 TRandom* z= new TRandom();


 TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,1400,1000); 


   

   double c;

 for (int i=0; i<nb_event; i++)
      {
	c=z->Gaus(mu1,sig1);
	histo_b_s->Fill(c);
	histo_CLS_b_s->Fill(-2 *log(Gaus(c,mu1,sig1)/Gaus(c,mu2,sig2)));


	c=z->Gaus(mu2,sig2);
	histo_b->Fill(c);
	histo_CLS_b->Fill(-2 *log(Gaus(c,mu1,sig1)/Gaus(c,mu2,sig2)));
      }



 histo_b->Draw();
 histo_b_s->Draw("same");
 histo_b_s->SetLineColor(kRed);
 c1->SaveAs("plotTest1.root");


 histo_CLS_b->Draw();
 histo_CLS_b_s->Draw("same");
 histo_CLS_b_s->SetLineColor(kRed);
 c1->SaveAs("plotTest2.root");

 value_dat_CLS=-2 * log(Gaus(value_dat,mu1,sig1)/Gaus(value_dat,mu2,sig2));

 Int_t min_b_s =histo_b_s->FindBin(value_dat);
 Int_t min_b =histo_b->FindBin(value_dat);
 Int_t min_CLS_b_s =histo_CLS_b_s->FindBin(value_dat_CLS);
 Int_t min_CLS_b =histo_CLS_b->FindBin(value_dat_CLS);
 
 cout<<value_dat_CLS<<endl;
 cout<<histo_b->Integral(0,min_b)/histo_b->Integral(0, nb_bin)<<endl;
cout<<histo_b_s->Integral(min_b_s, nb_bin)/histo_b_s->Integral(0, nb_bin)<<endl;
cout<<histo_CLS_b->Integral(0,min_CLS_b)/histo_CLS_b->Integral(0, nb_bin)<<endl;
cout<<histo_CLS_b_s->Integral(min_CLS_b_s, nb_bin)/histo_CLS_b_s->Integral(0, nb_bin)<<endl;


}

   /*

  
   // TGraph* gr_b = new TGraph(nb_event);
   // TGraph* gr_b_s = new TGraph(nb_event);

   vector<double> b;
   double c;
   double Q_b_s;
   double Q_b;
   double k=0;
   vector<double>  b_s;


	/*	histo_b_s->Fill(Q_b_s);
	histo_b->Fill(Q_b);
	
 //cout<<x->value()<<endl;
   }









 sort(b_s.begin() , b_s.end());
 sort(b.begin() , b.end());

 for (int i=0; i<nb_event; i++)
      {
	Q_b_s=-2 *log(Gaus(b_s[nb_event-i],mu1,sig1)/Gaus(b_s[nb_event-i],mu2,sig2));	
	gr_b_s->SetPoint(i,Q_b_s,Gaus(b_s[nb_event-i],mu1,sig1));
 	Q_b=-2 *log(Gaus(b[nb_event-i],mu1,sig1)/Gaus(b[nb_event-i],mu2,sig2));
	gr_b->SetPoint(i,Q_b,Gaus(b[nb_event-i],mu2,sig2));

   }




 gr_b_s->Draw();
 gr_b->Draw("same");
 c1->Update();
 c1->SaveAs("plotTest1.root");



 double vx;
 double vy;
 double ca=0;





vx=-1000;
 
 int i=1;
 while(value_dat_CLS > vx)
   { 
     gr_b_s->GetPoint(i,vx,vy);
     i++;
   }
 ca=i;
 //double ca2=0;
 cout<<value_dat_CLS<<endl;
 cout<<ca<<endl;
  gr_b_s->GetPoint(ca+1,vx,vy);
 cout<<vx<<endl;
 /*gr_b_s->GetPoint(ca,vx,vy);
 cout<<vx<<endl;
 gr_b_s->GetPoint(ca+10,vx,vy);
 cout<<vx<<endl;
 double tot_b_s;
 tot_b_s= gr_b_s->Integral(1,nb_event-1);
 Inte_b_s= gr_b_s->Integral(ca,nb_event-1);
 
cout<<Inte_b_s/tot_b_s<<endl;

vx=-1000;
 i=1;
 while(value_dat_CLS > vx)
   { 
     gr_b->GetPoint(i,vx,vy);
     i++;
   }

 ca=i;
 
 // cout<<value_dat_CLS<<endl;
 cout<<ca<<endl;
 /* gr_b_s->GetPoint(ca-10,vx,vy);
    cout<<vx<<endl;
 gr_b->GetPoint(ca+1,vx,vy);
 cout<<vx<<endl;
 gr_b_s->GetPoint(ca+10,vx,vy);
 cout<<vx<<endl;
 double tot_b;
 tot_b= gr_b->Integral(1,nb_event-1);
 Inte_b= gr_b->Integral(ca,nb_event-1);

 cout<<Inte_b/tot_b<<endl;
*/

