#include "TMath.h"
#include "TF1.h"
#include <iostream>
#include <algorithm>

using namespace std;

//--------------------------------------//
// The goal of the macro is to compute
// a p-value for 2 sets of measurement
// taken with different luminosities
// You can define many independant regions
// (signal or control)
// The code is based on MC toys
// It only take into account Poisson distribution
// No systematics are involved
//--------------------------------------//

bool myfunction (int i,int j) { return (i<j); }

void bins_chi2_n(){

//-------------------------------------------------
//Define here the luminosities for period 1 and 2
//-------------------------------------------------
float l1 = 3.99;
float l2 = 3.66;
float ltot = l1+l2;


//-------------------------------------------------
//Define here the number of boxes (indepednant) and the yields
//-------------------------------------------------
//observed 27 - 36 | CR 88 - 54
///*
int nbins = 2;
float n1[] = {27,88};
float n2[] = {36,54};
//*/

    //observed 31 - 13 | CR 59 - 43
/*
int nbins = 2;
float n1[] = {31,59};
float n2[] = {13,43};
*/
float *ntot = new float[nbins];
float *pred1 = new float[nbins];
float *pred2 = new float[nbins];
for(int i=0;i<nbins;i++){
    ntot[i] = n1[i]+n2[i];
    pred1[i] = ntot[i]/ltot*l1;
    pred2[i] = ntot[i]/ltot*l2;
}


double prob = 1;
//max should be larger that the maximum yield
int max = 200;
for(int i=0;i<nbins;i++){
	if(max<n1[i]) max = n1[i]*2;
	if(max<n2[i]) max = n2[i]*2;

}
double chi2=0;
double like0 = 1;
for(int i=0;i<nbins;i++){
    TF1 f1("poiss1","TMath::PoissonI(x,[1])",0,max);
    TF1 f2("poiss1","TMath::PoissonI(x,[1])",0,max);
    f1.SetParameter(1,pred1[i]);
    f2.SetParameter(1,pred2[i]);
    if(f1.Integral(n1[i],max)<0.5)
        prob*=f1.Integral(n1[i],max);
    else
        prob*=f1.Integral(0,n1[i]);
    if(f2.Integral(n2[i],max)<0.5)
        prob*=f2.Integral(n2[i],max);
    else
        prob*=f2.Integral(0,n2[i]);

    //likelihood
    
    like0*=TMath::PoissonI(n1[i],pred1[i]);
    like0*=TMath::PoissonI(n2[i],pred2[i]);
    //chi2-test
    
    chi2+=pow((n1[i]-pred1[i]),2)/pred1[i];
    chi2+=pow((n2[i]-pred2[i]),2)/pred2[i];
    //prob *= f1.Integral(n1[i],max)*f2.Integral(n2[i],max);
}
cout<<"prob = "<<prob<<endl;
cout<<"likelihood = "<<like0<<endl;
cout<<"chi2 = "<<chi2<<" "<<TMath::Prob(chi2,1)<<endl;

//MC toys
TRandom3 rand;
int ntoys = 1000000;
vector<double> vprob; 
for(int i=0;i<ntoys;i++){
    double p = 1;
    double like = 1;
    for(int b=0;b<nbins;b++){
          float v1 = rand.Poisson(pred1[b]);
          float v2 = rand.Poisson(pred2[b]);
          ///Compute prob
          like*=TMath::PoissonI(v1,pred1[b]);
          like*=TMath::PoissonI(v2,pred2[b]);
             
    }
    //vprob.push_back(prob);
    //cout<<like<<endl;
    vprob.push_back(like);
}
//sort the vector
//std::sort(vprob.begin(),vprob.end(),myfunction);
//compute prob
int count = 0;
for(unsigned int x=0;x<vprob.size();x++){
    if(like0<vprob[x]) count++;
       // break;
}
cout<<count<<endl;
double prob_toy = double(count)/vprob.size();
cout<<"prob from MC toys = "<<prob_toy<<endl;


}
