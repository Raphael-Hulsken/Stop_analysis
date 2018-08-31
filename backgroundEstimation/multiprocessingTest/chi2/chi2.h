#ifndef Chi2_h
#define Chi2_h

#include <vector>
#include <exception>
#include <cstring>
#include <algorithm>
#include "TFitter.h"
#include "TLorentzVector.h"
#include "Math/LorentzVector.h"

using namespace std;

static const float PDG_TOP_MASS = 173.5;
static const float PDG_W_MASS = 80.385;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
;
double fc2 (double c1, double m12, double m22, double m02, bool verbose = false);
double fchi2 (void *num, void* denom, uint32_t size);
//void minuitFunction(int&, double* , double &result, double par[], int);
void minuitFunction(int&, double* , double &result, double par[], int);
double calculateChi2(vector<LorentzVector>& jets, vector<float>& sigma_jets, vector<bool>& btag);
vector<double> calculateChi2ForWMass(vector<TLorentzVector>& jets, vector<float>& sigma_jets,int njets_max, bool new_formula);
vector<double> calculateChi2ForWdR(vector<TLorentzVector>& jets, vector<float>& sigma_jets,int njets_max, bool new_formula);
vector< vector<double> > calculateChi2ForTopMass(vector<TLorentzVector>& jets, vector<float>& sigma_jets, int njets_max);
vector< vector<double> > calculateChi2ForTopdR(vector<TLorentzVector>& jets, vector<float>& sigma_jets, int njets_max);

double Chi2(vector<float> jets_pt, vector<float> jets_eta, vector<float> jets_phi, vector<float> jets_m,int method, int njets_max, bool new_formula);

//double Chi2Interface(int njets, double* jets_pt, double* jets_eta, double* jets_phi, double* jets_E, double* sigma, bool* jets_btag){
//
//
//extern "C" double Chi2Interface(int njets, double* jets_pt, double* jets_eta, double* jets_phi, double* jets_E, bool* jets_btag);

#endif
