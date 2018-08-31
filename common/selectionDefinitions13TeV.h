#define MT_CUT    150
#define MET_CUT   250
#define MT_PSCUT    0
#define MET_PSCUT   50
#define NJET_CUT 4
#define NBJET_CUT 1
#define NLEP_CUT  1
#define DPHI_CUT 0.8
#define MET_CUTLL 50
#define MTW2_CUT 200
#define MLB_CUT 175
#define MET_BOUND_250 250
#define MET_BOUND_300 300
#define MET_BOUND_325 325
#define MET_BOUND_350 350
#define MET_BOUND_375 375
#define MET_BOUND_400 400
#define MET_BOUND_450 450
#define MET_BOUND_500 500

#define MET_CR 50

// Not sure that it is a good idea to include this here,
// since one often needs to use a modified format because
// of skimming or use of tiny tuples...

//#include "Reader.h"
//#include "Reader_final.h"  // Has the extended BDT info defined
#include "Reader_CommonFormat.h"  // Has the extended BDT info defined
#include "CrossSection.h" // used to compute the weight

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctgmath>

using namespace std;

// NB : When you call any of the following functions,
// these three variables must be filled with the current
// event being taken care of, the name of the sample and
// the sample type (background, signal or data)

babyEvent myEvent;
string sampleName;
string sampleType;



// MT cuts definitions
// ###################

bool goesInMTpeak()     { if ((myEvent.mt_met_lep > 30) && (myEvent.mt_met_lep < 80)) return true; else return false; }
bool goesInMTtail()     { if (myEvent.mt_met_lep > MT_CUT)                    return true; else return false; }
bool goesInMTinverted() { if (myEvent.mt_met_lep < MT_CUT)                    return true; else return false; }

// Control region definitions
// ##########################

bool goesInPreVetoSelection()
{

    if (myEvent.pfmet < MET_PSCUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;

    return true;
}

bool goesInPreVetoSelectionMTtail()     { return (goesInPreVetoSelection() && goesInMTtail());     }
bool goesInPreVetoSelectionMTpeak()     { return (goesInPreVetoSelection() && goesInMTpeak());     }
bool goesInPreVetoSelectionMTinverted() { return (goesInPreVetoSelection() && goesInMTinverted()); }


bool goesInPreselectionNoVeto()
{
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;

    return true;
}

bool goesInPreselectionNoVetoNoMetCut()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;

    return true;
}


bool goesInPreselection()
{
    
    if (myEvent.pfmet < MET_PSCUT) return false;
    //cout<<"1"<<endl;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    //cout<<"2"<<endl;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    //cout<<"3"<<endl;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    //cout<<"4"<<endl;
    
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    //cout<<"5"<<endl;

    return true;
}

bool goesInBaselineSearchSR() {
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < DPHI_CUT) return false; 
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return (goesInPreselection() && goesInMTtail() );
}
bool goesInBaselineSearchSR2b() { return (goesInPreselection() && goesInMTtail() && myEvent.ngoodbtags >=2 );}
//bool goesInLargeDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.MT2W > 200 && myEvent.dphi_Wlep> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}
//bool goesInSmallDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.dphi_Wlep> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}

/*
bool goesInLargeDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.MT2W > 200 && myEvent.minDPhi_jmet> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}
bool goesInSmallDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.minDPhi_jmet> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}
bool goesInSmallDMSR300() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.minDPhi_jmet> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 300);}
bool goesInSmallDMSR2b() { return (goesInSmallDMSR() && myEvent.ngoodbtags>=2);}
*/
//background

bool stBackground()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 2)  return false;
    if (myEvent.ngoodjets != 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET50()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1 && myEvent.ngoodjets != 2 )  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET80()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1 && myEvent.ngoodjets != 2 )  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 80) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501j()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt20()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 20)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt60()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 60)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET502jPt20()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 20)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET502jPt60()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 60)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt20b1()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if (myEvent.pfmet > 250) return false;
    if(im>70 && im<120) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 20)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt20b2()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet <= 250) return false;
    if (myEvent.pfmet > 350) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 20)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt20b3()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet <= 350) return false;
    if (myEvent.pfmet > 450) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 20)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt20b4()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet <= 450) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 20)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt40()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 40)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET502jPt40()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 40)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt40b1()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if (myEvent.pfmet > 250) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 40)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt40b2()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet <= 250) return false;
    if (myEvent.pfmet > 350) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 40)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt40b3()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet <= 350) return false;
    if (myEvent.pfmet > 450) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 40)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET501jPt40b4()
{

    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet <= 450) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;

            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            if(sum.Pt() > 40)
               return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET801j()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 1)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 80) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET502j()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool stBackgroundMET802j()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    if (myEvent.ngoodleps != 2) return false;
    if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets != 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 80) return false;
    if(im>70 && im<100) return false;
    if(im<20) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}


bool stBackgroundTest()
{
    float im = 2*myEvent.lep1_pt*myEvent.lep2_pt*(cosh(myEvent.lep1_eta - myEvent.lep2_eta) - cos(myEvent.lep1_phi - myEvent.lep2_phi));
    //if (myEvent.ngoodleps != 2) return false;
    //if (myEvent.ngoodbtags != 1)  return false;
    if (myEvent.ngoodjets < 2)  return false;
    //if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < 50) return false;
    //if(im>60 && im<120) return false;
    
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}




bool CR0b_presel_nojetreq()
{
    if (myEvent.pfmet < MET_CR) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodbtags > 0 )  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}


bool CR0b_presel()
{
    if (myEvent.pfmet < MET_CR) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    //if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodjets < 3)  return false;
    //WARNING: REMOVED TEMPORARILY
    //if (myEvent.ngoodbtags > 0 )  return false;
    //WARNING: ADDED TEMPORARILY
    if (myEvent.mt_met_lep > 100) return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool CR_MET250_lowMT(){
    if (myEvent.pfmet < 250) return false;
    if (myEvent.ngoodleps != 1) return false;
    if (myEvent.ngoodjets < 4)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    if (myEvent.mt_met_lep<100) return false;
    return true;
}
bool CR_MET250_lowMT_3j(){
    if (myEvent.pfmet < 250) return false;
    if (myEvent.ngoodleps != 1) return false;
    if (myEvent.ngoodjets < 3)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    if (myEvent.mt_met_lep<100) return false;
    return true;
}

bool CR0b_MET250_lowMT(){
	return (CR_MET250_lowMT() && myEvent.ngoodbtags == 0 );
}
bool CR0b_MET250_lowMT_3j(){
	return (CR_MET250_lowMT_3j() && myEvent.ngoodbtags == 0 );
}


bool CR0b_presel_MET50(){ return (CR0b_presel() && myEvent.pfmet>50 && myEvent.pfmet<100); }
bool CR0b_presel_MET100(){ return (CR0b_presel() && myEvent.pfmet>100 && myEvent.pfmet<150); }
bool CR0b_presel_MET150(){ return (CR0b_presel() && myEvent.pfmet>150 && myEvent.pfmet<200); }
bool CR0b_presel_MET200(){ return (CR0b_presel() && myEvent.pfmet>200 && myEvent.pfmet<250); }
bool CR0b_presel_MET250(){ return (CR0b_presel() && myEvent.pfmet>250 && myEvent.pfmet<300); }
bool CR0b_presel_MET300(){ return (CR0b_presel() && myEvent.pfmet>300); }

bool CR0b_presel_MT2Wtail()
{
	return (CR0b_presel() && myEvent.MT2W>200);
}

bool CR0b_presel_MT2W200(){return (CR0b_presel() && myEvent.MT2W>200 && myEvent.MT2W<250);}
bool CR0b_presel_MT2W250(){return (CR0b_presel() && myEvent.MT2W>250);}


bool CR0b_presel_MTpeak() { return CR0b_presel() && (myEvent.mt_met_lep> 50 && myEvent.mt_met_lep < 80) ; }
bool CR0b_presel_MTtail() { return CR0b_presel() && (myEvent.mt_met_lep > 100) ; }
bool CR0b_presel_MTtail_80() { return CR0b_presel() && (myEvent.mt_met_lep > 80) ; }
bool CR0b_presel_MTtail_90() { return CR0b_presel() && (myEvent.mt_met_lep > 90) ; }
bool CR0b_presel_MTtail_100() { return CR0b_presel() && (myEvent.mt_met_lep > 100) ; }
bool CR0b_presel_MTtail_110() { return CR0b_presel() && (myEvent.mt_met_lep > 110) ; }
bool CR0b_presel_MTtail_120() { return CR0b_presel() && (myEvent.mt_met_lep > 120) ; }
bool CR0b_presel_MTtail_130() { return CR0b_presel() && (myEvent.mt_met_lep > 120) ; }
bool CR0b_presel_MTtail_140() { return CR0b_presel() && (myEvent.mt_met_lep > 130) ; }
bool CR0b_presel_MTtail_150() { return CR0b_presel() && (myEvent.mt_met_lep > 140) ; }

bool CR0b_presel_MTtail_80_ex() { return CR0b_presel() && (myEvent.mt_met_lep > 80 && myEvent.mt_met_lep < 90) ; }
bool CR0b_presel_MTtail_90_ex() { return CR0b_presel() && (myEvent.mt_met_lep > 90 && myEvent.mt_met_lep < 100) ; }
bool CR0b_presel_MTtail_100_ex() { return CR0b_presel() && (myEvent.mt_met_lep > 100 &&  myEvent.mt_met_lep < 110) ; }
bool CR0b_presel_MTtail_110_ex() { return CR0b_presel() && (myEvent.mt_met_lep > 110 &&  myEvent.mt_met_lep < 120) ; }
bool CR0b_presel_MTtail_120_ex() { return CR0b_presel() && (myEvent.mt_met_lep > 120 &&  myEvent.mt_met_lep < 130) ; }
bool CR0b_presel_MTtail_130_ex() { return CR0b_presel() && (myEvent.mt_met_lep > 120 &&  myEvent.mt_met_lep < 140) ; }

bool CR0b_presel_2j_MTpeak() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 2 && (myEvent.mt_met_lep> 50 && myEvent.mt_met_lep < 80) ; }
bool CR0b_presel_2j_MTtail() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 2 && (myEvent.mt_met_lep > 100) ; }
bool CR0b_presel_3j_MTpeak() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 3 && (myEvent.mt_met_lep> 50 && myEvent.mt_met_lep < 80) ; }
bool CR0b_presel_3j_MTtail() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 3 && (myEvent.mt_met_lep > 100) ; }
bool CR0b_presel_4j_MTpeak() { return CR0b_presel_nojetreq() && myEvent.ngoodjets >= 4 && (myEvent.mt_met_lep> 50 && myEvent.mt_met_lep < 80) ; }
bool CR0b_presel_4j_MTtail() { return CR0b_presel_nojetreq() && myEvent.ngoodjets >= 4 && (myEvent.mt_met_lep > 100) ; }


bool CR0b_presel_2j() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 2  ; }
bool CR0b_presel_3j() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 3 ; }
bool CR0b_presel_4j() { return CR0b_presel_nojetreq() && myEvent.ngoodjets >= 4 ; }

bool CR0b_presel_2j_MET100() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 2  && myEvent.pfmet > 100; }
bool CR0b_presel_3j_MET100() { return CR0b_presel_nojetreq() && myEvent.ngoodjets == 3  && myEvent.pfmet > 100; }
bool CR0b_presel_4j_MET100() { return CR0b_presel_nojetreq() && myEvent.ngoodjets >= 4  && myEvent.pfmet > 100; }

// Categorization of Events
// ########################
bool boostedSearch()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    //if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < MET_CUT) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearch()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchSpec3()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchSpec4()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < 4)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchSpec4250()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < 4)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}
bool defaultSearchMETCut()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchMTCut()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchMT2WCut()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsExcMET()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsExcMETLooseMT()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_PSCUT) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsExcMT()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}


bool defaultSearchAllCuts()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}


bool defaultSearchAllCutsMS()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.met_sig < 88) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsMSToMET()
{
    Double_t met_sigTomet = myEvent.met_sig / myEvent.pfmet;
   
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (met_sigTomet < 0.35) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsLooseMTMS()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.mt_met_lep < MT_PSCUT) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.met_sig < 144) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsMTS()
{
    Double_t deltaPhi = myEvent.lep1_phi - myEvent.pfmet_phi;
    Double_t MT_sig = sqrt(2 * myEvent.lep1_pt * myEvent.met_sig * (1 - cos(deltaPhi) ));

    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (MT_sig < 96) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsMTSToMT()
{
    Double_t deltaPhi = myEvent.lep1_phi - myEvent.pfmet_phi;
    Double_t MT_sig = sqrt(2 * myEvent.lep1_pt * myEvent.met_sig * (1 - cos(deltaPhi) ));
    Double_t MT_sigToMT = MT_sig / myEvent.mt_met_lep;

    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (MT_sigToMT < 0.75) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsZ()
{
    Double_t met_sigTomet = myEvent.met_sig / myEvent.pfmet;
    Double_t Z = (myEvent.mt_met_lep - 80.385)/(myEvent.mt_met_lep*met_sigTomet);    

    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (Z < 1.2) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsLooser()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W < 150)  return false;
    if (myEvent.mt_met_lep < 100) return false;
    if (myEvent.pfmet < 200) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool defaultSearchAllCutsFakes()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin1()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex > 5)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false; 
    return true;
}

bool preselectionPUBin2()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 5)  return false;
    if (myEvent.nvertex > 10)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin3()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 10)  return false;
    if (myEvent.nvertex > 15)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin4()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 15)  return false;
    if (myEvent.nvertex > 20)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin5()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 20)  return false;
    if (myEvent.nvertex > 25)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin6()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 25)  return false;
    if (myEvent.nvertex > 30)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin7()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 30)  return false;
    if (myEvent.nvertex > 35)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool preselectionPUBin8()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.nvertex <= 35)  return false;
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;
    return true;
}

bool noCuts()
{
   return true;
}

bool defaultSearchAllCutsLowMT2W()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool lowMT2WboostedSearch()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.mt_met_lep < MT_CUT) return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.dphi_ak4pfjets_met < 0.8) return false;
    if (myEvent.pfmet < MET_CUT) return false;
     
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool lowMT2WDefaultBin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    return lowMT2WboostedSearch();
}

bool lowMT2WDefaultBin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return lowMT2WboostedSearch();
}

bool lowMT2WDefaultBin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    return lowMT2WboostedSearch();
}

bool lowMT2WDefaultBin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    return lowMT2WboostedSearch();
}

bool lowMT2WDefaultBin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    return lowMT2WboostedSearch();
}

bool DefaultBin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_325) return false;

    return boostedSearch();
}

bool DefaultBin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_325) return false;

    return boostedSearch();
}

bool DefaultBin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return boostedSearch();
}

bool DefaultBin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_450) return false;

    return boostedSearch();
}

bool DefaultBin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_450) return false;

    return boostedSearch();
}

bool DefaultBin6()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return boostedSearch();
}

bool DefaultBin7()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    return boostedSearch();
}


bool Default2Bin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    return boostedSearch();
}

bool Default2Bin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    return boostedSearch();
}

bool Default2Bin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    return boostedSearch();
}

bool Default3Bin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    return boostedSearch();
}

bool Default3Bin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    return boostedSearch();
}

bool Default3Bin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;

    return boostedSearch();
}

bool NoAk8JetsBin1()
{
    //if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool NoAk8JetsBin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool NoAk8JetsBin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool NoAk8JetsBin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool NoAk8JetsBin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool NoAk8()
{
    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_corrpruned_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_corrpruned_mass.at(idx) > 60 && myEvent.ak8pfjets_corrpruned_mass.at(idx) < 100 && myEvent.ak8pfjets_pt.at(idx) > 250 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
             return false;
    } 
    return true;
}

bool Ak8(uint8_t nr)
{
    uint32_t ak8ts = 0;
    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_corrpruned_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_corrpruned_mass.at(idx) > 60 && myEvent.ak8pfjets_corrpruned_mass.at(idx) < 100 && myEvent.ak8pfjets_pt.at(idx) > 250 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
         {
             ak8ts++;
         }
         if(idx == (myEvent.ak8pfjets_corrpruned_mass.size() -1))
         {
             if(ak8ts == 0)
                 return false;
             else
                 break;
         }
    }    

    if (nr==1 && ak8ts<nr) 
        return false;
    else if(nr >1 && ak8ts<nr)
        return false;

    return true;
}


bool lowDmNoAk8JetsBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;


    return (boostedSearch() && NoAk8());
}

bool lowDmOnePlusAk8JetBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_corrpruned_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;

    return (boostedSearch() && Ak8(1));
}

bool NoAk8JetsBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    
    return (boostedSearch() && NoAk8());
}

bool OnePlusAk8JetBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_corrpruned_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    return (boostedSearch() && Ak8(1));
}

bool threeJetsNoAk8JetsBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    
    return (boostedSearch() && NoAk8());
}

bool threeJetsOnePlusAk8JetBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_corrpruned_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;

    return (boostedSearch() && Ak8(1));
}


bool doLeptonAndAk10NotOverlap( float jet_pt, float jet_eta, float jet_phi, float jet_mass)
{

        TLorentzVector l1;
        TLorentzVector j1;


        l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
        j1.SetPtEtaPhiM(jet_pt, jet_eta, jet_phi, jet_mass);

        Double_t dR = l1.DeltaR(j1);

        bool result = dR > 1.0 ? true : false;

        return result;

}


bool NoAk10()
{
    for(uint32_t idx = 0; idx < myEvent.ak10pfjets_pruned_mass.size(); idx++)
    {
         bool noLeptonOverlap = doLeptonAndAk10NotOverlap( myEvent.ak10pfjets_pt.at(idx), myEvent.ak10pfjets_eta.at(idx), myEvent.ak10pfjets_phi.at(idx), myEvent.ak10pfjets_pruned_mass.at(idx));
         if(myEvent.ak10pfjets_pruned_mass.at(idx) > 60 && myEvent.ak10pfjets_pruned_mass.at(idx) < 100 && myEvent.ak10pfjets_pt.at(idx) > 200 && (myEvent.ak10pfjets_tau2.at(idx) / myEvent.ak10pfjets_tau1.at(idx)) < 0.5 && noLeptonOverlap)
             return false;
    } 
    return true;
}

bool Ak10(uint8_t nr)
{
    uint32_t ak10ts = 0;
    for(uint32_t idx = 0; idx < myEvent.ak10pfjets_pruned_mass.size(); idx++)
    {
         bool noLeptonOverlap = doLeptonAndAk10NotOverlap( myEvent.ak10pfjets_pt.at(idx), myEvent.ak10pfjets_eta.at(idx), myEvent.ak10pfjets_phi.at(idx), myEvent.ak10pfjets_pruned_mass.at(idx));
         if(myEvent.ak10pfjets_pruned_mass.at(idx) > 60 && myEvent.ak10pfjets_pruned_mass.at(idx) < 100 && myEvent.ak10pfjets_pt.at(idx) > 200 && (myEvent.ak10pfjets_tau2.at(idx) / myEvent.ak10pfjets_tau1.at(idx)) < 0.5 && noLeptonOverlap)
         {
             ak10ts++;
         }
         if(idx == (myEvent.ak10pfjets_pruned_mass.size() -1))
         {
             if(ak10ts == 0)
                 return false;
             else
                 break;
         }
    }    

    if (nr==1 && ak10ts<nr) 
        return false;
    else if(nr>1 && ak10ts<nr)
        return false;

    return true;
}

bool lowDmNoAk10JetsBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;


    return (boostedSearch() && NoAk10());
}

bool lowDmOnePlusAk10JetBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.ak10pfjets_pruned_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;

    return (boostedSearch() && Ak10(1));
}


bool NoAk10JetsBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    return (boostedSearch() && NoAk10());
}

bool OnePlusAk10JetBin()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ak10pfjets_pruned_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    
    return (boostedSearch() && Ak10(1));
}


bool threeJetsNoAk10JetsBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    
    return (boostedSearch() && NoAk10());
}

bool threeJetsOnePlusAk10JetBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ak10pfjets_pruned_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;

    return (boostedSearch() && Ak10(1));
}


bool AtLeastOneAk8JetBin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool AtLeastOneAk8JetBin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool AtLeastOneAk8JetBin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool AtLeastOneAk8JetBin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool AtLeastOneAk8JetBin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool AtLeastOneAk8JetNoCutBin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;
 
    return boostedSearch();
}

bool AtLeastOneAk8JetNoCutBin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return boostedSearch();
}

bool AtLeastOneAk8JetNoCutBin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    return boostedSearch();
}

bool AtLeastOneAk8JetNoCutBin4()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    return boostedSearch();
}

bool AtLeastOneAk8JetNoCutBin5()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    return boostedSearch();
}

bool NoAk10Jets()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak10pfjets_mass.size() != 0) return false; //@MJ@ TODO uncomment this one day

    return boostedSearch();
}

bool AtLeastOneAk10Jet()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak10pfjets_mass.size() == 0) return false; //@MJ@ TODO uncomment this one day

    return boostedSearch();
}

bool ThreeJetsDefaultBin1()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    return boostedSearch();
}

bool ThreeJetsDefaultBin2()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return boostedSearch();
}

bool ThreeJetsDefaultBin3()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    return boostedSearch();
}

bool ThreeJetsDefaultBin4()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    return boostedSearch();
}

bool ThreeJetsDefaultBin5()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    return boostedSearch();
}

bool ThreeJetsDefaultMiasBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    return boostedSearch();
}

bool lowMT2WThreeJetsDefaultMiasBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    return lowMT2WboostedSearch();
}

bool ThreeJetsInEventNoAk8Bin1()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool ThreeJetsInEventNoAk8Bin2()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool ThreeJetsInEventNoAk8Bin3()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool ThreeJetsInEventNoAk8Bin4()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 

    return boostedSearch();
}

bool ThreeJetsInEventNoAk8Bin5()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 
 
    return boostedSearch();
}

bool ThreeJetsInEventNoAk8MiasBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 
 
    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk8Bin1()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_300) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk8Bin2()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_300) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk8Bin3()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_400) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk8Bin4()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_400) return false;
    if (myEvent.pfmet >= MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk8Bin5()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_500) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk8MiasBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             break;
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool ThreeJetsInEventOneOfThemAk10()
{
    if (myEvent.ngoodjets < 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak10pfjets_mass.size() == 0) return false; //@MJ@ TODO uncomment this one day

    return boostedSearch();
}

bool Default2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;

    return boostedSearch();
}

bool oDefault2Bin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return boostedSearch();
}

bool oDefault2Bin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_450) return false;
    
    return boostedSearch();
}

bool oDefault2Bin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_450) return false;
    
    return boostedSearch();
}

bool Default2OneAk8Bin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
         {
            if ((myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
                 break;
         }
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    } 

    return boostedSearch();
}

bool Default2OneAk8Bin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_450) return false;
 
    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
         {
            if ((myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
                 break;
         }
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    } 
   
    return boostedSearch();
}

bool Default2OneAk8Bin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_450) return false;
 
    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
         {
            if ((myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
                 break;
         }
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    } 
   
    return boostedSearch();
}

bool Default2NoAk8Bin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
             return false;
    } 

    return boostedSearch();
}

bool Default2NoAk8Bin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;
    if (myEvent.pfmet >= MET_BOUND_450) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
             return false;
    } 
    
    return boostedSearch();
}

bool Default2NoAk8Bin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_450) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
             return false;
    } 
    
    return boostedSearch();
}

bool NewSR0ak8Bin2WOCuts()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             return false;
    } 
    return boostedSearch();
}

bool NewSR1andMoreak8Bin3WOCuts()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;   

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
         {
                 break;
         }
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    } 
    return boostedSearch();
}

bool NewSRDefaultForLowMETBin1()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250) return false;
    if (myEvent.pfmet >= MET_BOUND_350) return false;

    return boostedSearch();
}

bool NewSR0ak8Bin2()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
             return false;
    } 
    return boostedSearch();
}

bool NewSR1andMoreak8Bin3()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;   

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
         {
            if ((myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
                 break;
         }
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    } 
    return boostedSearch();
}

bool NewThreeJetsInEventOneOfThemAk8MiasBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
         {
            if ((myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
                 break;
         }
         if(idx == (myEvent.ak8pfjets_mass.size() -1))
             return false;
    }    

    return boostedSearch();
}

bool NewThreeJetsInEventNoAk8MiasBin()
{
    if (myEvent.ngoodjets != 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak8pfjets_mass.size() != 0) return false;
    if (myEvent.pfmet < MET_BOUND_350) return false;

    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200 && (myEvent.ak8pfjets_tau2.at(idx) / myEvent.ak8pfjets_tau1.at(idx)) < 0.5)
             return false;
    } 
 
    return boostedSearch();
}

// New signal regions
// ###################
/*
bool SRbin1(ofstream &outputFile)
{
    static uint32_t j;
    j++;

    if(j == 1)
    {
        outputFile << "Signal region bin 1" << endl;
        outputFile << "----------------------------" << endl;
        outputFile << "myEvent.ngoodleps == " << NLEP_CUT << endl;
        outputFile << "myEvent.ngoodjets > " << NJET_CUT << endl;
        outputFile << "myEvent.ngoodbtags > " << NBJET_CUT << endl;
        outputFile << "myEvent.Mlb < " << MLB_CUT << endl;
        outputFile << "myEvent.MT2W < " << MTW2_CUT << endl;
        outputFile <<  MET_BOUND_250 <=" < myEvent.pfmet <= MET_BOUND_325 << endl;
    } 
    
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250)  return false;
    if (myEvent.pfmet > MET_BOUND_325)  return false;

    return true;
}
*/
bool SRbin1()
{

    
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250)  return false;
    if (myEvent.pfmet > MET_BOUND_325)  return false;

    return true;
}

bool SRbin2()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet <= MET_BOUND_325)  return false;

    return true;
}

bool SRbin3()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250)  return false;
    if (myEvent.pfmet > MET_BOUND_375)  return false;

    return true;
}

bool SRbin4()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet <= MET_BOUND_375)  return false;

    return true;
}

bool SRbin5()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.pfmet <= MET_BOUND_250)  return false;

    return true;
}

bool SRbin6()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet < MET_BOUND_250)  return false;
    if (myEvent.pfmet > MET_BOUND_300)  return false;

    return true;
}

bool SRbin7()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet <= MET_BOUND_300)  return false;
    if (myEvent.pfmet > MET_BOUND_400)  return false;

    return true;
}

bool SRbin8()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.pfmet <= MET_BOUND_400)  return false;

    return true;
}

bool goesInPreselectionNoBVeto()
{
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    //if ((!myEvent.isolatedTrackVeto) || (!myEvent.tauVeto)) return false;

    return true;
}


bool goesInPreselectionMTtailNoBeto(){return (goesInPreselectionNoBVeto() && goesInMTtail());}
bool goesInPreselectionMTtail()     { return (goesInPreselection() && goesInMTtail());     }
bool goesInPreselectionMTpeak()     { return (goesInPreselection() && goesInMTpeak());     }
bool goesInPreselectionMTinverted() { return (goesInPreselection() && goesInMTinverted()); }

bool goesIn0BtagControlRegion()
{

    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags != 0 )  return false;
    //if ((!myEvent.isolatedTrackVeto) || (!myEvent.tauVeto)) return false;

    return true;
}

bool goesIn0BtagControlRegionMTtail()     { return (goesIn0BtagControlRegion() && goesInMTtail());     }
bool goesIn0BtagControlRegionMTpeak()     { return (goesIn0BtagControlRegion() && goesInMTpeak());     }
bool goesIn0BtagControlRegionMTinverted() { return (goesIn0BtagControlRegion() && goesInMTinverted()); }

bool goesInDileptonControlRegion(short int nJetCut = -1)
{
    if (myEvent.ngoodleps != 2) return false;

    if (nJetCut == -1) { if (myEvent.ngoodjets < 1) return false; }
    else if (nJetCut == 2)  { if ((myEvent.ngoodjets != 1) && (myEvent.ngoodjets != 2)) return false; }
    else if (nJetCut == 3)  { if (myEvent.ngoodjets != 3) return false; }
    else if (nJetCut == 4)  { if (myEvent.ngoodjets <  4) return false; }

    if (myEvent.ngoodbtags < NBJET_CUT)   return false;
    if (myEvent.pfmet   < MET_CUTLL)   return false;

    // Remove same-sign events
    if ((myEvent.lep1_pdgid < 0) && (myEvent.lep2_pdgid < 0)) return false;
    if ((myEvent.lep1_pdgid > 0) && (myEvent.lep2_pdgid > 0)) return false;

    // Remove Z mass peak
    if (fabs((myEvent.leadingLepton + myEvent.secondLepton).M() - 91) < 15) return false;

    return true;
}

bool goesInVetosControlRegion()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.pfmet   < MET_CUT) return false;

    // Apply reversed vetos
    //if ((myEvent.isolatedTrackVeto) && (myEvent.tauVeto)) return false;

    return true;
}

bool goesInDileptonControlRegionMTtail()     { return (goesInDileptonControlRegion() && goesInMTtail());     }
bool goesInDileptonControlRegionMTpeak()     { return (goesInDileptonControlRegion() && goesInMTpeak());     }
bool goesInDileptonControlRegionMTinverted() { return (goesInDileptonControlRegion() && goesInMTinverted()); }

bool goesInVetoControlRegionMTtail()     { return (goesInVetosControlRegion() && goesInMTtail());     }
bool goesInVetoControlRegionMTpeak()     { return (goesInVetosControlRegion() && goesInMTpeak());     }
bool goesInVetoControlRegionMTinverted() { return (goesInVetosControlRegion() && goesInMTinverted()); }

// Single-lepton channels definitions
// ##################################

bool goesInSingleElecChannel()
{
    // Keep only events with ngoodleps == 1
    if (myEvent.ngoodleps != 1) return false;
    //if (!myEvent.HLT_SingleE) return false;
    // For data, keep only events from SingleElec dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "SingleElec") || (!myEvent.triggerElec)) return false;
    }*/

    // Remove electrons with pT < 30 GeV
    //if (myEvent.leadingLepton.Pt() < 30)  return false;

    // Keep only events with an electron as leading lepton
    return (abs(myEvent.lep1_pdgid) == 11);
}

bool goesInSingleElecBarrelChannel()
{
    // Keep only events with ngoodleps == 1
    if (myEvent.ngoodleps != 1) return false;
    if (myEvent.lep1_eta>=1.479) return false;
    //if (!myEvent.HLT_SingleE) return false;
    // For data, keep only events from SingleElec dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "SingleElec") || (!myEvent.triggerElec)) return false;
    }*/

    // Remove electrons with pT < 30 GeV
    //if (myEvent.leadingLepton.Pt() < 30)  return false;

    // Keep only events with an electron as leading lepton
    return (abs(myEvent.lep1_pdgid) == 11);
}

bool goesInSingleElecECChannel()
{
    // Keep only events with ngoodleps == 1
    if (myEvent.ngoodleps != 1) return false;
    if (myEvent.lep1_eta<1.479) return false;
    //if (!myEvent.HLT_SingleE) return false;
    // For data, keep only events from SingleElec dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "SingleElec") || (!myEvent.triggerElec)) return false;
    }*/

    // Remove electrons with pT < 30 GeV
    //if (myEvent.leadingLepton.Pt() < 30)  return false;

    // Keep only events with an electron as leading lepton
    return (abs(myEvent.lep1_pdgid) == 11);
}


bool goesInSingleMuonChannel()
{
    // Keep only events with ngoodleps == 1
    if (myEvent.ngoodleps != 1) return false;
    //if (!myEvent.HLT_SingleMu) return false;
    // For data, keep only events from SingleMuon dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "SingleMuon") || ((!myEvent.triggerMuon) && (!myEvent.xtriggerMuon))) return false;

        // Take care of the splitting due to x-trigger
        if ((myEvent.leadingLepton.Pt() >= 26) && (!myEvent.triggerMuon))  return false;
        if ((myEvent.leadingLepton.Pt() <  26) && (!myEvent.xtriggerMuon)) return false;
    }
    */

    // Keep only events with a muon as leading lepton
    return (abs(myEvent.lep1_pdgid) == 13);
}


bool goesInSingleLeptonChannel()
{
    // Single-lepton channel is the union of e-channel + mu-channel
    bool test = (goesInSingleElecChannel()
              || goesInSingleMuonChannel());
    return test;
}



// Double-lepton channels definitions
// ##################################

bool goesInDoubleElecChannel()
{
    // Keep only events with ngoodleps == 2
    if (myEvent.ngoodleps != 2) return false;
    // For data, keep only events from DoubleElec dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "DoubleElec") || (!myEvent.triggerDoubleElec)) return false;
    }
    */
    // Keep only events with two electrons
    return ((abs(myEvent.lep1_pdgid) == 11)
         && (abs(myEvent.lep2_pdgid)  == 11));
}

bool goesInDoubleMuonChannel()
{
    // Keep only events with ngoodleps == 2
    if (myEvent.ngoodleps != 2) return false;
    // For data, keep only events from DoubleMuon dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "DoubleMuon") || (!myEvent.triggerDoubleMuon)) return false;
    }
    */
    // Keep only events with two muons
    return ((abs(myEvent.lep1_pdgid) == 13)
         && (abs(myEvent.lep2_pdgid)  == 13));
}

bool goesInMuonElecChannel()
{
    // Keep only events with ngoodleps == 2
    if (myEvent.ngoodleps != 2) return false;
    // For data, keep only events from SingleMuon channel that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "MuEl") || (!myEvent.triggerMuonElec)) return false;
    }
    */
    // Keep only events with an electron and a muon
    return   (((abs(myEvent.lep1_pdgid) == 13)
            && (abs(myEvent.lep2_pdgid)  == 11))
       ||     ((abs(myEvent.lep1_pdgid) == 11)
            && (abs(myEvent.lep2_pdgid)  == 13)));
}



bool goesInDoubleLeptonChannel()
{
    // Double-lepton channel is the union of ee, mumu and emu channels
    bool test = (goesInDoubleElecChannel()
              || goesInDoubleMuonChannel()
              ||   goesInMuonElecChannel()  );
    return test;
}



// Weighting information
// #####################

float getLumi()
{
         if (goesInSingleElecChannel())  return 19508.0;
    else if (goesInSingleMuonChannel())  return 19514.0;
    else if (goesInMuonElecChannel  ())  return 19513.0;
    else if (goesInDoubleMuonChannel())  return 19504.0;
    else if (goesInDoubleElecChannel())  return 19517.0;
    else                                 return 0.0;
}

float getWeight(const string& dataset, int nofEventsInTree, float sf_fracEvent = 1) //sf_fracEvent = used if we read only a part of the events
{
    float weight = 1.0;
    if (sampleType == "data") return weight;

    // Get the lumi
    //float lumi = getLumi();
    float lumi = 10000.;

    // Normalize to cross section times lumi
    int nofInitEvent = myEvent.totalNumberOfInitialEvent;
    /*
    cout<<"nofIniEvent = "<<nofInitEvent<<endl;
    cout<<"xs = "<<CrossSection(dataset)<<endl;
    cout<<"lumi = "<<lumi<<" sf = "<<sf_fracEvent<<endl;
    */
    //nofEventsInTree can be different from nofIniEvent due to skimming and preselection
    if(nofInitEvent!=0)
    	//weight = CrossSection(dataset) * lumi * sf_fracEvent *  nofEventsInTree/ nofInitEvent;
    	weight = (float) CrossSection(dataset) * lumi / sf_fracEvent / nofInitEvent;
    else 
    	weight = 0;
    
    // Normalize to cross section times lumi
    //weight *= myEvent.weightCrossSection * lumi;

    /*
    // Apply trigger efficiency weights for singleLepton channels
    if (myEvent.ngoodleps == 1)
    {
        weight *= myEvent.weightTriggerEfficiency;
    }
    // Apply trigger efficiency weights for doubleLepton channels
    else if (myEvent.ngoodleps == 2)
    {
        if (goesInDoubleElecChannel()) weight *= 0.96;
        if (goesInDoubleMuonChannel()) weight *= 0.88;
        if (goesInMuonElecChannel())   weight *= 0.93;
    }

    // Apply lepton ID efficiency and isolation scale factor to singleLepton channels
    if (myEvent.ngoodleps == 1)
    {
        weight *= myEvent.leadingLeptonIdEfficiency * myEvent.leadingLeptonIsoScaleFactor;

    }
    // TODO not sure about this, to be confirmed
    else if (myEvent.ngoodleps == 2)
    {
        weight *= myEvent.leadingLeptonIdEfficiency * myEvent.leadingLeptonIsoScaleFactor;
        weight *= myEvent.secondLeptonIdEfficiency  * myEvent.secondLeptonIsoScaleFactor;
    }

    // Apply pile-up weight
    // TODO : Do we confirm we'll use also this for signal ?
        weight *= myEvent.weightPileUp;


    if (sampleType == "signal")
    {
        // For signal, apply ISR reweighting
        weight *= myEvent.weightISRmodeling;

        // Check if event has been used in BDT training
        if (myEvent.isUsedInBDT == 0) weight *= 2;
        else                          weight  = 0;
    }


    // For ttbar only, apply topPt reweighting
    if (sampleName.find("ttbar") != string::npos)
        weight *= myEvent.weightTopPt;

    */
    return weight;
}


