// ############################################################
// # Usage                                                    #
// #                                                          #
// # - Include this header in your code                       #
// # - Create an instance of babyEvent, for instance event :  #
// #      babyEvent event;                                    #
// # - Open your tree, for example :                          #
// #      TTree* theTree = f.Get("babyTuple");                #
// # - Init branches by calling :                             #
// #      InitializeBranchesForReading(theTree,&event);       #
// # - To read the i-ish event, call :                        #
// #      ReadEvent(theTree,i,&event);                        #
// # - Get the value of your branch by acessing               #
// #      event.branchName;                                   #
// ############################################################
 
#ifndef babyFormat
#define babyFormat

#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "TMath.h"
#include "TVector3.h"

//#define USE_GEN_INFO

// ##########################
// #  Baby event structure  #
// ##########################
 
typedef struct
{
    bool hasGenInfo;
    //int           numberOfSelectedLeptons = -13;
    vector<bool>*   ak4pfjets_passMEDbtag;
    vector<float>* ak4pfjets_mva;
    vector<float>* ak4pfjets_loose_puid;

   float ak4pfjets_pt;
   float ak4pfjets_pt2;
  float ak4pfjets_pt3;
  float ak4pfjets_pt4;
  float ak4pfjets_pt_B;
  float ak4pfjets_leadJet_identification;
    vector<float>* ak4pfjets_phi;
    vector<float>* ak4pfjets_eta;
    vector<float>* ak4pfjets_mass;
    vector<int>*   ak4pfjets_partonFlavour;
    vector<int>*   ak4pfjets_hadronFlavour;
    vector<bool>*  ak4pfjets_loose_pfid;
    vector<float>* ak4pfjets_CSV;
    //vector<float> ak4pfjets_qgtag;
    float         ak4pfjets_rho = -13;
    float         dphi_ak4pfjets_met = -13;
    float         dphi_ak4pfjets_met_rl = -13;
    //float         jetsPt = -13;
    float         pfmet_phi = -13;
    float         pfmet_phi_rl = -13;
    //float         met_sig = -13;
    //vector<float> ak8pfjets_pt;
    //vector<float> ak10pfjets_pt;
    int          HLT_SingleEl;
    int          HLT_SingleMu;
    int          HLT_MET;
    int          HLT_MET100_MHT100;
    int          HLT_MET110_MHT110;
    int          HLT_MET120_MHT120;
    int          HLT_DiEl;
    int          HLT_DiMu;
    int          HLT_MuE;
    bool         trigger = 0;
    bool         json_ok=0;
    int          cutflow=-1;
    int          cutflow2=-1;
    int          cutflow3=-1;
    //vector<float> ak8pfjets_phi;
    //vector<float> ak10pfjets_phi;
    bool         lep1_passMediumID = -13;
    bool         lep1_passTightID = -13;
    bool         lep1_passLooseID = -13;
    bool         lep2_mass = -13;
    UInt_t           runId = -13;
    //float         pv_ndof = -13;
    //int           numberOfGeneratedLeptons = -13;
    float         HT = -13;
    float         lep2_passMediumID = -13;
    float         lep2_passLooseID = -13;
    float         lep2_passTightID = -13;
    float         lep2_phi = -13;
    //float         secondLeptonPt = -13;
    float         lep1_pt = -13;
    UInt_t           eventId = -13;
    //int           numberOfSelectedJets = -13;
    //int           pv_isFake = -13;
    float         lep2_pt = -13;
    float         leadingLeptonPt = -13;
    //vector<float> ak8pfjets_eta;
    //vector<float> ak10pfjets_eta;
    float         jetsEta = -13;
    float         leadingLeptonEta = -13;
    //float         Mjjj = -13;
    //vector<float> ak8pfjets_mass;
    //vector<float> ak10pfjets_mass;
    //vector<float> ak8pfjets_pruned_mass;
    //vector<float> ak10pfjets_pruned_mass;
    //vector<float> ak8pfjets_corrpruned_mass;
    //vector<float> ak8pfjets_trimmed_mass;
    //vector<float> ak10pfjets_trimmed_mass;
    //vector<float> ak8pfjets_softdrop_mass;
    //vector<float> ak10pfjets_softdrop_mass;
    //vector<int>   ak8pfjets_nSubJets;
    //vector<int>   ak10pfjets_nSubJets;
    //vector<float> ak8pfjets_tau1;
    //vector<float> ak10pfjets_tau1;
    //vector<float> ak8pfjets_tau2;
    //vector<float> ak10pfjets_tau2;
    int           ngoodbtags = -13;
    int           ntightbtags = -13;
    float         lep2_eta = -13;
    float         crossSection = -13;
    float         kfactor = -13;
    float         mt_met_lep = -13;
    float         mt_met_lep_rl = -13;
    //float         M3b = -13;
    int           lep1_pdgid = -13;
    int           lep1_mc_motherid = -13;
  int           lep1_mc_gmotherid = -13;
    //float         jetsCSVv2 = -13;
    float         scale1fb = -13;
    float         MT = -13;
    //float         jetsPhi = -13;
    int           ngoodjets = -13;
    bool         lep1_passVeto = -13;
    float         pfmet = -13;
    float         pfmet_rl = -13;
    float         lep1_MiniIso = -13;
    float         lep1_relIso = -13;
    UInt_t           lumiId = -13;
    float         lep1_phi = -13;
    int           lep2_pdgid = -13;
    //float         pv_z = -13;
    float         Mlb = -13;
    float         selectionCode = -13;
    float         jetsPUid = -13;
    //float         Mlb_leadb = -13;
    //float         leadingLeptonPhi = -13;
    //int           numberOfBTaggedJets = -13;
    //int           genlepsfromtop = -13;
    int          is1lepFromW = -13;
    int          is1lepFromTop = -13;
    ///int           secondLeptonId = -13;
    //float         jetsCSV = -13;
    float         weight_PU = -13;
    float         weight_PUup = -13;
    float         weight_PUdown = -13;
    float         weight_btagsf = -13;
    float         weight_btagsf_heavy_DN = -13;
    float         weight_btagsf_heavy_UP = -13;
    float         weight_btagsf_light_DN = -13;
    float         weight_btagsf_light_UP = -13;
    float         weight_tightbtagsf = -13;
    float         weight_tightbtagsf_heavy_DN = -13;
    float         weight_tightbtagsf_heavy_UP = -13;
    float         weight_tightbtagsf_light_DN = -13;
    float         weight_tightbtagsf_light_UP = -13;
    float         weight_lepSF = -13;
    float         weight_lepSF_down = -13;
    float         weight_lepSF_up = -13;
    float         weight_vetoLepSF = -13;
    float         weight_lepSF_fastSim = -13;
    float         weight_ISR = -13;
    float         weight_ISRdown = -13;
    float         weight_ISRup = -13;
    float         weight_ISRnjets = -13;
    float         weight_ISRnjets_UP = -13;
    float         weight_ISRnjets_DN = -13;
    float         pdf_down_weight = -13;
    float         pdf_up_weight = -13;
    float         weight_q2_down = -13;
    float         weight_q2_up = -13;
    float         weight_alphas_down = -13;
    float         weight_alphas_up = -13;
    float         nupt=-13;
    float         lep1_mass = -13;
    float         lep2_passVeto = -13;
    float         lep1_eta = -13;
    float         lep2_MiniIso = -13;
    float         lep2_relIso = -13;
    //float         pv_rho = -13;
    float         MT2W = -13;
    float         MT2W_rl = -13;
    int          isZtoNuNu = -13;
    //float         ETmiss = -13;
    //int           numberOfSelectedElectrons = -13;
    int           ngoodleps = -13;
    int           nlooseleps = -13;
    int           nvetoleps = -13;
    float         topness = -13;
    float         topnessMod = -13;
    float         topnessMod_rl = -13;
    //float         chi2 = -13;
    //float         ETmissPhi = -13;
    //int64_t      totalNumberOfInitialEvent = -13;
    //int           numberOfSelectedMuons = -13;
    float         hadronic_top_chi2 = -13;
    int           nvertex = -13;
    //int           puIntime = -13;
    float           puTrue = -13;
    bool           is_data = 0;
    float           filt_met = 0;
    bool          filt_fastsimjets;
    float          filt_badChargedCandidateFilter = 0;
    bool           filt_jetWithBadMuon = 0;
    bool           filt_pfovercalomet = 0;
    float           filt_badMuonFilter = 0;
    bool           filt_duplicatemuons = 0;
    bool           filt_badmuons = 0;
    bool           filt_nobadmuons = 0;
    int           is0lep = -13;
    int           is1lep = -13;
    int           is2lep = -13;
    float         mass_lsp = -13;
    float         mass_chargino = -13;
    float         mass_stop = -13;
    uint32_t      nentries = 0;
    vector<int>*   genqs_id;
    vector<int>*   genqs_motherid;
    vector<int>*   genqs_motheridx;
    vector<int>*   genqs_status;
    float         top_pt = -13;
    float         atop_pt = -13;
    float         lep1_dphiMET = -13;
    float         lep1_dphiMET_rl = -13;
    float         lep2_dphiMET = -13;
    float         lep2_dphiMET_rl = -13;
  // float         ak4pfjets_CSV=-13;

    int            ngoodtaus = 0;
    vector<bool>*  tau_passID ;
    vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >* p4tau_p4;
    vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > tau_p4;

    uint32_t      nEvents = 0;
    //uint32_t      nthentry = 0;

    //@MJ@ TODO fill these
    /*float         pfmet_rl_jup = -13;
    float         pfmet_rl_jdown = -13;
    float         MT2W_rl_jup = -13; //dphi?
    float         MT2W_rl_jdown = -13;
    float         mindphi_met_j1_j2_rl_jup = -13;
    float         mindphi_met_j1_j2_rl_jdown = -13;
    float         mt_met_lep_rl_jup = -13;
    float         mt_met_lep_rl_jdown = -13;
    float         topnessMod_rl_jup = -13;
    float         topnessMod_rl_jdown = -13;
    float         jup_ngoodjets = -13;
    float         jup_ngoodbtags = -13;
    float         jdown_ngoodjets = -13;
    float         jdown_ngoodbtags = -13;
*/
    vector<float> wNormalization;


   
    /*vector<float> gen_neutralino_m;
    vector<float> gen_stop_m;
    vector<float>* pointerForgen_neutralino_m;
    vector<float>* pointerForgen_stop_m;
    vector<float> genTops_pt;
    vector<int> genTops_pdgid;
    vector<float>* pointerForgenTops_pt;
    vector<int>* pointerForgenTops_pdgid;*/

  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep1_p4;
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >* p4lep1_p4 = NULL;
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep2_p4;
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >* p4lep2_p4 = NULL;
  vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genqs_p4;
  vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >* p4genqs_p4 = NULL;
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > ak4pfjets_leadbtag_p4;
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >* p4ak4pfjets_leadbtag_p4 = NULL;
 vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >  ak4pfjets_p4;
 vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >* p4ak4pfjets_p4 = NULL;

    bool PassTrackVeto = -13;
    bool PassTauVeto = -13;

    float metGen_pt = -13;
    float metGen_phi = -13;


    vector<float>* genweights;
    vector<string>* genweightsID;
    vector<int>* 	  gensusy_id;
    vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >* gensusy_p4;

    //float DeltaRlj = -13;
    //vector<float> ak8pfjets_nsubjettines; //WARNING! this is empty now
    //vector<float> ak10pfjets_nsubjettines; //WARNING! this is empty now

    #ifdef USE_NEW_VAR
    //float ST = -13;
    //float LP = -13;
    //float Meff = -13;
    //float MTdeco_Q = -13;
    //float DeltaPtbb = -13;
    //float DeltaPhibb = -13;
    //float DeltaRbb = -13;
    //float dphi_Wlep = -13;
    //float dR_lep_leadb = -13;
    float ak4_HT = -13;
    //float ak4_htssm = -13;
    //float ak4_htosm = -13;
    #endif

    vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genleps_p4;
    vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >* p4genleps_p4 = NULL;
    #ifdef USE_GEN_LEP
    vector<int>*  pointer_genleps_id;
    vector<int>*  pointer_genleps_motherid;
    vector<int>*  pointer_genleps_gmotherid;
    vector<bool>* pointer_genleps_isLastCopy;
    vector<bool>* pointer_genleps_fromHardProcessFinalState;
    vector<bool>* pointer_genleps_fromHardProcessDecayed;
    vector<int>*  pointer_genleps_gentaudecay;    
    vector<int>  genleps_id;
    vector<int>  genleps_motherid;
    vector<int>  genleps_gmotherid;
    vector<bool> genleps_isLastCopy;
    vector<bool> genleps_fromHardProcessFinalState;
    vector<bool> genleps_fromHardProcessDecayed;
    vector<int>  genleps_gentaudecay;    
    #endif

    //@MJ@ TODO gen info!
    #ifdef USE_GEN_LOSTLEPTON
    vector<float>* pointerForGenLostLeptons_pt;
    vector<float>* pointerForGenLostLeptons_eta;
    vector<int>* pointerForGenLostLeptons_pdgid;
    
    vector<float> genLostLeptons_pt;
    vector<float> genLostLeptons_eta;
    vector<int> genLostLeptons_pdgid;
    #endif

    #ifdef USE_GEN_INFO
    int gen_n;
    vector<float> gen_pt;
    vector<float> gen_eta;
    vector<float> gen_phi;
    vector<float> gen_m;
    vector<int> gen_status;
    vector<int> gen_id;
    vector<int> gen_daughter_n;
    vector<vector< int> > gen_daughter_index;
    
    vector<float>* pointerForgen_pt;
    vector<float>* pointerForgen_eta;
    vector<float>* pointerForgen_phi;
    vector<float>* pointerForgen_m;
    vector<int>* pointerForgen_status;
    vector<int>* pointerForgen_id;
    vector<int>* pointerForgen_daughter_n;
    vector<vector< int> >* pointerForgen_daughter_index;

    #endif
    
    #ifdef USE_GEN_INFO_EXT
    vector<int> gen_charge;
    vector<int> gen_index;
    vector<int> gen_mother_index;
    vector<int>* pointerForgen_charge;
    vector<int>* pointerForgen_mother_index;
    vector<int>* pointerForgen_index;
    
    #endif
    
    #ifdef USE_SKIMMING_VAR
    float minDPhi_jmet = -13; 
    float genMET = -13;
    //float Pz_ttbar = -13;
    //float Pt_ttbar = -13;
    //float E_ttbar = -13;
    //float M_ttbar = -13;
    //int ttbar_decay = -13;
    //int nofISR = -13;
    //int index_b = -13;
    //int index_bbar = -13;
    //bool matched;
    //vector<int> isr_id;
    //vector<TLorentzVector> isr_p4;
    //vector<TLorentzVector> neutrinos;
    //vector<TLorentzVector> neutralinos;
    
    //vector<int>* pointerForisr_id;
    //vector<TLorentzVector>* pointerForisr_p4;
    //vector<TLorentzVector>* pointerForneutrinos;
    //vector<TLorentzVector>* pointerForneutralinos;
   #endif

   // Intermediate pointers for special types
    // Yes, this shit is needed because ROOT is crap.
 
  //  vector<float>* pointerForak4pfjets_puid;
  //  vector<float>* pointerForak4pfjets_pt;
  /*  vector<float>* pointerForak8pfjets_pt;
      vector<float>* pointerForak10pfjets_pt;
    vector<float>* pointerForak4pfjets_phi;
    vector<float>* pointerForak8pfjets_phi;
     vector<float>* pointerForak10pfjets_phi;
     vector<float>* pointerForak4pfjets_eta;
    vector<float>* pointerForak8pfjets_eta;
      vector<float>* pointerForak10pfjets_eta;
      vector<float>* pointerForak4pfjets_mass;
      vector<float>* pointerForak8pfjets_mass;
      vector<float>* pointerForak8pfjets_tau1;
      vector<float>* pointerForak8pfjets_tau2;
    vector<float>* pointerForak8pfjets_pruned_mass;
    vector<float>* pointerForak8pfjets_corrpruned_mass;
    vector<float>* pointerForak8pfjets_trimmed_mass;
    vector<float>* pointerForak8pfjets_softdrop_mass;
    vector<int>*   pointerForak8pfjets_nSubJets;
    vector<float>* pointerForak10pfjets_mass;
    vector<float>* pointerForak10pfjets_tau1;
    vector<float>* pointerForak10pfjets_tau2;
    vector<float>* pointerForak10pfjets_pruned_mass;
    vector<float>* pointerForak10pfjets_trimmed_mass;
    vector<float>* pointerForak10pfjets_softdrop_mass;
    vector<int>*   pointerForak10pfjets_nSubJets;
    vector<bool>*  pointerForak4pfjets_loose_pfid;
    vector<float>* pointerForak4pfjets_CSV;
    vector<int>* pointerForak4pfjets_partonFlavour;
    vector<int>* pointerForak4pfjets_hadronFlavour;
    vector<float>* pointerForak4pfjets_qgtag;

    vector<string>* pointerForTriggerName;
    vector<bool>* pointerForTriggerPass;*/

   //Add content
   //TLorentzVector leadingLepton;
   //TLorentzVector secondLepton;

}
babyEvent;
 
// #############################
// #  Branches initialization  #
// #############################
 
void InitializeBranchesForReading(TTree* theTree, babyEvent* myEvent, TFile* f)
{
    myEvent->hasGenInfo = false;

    //fill weight normalization for data and MC
    
    TH1D* hWeight = NULL;
    hWeight = dynamic_cast<TH1D*>(f->Get("h_counter"));
    if (hWeight != NULL)
    {
	// cout << "weight hist " << hWeight << " get " <<  f->Get("h_counter")<< " file " << f << endl;
	myEvent->wNormalization.clear();
	myEvent->wNormalization.push_back(-13); //@MJ@ ATTENTION vector.at(0) random number to sy synchronized with binning numbering
	for(uint32_t n=0; n<hWeight->GetNbinsX(); n++)
	{
            myEvent->wNormalization.push_back(hWeight->GetBinContent(n+1));
	}
	// cout << "normalization size " << myEvent->wNormalization.size() << endl;
    }
    /* myEvent->pointerForTriggerName = NULL;
    myEvent->pointerForTriggerPass = NULL;
    myEvent->pointerForgen_neutralino_m = 0;
    myEvent->pointerForgen_stop_m = 0;
    myEvent->pointerForgenTops_pt = 0;
    myEvent->pointerForgenTops_pdgid = 0;
    */
    #ifdef USE_JETS
    // myEvent->pointerForak4pfjets_pt = 0;
    /*myEvent->pointerForak4pfjets_phi = 0;
    myEvent->pointerForak4pfjets_eta = 0;
    myEvent->pointerForak4pfjets_mass = 0;
    myEvent->pointerForak4pfjets_puid = 0;
    myEvent->pointerForak4pfjets_loose_pfid = 0;
    myEvent->pointerForak4pfjets_CSV = 0;
    myEvent->pointerForak4pfjets_partonFlavour = 0;
    myEvent->pointerForak4pfjets_hadronFlavour = 0;
    myEvent->pointerForak4pfjets_qgtag = 0;
    */#endif

    /*
    #ifdef USE_GEN_INFO_EXT
    myEvent->pointerForgen_charge = 0;
    myEvent->pointerForgen_index = 0;
    myEvent->pointerForgen_mother_index = 0;
    #endif 
    
    #ifdef USE_GEN_LOSTLEPTON
    myEvent->pointerForGenLostLeptons_pt = 0;
    myEvent->pointerForGenLostLeptons_eta = 0;
    myEvent->pointerForGenLostLeptons_pdgid = 0;
    #endif
*/

   //entries in tree
   myEvent->nentries = theTree->GetEntries();


    if(theTree->GetListOfBranches()->FindObject("xsec"))
        theTree->SetBranchAddress("xsec",            &(myEvent->crossSection));
    if(theTree->GetListOfBranches()->FindObject("kfactor"))
        theTree->SetBranchAddress("kfactor",            &(myEvent->kfactor));
    if(theTree->GetListOfBranches()->FindObject("scale1fb"))
        theTree->SetBranchAddress("scale1fb", &(myEvent->scale1fb));


    //@MJ@ TODO gen info - move
    if(theTree->GetListOfBranches()->FindObject("genweights"))
        theTree->SetBranchAddress("genweights", &(myEvent->genweights));
    if(theTree->GetListOfBranches()->FindObject("genweightsID"))
        theTree->SetBranchAddress("genweightsID", &(myEvent->genweightsID));
    if(theTree->GetListOfBranches()->FindObject("gensusy_id"))
        theTree->SetBranchAddress("gensusy_id", &(myEvent->gensusy_id));
    if(theTree->GetListOfBranches()->FindObject("gensusy_p4"))
        theTree->SetBranchAddress("gensusy_p4", &(myEvent->gensusy_p4));


    // multiplicities
    if(theTree->GetListOfBranches()->FindObject("ngoodleps"))
        theTree->SetBranchAddress("ngoodleps",               &(myEvent->ngoodleps));
    if(theTree->GetListOfBranches()->FindObject("nlooseleps"))
        theTree->SetBranchAddress("nlooseleps",               &(myEvent->nlooseleps));
    if(theTree->GetListOfBranches()->FindObject("nvetoleps"))
        theTree->SetBranchAddress("nvetoleps",               &(myEvent->nvetoleps));
    if(theTree->GetListOfBranches()->FindObject("is_data"))
        theTree->SetBranchAddress("is_data",               &(myEvent->is_data));
    if(theTree->GetListOfBranches()->FindObject("is0lep"))
        theTree->SetBranchAddress("is0lep",               &(myEvent->is0lep));
    if(theTree->GetListOfBranches()->FindObject("is1lep"))
        theTree->SetBranchAddress("is1lep",               &(myEvent->is1lep));
    if(theTree->GetListOfBranches()->FindObject("is2lep"))
        theTree->SetBranchAddress("is2lep",               &(myEvent->is2lep));
    if(theTree->GetListOfBranches()->FindObject("is1lepFromW"))
        theTree->SetBranchAddress("is1lepFromW",          &(myEvent->is1lepFromW));
    if(theTree->GetListOfBranches()->FindObject("is1lepFromTop"))
        theTree->SetBranchAddress("is1lepFromTop",          &(myEvent->is1lepFromTop));
    if(theTree->GetListOfBranches()->FindObject("isZtoNuNu"))
        theTree->SetBranchAddress("isZtoNuNu",                &(myEvent->isZtoNuNu));
    if(theTree->GetListOfBranches()->FindObject("nEvents"))
        theTree->SetBranchAddress("nEvents",                &(myEvent->nEvents));
    if(theTree->GetListOfBranches()->FindObject("filt_met"))
        theTree->SetBranchAddress("filt_met",               &(myEvent->filt_met));
    if(theTree->GetListOfBranches()->FindObject("filt_badChargedCandidateFilter"))
    if(theTree->GetListOfBranches()->FindObject("filt_fastsimjets"))
        theTree->SetBranchAddress("filt_fastsimjets",               &(myEvent->filt_fastsimjets));
        theTree->SetBranchAddress("filt_badChargedCandidateFilter",               &(myEvent->filt_badChargedCandidateFilter));
    if(theTree->GetListOfBranches()->FindObject("filt_jetWithBadMuon"))
        theTree->SetBranchAddress("filt_jetWithBadMuon",               &(myEvent->filt_jetWithBadMuon));
    if(theTree->GetListOfBranches()->FindObject("filt_pfovercalomet"))
        theTree->SetBranchAddress("filt_pfovercalomet",               &(myEvent->filt_pfovercalomet));
    if(theTree->GetListOfBranches()->FindObject("filt_badMuonFilter"))
        theTree->SetBranchAddress("filt_badMuonFilter",               &(myEvent->filt_badMuonFilter));
    if(theTree->GetListOfBranches()->FindObject("filt_duplicatemuons"))
        theTree->SetBranchAddress("filt_duplicatemuons",               &(myEvent->filt_duplicatemuons));
    if(theTree->GetListOfBranches()->FindObject("filt_badmuons"))
        theTree->SetBranchAddress("filt_badmuons",               &(myEvent->filt_badmuons));
    if(theTree->GetListOfBranches()->FindObject("filt_nobadmuons"))
        theTree->SetBranchAddress("filt_nobadmuons",               &(myEvent->filt_nobadmuons));

    //trigger
    if(theTree->GetListOfBranches()->FindObject("HLT_SingleEl"))
        theTree->SetBranchAddress("HLT_SingleEl",               &(myEvent->HLT_SingleEl));
    if(theTree->GetListOfBranches()->FindObject("HLT_SingleMu"))
        theTree->SetBranchAddress("HLT_SingleMu",               &(myEvent->HLT_SingleMu));
    if(theTree->GetListOfBranches()->FindObject("HLT_MET"))
        theTree->SetBranchAddress("HLT_MET",               &(myEvent->HLT_MET));
    if(theTree->GetListOfBranches()->FindObject("HLT_MET100_MHT100"))
        theTree->SetBranchAddress("HLT_MET100_MHT100",               &(myEvent->HLT_MET100_MHT100));
    if(theTree->GetListOfBranches()->FindObject("HLT_MET110_MHT110"))
        theTree->SetBranchAddress("HLT_MET110_MHT110",               &(myEvent->HLT_MET110_MHT110));
    if(theTree->GetListOfBranches()->FindObject("HLT_MET120_MHT120"))
        theTree->SetBranchAddress("HLT_MET120_MHT120",               &(myEvent->HLT_MET120_MHT120));
    if(theTree->GetListOfBranches()->FindObject("HLT_DiEl"))
        theTree->SetBranchAddress("HLT_DiEl",               &(myEvent->HLT_DiEl));
    if(theTree->GetListOfBranches()->FindObject("HLT_DiMu"))
        theTree->SetBranchAddress("HLT_DiMu",               &(myEvent->HLT_DiMu));
    if(theTree->GetListOfBranches()->FindObject("HLT_MuE"))
        theTree->SetBranchAddress("HLT_MuE",               &(myEvent->HLT_MuE));
    //theTree->SetBranchAddress("trigger_name",            &(myEvent->pointerForTriggerName));
    //theTree->SetBranchAddress("trigger_pass",            &(myEvent->pointerForTriggerPass));
    //
    //signal related
    //theTree->SetBranchAddress("gen_neutralino_m", &(myEvent->pointerForgen_neutralino_m));
    //theTree->SetBranchAddress("gen_stop_m", &(myEvent->pointerForgen_stop_m));
    //gen top
    //theTree->SetBranchAddress("genTops_pt", &(myEvent->pointerForgenTops_pt));
    //theTree->SetBranchAddress("genTops_pdgid", &(myEvent->pointerForgenTops_pdgid));
    
    //vetos
    if(theTree->GetListOfBranches()->FindObject("PassTrackVeto"))
        theTree->SetBranchAddress("PassTrackVeto",            &(myEvent->PassTrackVeto));
    if(theTree->GetListOfBranches()->FindObject("PassTauVeto"))
        theTree->SetBranchAddress("PassTauVeto",             &(myEvent->PassTauVeto));
    
    //gen met
    if(theTree->GetListOfBranches()->FindObject("genmet"))
        theTree->SetBranchAddress("genmet",              &(myEvent->metGen_pt));
    if(theTree->GetListOfBranches()->FindObject("genmet_phi"))
        theTree->SetBranchAddress("genmet_phi",             &(myEvent->metGen_phi));
   
 
    if(theTree->GetListOfBranches()->FindObject("lep1_pdgid"))
        theTree->SetBranchAddress("lep1_pdgid",              &(myEvent->lep1_pdgid));
    if(theTree->GetListOfBranches()->FindObject("lep1_mc_motherid"))
        theTree->SetBranchAddress("lep1_mc_motherid",              &(myEvent->lep1_mc_motherid));
    if(theTree->GetListOfBranches()->FindObject("lep1_mc_gmotherid"))
      theTree->SetBranchAddress("lep1_mc_gmotherid",              &(myEvent->lep1_mc_gmotherid));
    #ifdef USE_LEP1
    if(theTree->GetListOfBranches()->FindObject("lep1_p4"))
    {
        //myEvent->lep1_p4 = NULL;
        theTree->SetBranchAddress("lep1_p4",              &(myEvent->p4lep1_p4));
        //myEvent->lep1_pt = myEvent->lep1_p4->Pt();
        //myEvent->lep1_eta = myEvent->lep1_p4->Eta();
        //myEvent->lep1_phi = myEvent->lep1_p4->Phi();
        //myEvent->lep1_mass = myEvent->lep1_p4->M();//@MJ@ TODO something is going on
    }
//    #endif
//    #ifdef USE_LEP1
    if(theTree->GetListOfBranches()->FindObject("lep1_passVeto"))
        theTree->SetBranchAddress("lep1_passVeto",           &(myEvent->lep1_passVeto));
    if(theTree->GetListOfBranches()->FindObject("lep1_passLooseID"))
        theTree->SetBranchAddress("lep1_passLooseID",       &(myEvent->lep1_passLooseID));
    if(theTree->GetListOfBranches()->FindObject("lep1_passMediumID"))
        theTree->SetBranchAddress("lep1_passMediumID",       &(myEvent->lep1_passMediumID));
    if(theTree->GetListOfBranches()->FindObject("lep1_passTightID"))
        theTree->SetBranchAddress("lep1_passTightID",       &(myEvent->lep1_passTightID));
    if(theTree->GetListOfBranches()->FindObject("lep1_MiniIso"))
        theTree->SetBranchAddress("lep1_MiniIso",            &(myEvent->lep1_MiniIso));
    if(theTree->GetListOfBranches()->FindObject("lep1_relIso"))
        theTree->SetBranchAddress("lep1_relIso",            &(myEvent->lep1_relIso));
    #endif


    if(theTree->GetListOfBranches()->FindObject("lep2_pdgid"))
        theTree->SetBranchAddress("lep2_pdgid",              &(myEvent->lep2_pdgid));
    #ifdef USE_LEP2
    if(theTree->GetListOfBranches()->FindObject("lep2_p4"))
    {
        theTree->SetBranchAddress("lep2_p4",              &(myEvent->p4lep2_p4));
        //myEvent->lep2_pt = myEvent->lep2_p4->Pt();
        //myEvent->lep2_eta = myEvent->lep2_p4->Eta();
        //myEvent->lep2_phi = myEvent->lep2_p4->Phi();
        //myEvent->lep2_mass = myEvent->lep2_p4->M();//@MJ@ TODO something is going on
    }
    #endif
    #ifdef USE_LEP2_EXT
    if(theTree->GetListOfBranches()->FindObject("lep2_passVeto"))
        theTree->SetBranchAddress("lep2_passVeto",           &(myEvent->lep2_passVeto));
    if(theTree->GetListOfBranches()->FindObject("lep2_passLooseID"))
        theTree->SetBranchAddress("lep2_passLooseID",       &(myEvent->lep2_passLooseID));
    if(theTree->GetListOfBranches()->FindObject("lep2_passMediumID"))
        theTree->SetBranchAddress("lep2_passMediumID",       &(myEvent->lep2_passMediumID));
    if(theTree->GetListOfBranches()->FindObject("lep2_passTightID"))
        theTree->SetBranchAddress("lep2_passTightID",       &(myEvent->lep2_passTightID));
    if(theTree->GetListOfBranches()->FindObject("lep2_MiniIso"))
        theTree->SetBranchAddress("lep2_MiniIso",            &(myEvent->lep2_MiniIso));
    if(theTree->GetListOfBranches()->FindObject("lep2_relIso"))
        theTree->SetBranchAddress("lep2_relIso",            &(myEvent->lep2_relIso));
    #endif



    #ifdef USE_GEN_LEP
    if(theTree->GetListOfBranches()->FindObject("genleps_p4"))
    {
        theTree->SetBranchAddress("genleps_p4",              &(myEvent->p4genleps_p4));
    }
    if(theTree->GetListOfBranches()->FindObject("genleps_id"))
        theTree->SetBranchAddress("genleps_id",              &(myEvent->pointer_genleps_id));
    if(theTree->GetListOfBranches()->FindObject("genleps_motherid"))
        theTree->SetBranchAddress("genleps_motherid",              &(myEvent->pointer_genleps_motherid));
    if(theTree->GetListOfBranches()->FindObject("genleps_gmotherid"))
        theTree->SetBranchAddress("genleps_gmotherid",              &(myEvent->pointer_genleps_gmotherid));
    if(theTree->GetListOfBranches()->FindObject("genleps_isLastCopy"))
        theTree->SetBranchAddress("genleps_isLastCopy",              &(myEvent->pointer_genleps_isLastCopy));
    if(theTree->GetListOfBranches()->FindObject("genleps_fromHardProcessFinalState"))
        theTree->SetBranchAddress("genleps_fromHardProcessFinalState",              &(myEvent->pointer_genleps_fromHardProcessFinalState));
    if(theTree->GetListOfBranches()->FindObject("genleps_fromHardProcessDecayed"))
        theTree->SetBranchAddress("genleps_fromHardProcessDecayed",              &(myEvent->pointer_genleps_fromHardProcessDecayed));
    if(theTree->GetListOfBranches()->FindObject("genleps_gentaudecay"))
        theTree->SetBranchAddress("genleps_gentaudecay",              &(myEvent->pointer_genleps_gentaudecay));
    #endif

    #ifdef USE_JETS
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_passMEDbtag"))
        theTree->SetBranchAddress("ak4pfjets_passMEDbtag",           &(myEvent->ak4pfjets_passMEDbtag));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_mva"))
        theTree->SetBranchAddress("ak4pfjets_mva",           &(myEvent->ak4pfjets_mva));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_CSV"))
        theTree->SetBranchAddress("ak4pfjets_CSV",           &(myEvent->ak4pfjets_CSV));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_parton_flavor"))
        theTree->SetBranchAddress("ak4pfjets_parton_flavor", &(myEvent->ak4pfjets_partonFlavour));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_hadron_flavor"))
        theTree->SetBranchAddress("ak4pfjets_hadron_flavor", &(myEvent->ak4pfjets_hadronFlavour));
    //if(theTree->GetListOfBranches()->FindObject("ak4pfjets_loose_puid"))
    //    theTree->SetBranchAddress("ak4pfjets_loose_puid",          &(myEvent->ak4pfjets_loose_puid));
    //if(theTree->GetListOfBranches()->FindObject("ak4pfjets_loose_pfid"))
    //    theTree->SetBranchAddress("ak4pfjets_loose_pfid",    &(myEvent->ak4pfjets_loose_pfid));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_rho"))
        theTree->SetBranchAddress("ak4pfjets_rho",           &(myEvent->ak4pfjets_rho));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_leadbtag_p4"))
        theTree->SetBranchAddress("ak4pfjets_leadbtag_p4",           &(myEvent->p4ak4pfjets_leadbtag_p4));
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_p4"))
        theTree->SetBranchAddress("ak4pfjets_p4",           &(myEvent->p4ak4pfjets_p4));

    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_CSV"))
      theTree->SetBranchAddress("ak4pfjets_CSV",               &(myEvent->ak4pfjets_CSV));
    #endif

    #ifdef USE_TAU
    //taus
    if (theTree->GetListOfBranches()->FindObject("ngoodtaus"))
        theTree->SetBranchAddress("ngoodtaus",            &(myEvent->ngoodtaus));
    if (theTree->GetListOfBranches()->FindObject("tau_passID"))
        theTree->SetBranchAddress("tau_passID",            &(myEvent->tau_passID));
    if (theTree->GetListOfBranches()->FindObject("tau_p4"))
        theTree->SetBranchAddress("tau_p4",            &(myEvent->p4tau_p4));
    #endif


    //Fat jets, PV
    //
    /*#ifdef USE_AK8_JETS
    //std::cout << "reading bak8 info" << std::endl;
    theTree->SetBranchAddress("ak8pfjets_pt",            &(myEvent->pointerForak8pfjets_pt));
    theTree->SetBranchAddress("ak8pfjets_phi",           &(myEvent->pointerForak8pfjets_phi));
    theTree->SetBranchAddress("ak8pfjets_eta",           &(myEvent->pointerForak8pfjets_eta));
    theTree->SetBranchAddress("ak8pfjets_mass",          &(myEvent->pointerForak8pfjets_mass));
    theTree->SetBranchAddress("ak8pfjets_pruned_mass",   &(myEvent->pointerForak8pfjets_pruned_mass));
    theTree->SetBranchAddress("ak8pfjets_corrpruned_mass",   &(myEvent->pointerForak8pfjets_corrpruned_mass));
    theTree->SetBranchAddress("ak8pfjets_trimmed_mass",  &(myEvent->pointerForak8pfjets_trimmed_mass));
    theTree->SetBranchAddress("ak8pfjets_softdrop_mass", &(myEvent->pointerForak8pfjets_softdrop_mass));
    theTree->SetBranchAddress("ak8pfjets_nSubJets",      &(myEvent->pointerForak8pfjets_nSubJets));
    theTree->SetBranchAddress("ak8pfjets_tau1",          &(myEvent->pointerForak8pfjets_tau1));
    theTree->SetBranchAddress("ak8pfjets_tau2",          &(myEvent->pointerForak8pfjets_tau2));
    //std::cout << "reading bak8 info" << std::endl;
    #endif

    #ifdef USE_AK10_JETS
    //std::cout << "reading bak10 info" << std::endl;
    theTree->SetBranchAddress("ak10pfjets_pt",            &(myEvent->pointerForak10pfjets_pt));
    theTree->SetBranchAddress("ak10pfjets_phi",           &(myEvent->pointerForak10pfjets_phi));
    theTree->SetBranchAddress("ak10pfjets_eta",           &(myEvent->pointerForak10pfjets_eta));
    theTree->SetBranchAddress("ak10pfjets_mass",          &(myEvent->pointerForak10pfjets_mass));
    theTree->SetBranchAddress("ak10pfjets_pruned_mass",   &(myEvent->pointerForak10pfjets_pruned_mass));
    theTree->SetBranchAddress("ak10pfjets_trimmed_mass",  &(myEvent->pointerForak10pfjets_trimmed_mass));
    theTree->SetBranchAddress("ak10pfjets_softdrop_mass", &(myEvent->pointerForak10pfjets_softdrop_mass));
    theTree->SetBranchAddress("ak10pfjets_nSubJets",      &(myEvent->pointerForak10pfjets_nSubJets));
    theTree->SetBranchAddress("ak10pfjets_tau1",          &(myEvent->pointerForak10pfjets_tau1));
    theTree->SetBranchAddress("ak10pfjets_tau2",          &(myEvent->pointerForak10pfjets_tau2));
    //std::cout << "reading bak10 info" << std::endl;
    #endif

    #ifdef USE_PV
    theTree->SetBranchAddress("pv_ndof",                 &(myEvent->pv_ndof));
    theTree->SetBranchAddress("pv_isFake",               &(myEvent->pv_isFake));
    theTree->SetBranchAddress("pv_rho",                  &(myEvent->pv_rho));
    theTree->SetBranchAddress("pv_z",                    &(myEvent->pv_z));
    #endif */

    #ifdef USE_WEIGHTS
    if(theTree->GetListOfBranches()->FindObject("weight_PU"))
        theTree->SetBranchAddress("weight_PU",                  &(myEvent->weight_PU));
    if(theTree->GetListOfBranches()->FindObject("weight_PUup"))
        theTree->SetBranchAddress("weight_PUup",                  &(myEvent->weight_PUup));
    if(theTree->GetListOfBranches()->FindObject("weight_PUdown"))
        theTree->SetBranchAddress("weight_PUdown",                  &(myEvent->weight_PUdown));
    if(theTree->GetListOfBranches()->FindObject("genweights"))
        theTree->SetBranchAddress("genweights",                  &(myEvent->genweights));
    if(theTree->GetListOfBranches()->FindObject("weight_btagsf"))
        theTree->SetBranchAddress("weight_btagsf",                  &(myEvent->weight_btagsf));
    if(theTree->GetListOfBranches()->FindObject("weight_btagsf_heavy_DN"))
        theTree->SetBranchAddress("weight_btagsf_heavy_DN",                  &(myEvent->weight_btagsf_heavy_DN));
    if(theTree->GetListOfBranches()->FindObject("weight_btagsf_heavy_UP"))
        theTree->SetBranchAddress("weight_btagsf_heavy_UP",                  &(myEvent->weight_btagsf_heavy_UP));
    if(theTree->GetListOfBranches()->FindObject("weight_btagsf_light_DN"))
        theTree->SetBranchAddress("weight_btagsf_light_DN",                  &(myEvent->weight_btagsf_light_DN));
    if(theTree->GetListOfBranches()->FindObject("weight_btagsf_light_UP"))
        theTree->SetBranchAddress("weight_btagsf_light_UP",                  &(myEvent->weight_btagsf_light_UP));
    if(theTree->GetListOfBranches()->FindObject("weight_tightbtagsf"))
        theTree->SetBranchAddress("weight_tightbtagsf",                  &(myEvent->weight_tightbtagsf));
    if(theTree->GetListOfBranches()->FindObject("weight_tightbtagsf_heavy_DN"))
        theTree->SetBranchAddress("weight_tightbtagsf_heavy_DN",                  &(myEvent->weight_tightbtagsf_heavy_DN));
    if(theTree->GetListOfBranches()->FindObject("weight_tightbtagsf_heavy_UP"))
        theTree->SetBranchAddress("weight_tightbtagsf_heavy_UP",                  &(myEvent->weight_tightbtagsf_heavy_UP));
    if(theTree->GetListOfBranches()->FindObject("weight_tightbtagsf_light_DN"))
        theTree->SetBranchAddress("weight_tightbtagsf_light_DN",                  &(myEvent->weight_tightbtagsf_light_DN));
    if(theTree->GetListOfBranches()->FindObject("weight_tightbtagsf_light_UP"))
        theTree->SetBranchAddress("weight_tightbtagsf_light_UP",                  &(myEvent->weight_tightbtagsf_light_UP));
    if(theTree->GetListOfBranches()->FindObject("weight_lepSF"))
        theTree->SetBranchAddress("weight_lepSF",                  &(myEvent->weight_lepSF));
    if(theTree->GetListOfBranches()->FindObject("weight_lepSF_down"))
        theTree->SetBranchAddress("weight_lepSF_down",                  &(myEvent->weight_lepSF_down));
    if(theTree->GetListOfBranches()->FindObject("weight_lepSF_up"))
        theTree->SetBranchAddress("weight_lepSF_up",                  &(myEvent->weight_lepSF_up));
    if(theTree->GetListOfBranches()->FindObject("weight_vetoLepSF"))
        theTree->SetBranchAddress("weight_vetoLepSF",                  &(myEvent->weight_vetoLepSF));
    if(theTree->GetListOfBranches()->FindObject("weight_lepSF_fastSim"))
        theTree->SetBranchAddress("weight_lepSF_fastSim",                  &(myEvent->weight_lepSF_fastSim));
    if(theTree->GetListOfBranches()->FindObject("weight_ISR"))
        theTree->SetBranchAddress("weight_ISR",                  &(myEvent->weight_ISR));
    if(theTree->GetListOfBranches()->FindObject("weight_ISRdown"))
        theTree->SetBranchAddress("weight_ISRdown",                  &(myEvent->weight_ISRdown));
    if(theTree->GetListOfBranches()->FindObject("weight_ISRup"))
        theTree->SetBranchAddress("weight_ISRup",                  &(myEvent->weight_ISRup));
    if(theTree->GetListOfBranches()->FindObject("weight_ISRnjets"))
        theTree->SetBranchAddress("weight_ISRnjets",                  &(myEvent->weight_ISRnjets));
    if(theTree->GetListOfBranches()->FindObject("weight_ISRnjets_DN"))
        theTree->SetBranchAddress("weight_ISRnjets_DN",                  &(myEvent->weight_ISRnjets_DN));
    if(theTree->GetListOfBranches()->FindObject("weight_ISRnjets_UP"))
        theTree->SetBranchAddress("weight_ISRnjets_UP",                  &(myEvent->weight_ISRnjets_UP));
    if(theTree->GetListOfBranches()->FindObject("pdf_down_weight"))
        theTree->SetBranchAddress("pdf_down_weight",                  &(myEvent->pdf_down_weight));
    if(theTree->GetListOfBranches()->FindObject("pdf_up_weight"))
        theTree->SetBranchAddress("pdf_up_weight",                  &(myEvent->pdf_up_weight));
    #endif

    #ifdef USE_VAR_BASELINE
    //theTree->SetBranchAddress("met_sig",                 &(myEvent->met_sig));
    if(theTree->GetListOfBranches()->FindObject("pfmet"))
        theTree->SetBranchAddress("pfmet",                   &(myEvent->pfmet));
    if(theTree->GetListOfBranches()->FindObject("pfmet_rl"))
        theTree->SetBranchAddress("pfmet_rl",                   &(myEvent->pfmet_rl));
    if(theTree->GetListOfBranches()->FindObject("MT2W"))
        theTree->SetBranchAddress("MT2W",                    &(myEvent->MT2W));
    if(theTree->GetListOfBranches()->FindObject("MT2W_rl"))
        theTree->SetBranchAddress("MT2W_rl",                    &(myEvent->MT2W_rl));
    if(theTree->GetListOfBranches()->FindObject("mt_met_lep"))
        theTree->SetBranchAddress("mt_met_lep",              &(myEvent->mt_met_lep));
    if(theTree->GetListOfBranches()->FindObject("mt_met_lep_rl"))
        theTree->SetBranchAddress("mt_met_lep_rl",              &(myEvent->mt_met_lep_rl));
    if(theTree->GetListOfBranches()->FindObject("topnessMod"))
        theTree->SetBranchAddress("topnessMod",                 &(myEvent->topnessMod));
    if(theTree->GetListOfBranches()->FindObject("topnessMod_rl"))
        theTree->SetBranchAddress("topnessMod_rl",                 &(myEvent->topnessMod_rl));
    if(theTree->GetListOfBranches()->FindObject("ngoodjets"))
        theTree->SetBranchAddress("ngoodjets",               &(myEvent->ngoodjets));
    if(theTree->GetListOfBranches()->FindObject("ngoodbtags"))
        theTree->SetBranchAddress("ngoodbtags",              &(myEvent->ngoodbtags));
    if(theTree->GetListOfBranches()->FindObject("ntightbtags"))
        theTree->SetBranchAddress("ntightbtags",              &(myEvent->ntightbtags));
    if(theTree->GetListOfBranches()->FindObject("mindphi_met_j1_j2"))
        theTree->SetBranchAddress("mindphi_met_j1_j2",      &(myEvent->dphi_ak4pfjets_met));
    if(theTree->GetListOfBranches()->FindObject("mindphi_met_j1_j2_rl"))
        theTree->SetBranchAddress("mindphi_met_j1_j2_rl",      &(myEvent->dphi_ak4pfjets_met_rl));
    if(theTree->GetListOfBranches()->FindObject("lep1_dphiMET"))
        theTree->SetBranchAddress("lep1_dphiMET",      &(myEvent->lep1_dphiMET));
    if(theTree->GetListOfBranches()->FindObject("lep1_dphiMET_rl"))
        theTree->SetBranchAddress("lep1_dphiMET_rl",      &(myEvent->lep1_dphiMET_rl));
    if(theTree->GetListOfBranches()->FindObject("lep2_dphiMET"))
        theTree->SetBranchAddress("lep2_dphiMET",      &(myEvent->lep2_dphiMET));
    if(theTree->GetListOfBranches()->FindObject("lep2_dphiMET_rl"))
        theTree->SetBranchAddress("lep2_dphiMET_rl",      &(myEvent->lep2_dphiMET_rl));
    if(theTree->GetListOfBranches()->FindObject("pfmet_phi"))
        theTree->SetBranchAddress("pfmet_phi",               &(myEvent->pfmet_phi));
    if(theTree->GetListOfBranches()->FindObject("pfmet_phi_rl"))
        theTree->SetBranchAddress("pfmet_phi_rl",               &(myEvent->pfmet_phi_rl));
    if(theTree->GetListOfBranches()->FindObject("Mlb_closestb"))
        theTree->SetBranchAddress("Mlb_closestb",                     &(myEvent->Mlb));
    #endif


    //@MJ@ different define for j_up/down

    #ifdef USE_VAR_BASELINE_UP
    if(theTree->GetListOfBranches()->FindObject("pfmet_rl_jup"))
        theTree->SetBranchAddress("pfmet_rl_jup",                   &(myEvent->pfmet));
    if(theTree->GetListOfBranches()->FindObject("MT2W_rl_jup"))
        theTree->SetBranchAddress("MT2W_rl_jup",                    &(myEvent->MT2W));
    if(theTree->GetListOfBranches()->FindObject("mt_met_lep_rl_jup"))
        theTree->SetBranchAddress("mt_met_lep_rl_jup",              &(myEvent->mt_met_lep));
    if(theTree->GetListOfBranches()->FindObject("topnessMod_rl_jup"))
        theTree->SetBranchAddress("topnessMod_rl_jup",                 &(myEvent->topnessMod));
    if(theTree->GetListOfBranches()->FindObject("jup_ngoodjets"))
        theTree->SetBranchAddress("jup_ngoodjets",               &(myEvent->ngoodjets));
    if(theTree->GetListOfBranches()->FindObject("jup_ngoodbtags"))
        theTree->SetBranchAddress("jup_ngoodbtags",              &(myEvent->ngoodbtags));
    if(theTree->GetListOfBranches()->FindObject("jup_ntightbtags"))
        theTree->SetBranchAddress("jup_ntightbtags",              &(myEvent->ntightbtags));
    if(theTree->GetListOfBranches()->FindObject("mindphi_met_j1_j2_rl_jup"))
        theTree->SetBranchAddress("mindphi_met_j1_j2_rl_jup",      &(myEvent->dphi_ak4pfjets_met));
    if(theTree->GetListOfBranches()->FindObject("lep1_dphiMET_rl_jup"))
        theTree->SetBranchAddress("lep1_dphiMET_rl_jup",      &(myEvent->lep1_dphiMET));
    if(theTree->GetListOfBranches()->FindObject("pfmet_phi_rl_jup"))
        theTree->SetBranchAddress("pfmet_phi_rl_jup",               &(myEvent->pfmet_phi));
    if(theTree->GetListOfBranches()->FindObject("Mlb_closestb_jup"))
        theTree->SetBranchAddress("Mlb_closestb_jup",                     &(myEvent->Mlb));
    #endif

    #ifdef USE_VAR_BASELINE_DOWN
    if(theTree->GetListOfBranches()->FindObject("pfmet_rl_jdown"))
        theTree->SetBranchAddress("pfmet_rl_jdown",                   &(myEvent->pfmet));
    if(theTree->GetListOfBranches()->FindObject("MT2W_rl_jdown"))
        theTree->SetBranchAddress("MT2W_rl_jdown",                    &(myEvent->MT2W));
    if(theTree->GetListOfBranches()->FindObject("mt_met_lep_rl_jdown"))
        theTree->SetBranchAddress("mt_met_lep_rl_jdown",              &(myEvent->mt_met_lep));
    if(theTree->GetListOfBranches()->FindObject("topnessMod_rl_jdown"))
        theTree->SetBranchAddress("topnessMod_rl_jdown",                 &(myEvent->topnessMod));
    if(theTree->GetListOfBranches()->FindObject("jdown_ngoodjets"))
        theTree->SetBranchAddress("jdown_ngoodjets",               &(myEvent->ngoodjets));
    if(theTree->GetListOfBranches()->FindObject("jdown_ngoodbtags"))
        theTree->SetBranchAddress("jdown_ngoodbtags",              &(myEvent->ngoodbtags));
    if(theTree->GetListOfBranches()->FindObject("jdown_ntightbtags"))
        theTree->SetBranchAddress("jdown_ntightbtags",              &(myEvent->ntightbtags));
    if(theTree->GetListOfBranches()->FindObject("mindphi_met_j1_j2_rl_jdown"))
        theTree->SetBranchAddress("mindphi_met_j1_j2_rl_jdown",      &(myEvent->dphi_ak4pfjets_met));
    if(theTree->GetListOfBranches()->FindObject("lep1_dphiMET_rl_jdown"))
        theTree->SetBranchAddress("lep1_dphiMET_rl_jdown",      &(myEvent->lep1_dphiMET));
    if(theTree->GetListOfBranches()->FindObject("pfmet_phi_rl_jdown"))
        theTree->SetBranchAddress("pfmet_phi_rl_jdown",               &(myEvent->pfmet_phi));
    if(theTree->GetListOfBranches()->FindObject("Mlb_closestb_jdown"))
        theTree->SetBranchAddress("Mlb_closestb_jdown",                     &(myEvent->Mlb));
    #endif

    #ifdef USE_GLOBAL_VAR
    //theTree->SetBranchAddress("Mjjj",                    &(myEvent->Mjjj));
    if(theTree->GetListOfBranches()->FindObject("hadronic_top_chi2"))
        theTree->SetBranchAddress("hadronic_top_chi2",       &(myEvent->hadronic_top_chi2));
    if(theTree->GetListOfBranches()->FindObject("nvtxs"))
    {
        theTree->SetBranchAddress("nvtxs",                 &(myEvent->nvertex));
        //cout << " nvertex filled by "  << myEvent->nvertex << endl;
    }
    //theTree->SetBranchAddress("puIntime",                &(myEvent->puIntime));
    if(theTree->GetListOfBranches()->FindObject("pu_ntrue"))
        theTree->SetBranchAddress("pu_ntrue",                  &(myEvent->puTrue));
    if(theTree->GetListOfBranches()->FindObject("mass_lsp"))
        theTree->SetBranchAddress("mass_lsp",                  &(myEvent->mass_lsp));
    if(theTree->GetListOfBranches()->FindObject("mass_chargino"))
        theTree->SetBranchAddress("mass_chargino",                  &(myEvent->mass_chargino));
    if(theTree->GetListOfBranches()->FindObject("mass_stop"))
        theTree->SetBranchAddress("mass_stop",                  &(myEvent->mass_stop));
    //theTree->SetBranchAddress("Mlb_leadb",               &(myEvent->Mlb_leadb));
    if(theTree->GetListOfBranches()->FindObject("topness"))
        theTree->SetBranchAddress("topness",                 &(myEvent->topness));
    //theTree->SetBranchAddress("leadingLeptonId",         &(myEvent->leadingLeptonId));
    //theTree->SetBranchAddress("numberOfGeneratedLeptons",&(myEvent->numberOfGeneratedLeptons));
    if(theTree->GetListOfBranches()->FindObject("nupt"))
        theTree->SetBranchAddress("nupt",       &(myEvent->nupt));
    #endif
    //generated info
    //
    if(theTree->GetListOfBranches()->FindObject("genqs_p4"))
        theTree->SetBranchAddress("genqs_p4",                 &(myEvent->p4genqs_p4));
    if(theTree->GetListOfBranches()->FindObject("genqs_id"))
        theTree->SetBranchAddress("genqs_id",                 &(myEvent->genqs_id));
    if(theTree->GetListOfBranches()->FindObject("genqs_motherid"))
        theTree->SetBranchAddress("genqs_motherid",                 &(myEvent->genqs_motherid));
    if(theTree->GetListOfBranches()->FindObject("genqs_motheridx"))
        theTree->SetBranchAddress("genqs_motheridx",                 &(myEvent->genqs_motheridx));
    if(theTree->GetListOfBranches()->FindObject("genqs_status"))
        theTree->SetBranchAddress("genqs_status",                 &(myEvent->genqs_status));
    #ifdef USE_NEW_VAR
    //theTree->SetBranchAddress("ST"		,&(myEvent->ST));
    //theTree->SetBranchAddress("LP"		,&(myEvent->LP));
    //theTree->SetBranchAddress("Meff"		,&(myEvent->Meff));
    //theTree->SetBranchAddress("MTdeco_Q"	,&(myEvent->MTdeco_Q));
    //theTree->SetBranchAddress("DeltaPtbb"	,&(myEvent->DeltaPtbb));
    //theTree->SetBranchAddress("DeltaPhibb"	,&(myEvent->DeltaPhibb));
    //theTree->SetBranchAddress("DeltaRbb"	,&(myEvent->DeltaRbb));
    //theTree->SetBranchAddress("dphi_Wlep",	&(myEvent->dphi_Wlep));
    //theTree->SetBranchAddress("dR_lep_leadb",	&(myEvent->dR_lep_leadb));
    
    if(theTree->GetListOfBranches()->FindObject("ak4_HT"))
        theTree->SetBranchAddress("ak4_HT",		&(myEvent->ak4_HT));
    //theTree->SetBranchAddress("ak4_htssm",	&(myEvent->ak4_htssm));
    //theTree->SetBranchAddress("ak4_htosm",	&(myEvent->ak4_htosm));
    #endif

    /*#ifdef USE_GEN_LOSTLEPTON
    theTree->SetBranchAddress("genLostLeptons_pt",	&(myEvent->pointerForGenLostLeptons_pt));
    theTree->SetBranchAddress("genLostLeptons_eta",	&(myEvent->pointerForGenLostLeptons_eta));
    theTree->SetBranchAddress("genLostLeptons_pdgid",   &(myEvent->pointerForGenLostLeptons_pdgid));
    //cout<<"Branch adress done !"<<endl;
    #endif*/

    #ifdef USE_OLD_VAR
    if(theTree->GetListOfBranches()->FindObject("run"))
        theTree->SetBranchAddress("run",                   &(myEvent->runId));
    
    //theTree->SetBranchAddress("numberOfSelectedLeptons", &(myEvent->numberOfSelectedLeptons));
    //theTree->SetBranchAddress("numberOfSelectedElectrons", &(myEvent->numberOfSelectedElectrons));
    //theTree->SetBranchAddress("numberOfSelectedMuons",   &(myEvent->numberOfSelectedMuons));
    //theTree->SetBranchAddress("numberOfBTaggedJets",     &(myEvent->numberOfBTaggedJets));
    //theTree->SetBranchAddress("numberOfSelectedJets",    &(myEvent->numberOfSelectedJets));
    
    //theTree->SetBranchAddress("jetsPt",                  &(myEvent->jetsPt));
    //theTree->SetBranchAddress("jetsEta",                 &(myEvent->jetsEta));
    //theTree->SetBranchAddress("jetsPhi",                 &(myEvent->jetsPhi));
    //theTree->SetBranchAddress("jetsCSV",                 &(myEvent->jetsCSV));
    //theTree->SetBranchAddress("jetsCSVv2",               &(myEvent->jetsCSVv2));
    //theTree->SetBranchAddress("jetsPUid",                &(myEvent->jetsPUid));
    
    if(theTree->GetListOfBranches()->FindObject("evt"))
        theTree->SetBranchAddress("evt",                 &(myEvent->eventId));
    if(theTree->GetListOfBranches()->FindObject("ls"))
        theTree->SetBranchAddress("ls",                  &(myEvent->lumiId));
    //theTree->SetBranchAddress("selectionCode",           &(myEvent->selectionCode));
    
    if(theTree->GetListOfBranches()->FindObject("ak4_HT"))
        theTree->SetBranchAddress("ak4_HT",                      &(myEvent->HT));
    //theTree->SetBranchAddress("MT",                      &(myEvent->MT));
    //theTree->SetBranchAddress("ETmiss",                  &(myEvent->ETmiss));
    //theTree->SetBranchAddress("ETmissPhi",               &(myEvent->ETmissPhi));
    //theTree->SetBranchAddress("chi2",                    &(myEvent->chi2));
    //theTree->SetBranchAddress("M3b",                     &(myEvent->M3b));
    #endif

    #ifdef USE_EVT_VAR
    if(theTree->GetListOfBranches()->FindObject("run"))
        theTree->SetBranchAddress("run",                   &(myEvent->runId));
    if(theTree->GetListOfBranches()->FindObject("evt"))
        theTree->SetBranchAddress("evt",                 &(myEvent->eventId));
    if(theTree->GetListOfBranches()->FindObject("ls"))
        theTree->SetBranchAddress("ls",                  &(myEvent->lumiId));
    #endif
   
   //*/
    /*#ifdef USE_GEN_INFO
    myEvent->hasGenInfo = false;
    try{
    	if(theTree->GetListOfBranches()->FindObject("gen_pt")!=0){
    		myEvent->hasGenInfo = true;
    	}
    }
    catch(...){
    	myEvent->hasGenInfo = false;
    }
    //if(theTree->GetListOfBranches()->FindObject("gen_pt")!=0){
    //myEvent->hasGenInfo = true;
    if(myEvent->hasGenInfo){
    theTree->SetBranchAddress("gen_pt", &(myEvent->pointerForgen_pt));
    theTree->SetBranchAddress("gen_eta",  &(myEvent->pointerForgen_eta));
    theTree->SetBranchAddress("gen_phi",  &(myEvent->pointerForgen_phi));
    theTree->SetBranchAddress("gen_m", &(myEvent->pointerForgen_m));
    theTree->SetBranchAddress("gen_status", &(myEvent->pointerForgen_status));
    theTree->SetBranchAddress("gen_id",  &(myEvent->pointerForgen_id));
    theTree->SetBranchAddress("gen_daughter_n",  &(myEvent->pointerForgen_daughter_n));
    theTree->SetBranchAddress("gen_daughter_index", &(myEvent->pointerForgen_daughter_index));
    }
    #endif
    #ifdef USE_GEN_INFO_EXT
    theTree->SetBranchAddress("gen_mother_index",  &(myEvent->pointerForgen_mother_index));
    theTree->SetBranchAddress("gen_charge",  &(myEvent->pointerForgen_charge));
    theTree->SetBranchAddress("gen_index",  &(myEvent->pointerForgen_index));
    #endif*/
    
    
    #ifdef USE_SKIMMING_VAR
    //myEvent->pointerForisr_id = 0;
    //myEvent->pointerForisr_p4 = 0;
    //myEvent->pointerForneutrinos = 0;
    //myEvent->pointerForneutralinos = 0;
    if(theTree->GetListOfBranches()->FindObject("mindphi_met_j1_j2"))
        theTree->SetBranchAddress("mindphi_met_j1_j2", &myEvent->minDPhi_jmet); 
    if(theTree->GetListOfBranches()->FindObject("genmet"))
        theTree->SetBranchAddress("genmet", &myEvent->genMET); 
    //theTree->SetBranchAddress("Pz_ttbar", &myEvent->Pz_ttbar); 
    //theTree->SetBranchAddress("Pt_ttbar", &myEvent->Pt_ttbar); 
    //theTree->SetBranchAddress("E_ttbar", &myEvent->E_ttbar); 
    //theTree->SetBranchAddress("M_ttbar", &myEvent->M_ttbar); 
    //theTree->SetBranchAddress("ttbar_decay", &myEvent->ttbar_decay); 
    //theTree->SetBranchAddress("nofISR", &myEvent->nofISR); 
    //theTree->SetBranchAddress("index_b", &myEvent->index_b); 
    //theTree->SetBranchAddress("index_bbar", &myEvent->index_bbar); 
    //theTree->SetBranchAddress("matched", &myEvent->matched); 
    //theTree->SetBranchAddress("isr_id", &(myEvent->pointerForisr_id)); 
    //theTree->SetBranchAddress("isr_p4", &myEvent->pointerForisr_p4); 
    //theTree->SetBranchAddress("neutralinos", &myEvent->pointerForneutralinos); 
    //theTree->SetBranchAddress("neutrinos", &myEvent->pointerForneutrinos); 
    #endif
}



// ################################
// #  Function to read one event  #
// ################################
 
void ReadEvent(TTree* theTree, long int i, babyEvent* myEvent, TFile* f)
{
    theTree->GetEntry(i);
    //myEvent->nthentry = i;
    //
    // cout << "new evennt" << endl;

    if(myEvent->genweights->size()>109)
    {
        myEvent->weight_q2_down = myEvent->genweights->at(8)/myEvent->genweights->at(0);
        myEvent->weight_q2_up = myEvent->genweights->at(4)/myEvent->genweights->at(0);
        myEvent->weight_alphas_down = myEvent->genweights->at(110)/myEvent->genweights->at(0);
        myEvent->weight_alphas_up = myEvent->genweights->at(109)/myEvent->genweights->at(0);
    }
/*
    else
        cout << "gen weights not properly set, the size is: " << myEvent->genweights->size() << endl;
*/

    #ifdef USE_LEP1
    myEvent->lep1_p4 = *(myEvent->p4lep1_p4);
    myEvent->lep1_pt = myEvent->lep1_p4.Pt();
    myEvent->lep1_eta = myEvent->lep1_p4.Eta();
    myEvent->lep1_phi = myEvent->lep1_p4.Phi();
    myEvent->lep1_mass = sqrt(myEvent->lep1_p4.M2());
    #endif

    #ifdef USE_LEP2
    myEvent->lep2_p4 = *(myEvent->p4lep2_p4);
    myEvent->lep2_pt = myEvent->lep2_p4.Pt();
    myEvent->lep2_eta = myEvent->lep2_p4.Eta();
    myEvent->lep2_phi = myEvent->lep2_p4.Phi();
    myEvent->lep2_mass = sqrt(myEvent->lep2_p4.M2());
    #endif

    #ifdef USE_TAU
    myEvent->tau_p4= *(myEvent->p4tau_p4);
    #endif

     #ifdef USE_GEN_LEP
     myEvent->genleps_p4= *(myEvent->p4genleps_p4);
     myEvent->genleps_id= *(myEvent->pointer_genleps_id);
     myEvent->genleps_motherid= *(myEvent->pointer_genleps_motherid);
     myEvent->genleps_gmotherid= *(myEvent->pointer_genleps_gmotherid);
     myEvent->genleps_isLastCopy= *(myEvent->pointer_genleps_isLastCopy);
     myEvent->genleps_fromHardProcessFinalState= *(myEvent->pointer_genleps_fromHardProcessFinalState);
     myEvent->genleps_fromHardProcessDecayed= *(myEvent->pointer_genleps_fromHardProcessDecayed);
     myEvent->genleps_gentaudecay= *(myEvent->pointer_genleps_gentaudecay);    
      #endif

    //myEvent->lep1_dphiMET = abs(myEvent->lep1_phi - myEvent->pfmet_phi) < TMath::Pi()  ? abs(myEvent->lep1_phi - myEvent->pfmet_phi) : 2*TMath::Pi() - abs(myEvent->lep1_phi - myEvent->pfmet_phi) ; //@MJ@ TODO remove

    myEvent->genqs_p4 = *(myEvent->p4genqs_p4);
    for(uint32_t t = 0; t<myEvent->genqs_id->size(); t++)
    {
        if(myEvent->genqs_id->at(t) == 6)
        {
            //cout << "top " << myEvent->genqs_status->at(t) << " pt " << myEvent->genqs_p4.at(t).Pt() << endl;
            if(myEvent->genqs_status->at(t) == 22)
                myEvent->top_pt = myEvent->genqs_p4.at(t).Pt();
        }
        else if(myEvent->genqs_id->at(t) == -6)
        {
            //cout << "atop " << myEvent->genqs_status->at(t) << " pt " << myEvent->genqs_p4.at(t).Pt() << endl;
            if(myEvent->genqs_status->at(t) == 22)
                myEvent->atop_pt = myEvent->genqs_p4.at(t).Pt();
        }

    }
    //jets 4-vector

 

    #ifdef USE_JETS
    myEvent->ak4pfjets_leadbtag_p4 = *(myEvent->p4ak4pfjets_leadbtag_p4);
    myEvent->ak4pfjets_p4 = *(myEvent->p4ak4pfjets_p4);
    //  cout<< myEvent->ak4pfjets_p4.size()<<endl;
    // for (uint32_t t = 0; t<myEvent->ak4pfjets_p4.size(); t++)
    if (myEvent->ak4pfjets_p4.size()>0)
      {
	myEvent->ak4pfjets_pt =myEvent->ak4pfjets_p4.at(0).Pt();
	myEvent->ak4pfjets_leadJet_identification = myEvent->ak4pfjets_CSV->at(0);
	  
	if (myEvent->ak4pfjets_p4.size()>1)
	  {
	    myEvent->ak4pfjets_pt2 =myEvent->ak4pfjets_p4.at(1).Pt();
	    if (myEvent->ak4pfjets_p4.size()>2)
	      {
		myEvent->ak4pfjets_pt3 =myEvent->ak4pfjets_p4.at(2).Pt();
		if (myEvent->ak4pfjets_p4.size()>3)
		  {
		    myEvent->ak4pfjets_pt4 =myEvent->ak4pfjets_p4.at(3).Pt();
		    
		  }
		
	      }
	    
	    
	  }
  
      }
    #endif 

    /*  for (uint32_t t = 0; t<myEvent->ak4pfjets_p4.size(); t++)
      {
	if (myEvent->ak4pfjets_CSV.ak4pfjets_CSV < 0.8)
	  {
	    myEvent->ak4pfjets_pt_B =myEvent->ak4pfjets_p4.at(t).Pt();
	  }
	  }*/
    #endif 

    //fill weight normalization
    
    TH1D* hWeight = NULL;
    hWeight = dynamic_cast<TH1D*>(f->Get("h_counter"));
    if(hWeight==NULL)
    {
        TH3D* hSigWeight = NULL;
            hSigWeight = dynamic_cast<TH3D*>(f->Get("h_counterSMS"));
        if( hSigWeight == NULL)
            throw std::runtime_error("weights cannot be filled");

        TH2F* h2Dplane = dynamic_cast<TH2F*>(f->Get("histNEvts"));
        TH1D* projX = h2Dplane->ProjectionX();
        TH1D* projY = h2Dplane->ProjectionY();
        int binNX = projX->FindBin(myEvent->mass_stop);
        int binNY = projY->FindBin(myEvent->mass_lsp);
        //cout << "binNx " << binNX << " binNY " << binNY << endl;
	//cout << "weight hist " << hSigWeight << " get " <<  f->Get("h_counter")<< " file " << f << endl;
	myEvent->wNormalization.clear();
	myEvent->wNormalization.push_back(-13); //@MJ@ ATTENTION vector.at(0) random number to sy synchronized with binning numbering
	for(uint32_t n=0; n<hSigWeight->GetNbinsZ(); n++)
	{
            myEvent->wNormalization.push_back(hSigWeight->GetBinContent(binNX,binNY, n+1));
	}
	//cout << "normalization size " << myEvent->wNormalization.size() << "first element " << myEvent->wNormalization.at(1) << endl;
    }
    
    // Put actual content of special type branches where they should be...
    ///*

    /*
    myEvent->ak4pfjets_puid            = *(myEvent->pointerForak4pfjets_puid);
    cout<<myEvent->pointerForak4pfjets_pt<<endl;
    cout<<(*myEvent->pointerForak4pfjets_pt).size()<<endl;
    cout<<(*myEvent->pointerForak4pfjets_pt)[0]<<endl;
    vector<float> pts = *(myEvent->pointerForak4pfjets_pt);
*/    
//myEvent->pts = *(myEvent->pointerForak4pfjets_pt);
    //myEvent->ak4pfjets_pt.clear();
    //myEvent->ak4pfjets_pt = pts;
    //myEvent->ak4pfjets_phi              = *(myEvent->pointerForak4pfjets_pt);
    //myEvent->ak4pfjets_phi             = *(myEvent->pointerForak4pfjets_phi);
    //myEvent->ak4pfjets_eta             = *(myEvent->pointerForak4pfjets_eta);
    //myEvent->ak4pfjets_mass            = *(myEvent->pointerForak4pfjets_mass);
    //myEvent->ak4pfjets_loose_pfid      = *(myEvent->pointerForak4pfjets_loose_pfid);
    //myEvent->ak4pfjets_CSV             = *(myEvent->pointerForak4pfjets_CSV);
   
    //

    //myEvent->trigger_name        = *(myEvent->pointerForTriggerName);
    //myEvent->trigger_pass        = *(myEvent->pointerForTriggerPass);
    //protect against branches not found
    /*if(myEvent->pointerForgen_neutralino_m!=0 && myEvent->pointerForgen_stop_m!=0){
      myEvent->gen_neutralino_m	         =                      *(myEvent->pointerForgen_neutralino_m);
      myEvent->gen_stop_m			 =                      *(myEvent->pointerForgen_stop_m);
    }
    if(myEvent->pointerForgenTops_pt!=0 && myEvent->pointerForgenTops_pdgid!=0){
    	myEvent->genTops_pt = *(myEvent->pointerForgenTops_pt);
    	myEvent->genTops_pdgid = *(myEvent->pointerForgenTops_pdgid);
    }*/
    /*
    #ifdef USE_JETS
    //std::cout << "setting jet info" << std::endl;
    myEvent->jet_pt              = *(myEvent->pointerForak4pfjets_pt);
    myEvent->jet_phi             = *(myEvent->pointerForak4pfjets_phi);
    myEvent->jet_eta             = *(myEvent->pointerForak4pfjets_eta);
    myEvent->jet_mass            = *(myEvent->pointerForak4pfjets_mass);
    myEvent->jet_CSV             = *(myEvent->pointerForak4pfjets_CSV);
    myEvent->ak4pfjets_partonFlavour = *(myEvent->pointerForak4pfjets_partonFlavour);
    myEvent->ak4pfjets_hadronFlavour = *(myEvent->pointerForak4pfjets_hadronFlavour);
    myEvent->ak4pfjets_qgtag     = *(myEvent->pointerForak4pfjets_qgtag);
    myEvent->ak4pfjets_pt              = *(myEvent->pointerForak4pfjets_pt);
    myEvent->ak4pfjets_eta              = *(myEvent->pointerForak4pfjets_eta);
    myEvent->ak4pfjets_CSV              = *(myEvent->pointerForak4pfjets_CSV);
    myEvent->ak4pfjets_phi             = *(myEvent->pointerForak4pfjets_phi);
    #endif
    #ifdef USE_JETS_EXT
    myEvent->jet_puid      	= *(myEvent->pointerForak4pfjets_puid);
    myEvent->jet_loose_pfid      = *(myEvent->pointerForak4pfjets_loose_pfid);
    #endif
    #ifdef USE_AK8_JETS
    //std::cout << "setting ak8 jet info" << std::endl;
    myEvent->ak8pfjets_pt              = *(myEvent->pointerForak8pfjets_pt);
    myEvent->ak8pfjets_phi             = *(myEvent->pointerForak8pfjets_phi);
    myEvent->ak8pfjets_eta             = *(myEvent->pointerForak8pfjets_eta);
    myEvent->ak8pfjets_mass            = *(myEvent->pointerForak8pfjets_mass);
    myEvent->ak8pfjets_pruned_mass     = *(myEvent->pointerForak8pfjets_pruned_mass);
    myEvent->ak8pfjets_corrpruned_mass = *(myEvent->pointerForak8pfjets_corrpruned_mass);
    myEvent->ak8pfjets_trimmed_mass    = *(myEvent->pointerForak8pfjets_trimmed_mass);
    myEvent->ak8pfjets_softdrop_mass   = *(myEvent->pointerForak8pfjets_softdrop_mass);
    myEvent->ak8pfjets_nSubJets        = *(myEvent->pointerForak8pfjets_nSubJets);
    myEvent->ak8pfjets_tau1            = *(myEvent->pointerForak8pfjets_tau1);
    myEvent->ak8pfjets_tau2            = *(myEvent->pointerForak8pfjets_tau2);
    //std::cout << "setting ak8 jet info" << std::endl;
    #endif
    #ifdef USE_AK10_JETS
    //std::cout << "setting ak10 jet info" << std::endl;
    myEvent->ak10pfjets_pt              = *(myEvent->pointerForak10pfjets_pt);
    myEvent->ak10pfjets_phi             = *(myEvent->pointerForak10pfjets_phi);
    myEvent->ak10pfjets_eta             = *(myEvent->pointerForak10pfjets_eta);
    myEvent->ak10pfjets_mass            = *(myEvent->pointerForak10pfjets_mass);
    myEvent->ak10pfjets_pruned_mass     = *(myEvent->pointerForak10pfjets_pruned_mass);
    myEvent->ak10pfjets_trimmed_mass    = *(myEvent->pointerForak10pfjets_trimmed_mass);
    myEvent->ak10pfjets_softdrop_mass   = *(myEvent->pointerForak10pfjets_softdrop_mass);
    myEvent->ak10pfjets_nSubJets        = *(myEvent->pointerForak10pfjets_nSubJets);
    myEvent->ak10pfjets_tau1            = *(myEvent->pointerForak10pfjets_tau1);
    myEvent->ak10pfjets_tau2            = *(myEvent->pointerForak10pfjets_tau2);
    //std::cout << "setting ak10 jet info" << std::endl;
    #endif
    
    #ifdef USE_GEN_INFO
    //std::cout << "setting gen info" << std::endl;
    myEvent->gen_pt			 =			*(myEvent->pointerForgen_pt);            	 
    myEvent->gen_eta			 =                      *(myEvent->pointerForgen_eta);
    myEvent->gen_phi			 =                      *(myEvent->pointerForgen_phi);
    myEvent->gen_m			 =                      *(myEvent->pointerForgen_m);
    myEvent->gen_status			 =                      *(myEvent->pointerForgen_status);
    myEvent->gen_id			 =                      *(myEvent->pointerForgen_id);
    myEvent->gen_daughter_n			 =              *(myEvent->pointerForgen_daughter_n);
    myEvent->gen_daughter_index			 =      *(myEvent->pointerForgen_daughter_index);
    //std::cout << "setting gen info" << std::endl;*/
    //@MJ@ TODO try this!!
    /*if(myEvent->hasGenInfo){
	    myEvent->gen_pt			 =			*(myEvent->pointerForgen_pt);            	 
	    myEvent->gen_eta			 =                      *(myEvent->pointerForgen_eta);
 	    myEvent->gen_phi			 =                      *(myEvent->pointerForgen_phi);
  	    myEvent->gen_m			 =                      *(myEvent->pointerForgen_m);
   	    myEvent->gen_status			 =                      *(myEvent->pointerForgen_status);
   	    myEvent->gen_id			 =                      *(myEvent->pointerForgen_id);
 	    myEvent->gen_daughter_n			 =              *(myEvent->pointerForgen_daughter_n);
   	    myEvent->gen_daughter_index			 =      *(myEvent->pointerForgen_daughter_index);
    }*/
    /*#endif
    #ifdef USE_GEN_INFO_EXT
    myEvent->gen_charge			 =                      *(myEvent->pointerForgen_charge);
    myEvent->gen_index			 =                      *(myEvent->pointerForgen_index);
    myEvent->gen_mother_index			 =              *(myEvent->pointerForgen_mother_index);
    #endif
    
    #ifdef USE_GEN_LOSTLEPTON
    //protect against branches not found
    if(myEvent->pointerForGenLostLeptons_pt!=0 && myEvent->pointerForGenLostLeptons_eta!=0 && myEvent->pointerForGenLostLeptons_pdgid!=0){
      myEvent->genLostLeptons_pt = *myEvent->pointerForGenLostLeptons_pt;
      myEvent->genLostLeptons_eta = *myEvent->pointerForGenLostLeptons_eta;
      myEvent->genLostLeptons_pdgid = *myEvent->pointerForGenLostLeptons_pdgid;
    }
    #endif

    //Fill temporary info
    myEvent->leadingLepton.SetPtEtaPhiM(myEvent->lep1_pt, myEvent->lep1_eta, myEvent->lep1_phi, myEvent->lep1_mass);
    myEvent->secondLepton.SetPtEtaPhiM(myEvent->lep2_pt, myEvent->lep2_eta, myEvent->lep2_phi, myEvent->lep2_mass);
*/
}

/*void ReadEvent(TTree* theTree, uint64_t j, babyEvent* myEvent)
{
    //intermediate variables
    TLorentzVector ak4pfjets_p4;
    bool   one_ak4pfjets_passMEDbtag;
    float  one_ak4pfjets_mva;
    float  one_ak4pfjets_loose_puid;
    int    one_ak4pfjets_partonFlavour;
    int    one_ak4pfjets_hadronFlavour;
    bool   one_ak4pfjets_loose_pfid;
    float  one_ak4pfjets_CSV;

   
    theTree->GetEntry(j);
    #ifdef USE_JETS
    if(j==0)
    {
        myEvent->ak4pfjets_loose_puid.clear();
        myEvent->ak4pfjets_passMEDbtag.clear();
        myEvent->ak4pfjets_mva.clear();
        myEvent->ak4pfjets_pt.clear();
        myEvent->ak4pfjets_phi.clear();
        myEvent->ak4pfjets_eta.clear();
        myEvent->ak4pfjets_mass.clear();
        myEvent->ak4pfjets_partonFlavour.clear();
        myEvent->ak4pfjets_hadronFlavour.clear();
        myEvent->ak4pfjets_loose_pfid.clear();
        myEvent->ak4pfjets_CSV.clear();
    }

    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_p4"))
    {
        theTree->SetBranchAddress("ak4pfjets_p4",              &(ak4pfjets_p4));
        myEvent->ak4pfjets_p4ak4pfjets_pt.push_back(ak4pfjets_p4.Pt());
        myEvent->ak4pfjets_phi.push_back(ak4pfjets_p4.Phi());
        myEvent->pointerForak4pfjets_eta.push_back(ak4pfjets_p4.Eta());
        myEvent->ak4pfjets_mass.push_back(ak4pfjets_p4.M());

    } 

    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_passMEDbtag"))
    {
        theTree->SetBranchAddress("ak4pfjets_passMEDbtag",           &(one_ak4pfjets_passMEDbtag));
        myEvent->ak4pfjets_passMEDbtag.push_back(one_ak4pfjets_passMEDbtag);
    }
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_mva"))
    {
        theTree->SetBranchAddress("ak4pfjets_mva",           &(one_ak4pfjets_mva));
        myEvent->ak4pfjets_mva.push_back(one_ak4pfjets_mva);
    }
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_CSV"))
    {
        theTree->SetBranchAddress("ak4pfjets_CSV",           &(one_ak4pfjets_CSV));
        myEvent->ak4pfjets_CSV.push_back(one_ak4pfjets_CSV);
    }
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_parton_flavor"))
    {
        theTree->SetBranchAddress("ak4pfjets_parton_flavor", &(one_ak4pfjets_partonFlavour));
        myEvent->ak4pfjets_partonFlavour.push_back(one_ak4pfjets_partonFlavour);
    }
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_hadron_flavor"))
    {
        theTree->SetBranchAddress("ak4pfjets_hadron_flavor", &(one_ak4pfjets_hadronFlavour));
        myEvent->ak4pfjets_hadronFlavour.push_back(one_ak4pfjets_hadronFlavour);
    }
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_loose_puid"))
    {
        theTree->SetBranchAddress("ak4pfjets_loose_puid",          &(one_ak4pfjets_loose_puid));
        myEvent->ak4pfjets_loose_puid.push_back(one_ak4pfjets_loose_puid);
    }
    if(theTree->GetListOfBranches()->FindObject("ak4pfjets_loose_pfid"))
    {
        theTree->SetBranchAddress("ak4pfjets_loose_pfid",    &(one_ak4pfjets_loose_pfid));
        myEvent->ak4pfjets_loose_pfid.push_back(one_ak4pfjets_loose_pfid);
    }

    #endif
   
}*/
