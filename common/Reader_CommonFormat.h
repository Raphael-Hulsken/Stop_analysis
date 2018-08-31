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

//#define USE_GEN_INFO

// ##########################
// #  Baby event structure  #
// ##########################
 
typedef struct
{
    bool hasGenInfo;
    vector<float> jet_pt;
    vector<float> jet_eta;
    vector<float> jet_phi;
    vector<float> jet_mass;
    vector<float> jet_puid;
    vector<float> jet_CSV;
    vector<bool>  jet_loose_pfid;
    //vector<string> trigger_name;
    //vector<bool>   trigger_pass;
    float         lep1_d0;
    int           numberOfSelectedLeptons;
    float         lep2_d0;
    vector<float> ak4pfjets_puid;
    float         jetsPt;
    float         pfmet_phi;
    float         met_sig;
    vector<float> ak4pfjets_pt;
    vector<float> ak8pfjets_pt;
    vector<float> ak10pfjets_pt;
    bool          HLT_Ele25_eta2p1_WPLoose;
    bool          HLT_Ele27_eta2p1_WPLoose;
    bool          HLT_IsoMu20;
    bool          HLT_IsoMu22;
    bool          HLT_PFMET170;
    bool          HLT_PFMET100_PFMHT100_IDTight;
    vector<float> ak4pfjets_phi;
    vector<float> ak8pfjets_phi;
    vector<float> ak10pfjets_phi;
    float         lep1_passMediumID;
    float         lep2_dz;
    float         lep2_mass;
    float         lep1_dz;
    int           runId;
    float         pv_ndof;
    float         secondLeptonPhi;
    float         secondLeptonEta;
    float         leadingLeptonIso;
    int           numberOfGeneratedLeptons = 13;
    float         HT;
    float         lep2_passMediumID;
    float         lep2_phi;
    float         secondLeptonPt;
    float         lep1_pt;
    int           eventId;
    int           numberOfSelectedJets;
    int           leadingLeptonId = 0;
    int           pv_isFake;
    float         lep2_pt;
    float         leadingLeptonPt;
    vector<float> ak4pfjets_eta;
    vector<float> ak8pfjets_eta;
    vector<float> ak10pfjets_eta;
    float         jetsEta;
    float         leadingLeptonEta;
    float         Mjjj;
    vector<float> ak4pfjets_mass;
    vector<float> ak8pfjets_mass;
    vector<float> ak10pfjets_mass;
    vector<float> ak8pfjets_pruned_mass;
    vector<float> ak10pfjets_pruned_mass;
    vector<float> ak8pfjets_corrpruned_mass;
    vector<float> ak8pfjets_trimmed_mass;
    vector<float> ak10pfjets_trimmed_mass;
    vector<float> ak8pfjets_softdrop_mass;
    vector<float> ak10pfjets_softdrop_mass;
    vector<int>   ak8pfjets_nSubJets;
    vector<int>   ak10pfjets_nSubJets;
    vector<float> ak8pfjets_tau1;
    vector<float> ak10pfjets_tau1;
    vector<float> ak8pfjets_tau2;
    vector<float> ak10pfjets_tau2;
    vector<int>   ak4pfjets_partonFlavour;
    vector<int>   ak4pfjets_hadronFlavour;
    int           ngoodbtags;
    float         lep2_eta;
    float         crossSection = -13;
    float         mt_met_lep;
    float         M3b;
    int           lep1_pdgid;
    float         jetsCSVv2;
    float         scale1fb;
    float         MT;
    float         jetsPhi;
    vector<bool>  ak4pfjets_loose_pfid;
    int           ngoodjets;
    float         lep1_passVeto;
    float         pfmet = -13;
    float         lep1_MiniIso;
    int           lumiId;
    float         lep1_phi;
    int           lep2_pdgid;
    float         pv_z;
    float         ak4pfjets_rho;
    float         Mlb;
    float         selectionCode;
    float         jetsPUid;
    float         Mlb_leadb;
    float         leadingLeptonPhi;
    int           numberOfBTaggedJets;
    int           genlepsfromtop;
    int           secondLeptonId;
    float         jetsCSV;
    float         pu_weight;
    Double_t 	  mc_weight;
    float         lep1_mass;
    float         lep2_passVeto;
    float         lep1_eta;
    float         lep2_MiniIso;
    float         pv_rho;
    float         MT2W = -13;
    float         ETmiss;
    float         secondLeptonIso;
    int           numberOfSelectedElectrons;
    vector<float> ak4pfjets_CSV;
    vector<float> ak4pfjets_qgtag;
    int           ngoodleps;
    float         topness;
    float         dphi_ak4pfjets_met;
    float         chi2;
    float         ETmissPhi;
    //Int_t      totalNumberOfInitialEvent = -13;
    int64_t      totalNumberOfInitialEvent = -13;
    //double      totalNumberOfInitialEvent = -13;
    float         lep_sf;
    float         btag_sf;
    int           nvetoleps;
    int           numberOfSelectedMuons = -13;
    float         hadronic_top_chi2;
    int           nvertex = -13;
    int           puIntime = -13;
    int           puTrue = -13;
   
    vector<float> gen_neutralino_m;
    vector<float> gen_stop_m;
    vector<float>* pointerForgen_neutralino_m;
    vector<float>* pointerForgen_stop_m;
    vector<float> genTops_pt;
    vector<int> genTops_pdgid;
    vector<float>* pointerForgenTops_pt;
    vector<int>* pointerForgenTops_pdgid;


    bool PassTrackVeto;
    bool PassTauVeto;

    float metGen_pt;
    float metGen_phi;

    float DeltaRlj;
    vector<float> ak8pfjets_nsubjettines; //WARNING! this is empty now
    vector<float> ak10pfjets_nsubjettines; //WARNING! this is empty now

    #ifdef USE_NEW_VAR
    float ST;
    float LP;
    float Meff;
    float MTdeco_Q;
    float DeltaPtbb;
    float DeltaPhibb;
    float DeltaRbb;
    float dphi_Wlep;
    float dR_lep_leadb;
    float ak4_HT;
    float ak4_htssm;
    float ak4_htosm;
    #endif

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
    float minDPhi_jmet; 
    float genMET;
    float Pz_ttbar;
    float Pt_ttbar;
    float E_ttbar;
    float M_ttbar;
    int ttbar_decay;
    int nofISR;
    int index_b;
    int index_bbar;
    bool matched;
    vector<int> isr_id;
    vector<TLorentzVector> isr_p4;
    vector<TLorentzVector> neutrinos;
    vector<TLorentzVector> neutralinos;
    
    vector<int>* pointerForisr_id;
    vector<TLorentzVector>* pointerForisr_p4;
    vector<TLorentzVector>* pointerForneutrinos;
    vector<TLorentzVector>* pointerForneutralinos;
   #endif

   // Intermediate pointers for special types
    // Yes, this shit is needed because ROOT is crap.
 
    vector<float>* pointerForak4pfjets_puid;
    vector<float>* pointerForak4pfjets_pt;
    vector<float>* pointerForak8pfjets_pt;
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
    vector<bool>* pointerForTriggerPass;

   //Add content
   TLorentzVector leadingLepton;
   TLorentzVector secondLepton;

}
babyEvent;
 
// #############################
// #  Branches initialization  #
// #############################
 
void InitializeBranchesForReading(TTree* theTree, babyEvent* myEvent)
{
    myEvent->hasGenInfo = false;

    myEvent->pointerForTriggerName = NULL;
    myEvent->pointerForTriggerPass = NULL;
    myEvent->pointerForgen_neutralino_m = 0;
    myEvent->pointerForgen_stop_m = 0;
    myEvent->pointerForgenTops_pt = 0;
    myEvent->pointerForgenTops_pdgid = 0;

    #ifdef USE_JETS
    myEvent->pointerForak4pfjets_pt = 0;
    myEvent->pointerForak4pfjets_phi = 0;
    myEvent->pointerForak4pfjets_eta = 0;
    myEvent->pointerForak4pfjets_mass = 0;
    myEvent->pointerForak4pfjets_puid = 0;
    myEvent->pointerForak4pfjets_loose_pfid = 0;
    myEvent->pointerForak4pfjets_CSV = 0;
    myEvent->pointerForak4pfjets_partonFlavour = 0;
    myEvent->pointerForak4pfjets_hadronFlavour = 0;
    myEvent->pointerForak4pfjets_qgtag = 0;
    #endif

    #ifdef USE_AK8_JETS
    myEvent->pointerForak8pfjets_pt = 0;
    myEvent->pointerForak8pfjets_phi = 0;
    myEvent->pointerForak8pfjets_eta = 0;
    myEvent->pointerForak8pfjets_mass = 0;
    myEvent->pointerForak8pfjets_pruned_mass = 0;
    myEvent->pointerForak8pfjets_corrpruned_mass = 0;
    myEvent->pointerForak8pfjets_trimmed_mass = 0;
    myEvent->pointerForak8pfjets_softdrop_mass = 0;
    myEvent->pointerForak8pfjets_nSubJets = 0;
    myEvent->pointerForak8pfjets_tau1 = 0;
    myEvent->pointerForak8pfjets_tau2 = 0;
    #endif

    #ifdef USE_AK10_JETS
    myEvent->pointerForak10pfjets_pt = 0;
    myEvent->pointerForak10pfjets_phi = 0;
    myEvent->pointerForak10pfjets_eta = 0;
    myEvent->pointerForak10pfjets_mass = 0;
    myEvent->pointerForak10pfjets_pruned_mass = 0;
    myEvent->pointerForak10pfjets_trimmed_mass = 0;
    myEvent->pointerForak10pfjets_softdrop_mass = 0;
    myEvent->pointerForak10pfjets_nSubJets = 0;
    myEvent->pointerForak10pfjets_tau1 = 0;
    myEvent->pointerForak10pfjets_tau2 = 0;
    #endif
    
    #ifdef USE_GEN_INFO
    myEvent->pointerForgen_pt = 0;
    myEvent->pointerForgen_eta = 0;
    myEvent->pointerForgen_phi = 0;
    myEvent->pointerForgen_m = 0;
    myEvent->pointerForgen_status = 0;
    myEvent->pointerForgen_id = 0;
    myEvent->pointerForgen_daughter_n = 0;
    myEvent->pointerForgen_daughter_index = 0;
    #endif 
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

    //not in the common format but needed
    //theTree->SetBranchAddress("scale1fb",                &(myEvent->scale1fb));
    //theTree->SetBranchAddress("crossSection",            &(myEvent->crossSection));
    ////////////////
    //
    //used for weight
    //std::cout << "reading basic info" << std::endl;
    theTree->SetBranchAddress("totalNumberOfInitialEvent", &(myEvent->totalNumberOfInitialEvent));

    // multiplicities
    theTree->SetBranchAddress("nvetoleps",               &(myEvent->nvetoleps));
    theTree->SetBranchAddress("ngoodjets",               &(myEvent->ngoodjets));
    theTree->SetBranchAddress("ngoodbtags",              &(myEvent->ngoodbtags));
    theTree->SetBranchAddress("ngoodleps",               &(myEvent->ngoodleps));
    
    
    //gen info
    theTree->SetBranchAddress("genlepsfromtop",          &(myEvent->genlepsfromtop));
    
    //trigger
    theTree->SetBranchAddress("HLT_Ele25_eta2p1_WPLoose",            &(myEvent->HLT_Ele25_eta2p1_WPLoose));
    theTree->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose",            &(myEvent->HLT_Ele27_eta2p1_WPLoose));
    theTree->SetBranchAddress("HLT_IsoMu20",            &(myEvent->HLT_IsoMu20));
    theTree->SetBranchAddress("HLT_IsoMu22",            &(myEvent->HLT_IsoMu22));
    theTree->SetBranchAddress("HLT_PFMET170",            &(myEvent->HLT_PFMET170));
    theTree->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight",            &(myEvent->HLT_PFMET100_PFMHT100_IDTight));
    //theTree->SetBranchAddress("trigger_name",            &(myEvent->pointerForTriggerName));
    //theTree->SetBranchAddress("trigger_pass",            &(myEvent->pointerForTriggerPass));
    //
    //signal related
    theTree->SetBranchAddress("gen_neutralino_m", &(myEvent->pointerForgen_neutralino_m));
    theTree->SetBranchAddress("gen_stop_m", &(myEvent->pointerForgen_stop_m));
    //gen top
    theTree->SetBranchAddress("genTops_pt", &(myEvent->pointerForgenTops_pt));
    theTree->SetBranchAddress("genTops_pdgid", &(myEvent->pointerForgenTops_pdgid));
    
    //vetos
    theTree->SetBranchAddress("PassTrackVeto",            &(myEvent->PassTrackVeto));
    theTree->SetBranchAddress("PassTauVeto",             &(myEvent->PassTauVeto));
    
    //gen met
    theTree->SetBranchAddress("metGen_pt",              &(myEvent->metGen_pt));
    theTree->SetBranchAddress("metGen_phi",             &(myEvent->metGen_phi));
    
    //needed to define the channel !
    theTree->SetBranchAddress("lep1_pdgid",              &(myEvent->lep1_pdgid));
    //std::cout << "reading basic info" << std::endl;
    #ifdef USE_LEP1
    theTree->SetBranchAddress("lep1_pt",                 &(myEvent->lep1_pt));
    theTree->SetBranchAddress("lep1_eta",                &(myEvent->lep1_eta));
    theTree->SetBranchAddress("lep1_phi",                &(myEvent->lep1_phi));
    theTree->SetBranchAddress("lep1_mass",               &(myEvent->lep1_mass));
    #endif
    #ifdef USE_LEP1_EXT
    theTree->SetBranchAddress("lep1_d0",                 &(myEvent->lep1_d0));
    theTree->SetBranchAddress("lep1_dz",                 &(myEvent->lep1_dz));
    theTree->SetBranchAddress("lep1_passVeto",           &(myEvent->lep1_passVeto));
    theTree->SetBranchAddress("lep1_passMediumID",       &(myEvent->lep1_passMediumID));
    theTree->SetBranchAddress("lep1_MiniIso",            &(myEvent->lep1_MiniIso));
    #endif
    #ifdef USE_LEP2
    theTree->SetBranchAddress("lep2_pt",                 &(myEvent->lep2_pt));
    theTree->SetBranchAddress("lep2_eta",                &(myEvent->lep2_eta));
    theTree->SetBranchAddress("lep2_phi",                &(myEvent->lep2_phi));
    theTree->SetBranchAddress("lep2_mass",               &(myEvent->lep2_mass));
    theTree->SetBranchAddress("lep2_pdgid",              &(myEvent->lep2_pdgid));
    #endif
    #ifdef USE_LEP2_EXT
    theTree->SetBranchAddress("lep2_d0",                 &(myEvent->lep2_d0));
    theTree->SetBranchAddress("lep2_dz",                 &(myEvent->lep2_dz));
    theTree->SetBranchAddress("lep2_passVeto",           &(myEvent->lep2_passVeto));
    theTree->SetBranchAddress("lep2_passMediumID",       &(myEvent->lep2_passMediumID));
    theTree->SetBranchAddress("lep2_MiniIso",            &(myEvent->lep2_MiniIso));
    #endif

    #ifdef USE_JETS
    theTree->SetBranchAddress("ak4pfjets_pt",            &(myEvent->pointerForak4pfjets_pt));
    theTree->SetBranchAddress("ak4pfjets_phi",           &(myEvent->pointerForak4pfjets_phi));
    theTree->SetBranchAddress("ak4pfjets_eta",           &(myEvent->pointerForak4pfjets_eta));
    theTree->SetBranchAddress("ak4pfjets_mass",          &(myEvent->pointerForak4pfjets_mass));
    theTree->SetBranchAddress("ak4pfjets_CSV",           &(myEvent->pointerForak4pfjets_CSV));
    theTree->SetBranchAddress("ak4pfjets_partonFlavour", &(myEvent->pointerForak4pfjets_partonFlavour));
    theTree->SetBranchAddress("ak4pfjets_hadronFlavour", &(myEvent->pointerForak4pfjets_hadronFlavour));
    theTree->SetBranchAddress("ak4pfjets_qgtag",           &(myEvent->pointerForak4pfjets_qgtag));
    #endif

    #ifdef USE_JETS_EXT
    theTree->SetBranchAddress("ak4pfjets_puid",          &(myEvent->pointerForak4pfjets_puid));
    theTree->SetBranchAddress("ak4pfjets_loose_pfid",    &(myEvent->pointerForak4pfjets_loose_pfid));
    theTree->SetBranchAddress("ak4pfjets_rho",           &(myEvent->ak4pfjets_rho));
    theTree->SetBranchAddress("dphi_ak4pfjets_met",      &(myEvent->dphi_ak4pfjets_met));
    #endif

    #ifdef USE_AK8_JETS
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
    #endif 

    #ifdef USE_WEIGHTS
    theTree->SetBranchAddress("lep_sf",                  &(myEvent->lep_sf));
    theTree->SetBranchAddress("btag_sf",                 &(myEvent->btag_sf));
    theTree->SetBranchAddress("pu_weight",               &(myEvent->pu_weight));
    theTree->SetBranchAddress("mc_weight",               &(myEvent->mc_weight));
    #endif

    #ifdef USE_VAR_BASELINE
    theTree->SetBranchAddress("pfmet_phi",               &(myEvent->pfmet_phi));
    theTree->SetBranchAddress("met_sig",                 &(myEvent->met_sig));
    theTree->SetBranchAddress("pfmet",                   &(myEvent->pfmet));
    theTree->SetBranchAddress("MT2W",                    &(myEvent->MT2W));
    theTree->SetBranchAddress("mt_met_lep",              &(myEvent->mt_met_lep));
    theTree->SetBranchAddress("hadronic_top_chi2",       &(myEvent->hadronic_top_chi2));
    theTree->SetBranchAddress("nvertex",                 &(myEvent->nvertex));
    theTree->SetBranchAddress("puIntime",                &(myEvent->puIntime));
    theTree->SetBranchAddress("puTrue",                  &(myEvent->puTrue));
    #endif

    #ifdef USE_GLOBAL_VAR
    theTree->SetBranchAddress("Mjjj",                    &(myEvent->Mjjj));
    theTree->SetBranchAddress("Mlb",                     &(myEvent->Mlb));
    theTree->SetBranchAddress("Mlb_leadb",               &(myEvent->Mlb_leadb));
    theTree->SetBranchAddress("topness",                 &(myEvent->topness));
    theTree->SetBranchAddress("leadingLeptonId",         &(myEvent->leadingLeptonId));
    theTree->SetBranchAddress("numberOfGeneratedLeptons",&(myEvent->numberOfGeneratedLeptons));
    #endif
    
    #ifdef USE_NEW_VAR
    theTree->SetBranchAddress("ST"		,&(myEvent->ST));
    theTree->SetBranchAddress("LP"		,&(myEvent->LP));
    theTree->SetBranchAddress("Meff"		,&(myEvent->Meff));
    theTree->SetBranchAddress("MTdeco_Q"	,&(myEvent->MTdeco_Q));
    theTree->SetBranchAddress("DeltaPtbb"	,&(myEvent->DeltaPtbb));
    theTree->SetBranchAddress("DeltaPhibb"	,&(myEvent->DeltaPhibb));
    theTree->SetBranchAddress("DeltaRbb"	,&(myEvent->DeltaRbb));
    theTree->SetBranchAddress("dphi_Wlep",	&(myEvent->dphi_Wlep));
    theTree->SetBranchAddress("dR_lep_leadb",	&(myEvent->dR_lep_leadb));
    theTree->SetBranchAddress("ak4_HT",		&(myEvent->ak4_HT));
    theTree->SetBranchAddress("ak4_htssm",	&(myEvent->ak4_htssm));
    theTree->SetBranchAddress("ak4_htosm",	&(myEvent->ak4_htosm));
    #endif

    #ifdef USE_GEN_LOSTLEPTON
    theTree->SetBranchAddress("genLostLeptons_pt",	&(myEvent->pointerForGenLostLeptons_pt));
    theTree->SetBranchAddress("genLostLeptons_eta",	&(myEvent->pointerForGenLostLeptons_eta));
    theTree->SetBranchAddress("genLostLeptons_pdgid",   &(myEvent->pointerForGenLostLeptons_pdgid));
    //cout<<"Branch adress done !"<<endl;
    #endif

    //old framework
    #ifdef USE_OLD_VAR
    theTree->SetBranchAddress("leadingLeptonPt",         &(myEvent->leadingLeptonPt));
    theTree->SetBranchAddress("leadingLeptonEta",        &(myEvent->leadingLeptonEta));
    theTree->SetBranchAddress("leadingLeptonPhi",        &(myEvent->leadingLeptonPhi));
    theTree->SetBranchAddress("leadingLeptonId",         &(myEvent->leadingLeptonId));
    theTree->SetBranchAddress("leadingLeptonIso",        &(myEvent->leadingLeptonIso));
    
    theTree->SetBranchAddress("secondLeptonPt",          &(myEvent->secondLeptonPt));
    theTree->SetBranchAddress("secondLeptonPhi",         &(myEvent->secondLeptonPhi));
    theTree->SetBranchAddress("secondLeptonEta",         &(myEvent->secondLeptonEta));
    theTree->SetBranchAddress("secondLeptonId",          &(myEvent->secondLeptonId));
    theTree->SetBranchAddress("secondLeptonIso",         &(myEvent->secondLeptonIso));
    
    theTree->SetBranchAddress("runId",                   &(myEvent->runId));
    
    theTree->SetBranchAddress("numberOfSelectedLeptons", &(myEvent->numberOfSelectedLeptons));
    theTree->SetBranchAddress("numberOfSelectedElectrons", &(myEvent->numberOfSelectedElectrons));
    theTree->SetBranchAddress("numberOfSelectedMuons",   &(myEvent->numberOfSelectedMuons));
    theTree->SetBranchAddress("numberOfBTaggedJets",     &(myEvent->numberOfBTaggedJets));
    theTree->SetBranchAddress("numberOfSelectedJets",    &(myEvent->numberOfSelectedJets));
    
    theTree->SetBranchAddress("jetsPt",                  &(myEvent->jetsPt));
    theTree->SetBranchAddress("jetsEta",                 &(myEvent->jetsEta));
    theTree->SetBranchAddress("jetsPhi",                 &(myEvent->jetsPhi));
    theTree->SetBranchAddress("jetsCSV",                 &(myEvent->jetsCSV));
    theTree->SetBranchAddress("jetsCSVv2",               &(myEvent->jetsCSVv2));
    theTree->SetBranchAddress("jetsPUid",                &(myEvent->jetsPUid));
    
    theTree->SetBranchAddress("eventId",                 &(myEvent->eventId));
    theTree->SetBranchAddress("lumiId",                  &(myEvent->lumiId));
    theTree->SetBranchAddress("selectionCode",           &(myEvent->selectionCode));
    
    theTree->SetBranchAddress("HT",                      &(myEvent->HT));
    theTree->SetBranchAddress("MT",                      &(myEvent->MT));
    theTree->SetBranchAddress("ETmiss",                  &(myEvent->ETmiss));
    theTree->SetBranchAddress("ETmissPhi",               &(myEvent->ETmissPhi));
    theTree->SetBranchAddress("chi2",                    &(myEvent->chi2));
    theTree->SetBranchAddress("M3b",                     &(myEvent->M3b));
    #endif
   
   //*/
    #ifdef USE_GEN_INFO
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
    #endif
    
    
    #ifdef USE_SKIMMING_VAR
    myEvent->pointerForisr_id = 0;
    myEvent->pointerForisr_p4 = 0;
    myEvent->pointerForneutrinos = 0;
    myEvent->pointerForneutralinos = 0;
    theTree->SetBranchAddress("minDPhi_jmet", &myEvent->minDPhi_jmet); 
    theTree->SetBranchAddress("genMET", &myEvent->genMET); 
    theTree->SetBranchAddress("Pz_ttbar", &myEvent->Pz_ttbar); 
    theTree->SetBranchAddress("Pt_ttbar", &myEvent->Pt_ttbar); 
    theTree->SetBranchAddress("E_ttbar", &myEvent->E_ttbar); 
    theTree->SetBranchAddress("M_ttbar", &myEvent->M_ttbar); 
    theTree->SetBranchAddress("ttbar_decay", &myEvent->ttbar_decay); 
    theTree->SetBranchAddress("nofISR", &myEvent->nofISR); 
    theTree->SetBranchAddress("index_b", &myEvent->index_b); 
    theTree->SetBranchAddress("index_bbar", &myEvent->index_bbar); 
    theTree->SetBranchAddress("matched", &myEvent->matched); 
    theTree->SetBranchAddress("isr_id", &(myEvent->pointerForisr_id)); 
    theTree->SetBranchAddress("isr_p4", &myEvent->pointerForisr_p4); 
    theTree->SetBranchAddress("neutralinos", &myEvent->pointerForneutralinos); 
    theTree->SetBranchAddress("neutrinos", &myEvent->pointerForneutrinos); 
    #endif
}



// ################################
// #  Function to read one event  #
// ################################
 
void ReadEvent(TTree* theTree, long int i, babyEvent* myEvent)
{
    theTree->GetEntry(i);

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
    if(myEvent->pointerForgen_neutralino_m!=0 && myEvent->pointerForgen_stop_m!=0){
      myEvent->gen_neutralino_m	         =                      *(myEvent->pointerForgen_neutralino_m);
      myEvent->gen_stop_m			 =                      *(myEvent->pointerForgen_stop_m);
    }
    if(myEvent->pointerForgenTops_pt!=0 && myEvent->pointerForgenTops_pdgid!=0){
    	myEvent->genTops_pt = *(myEvent->pointerForgenTops_pt);
    	myEvent->genTops_pdgid = *(myEvent->pointerForgenTops_pdgid);
    }
    
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
    //*/
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
    //std::cout << "setting gen info" << std::endl;
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
    #endif
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

}
 
#endif
