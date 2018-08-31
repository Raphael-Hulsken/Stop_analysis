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
    vector<float> jet_pt;
    vector<float> jet_eta;
    vector<float> jet_phi;
    vector<float> jet_mass;
    vector<float> jet_puid;
    vector<float> jet_CSV;
    vector<bool> jet_loose_pfid;
    float         dphi_Wlep;
    float         ak4_htssm;
    float         lep1_d0;
    float         dR_lep_leadb;
    int           numberOfSelectedLeptons;
    float         lep2_d0;
    vector<float> ak4pfjets_puid;
    float         jetsPt;
    float         pfmet_phi;
    vector<float> ak4pfjets_pt;
    bool          HLT_SingleE;
    float         ak4_htosm;
    float         ak4_HT;
    vector<float> ak4pfjets_phi;
    float         lep1_passMediumID;
    float         lep2_dz;
    float         lep2_mass;
    float         lep1_dz;
    int           runId;
    float         pv_ndof;
    float         secondLeptonPhi;
    float         secondLeptonEta;
    float         leadingLeptonIso;
    int           numberOfGeneratedLeptons;
    float         HT;
    float         lep2_passMediumID;
    float         lep2_phi;
    float         secondLeptonPt;
    float         lep1_pt;
    int           eventId;
    int           numberOfSelectedJets;
    int           leadingLeptonId;
    int           pv_isFake;
    float         lep2_pt;
    float         leadingLeptonPt;
    vector<float> ak4pfjets_eta;
    float         jetsEta;
    float         leadingLeptonEta;
    float         Mjjj;
    vector<float> ak4pfjets_mass;
    int           ngoodbtags;
    float         lep2_eta;
    float         crossSection;
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
    float         pfmet;
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
    double 	  mc_weight;
    float         lep1_mass;
    float         lep2_passVeto;
    float         lep1_eta;
    float         lep2_MiniIso;
    float         pv_rho;
    float         MT2W;
    float         ETmiss;
    float         secondLeptonIso;
    int           numberOfSelectedElectrons;
    bool          HLT_SingleMu;
    vector<float> ak4pfjets_CSV;
    int           ngoodleps;
    float         topness;
    float         dphi_ak4pfjets_met;
    float         chi2;
    float         ETmissPhi;
    int           totalNumberOfInitialEvent;
    float         lep_sf;
    float         btag_sf;
    int           nvetoleps;
    int           numberOfSelectedMuons;
    float         hadronic_top_chi2;
   
    bool PassTrackVeto;
    bool PassTauVeto;

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
    vector<float>* pointerForak4pfjets_phi;
    vector<float>* pointerForak4pfjets_eta;
    vector<float>* pointerForak4pfjets_mass;
    vector<bool>*  pointerForak4pfjets_loose_pfid;
    vector<float>* pointerForak4pfjets_CSV;

   


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
    #ifdef USE_JETS
    myEvent->pointerForak4pfjets_pt = 0;
    myEvent->pointerForak4pfjets_phi = 0;
    myEvent->pointerForak4pfjets_eta = 0;
    myEvent->pointerForak4pfjets_mass = 0;
    myEvent->pointerForak4pfjets_puid = 0;
    myEvent->pointerForak4pfjets_loose_pfid = 0;
    myEvent->pointerForak4pfjets_CSV = 0;
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
    
    //not in the common format but needed
    //theTree->SetBranchAddress("scale1fb",                &(myEvent->scale1fb));
    //theTree->SetBranchAddress("crossSection",            &(myEvent->crossSection));
    ////////////////
    //
    //used for weight
    theTree->SetBranchAddress("totalNumberOfInitialEvent", &(myEvent->totalNumberOfInitialEvent));

    // multiplicities
    theTree->SetBranchAddress("nvetoleps",               &(myEvent->nvetoleps));
    theTree->SetBranchAddress("ngoodjets",               &(myEvent->ngoodjets));
    theTree->SetBranchAddress("ngoodbtags",              &(myEvent->ngoodbtags));
    theTree->SetBranchAddress("ngoodleps",               &(myEvent->ngoodleps));
    
    
    //gen info
    theTree->SetBranchAddress("genlepsfromtop",          &(myEvent->genlepsfromtop));
    
    //trigger
    theTree->SetBranchAddress("HLT_SingleMu",            &(myEvent->HLT_SingleMu));
    theTree->SetBranchAddress("HLT_SingleE",             &(myEvent->HLT_SingleE));
    
    //vetos
    theTree->SetBranchAddress("PassTrackVeto",            &(myEvent->PassTrackVeto));
    theTree->SetBranchAddress("PassTauVeto",             &(myEvent->PassTauVeto));
    
    //needed to define the channel !
    theTree->SetBranchAddress("lep1_pdgid",              &(myEvent->lep1_pdgid));
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
    #endif
    #ifdef USE_JETS_EXT
    theTree->SetBranchAddress("ak4pfjets_puid",          &(myEvent->pointerForak4pfjets_puid));
    theTree->SetBranchAddress("ak4pfjets_loose_pfid",    &(myEvent->pointerForak4pfjets_loose_pfid));
    theTree->SetBranchAddress("ak4pfjets_rho",           &(myEvent->ak4pfjets_rho));
    theTree->SetBranchAddress("dphi_ak4pfjets_met",      &(myEvent->dphi_ak4pfjets_met));
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

    cout<<"there anyway"<<endl;
    #ifdef USE_VAR_BASELINE
    theTree->SetBranchAddress("pfmet_phi",               &(myEvent->pfmet_phi));
    theTree->SetBranchAddress("pfmet",                   &(myEvent->pfmet));
    theTree->SetBranchAddress("MT2W",                    &(myEvent->MT2W));
    cout<<"load MT vari"<<endl;
    theTree->SetBranchAddress("mt_met_lep",              &(myEvent->mt_met_lep));
    theTree->SetBranchAddress("hadronic_top_chi2",       &(myEvent->hadronic_top_chi2));
    theTree->SetBranchAddress("dphi_Wlep",               &(myEvent->dphi_Wlep));
    #endif

    #ifdef USE_GLOBAL_VAR
    theTree->SetBranchAddress("ak4_htssm",               &(myEvent->ak4_htssm));
    theTree->SetBranchAddress("ak4_htosm",               &(myEvent->ak4_htosm));
    theTree->SetBranchAddress("ak4_HT",                  &(myEvent->ak4_HT));
    theTree->SetBranchAddress("dR_lep_leadb",            &(myEvent->dR_lep_leadb));
    theTree->SetBranchAddress("Mjjj",                    &(myEvent->Mjjj));
    theTree->SetBranchAddress("Mlb",                     &(myEvent->Mlb));
    theTree->SetBranchAddress("Mlb_leadb",               &(myEvent->Mlb_leadb));
    theTree->SetBranchAddress("topness",                 &(myEvent->topness));
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
    
    theTree->SetBranchAddress("numberOfGeneratedLeptons", &(myEvent->numberOfGeneratedLeptons));
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
    theTree->SetBranchAddress("gen_pt", &(myEvent->pointerForgen_pt));
    theTree->SetBranchAddress("gen_eta",  &(myEvent->pointerForgen_eta));
    theTree->SetBranchAddress("gen_phi",  &(myEvent->pointerForgen_phi));
    theTree->SetBranchAddress("gen_m", &(myEvent->pointerForgen_m));
    theTree->SetBranchAddress("gen_status", &(myEvent->pointerForgen_status));
    theTree->SetBranchAddress("gen_id",  &(myEvent->pointerForgen_id));
    theTree->SetBranchAddress("gen_daughter_n",  &(myEvent->pointerForgen_daughter_n));
    theTree->SetBranchAddress("gen_daughter_index", &(myEvent->pointerForgen_daughter_index));
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
    myEvent->pts = *(myEvent->pointerForak4pfjets_pt);
    myEvent->ak4pfjets_pt.clear();
    myEvent->ak4pfjets_pt = pts;
    myEvent->ak4pfjets_phi              = *(myEvent->pointerForak4pfjets_pt);
    myEvent->ak4pfjets_phi             = *(myEvent->pointerForak4pfjets_phi);
    myEvent->ak4pfjets_eta             = *(myEvent->pointerForak4pfjets_eta);
    //myEvent->ak4pfjets_mass            = *(myEvent->pointerForak4pfjets_mass);
    //myEvent->ak4pfjets_loose_pfid      = *(myEvent->pointerForak4pfjets_loose_pfid);
    myEvent->ak4pfjets_CSV             = *(myEvent->pointerForak4pfjets_CSV);
    */
    //
    #ifdef USE_JETS
    myEvent->jet_pt              = *(myEvent->pointerForak4pfjets_pt);
    myEvent->jet_phi             = *(myEvent->pointerForak4pfjets_phi);
    myEvent->jet_eta             = *(myEvent->pointerForak4pfjets_eta);
    myEvent->jet_mass            = *(myEvent->pointerForak4pfjets_mass);
    myEvent->jet_CSV             = *(myEvent->pointerForak4pfjets_CSV);
    #endif
    #ifdef USE_JETS_EXT
    myEvent->jet_puid      	= *(myEvent->pointerForak4pfjets_puid);
    myEvent->jet_loose_pfid      = *(myEvent->pointerForak4pfjets_loose_pfid);
    #endif
    //*/
    #ifdef USE_GEN_INFO
    myEvent->gen_pt			 =			*(myEvent->pointerForgen_pt);            	 
    myEvent->gen_eta			 =                      *(myEvent->pointerForgen_eta);
    myEvent->gen_phi			 =                      *(myEvent->pointerForgen_phi);
    myEvent->gen_m			 =                      *(myEvent->pointerForgen_m);
    myEvent->gen_status			 =                      *(myEvent->pointerForgen_status);
    myEvent->gen_id			 =                      *(myEvent->pointerForgen_id);
    myEvent->gen_daughter_n			 =              *(myEvent->pointerForgen_daughter_n);
    myEvent->gen_daughter_index			 =      *(myEvent->pointerForgen_daughter_index);
    #endif
    #ifdef USE_GEN_INFO_EXT
    myEvent->gen_charge			 =                      *(myEvent->pointerForgen_charge);
    myEvent->gen_index			 =                      *(myEvent->pointerForgen_index);
    myEvent->gen_mother_index			 =              *(myEvent->pointerForgen_mother_index);
    #endif

    //Fill temporary info
    myEvent->leadingLepton.SetPtEtaPhiM(myEvent->lep1_pt, myEvent->lep1_eta, myEvent->lep1_phi, myEvent->lep1_mass);
    myEvent->secondLepton.SetPtEtaPhiM(myEvent->lep2_pt, myEvent->lep2_eta, myEvent->lep2_phi, myEvent->lep2_mass);

}
 
#endif
