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

// ##########################
// #  Baby event structure  #
// ##########################

typedef struct
{
    float                dphi_Wlep;
    float                ak4_htssm;
    float                lep1_d0;
    float                dR_lep_leadb;
    int                  numberOfSelectedLeptons;
    float                lep2_d0;
    int                  ls;
    //vector<float>        ak4pfjets_puid;
    float                jetsPt;
    float                pfmet_phi;
    //vector<float>        ak4pfjets_pt;
    bool                 HLT_SingleE;
    float                ak4_htosm;
    float                ak4_HT;
    //vector<float>        ak4pfjets_phi;
    float                lep1_passMediumID;
    float                lep2_dz;
    //vector<vector<int> > gen_daughter_index;
    float                lep2_mass;
    float                lep1_dz;
    int                  runId;
    float                pv_ndof;
    float                secondLeptonPhi;
    float                secondLeptonEta;
    //vector<float>        gen_phi;
    int                  numberOfGeneratedLeptons;
    float                HT;
    float                lep2_passMediumID;
    float                jetsCSV;
    bool                 PassTauVeto;
    float                secondLeptonPt;
    float                lep1_pt;
    int                  event;
    int                  eventId;
    int                  numberOfSelectedJets;
    int                  leadingLeptonId;
    int                  pv_isFake;
    //vector<int>          gen_status;
    float                lep2_pt;
    float                leadingLeptonPt;
    //vector<float>        ak4pfjets_eta;
    float                leadingLeptonEta;
    //vector<int>          gen_charge;
    int                  run;
    float                Mjjj;
    //vector<float>        ak4pfjets_mass;
    int                  ngoodbtags;
    float                lep2_eta;
    float                crossSection;
    //vector<float>        gen_pt;
    float                mt_met_lep;
    float                M3b;
    int                  lep1_pdgid;
    bool                 PassTrackVeto;
    float                jetsCSVv2;
    float                scale1fb;
    float                MT;
    float                jetsPhi;
    //vector<bool>         ak4pfjets_loose_pfid;
    float                lep2_phi;
    int                  ngoodjets;
    float                lep1_passVeto;
    float                pfmet;
    float                lep1_MiniIso;
    int                  lumiId;
    float                lep1_phi;
    int                  lep2_pdgid;
    float                pv_z;
    //vector<int>          gen_daughter_n;
    //vector<int>          gen_mother_index;
    float                ak4pfjets_rho;
    float                Mlb;
    float                selectionCode;
    float                jetsPUid;
    float                Mlb_leadb;
    float                leadingLeptonPhi;
    //vector<int>          gen_id;
    //vector<float>        gen_eta;
    int                  numberOfBTaggedJets;
    int                  genlepsfromtop;
    int                  secondLeptonId;
    float                jetsEta;
    float                pu_weight;
    float                lep1_mass;
    float                lep2_passVeto;
    float                lep1_eta;
    float                lep2_MiniIso;
    float                pv_rho;
    float                MT2W;
    float                ETmiss;
    float                secondLeptonIso;
    int                  numberOfSelectedElectrons;
    //vector<int>          gen_index;
    bool                 HLT_SingleMu;
    //vector<float>        ak4pfjets_CSV;
    int                  ngoodleps;
    float                topness;
    float                dphi_ak4pfjets_met;
    float                chi2;
    float                ETmissPhi;
    float                hadronic_top_chi2;
    int                  totalNumberOfInitialEvent;
    float                lep_sf;
    int                  gen_n;
    //vector<float>        gen_m;
    float                btag_sf;
    int                  nvetoleps;
    int                  numberOfSelectedMuons;
    float                leadingLeptonIso;

    // Intermediate pointers for special types
    // Yes, this shit is needed because ROOT is crap.

    /*
    vector<float>*        pointerForak4pfjets_puid;
    vector<float>*        pointerForak4pfjets_pt;
    vector<float>*        pointerForak4pfjets_phi;
    vector<vector<int> >* pointerForgen_daughter_index;
    vector<float>*        pointerForgen_phi;
    vector<int>*          pointerForgen_status;
    vector<float>*        pointerForak4pfjets_eta;
    vector<int>*          pointerForgen_charge;
    vector<float>*        pointerForak4pfjets_mass;
    vector<float>*        pointerForgen_pt;
    vector<bool>*         pointerForak4pfjets_loose_pfid;
    vector<int>*          pointerForgen_daughter_n;
    vector<int>*          pointerForgen_mother_index;
    vector<int>*          pointerForgen_id;
    vector<float>*        pointerForgen_eta;
    vector<int>*          pointerForgen_index;
    vector<float>*        pointerForak4pfjets_CSV;
    vector<float>*        pointerForgen_m;
    */
}
babyEvent;

// #############################
// #  Branches initialization  #
// #############################

void InitializeBranchesForReading(TTree* theTree, babyEvent* myEvent)
{
    /*
    myEvent->pointerForak4pfjets_puid = 0;
    myEvent->pointerForak4pfjets_pt = 0;
    myEvent->pointerForak4pfjets_phi = 0;
    myEvent->pointerForgen_daughter_index = 0;
    myEvent->pointerForgen_phi = 0;
    myEvent->pointerForgen_status = 0;
    myEvent->pointerForak4pfjets_eta = 0;
    myEvent->pointerForgen_charge = 0;
    myEvent->pointerForak4pfjets_mass = 0;
    myEvent->pointerForgen_pt = 0;
    myEvent->pointerForak4pfjets_loose_pfid = 0;
    myEvent->pointerForgen_daughter_n = 0;
    myEvent->pointerForgen_mother_index = 0;
    myEvent->pointerForgen_id = 0;
    myEvent->pointerForgen_eta = 0;
    myEvent->pointerForgen_index = 0;
    myEvent->pointerForak4pfjets_CSV = 0;
    myEvent->pointerForgen_m = 0;
    */

    theTree->SetBranchAddress("dphi_Wlep",               &(myEvent->dphi_Wlep));
    theTree->SetBranchAddress("ak4_htssm",               &(myEvent->ak4_htssm));
    theTree->SetBranchAddress("lep1_d0",                 &(myEvent->lep1_d0));
    theTree->SetBranchAddress("dR_lep_leadb",            &(myEvent->dR_lep_leadb));
    theTree->SetBranchAddress("numberOfSelectedLeptons", &(myEvent->numberOfSelectedLeptons));
    theTree->SetBranchAddress("lep2_d0",                 &(myEvent->lep2_d0));
    theTree->SetBranchAddress("ls",                      &(myEvent->ls));
    //theTree->SetBranchAddress("ak4pfjets_puid",          &(myEvent->pointerForak4pfjets_puid));
    theTree->SetBranchAddress("jetsPt",                  &(myEvent->jetsPt));
    theTree->SetBranchAddress("pfmet_phi",               &(myEvent->pfmet_phi));
    //theTree->SetBranchAddress("ak4pfjets_pt",            &(myEvent->pointerForak4pfjets_pt));
    theTree->SetBranchAddress("HLT_SingleE",             &(myEvent->HLT_SingleE));
    theTree->SetBranchAddress("ak4_htosm",               &(myEvent->ak4_htosm));
    theTree->SetBranchAddress("ak4_HT",                  &(myEvent->ak4_HT));
    //theTree->SetBranchAddress("ak4pfjets_phi",           &(myEvent->pointerForak4pfjets_phi));
    theTree->SetBranchAddress("lep1_passMediumID",       &(myEvent->lep1_passMediumID));
    theTree->SetBranchAddress("lep2_dz",                 &(myEvent->lep2_dz));
    //theTree->SetBranchAddress("gen_daughter_index",      &(myEvent->pointerForgen_daughter_index));
    theTree->SetBranchAddress("lep2_mass",               &(myEvent->lep2_mass));
    theTree->SetBranchAddress("lep1_dz",                 &(myEvent->lep1_dz));
    theTree->SetBranchAddress("runId",                   &(myEvent->runId));
    theTree->SetBranchAddress("pv_ndof",                 &(myEvent->pv_ndof));
    theTree->SetBranchAddress("secondLeptonPhi",         &(myEvent->secondLeptonPhi));
    theTree->SetBranchAddress("secondLeptonEta",         &(myEvent->secondLeptonEta));
    //theTree->SetBranchAddress("gen_phi",                 &(myEvent->pointerForgen_phi));
    theTree->SetBranchAddress("numberOfGeneratedLeptons", &(myEvent->numberOfGeneratedLeptons));
    theTree->SetBranchAddress("HT",                      &(myEvent->HT));
    theTree->SetBranchAddress("lep2_passMediumID",       &(myEvent->lep2_passMediumID));
    theTree->SetBranchAddress("jetsCSV",                 &(myEvent->jetsCSV));
    theTree->SetBranchAddress("PassTauVeto",             &(myEvent->PassTauVeto));
    theTree->SetBranchAddress("secondLeptonPt",          &(myEvent->secondLeptonPt));
    theTree->SetBranchAddress("lep1_pt",                 &(myEvent->lep1_pt));
    theTree->SetBranchAddress("event",                   &(myEvent->event));
    theTree->SetBranchAddress("eventId",                 &(myEvent->eventId));
    theTree->SetBranchAddress("numberOfSelectedJets",    &(myEvent->numberOfSelectedJets));
    theTree->SetBranchAddress("leadingLeptonId",         &(myEvent->leadingLeptonId));
    theTree->SetBranchAddress("pv_isFake",               &(myEvent->pv_isFake));
    //theTree->SetBranchAddress("gen_status",              &(myEvent->pointerForgen_status));
    theTree->SetBranchAddress("lep2_pt",                 &(myEvent->lep2_pt));
    theTree->SetBranchAddress("leadingLeptonPt",         &(myEvent->leadingLeptonPt));
    //theTree->SetBranchAddress("ak4pfjets_eta",           &(myEvent->pointerForak4pfjets_eta));
    theTree->SetBranchAddress("leadingLeptonEta",        &(myEvent->leadingLeptonEta));
    //theTree->SetBranchAddress("gen_charge",              &(myEvent->pointerForgen_charge));
    theTree->SetBranchAddress("run",                     &(myEvent->run));
    theTree->SetBranchAddress("Mjjj",                    &(myEvent->Mjjj));
    //theTree->SetBranchAddress("ak4pfjets_mass",          &(myEvent->pointerForak4pfjets_mass));
    theTree->SetBranchAddress("ngoodbtags",              &(myEvent->ngoodbtags));
    theTree->SetBranchAddress("lep2_eta",                &(myEvent->lep2_eta));
    theTree->SetBranchAddress("crossSection",            &(myEvent->crossSection));
    //theTree->SetBranchAddress("gen_pt",                  &(myEvent->pointerForgen_pt));
    theTree->SetBranchAddress("mt_met_lep",              &(myEvent->mt_met_lep));
    theTree->SetBranchAddress("M3b",                     &(myEvent->M3b));
    theTree->SetBranchAddress("lep1_pdgid",              &(myEvent->lep1_pdgid));
    theTree->SetBranchAddress("PassTrackVeto",           &(myEvent->PassTrackVeto));
    theTree->SetBranchAddress("jetsCSVv2",               &(myEvent->jetsCSVv2));
    theTree->SetBranchAddress("scale1fb",                &(myEvent->scale1fb));
    theTree->SetBranchAddress("MT",                      &(myEvent->MT));
    theTree->SetBranchAddress("jetsPhi",                 &(myEvent->jetsPhi));
    //theTree->SetBranchAddress("ak4pfjets_loose_pfid",    &(myEvent->pointerForak4pfjets_loose_pfid));
    theTree->SetBranchAddress("lep2_phi",                &(myEvent->lep2_phi));
    theTree->SetBranchAddress("ngoodjets",               &(myEvent->ngoodjets));
    theTree->SetBranchAddress("lep1_passVeto",           &(myEvent->lep1_passVeto));
    theTree->SetBranchAddress("pfmet",                   &(myEvent->pfmet));
    theTree->SetBranchAddress("lep1_MiniIso",            &(myEvent->lep1_MiniIso));
    theTree->SetBranchAddress("lumiId",                  &(myEvent->lumiId));
    theTree->SetBranchAddress("lep1_phi",                &(myEvent->lep1_phi));
    theTree->SetBranchAddress("lep2_pdgid",              &(myEvent->lep2_pdgid));
    theTree->SetBranchAddress("pv_z",                    &(myEvent->pv_z));
    //theTree->SetBranchAddress("gen_daughter_n",          &(myEvent->pointerForgen_daughter_n));
    //theTree->SetBranchAddress("gen_mother_index",        &(myEvent->pointerForgen_mother_index));
    theTree->SetBranchAddress("ak4pfjets_rho",           &(myEvent->ak4pfjets_rho));
    theTree->SetBranchAddress("Mlb",                     &(myEvent->Mlb));
    theTree->SetBranchAddress("selectionCode",           &(myEvent->selectionCode));
    theTree->SetBranchAddress("jetsPUid",                &(myEvent->jetsPUid));
    theTree->SetBranchAddress("Mlb_leadb",               &(myEvent->Mlb_leadb));
    theTree->SetBranchAddress("leadingLeptonPhi",        &(myEvent->leadingLeptonPhi));
    //theTree->SetBranchAddress("gen_id",                  &(myEvent->pointerForgen_id));
    //theTree->SetBranchAddress("gen_eta",                 &(myEvent->pointerForgen_eta));
    theTree->SetBranchAddress("numberOfBTaggedJets",     &(myEvent->numberOfBTaggedJets));
    theTree->SetBranchAddress("genlepsfromtop",          &(myEvent->genlepsfromtop));
    theTree->SetBranchAddress("secondLeptonId",          &(myEvent->secondLeptonId));
    theTree->SetBranchAddress("jetsEta",                 &(myEvent->jetsEta));
    theTree->SetBranchAddress("pu_weight",               &(myEvent->pu_weight));
    theTree->SetBranchAddress("lep1_mass",               &(myEvent->lep1_mass));
    theTree->SetBranchAddress("lep2_passVeto",           &(myEvent->lep2_passVeto));
    theTree->SetBranchAddress("lep1_eta",                &(myEvent->lep1_eta));
    theTree->SetBranchAddress("lep2_MiniIso",            &(myEvent->lep2_MiniIso));
    theTree->SetBranchAddress("pv_rho",                  &(myEvent->pv_rho));
    theTree->SetBranchAddress("MT2W",                    &(myEvent->MT2W));
    theTree->SetBranchAddress("ETmiss",                  &(myEvent->ETmiss));
    theTree->SetBranchAddress("secondLeptonIso",         &(myEvent->secondLeptonIso));
    theTree->SetBranchAddress("numberOfSelectedElectrons", &(myEvent->numberOfSelectedElectrons));
    //theTree->SetBranchAddress("gen_index",               &(myEvent->pointerForgen_index));
    theTree->SetBranchAddress("HLT_SingleMu",            &(myEvent->HLT_SingleMu));
    //theTree->SetBranchAddress("ak4pfjets_CSV",           &(myEvent->pointerForak4pfjets_CSV));
    theTree->SetBranchAddress("ngoodleps",               &(myEvent->ngoodleps));
    theTree->SetBranchAddress("topness",                 &(myEvent->topness));
    theTree->SetBranchAddress("dphi_ak4pfjets_met",      &(myEvent->dphi_ak4pfjets_met));
    theTree->SetBranchAddress("chi2",                    &(myEvent->chi2));
    theTree->SetBranchAddress("ETmissPhi",               &(myEvent->ETmissPhi));
    theTree->SetBranchAddress("hadronic_top_chi2",       &(myEvent->hadronic_top_chi2));
    theTree->SetBranchAddress("totalNumberOfInitialEvent", &(myEvent->totalNumberOfInitialEvent));
    theTree->SetBranchAddress("lep_sf",                  &(myEvent->lep_sf));
    theTree->SetBranchAddress("gen_n",                   &(myEvent->gen_n));
    //theTree->SetBranchAddress("gen_m",                   &(myEvent->pointerForgen_m));
    theTree->SetBranchAddress("btag_sf",                 &(myEvent->btag_sf));
    theTree->SetBranchAddress("nvetoleps",               &(myEvent->nvetoleps));
    theTree->SetBranchAddress("numberOfSelectedMuons",   &(myEvent->numberOfSelectedMuons));
    theTree->SetBranchAddress("leadingLeptonIso",        &(myEvent->leadingLeptonIso));
}

// ################################
// #  Function to read one event  #
// ################################

void ReadEvent(TTree* theTree, long int i, babyEvent* myEvent)
{

    theTree->GetEntry(i);
/*
    // Put actual content of special type branches where they should be...
    myEvent->ak4pfjets_puid            = *(myEvent->pointerForak4pfjets_puid);
    myEvent->ak4pfjets_pt              = *(myEvent->pointerForak4pfjets_pt);
    myEvent->ak4pfjets_phi             = *(myEvent->pointerForak4pfjets_phi);
    myEvent->gen_daughter_index        = *(myEvent->pointerForgen_daughter_index);
    myEvent->gen_phi                   = *(myEvent->pointerForgen_phi);
    myEvent->gen_status                = *(myEvent->pointerForgen_status);
    myEvent->ak4pfjets_eta             = *(myEvent->pointerForak4pfjets_eta);
    myEvent->gen_charge                = *(myEvent->pointerForgen_charge);
    myEvent->ak4pfjets_mass            = *(myEvent->pointerForak4pfjets_mass);
    myEvent->gen_pt                    = *(myEvent->pointerForgen_pt);
    myEvent->ak4pfjets_loose_pfid      = *(myEvent->pointerForak4pfjets_loose_pfid);
    myEvent->gen_daughter_n            = *(myEvent->pointerForgen_daughter_n);
    myEvent->gen_mother_index          = *(myEvent->pointerForgen_mother_index);
    myEvent->gen_id                    = *(myEvent->pointerForgen_id);
    myEvent->gen_eta                   = *(myEvent->pointerForgen_eta);
    myEvent->gen_index                 = *(myEvent->pointerForgen_index);
    myEvent->ak4pfjets_CSV             = *(myEvent->pointerForak4pfjets_CSV);
    myEvent->gen_m                     = *(myEvent->pointerForgen_m);
    */
}

#endif
