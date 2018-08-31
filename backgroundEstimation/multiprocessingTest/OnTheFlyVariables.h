#ifndef ON_THE_FLY_VARIABLES
#define ON_THE_FLY_VARIABLES
#include <exception>

struct OnTheFlyVariables
{
public:
   OnTheFlyVariables()
    :  
      QGtag_leadJet(-1),
      m_genWPt(0),
      m_genTopPt(0),
      m_genBarTopPt(0),
      m_genWdR(0),
      m_genTopdR(0),
      m_genBarTopdR(0),
      m_mumIdx(-1),
      m_recEffDef1(-1),
      m_recEffDef2(-1),
      m_genEffDef1(-1),
      m_genEffDef2(-1),
      m_matchedFakes(-1),
      m_matchedAndTaggedFakes(-1),      
      m_recEffDef1Counter(0),
      m_recEffDef2Counter(0),
      m_genEffDef1Counter(0),
      m_genEffDef2Counter(0),
      m_matchedFakesCounter(0),
      m_matchedAndTaggedFakesCounter(0)
     {
     }

public:
    float someCoolVariable;

    int SRbins;
    float QGtag_leadJet;

public:
    float m_genWPt;
    float m_genTopPt;
    float m_genBarTopPt;
    
    float m_genWdR;
    float m_genTopdR;
    float m_genBarTopdR;

    float m_recWPt;
    float m_mindR;

    //float m_recTopPt; //@MJ@ TODO not used

    int m_mumIdx;
    
    vector<int> m_daughtersIdx;

public:
    float m_ak8WMass;
    float m_ak8WPrunedMass;
    float m_ak8WTrimmedMass;
    float m_ak8WSoftDropMass;
    float m_ak8recoWPt;
    float m_ak8WSubJets;
    float m_ak8MET;
    float m_ak8WGenPt;

    float m_ak8WFakeMass;
    float m_ak8WFakePrunedMass;
    float m_ak8WFakeTrimmedMass;
    float m_ak8WFakeSoftDropMass;
    float m_ak8recoWFakePt;
    float m_ak8WFakeSubJets;
    float m_ak8FakeMET;
    float m_ak8WFakeGenPt;

public:
    int m_recEffDef1;
    int m_recEffDef2;
    int m_genEffDef1;
    int m_genEffDef2;
    int m_matchedFakes;
    int m_matchedAndTaggedFakes;
    
    float m_subjetinessReals;
    float m_subjetinessFakes;

    uint32_t m_recEffDef1Counter;
    uint32_t m_recEffDef2Counter;
    uint32_t m_genEffDef1Counter;
    uint32_t m_genEffDef2Counter;
    uint32_t m_matchedFakesCounter;
    uint32_t m_matchedAndTaggedFakesCounter;

public:
   
   float m_genRecoDiff;
   int m_lepMumId;
   float m_lepMumPt;
   int m_lepMumMumId;
   float m_lepMumMumPt;
   float m_lepPt;
   float m_neutrinoPt;
   float m_topPt; 


   float m_lepDeta;
   float m_lepDphi;
   Double_t m_lepDr;
   Double_t m_selDr;

   Double_t m_ljDphi;
   Double_t m_ljDr;

   Double_t m_systPt;

  float m_fakeLep;
  
  float m_ak10WPrunedMass;
  float m_ak10FakeWPrunedMass;
  float m_genW4Ak10;
  float m_recW4Ak10;
  float m_recMatchedW4Ak10;
  float m_recMatchedFakeW4Ak10;
  float m_subjetinessRealsAk10;
};

OnTheFlyVariables onTheFlyVariables;

void ComputeOnTheFlyVariables()
{
    onTheFlyVariables.someCoolVariable = 3.14;


    //-- Compute SRbin variables (MT2W/MET binning)
    onTheFlyVariables.SRbins = 0;

    if(myEvent.MT2W<200 && myEvent.MT2W>=0){
    	if(myEvent.pfmet>=250 && myEvent.pfmet<300) onTheFlyVariables.SRbins = -1;
	if(myEvent.pfmet>=300 && myEvent.pfmet<350) onTheFlyVariables.SRbins = -2;
    	if(myEvent.pfmet>=350 && myEvent.pfmet<400) onTheFlyVariables.SRbins = -3;
    	if(myEvent.pfmet>=400) onTheFlyVariables.SRbins = -4;
    }
    if(myEvent.MT2W>=200){
    	if(myEvent.pfmet>=250 && myEvent.pfmet<300) onTheFlyVariables.SRbins = 1;
	if(myEvent.pfmet>=300 && myEvent.pfmet<350) onTheFlyVariables.SRbins = 2;
    	if(myEvent.pfmet>=350 && myEvent.pfmet<400) onTheFlyVariables.SRbins = 3;
    	if(myEvent.pfmet>=400 && myEvent.pfmet<450) onTheFlyVariables.SRbins = 4;
    	if(myEvent.pfmet>=450 ) onTheFlyVariables.SRbins = 5;
    }
    //-- End SRbin computation

    onTheFlyVariables.QGtag_leadJet = -999;
    if(myEvent.ak4pfjets_qgtag.size()>0){
	onTheFlyVariables.QGtag_leadJet = myEvent.ak4pfjets_qgtag[0];
    }
}

#ifdef USE_GEN_INFO

//find mother particle with given index
//looking for particle which is decaying to two particles, which are not lepton
void findMotherAndDaughters(int idx)
{

    //initial clean-up
    onTheFlyVariables.m_mumIdx = -1;
    
    static uint32_t evt = 0;
    evt++;
    for(uint32_t k = 0; k < myEvent.gen_id.size(); k++)
    {
      if(evt == 1)
      {
         //cout << "gen particle id: " << myEvent.gen_id.at(k) << ", its index: " << myEvent.gen_index.at(k) << ", mother index: " << myEvent.gen_mother_index.at(k) << ", its status: " << myEvent.gen_status.at(k) << endl;
        if(abs(myEvent.gen_id.at(k)) == 24 && myEvent.gen_daughter_index.at(k).size() == 2)
        {
           int m = myEvent.gen_index.at(k);
           int d1  = myEvent.gen_daughter_index.at(k).at(0);
           int d2  = myEvent.gen_daughter_index.at(k).at(1);
           cout << "daughter1: " << d1 << endl;
           cout << "daughter2: " << d2 << endl;
    
           //mother propeties
           cout << "mother pt: " << myEvent.gen_pt.at(m) << ", mother eta: " << myEvent.gen_eta.at(m) << ", mother phi: " << myEvent.gen_phi.at(m) << endl;
           cout << "d1 pt: " << myEvent.gen_pt.at(d1) << ", d1 eta: " << myEvent.gen_eta.at(d1) << ", d1 phi: " << myEvent.gen_phi.at(d1) << endl;
           cout << "d2 pt: " << myEvent.gen_pt.at(d2) << ", d2 eta: " << myEvent.gen_eta.at(d2) << ", d2 phi: " << myEvent.gen_phi.at(d2) << endl;
          
        }
      }
    }
    
    //loop over particles and find mum candidates
    for(uint32_t i = 0; i < myEvent.gen_id.size(); i++)
    {


       if((onTheFlyVariables.m_mumIdx == -1) && (abs(myEvent.gen_id.at(i)) == abs(idx)) && myEvent.gen_daughter_index.at(i).size() == 2) // find W or t = mother particles
       {
           int dp1  = myEvent.gen_daughter_index.at(i).at(0);
           int dp2  = myEvent.gen_daughter_index.at(i).at(1);
           if(! ( (11 <= myEvent.gen_id.at(dp1) && myEvent.gen_id.at(dp1) <= 18) || (11 <= myEvent.gen_id.at(dp2) && myEvent.gen_id.at(dp2)<= 18) )) //@EC@ if not leptonic decay investigate other options
           {	
                onTheFlyVariables.m_mumIdx = myEvent.gen_index.at(i);
                onTheFlyVariables.m_daughtersIdx.push_back(dp1);
                onTheFlyVariables.m_daughtersIdx.push_back(dp2); 
           }
       }

       if(onTheFlyVariables.m_mumIdx == -1) 
           continue;
       if(abs(myEvent.gen_id.at(onTheFlyVariables.m_mumIdx)) != abs(idx))
       {
	   onTheFlyVariables.m_daughtersIdx.clear();
	   onTheFlyVariables.m_mumIdx = -1;
           continue;
       }


       //@EC@ there are the other options, I am using W->qq' now, if the mother particle, which was selected before does no fulfill any constion, the procesure is redone
       //W decaying to two quarks found
       if(idx == 24  && abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0))) <= 6 && abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1))) <= 6)
       {
               //cout << "W found" << endl;
               break;
       }
      //t decaying to W+ and sth found 
       else if (idx == 6 && (myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0)) == 24 || myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1)) == 24))
       {
               //cout << "top found" << endl;
               break;
       }
       
      //bar t decaying to W- and sth found 
       else if (idx == -6 && (myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0)) == -24 || myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1)) == -24))
       {
               //cout << "bar top found" << endl;
               break;
       }
       //this was not correct mother particle
       else
       {
               //cout << "mother condition not fulfilled" << endl;
	       onTheFlyVariables.m_daughtersIdx.clear();
	       onTheFlyVariables.m_mumIdx = -1;
	       
       }
    }
}
#endif


      //@MJ@ TODO what about ak10 jets?
      void fillRealAndFakes(float deltaR, uint32_t i)
      {

           //initialize vars for efficiencies to zero
           onTheFlyVariables.m_recEffDef1 = -1;
           onTheFlyVariables.m_recEffDef2 = -1;
           onTheFlyVariables.m_genEffDef1 = -1;
           onTheFlyVariables.m_genEffDef2 = -1;
           onTheFlyVariables.m_matchedFakes = -1;
           onTheFlyVariables.m_matchedAndTaggedFakes = -1; 
           onTheFlyVariables.m_subjetinessReals = -1;
           onTheFlyVariables.m_subjetinessFakes = -1;
    
           if(myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx) > 250)
           {
               onTheFlyVariables.m_genEffDef2 = 1;
           }
           //Ws
           if(abs(deltaR) < 0.1 && myEvent.ak8pfjets_pt.at(i) > 200) //@MJ@ TODO put the correct number
           {
               onTheFlyVariables.m_subjetinessReals = myEvent.ak8pfjets_tau2.at(i) / myEvent.ak8pfjets_tau1.at(i);
               
               //definition one efficiency
               if(myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx) > 250)
               {
                   onTheFlyVariables.m_genEffDef1 = 1;
                   if(myEvent.ak8pfjets_mass.at(i) > 80  && (myEvent.ak8pfjets_tau2.at(i) / myEvent.ak8pfjets_tau1.at(i)) < 0.5)
                   {
                       onTheFlyVariables.m_recEffDef1 = 1;
                       onTheFlyVariables.m_recEffDef2 = 1;
                   }
               }
               
               
               //@MJ@ TODO count subjets?
               onTheFlyVariables.m_ak8WMass = myEvent.ak8pfjets_mass.at(i);
               onTheFlyVariables.m_ak8WPrunedMass = myEvent.ak8pfjets_pruned_mass.at(i);
               onTheFlyVariables.m_ak8WTrimmedMass = myEvent.ak8pfjets_trimmed_mass.at(i);
               onTheFlyVariables.m_ak8WSoftDropMass = myEvent.ak8pfjets_softdrop_mass.at(i);
               onTheFlyVariables.m_ak8recoWPt = myEvent.ak8pfjets_pt.at(i);
               onTheFlyVariables.m_ak8WSubJets = myEvent.ak8pfjets_nSubJets.at(i); //for the event
               onTheFlyVariables.m_ak8WGenPt = myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx);
               onTheFlyVariables.m_ak8MET = myEvent.pfmet; //for the event
           
               onTheFlyVariables.m_ak8WFakeMass = -1;
               onTheFlyVariables.m_ak8WFakePrunedMass = -1;
               onTheFlyVariables.m_ak8WFakeTrimmedMass = -1;
               onTheFlyVariables.m_ak8WFakeSoftDropMass = -1;
               onTheFlyVariables.m_ak8recoWFakePt = -1;
               onTheFlyVariables.m_ak8WFakeSubJets = -1; //for the event
               onTheFlyVariables.m_ak8WFakeGenPt = -1;
               onTheFlyVariables.m_ak8FakeMET = -1; //for the event
           }
           //fakes
           else if(abs(deltaR) > 0.8 && myEvent.ak8pfjets_pt.at(i) > 200)
           {
               onTheFlyVariables.m_subjetinessFakes = myEvent.ak8pfjets_tau2.at(i) / myEvent.ak8pfjets_tau1.at(i);
               //definition of fake rate
               onTheFlyVariables.m_matchedFakes = 1;
               if(myEvent.ak8pfjets_mass.at(i) > 80  && (myEvent.ak8pfjets_tau2.at(i) / myEvent.ak8pfjets_tau1.at(i)) < 0.5)
               {
                   onTheFlyVariables.m_matchedAndTaggedFakes = 1;
               }
 
               onTheFlyVariables.m_ak8WFakeMass = myEvent.ak8pfjets_mass.at(i);
               onTheFlyVariables.m_ak8WFakePrunedMass = myEvent.ak8pfjets_pruned_mass.at(i);
               onTheFlyVariables.m_ak8WFakeTrimmedMass = myEvent.ak8pfjets_trimmed_mass.at(i);
               onTheFlyVariables.m_ak8WFakeSoftDropMass = myEvent.ak8pfjets_softdrop_mass.at(i);
               onTheFlyVariables.m_ak8recoWFakePt = myEvent.ak8pfjets_pt.at(i);
               onTheFlyVariables.m_ak8WFakeSubJets = myEvent.ak8pfjets_nSubJets.at(i); //for the event
               onTheFlyVariables.m_ak8WFakeGenPt = myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx);
               onTheFlyVariables.m_ak8FakeMET = myEvent.pfmet; //for the event
           
               onTheFlyVariables.m_ak8WMass = -1;
               onTheFlyVariables.m_ak8WPrunedMass = -1;
               onTheFlyVariables.m_ak8WTrimmedMass = -1;
               onTheFlyVariables.m_ak8WSoftDropMass = -1;
               onTheFlyVariables.m_ak8recoWPt = -1;
               onTheFlyVariables.m_ak8WSubJets = -1; //for the event
               onTheFlyVariables.m_ak8WGenPt = -1;
               onTheFlyVariables.m_ak8MET = -1; //for the event
           }
           //nothing found
           else
           {
               //cout << "W boosted jet not found" << endl;
           }
       }

      //ak10
      void fillRealAndFakesAk10(float deltaR, uint32_t i)
      {
         
           //initialize vars for efficiencies to zero
           onTheFlyVariables.m_genW4Ak10 = -1;
           onTheFlyVariables.m_recW4Ak10 = -1;
           onTheFlyVariables.m_recMatchedW4Ak10 = -1;
           onTheFlyVariables.m_recMatchedFakeW4Ak10 = -1;
           onTheFlyVariables.m_subjetinessRealsAk10 = -1;
    
           //gen W
           if(myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx) > 250)
           {
               onTheFlyVariables.m_genW4Ak10 = 1;
           }
           
           onTheFlyVariables.m_subjetinessRealsAk10 = myEvent.ak10pfjets_tau2.at(i) / myEvent.ak10pfjets_tau1.at(i);
           if(myEvent.ak10pfjets_pruned_mass.at(i) > 60 && myEvent.ak10pfjets_pruned_mass.at(i) < 100 &&  myEvent.ak10pfjets_pt.at(i) > 200 && onTheFlyVariables.m_subjetinessRealsAk10 < 0.5)
           {
               onTheFlyVariables.m_recW4Ak10 = 1;
               //Ws
               if(abs(deltaR) < 0.1 && myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx) > 250)
               {
                   //real
                   onTheFlyVariables.m_recMatchedW4Ak10 = 1;
                   onTheFlyVariables.m_ak10WPrunedMass = myEvent.ak10pfjets_pruned_mass.at(i);
               }
               else if(abs(deltaR) > 0.1 && myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx) > 250)
               {
                   //fake
                   onTheFlyVariables.m_recMatchedFakeW4Ak10 = 1;
                   onTheFlyVariables.m_ak10FakeWPrunedMass = myEvent.ak10pfjets_pruned_mass.at(i);
               }
               
           }
           //nothing found
           else
           {
               //cout << "W boosted jet not found" << endl;
           }
       }
#ifdef USE_GEN_INFO

void findGenParticleProps(int idx, float* genPt, float* gendR)
{ 
    //do initial clean up
    onTheFlyVariables.m_mumIdx = -1;
    *genPt = -1;
    *gendR = -1;
    onTheFlyVariables.m_daughtersIdx.clear();
    
    //loop over the generated paricles, find mother and its daughters
    findMotherAndDaughters(idx);

    //no mother found
    if(onTheFlyVariables.m_mumIdx == -1)
    {
        //cout << "no mother particle was found! " << "mum Id is: " << onTheFlyVariables.m_mumIdx << endl;
    }

    //fill the kinematic variables
    if(onTheFlyVariables.m_mumIdx != -1)
    {
        *genPt = myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx);
    
        float deta = myEvent.gen_eta.at(onTheFlyVariables.m_daughtersIdx.at(0)) - myEvent.gen_eta.at(onTheFlyVariables.m_daughtersIdx.at(1));
        float dphi = myEvent.gen_phi.at(onTheFlyVariables.m_daughtersIdx.at(0)) - myEvent.gen_phi.at(onTheFlyVariables.m_daughtersIdx.at(1));
        *gendR =  TMath::Sqrt( deta*deta+dphi*dphi );
    }
    else
    {
        *genPt = -1;
        *gendR = -1;
    }
    
    //final cleanup
    onTheFlyVariables.m_daughtersIdx.clear();
    onTheFlyVariables.m_mumIdx = -1; 
}

void minJetdRToGenParticle(int idx, float* mindR)
{

    onTheFlyVariables.m_mumIdx = -1;
    *mindR = -130;
    onTheFlyVariables.m_daughtersIdx.clear();
    
    //loop over the generated paricles, find mother and its daughters
    findMotherAndDaughters(idx);

    //fill the kinematic variables
    if(onTheFlyVariables.m_mumIdx != -1 && myEvent.gen_eta.at(onTheFlyVariables.m_mumIdx) > 250)
    { 
        for(uint32_t j = 0; j < myEvent.ak8pfjets_mass.size(); j++)
        {  
            if(myEvent.ak8pfjets_pt.at(j) < 200)
                continue;  
            float deta = myEvent.gen_eta.at(onTheFlyVariables.m_mumIdx) - myEvent.ak8pfjets_eta.at(j);
            float dphi = myEvent.gen_phi.at(onTheFlyVariables.m_daughtersIdx.at(0)) - myEvent.ak8pfjets_phi.at(j);
            float dR =  TMath::Sqrt( deta*deta+dphi*dphi );
            
            if(abs(dR) < abs(*mindR))
                *mindR = dR;
        }
    }
}
#endif


#ifdef USE_GEN_INFO

//count efficiency of particle tagging
void countEfficiency(int idx, string cutName = "", float lowCut = -1, float upCut = -1)
{
  //initial clean up
  onTheFlyVariables.m_daughtersIdx.clear();
  onTheFlyVariables.m_ak8WMass = -1;
  onTheFlyVariables.m_ak8WPrunedMass = -1;
  onTheFlyVariables.m_ak8WTrimmedMass = -1;
  onTheFlyVariables.m_ak8WSoftDropMass = -1;
  onTheFlyVariables.m_ak8recoWPt = -1;
  onTheFlyVariables.m_ak8WSubJets = -1; //for the event
  onTheFlyVariables.m_ak8MET = -1;
  onTheFlyVariables.m_recWPt = -1;
  onTheFlyVariables.m_ak8WGenPt = -1;
 
  onTheFlyVariables.m_ak8WFakeMass = -1;
  onTheFlyVariables.m_ak8WFakePrunedMass = -1;
  onTheFlyVariables.m_ak8WFakeTrimmedMass = -1;
  onTheFlyVariables.m_ak8WFakeSoftDropMass = -1;
  onTheFlyVariables.m_ak8recoWFakePt = -1;
  onTheFlyVariables.m_ak8WFakeSubJets = -1; //for the event
  onTheFlyVariables.m_ak8FakeMET = -1;
  onTheFlyVariables.m_ak8WFakeGenPt = -1;
  
   //choose on which variable cut
   vector<float> masses(myEvent.ak8pfjets_mass.size(), -1.0);
   if(cutName == "raw")
       masses = myEvent.ak8pfjets_mass;
   else if(cutName == "pruned")
       masses = myEvent.ak8pfjets_pruned_mass;
   else if(cutName == "trimmed")
       masses = myEvent.ak8pfjets_trimmed_mass;
   else if(cutName == "softdrop")
       masses = myEvent.ak8pfjets_softdrop_mass;
   else
   {
       //cout << "no cut" << endl;
   }

  uint32_t nrOfJets = myEvent.ak8pfjets_mass.size();
  
  for(uint32_t i = 0; i < nrOfJets; i++ )
  {

       //@MJ@ TODO add another info about W tagged jet in order to discriminate fakes rate
       // e.g. sufficient mass of jet, number of subjets, N-subjettiness, PUPPI
       // find mother particle
       findMotherAndDaughters(idx);
       float deta = -13;
       float dphi = -13;
       float dR = -13;
      
       //count dR
       if(onTheFlyVariables.m_mumIdx != -1)
       {
           deta = myEvent.ak8pfjets_eta.at(i) - myEvent.gen_eta.at(onTheFlyVariables.m_mumIdx);
           dphi = myEvent.ak8pfjets_phi.at(i) - myEvent.gen_phi.at(onTheFlyVariables.m_mumIdx);
           dR = TMath::Sqrt( deta*deta+dphi*dphi );
       }
       else
          continue;

       //perform cuts and filling fakes/reals
       if(cutName == "")
       {
           fillRealAndFakes(dR, i);
           if(onTheFlyVariables.m_ak8WMass == myEvent.ak8pfjets_mass.at(i)) //end the loop when W jet found //@MJ@ TODO find better solution
               break;
       }
       else if(cutName != "" && lowCut != -1 && upCut != -1)
       {
           if(masses.at(i) > lowCut && masses.at(i) < lowCut)
               fillRealAndFakes(dR, i);
           if(onTheFlyVariables.m_ak8WMass == myEvent.ak8pfjets_mass.at(i)) //end the loop when W jet found
               break;
       }
       else if(cutName != "" && (lowCut != -1 || upCut != -1))
       {
           if(lowCut != -1 && masses.at(i) > lowCut)
               fillRealAndFakes(dR, i);
           if(upCut != -1 && masses.at(i) < upCut)
               fillRealAndFakes(dR, i);
           if(onTheFlyVariables.m_ak8WMass == myEvent.ak8pfjets_mass.at(i))  //end the loop when W jet found
               break;
       }
      else
      {
          cout << "cuts were defined in wrog way!" << endl;
      }
  }
  
}

//count efficiency of particle tagging AK10
void countEffAndFRAK10(int idx)  //@EC@ this is the method  you can use for matching gen particle to jets (or at least part of it)
{
  //initial clean up
  onTheFlyVariables.m_daughtersIdx.clear();

  onTheFlyVariables.m_ak10WPrunedMass = -1;
  onTheFlyVariables.m_ak10FakeWPrunedMass = -1;
  
   //choose on which variable cut
  vector<float> masses = myEvent.ak10pfjets_pruned_mass;

  uint32_t nrOfJets = myEvent.ak10pfjets_pruned_mass.size();
  
  for(uint32_t i = 0; i < nrOfJets; i++ )
  {
       //clear
       onTheFlyVariables.m_mumIdx = -1;

       // find mother particle
       findMotherAndDaughters(idx); //@EC@ checks if there is generated particle with given index idx, if you look in the code of this method I also specify the decays of this particle. I do not know, what you want to check, but maybe you will have to implement your own part (I did some develpement for W and top). If it findes the particle specified, in onTheFlyVariables.m_mumIdx its index is filled
       TLorentzVector j1;
       TLorentzVector g2;
       double dR = -13;
      
       //count dR
       if(onTheFlyVariables.m_mumIdx != -1) //@EC@ if you found gen particle compute dR
       {
           j1.SetPtEtaPhiM(myEvent.ak10pfjets_pt.at(i), myEvent.ak10pfjets_eta.at(i), myEvent.ak10pfjets_phi.at(i), myEvent.ak10pfjets_pruned_mass.at(i));
           g2.SetPtEtaPhiM(myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx), myEvent.gen_eta.at(onTheFlyVariables.m_mumIdx), myEvent.gen_phi.at(onTheFlyVariables.m_mumIdx), myEvent.gen_m.at(onTheFlyVariables.m_mumIdx));
           dR = j1.DeltaR(g2);
       }
       else
          continue;

       fillRealAndFakesAk10(dR, i); //@EC@ here you fill, real fakes, number of generated (for purity, efficiency, fake rate)
       if(onTheFlyVariables.m_ak10WPrunedMass == myEvent.ak10pfjets_pruned_mass.at(i)) //@EC@ end the loop when matched W found
           break;
  }
  
}

float calculateDPhi()
{
    float dPhi1 = myEvent.jet_phi.at(0) - myEvent.pfmet_phi;
    float dPhi2 = myEvent.jet_phi.at(1) - myEvent.pfmet_phi;
    
    return (abs(dPhi1) < abs(dPhi2)) ? dPhi1 : dPhi2;
}

void findLeptonOrigin()
{
    onTheFlyVariables.m_genRecoDiff = -1;
    onTheFlyVariables.m_lepMumId = -1;
    onTheFlyVariables.m_lepMumPt = -1;
    onTheFlyVariables.m_lepMumMumId = -1;
    onTheFlyVariables.m_lepMumMumPt = -1;
    onTheFlyVariables.m_lepPt = -1;
    onTheFlyVariables.m_neutrinoPt = -1;

    //onTheFlyVariables.m_genRecoDiff = myEvent.pfmet - myEvent.genMET;

    for(uint32_t p = 0; p < myEvent.gen_index.size(); p++)
    {
        uint32_t idx = myEvent.gen_index.at(p);
        //cout << "in here 1" << endl;
        if(abs(myEvent.gen_id.at(idx)) == 11 || abs(myEvent.gen_id.at(idx)) == 13 || abs(myEvent.gen_id.at(idx)) == 15 )
        {
      //  cout << "in here 2" << endl;
            uint32_t motherIdx = myEvent.gen_mother_index.at(idx);
            onTheFlyVariables.m_lepMumId = myEvent.gen_id.at(motherIdx);
            onTheFlyVariables.m_lepMumPt = myEvent.gen_pt.at(motherIdx);
            uint32_t motherMotherIdx = myEvent.gen_mother_index.at(motherIdx);
            onTheFlyVariables.m_lepMumMumId = myEvent.gen_id.at(motherMotherIdx);
            onTheFlyVariables.m_lepMumMumPt = myEvent.gen_pt.at(motherMotherIdx);
            
            if(abs(onTheFlyVariables.m_lepMumId) == 24)
                onTheFlyVariables.m_topPt = myEvent.gen_pt.at(motherMotherIdx);
            //@MJ@ TODO selected lepton?
            
            onTheFlyVariables.m_lepPt = myEvent.gen_pt.at(idx);
            //onTheFlyVariables.m_lepPt = myEvent.lep1_pt;
            //cout << "gen lep pt " << myEvent.gen_pt.at(idx) << " lep1 pt " << myEvent.lep1_pt << endl;

       // cout << "in here 3" << endl;
            if(abs(onTheFlyVariables.m_lepMumId) == 24)
            {
       // cout << "in here 4" << endl;
                for(uint32_t d = 0; d < myEvent.gen_daughter_index.at(idx).size(); d++)
                {
       // cout << "in here 5" << endl;
                    uint32_t daughterIdx = myEvent.gen_daughter_index.at(idx).at(d);
       // cout << "in here 5.5" << endl;
                    if(abs(myEvent.gen_id.at(daughterIdx) == 12) || abs(myEvent.gen_id.at(daughterIdx) == 14) || abs(myEvent.gen_id.at(daughterIdx) == 16))
                    {
        //cout << "in here 6" << endl;
                        onTheFlyVariables.m_neutrinoPt = myEvent.gen_pt.at(daughterIdx); 
                        break;
                    }
                }
            }
            break;
        }
    }
}

void closeLepAndJet()
{
        onTheFlyVariables.m_ljDphi = -13;
        onTheFlyVariables.m_ljDr =   -13;
        
        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;

        if(myEvent.jet_eta.size() == 1 && myEvent.ngoodleps == 2)
        {

            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));

            Double_t dR11 = l1.DeltaR(j1);
            Double_t dR21 = l2.DeltaR(j1);
        
            onTheFlyVariables.m_ljDr = dR11 < dR21 ? dR11 : dR21;

           Double_t dPhi1 = l1.DeltaPhi(j1);
           Double_t dPhi2 = l2.DeltaPhi(j1);

           onTheFlyVariables.m_ljDphi = abs(dPhi1) < abs(dPhi2) ? dPhi1 : dPhi2;

        }


}


void countLepDeltas()
{
        //cout << "in countLepDelats" << endl;
        onTheFlyVariables.m_lepDeta = -13;
        onTheFlyVariables.m_lepDphi = -13;
        onTheFlyVariables.m_lepDr =   -13;
        onTheFlyVariables.m_selDr =   -13;

        onTheFlyVariables.m_lepDeta = myEvent.lep1_eta - myEvent.lep2_eta;
        onTheFlyVariables.m_lepDphi = myEvent.lep1_phi - myEvent.lep2_phi;
        onTheFlyVariables.m_lepDr =  TMath::Sqrt( (onTheFlyVariables.m_lepDeta*onTheFlyVariables.m_lepDeta)+(onTheFlyVariables.m_lepDphi*onTheFlyVariables.m_lepDphi) );

        vector<uint32_t> bId;
        bId.clear();



        //cout << "in here 1" << endl;
        for(uint32_t j = 0; j < myEvent.ak4pfjets_partonFlavour.size(); j++)
        {
            if(abs(myEvent.ak4pfjets_partonFlavour.at(j)) == 5)
        //cout << "in here 2" << endl;
            {
               bId.push_back(j);        
            }
        }
        if(bId.size() == 2)

        if(myEvent.jet_eta.size() == 2 && myEvent.ngoodleps == 2)
        {
        //cout << "in here 3" << endl;
        /*    float dR11 = (myEvent.lep1_eta - myEvent.ak4pfjets_eta.at(bId.at(0)))*(myEvent.lep1_eta - myEvent.ak4pfjets_eta.at(bId.at(0))) +(myEvent.lep1_phi - myEvent.ak4pfjets_phi.at(bId.at(0)))*(myEvent.lep1_phi - myEvent.ak4pfjets_phi.at(bId.at(0)));
            float dR12 = (myEvent.lep1_eta - myEvent.ak4pfjets_eta.at(bId.at(1)))*(myEvent.lep1_eta - myEvent.ak4pfjets_eta.at(bId.at(1))) +(myEvent.lep1_phi - myEvent.ak4pfjets_phi.at(bId.at(1)))*(myEvent.lep1_phi - myEvent.ak4pfjets_phi.at(bId.at(1)));
            float dR22 = (myEvent.lep2_eta - myEvent.ak4pfjets_eta.at(bId.at(1)))*(myEvent.lep2_eta - myEvent.ak4pfjets_eta.at(bId.at(1))) +(myEvent.lep2_phi - myEvent.ak4pfjets_phi.at(bId.at(1)))*(myEvent.lep2_phi - myEvent.ak4pfjets_phi.at(bId.at(1)));
            float dR21 = (myEvent.lep2_eta - myEvent.ak4pfjets_eta.at(bId.at(0)))*(myEvent.lep2_eta - myEvent.ak4pfjets_eta.at(bId.at(0))) +(myEvent.lep2_phi - myEvent.ak4pfjets_phi.at(bId.at(0)))*(myEvent.lep2_phi - myEvent.ak4pfjets_phi.at(bId.at(0)));
*/

        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector j2;


            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            j2.SetPtEtaPhiM(myEvent.jet_pt.at(1), myEvent.jet_eta.at(1), myEvent.jet_phi.at(1), myEvent.jet_mass.at(1));

            //@MJ@ TODO recompute with tlorentz vector def constr and then set
            Double_t dR11 = l1.DeltaR(j1);
            //cout  << "here 1" << endl;
            Double_t dR12 = l1.DeltaR(j2);
            //cout  << "here 2" << endl;
            Double_t dR21 = l2.DeltaR(j1);
            //cout  << "here 3" << endl;
            Double_t dR22 = l2.DeltaR(j2);
            //cout  << "here 4" << endl;
            


            Double_t win1 = abs(dR11) < abs(dR12) ? dR11 : dR12;
            Double_t win2 = abs(win1) < abs(dR21) ? win1 : dR21;
            Double_t win3 = abs(win2) < abs(dR22) ? win2 : dR22;

            if(win3 == dR11)
                onTheFlyVariables.m_selDr = dR22;
            else if(win3 == dR12)
                onTheFlyVariables.m_selDr = dR21;
            else if(win3 == dR21)
                onTheFlyVariables.m_selDr = dR12;
            else if(win3 == dR22)
                onTheFlyVariables.m_selDr = dR11;
            else
                cout << "something wrong in countLeptonDelates method" << endl;
       
        }

}

void countAllPt()
{
        TLorentzVector l1;
        TLorentzVector l2;
        TLorentzVector j1;
        TLorentzVector met;
        TLorentzVector sum;
       
        if(myEvent.jet_eta.size() == 1 && myEvent.ngoodleps == 2)
        {           
 
            // @MJ@ TODO eta
            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            l2.SetPtEtaPhiM(myEvent.lep2_pt, myEvent.lep2_eta, myEvent.lep2_phi, myEvent.lep2_mass);
            j1.SetPtEtaPhiM(myEvent.jet_pt.at(0), myEvent.jet_eta.at(0), myEvent.jet_phi.at(0), myEvent.jet_mass.at(0));
            met.SetPtEtaPhiE(myEvent.pfmet, 0, myEvent.pfmet_phi, myEvent.pfmet);

            sum = l1 + l2 + j1 + met;
            onTheFlyVariables.m_systPt = sum.Pt();

        }

}

uint32_t nrOfWTaggs() //@MJ@ TODO is this sufficient
{
    uint32_t w = 0;
    for(uint32_t idx = 0; idx < myEvent.ak8pfjets_mass.size(); idx++)
    {
         if(myEvent.ak8pfjets_mass.at(idx) > 80 && myEvent.ak8pfjets_pt.at(idx) > 200)
             w++;
    }
    return w;
}

void misidentifiedLepton()
{
    uint32_t leps = 0;
    uint32_t idx = 3333;
    uint32_t i = 4444;
    uint32_t j = 5555;
    onTheFlyVariables.m_fakeLep = -13;
    onTheFlyVariables.m_lepDr = 0;
    for(uint32_t p = 0; p < myEvent.gen_index.size(); p++)
    {
        i = myEvent.gen_index.at(p);
        j = myEvent.gen_mother_index.at(i);
        if((abs(myEvent.gen_id.at(i)) == 11 || abs(myEvent.gen_id.at(i)) == 13 || abs(myEvent.gen_id.at(i)) == 15) && (abs(myEvent.gen_id.at(j)) == 24) && (myEvent.gen_status.at(i) <=3)) //sth fishy in here, statuses are weit\rd, much more leps than nr of Gen leps
        {
           idx = myEvent.gen_index.at(p);
           //cout << "motherIdx" <<myEvent.gen_mother_index.at(i) << endl;
           //cout << "motherId" <<myEvent.gen_id.at(j) << endl;
           //cout << "index " << i << "i " << p << endl;
           //cout << "status" << myEvent.gen_status.at(i) << endl;
           //break;
           leps++;
        }
   }
 
   //if(leps != myEvent.numberOfGeneratedLeptons) cout << "this is nonsense" << endl;
   if(myEvent.ngoodleps == 1 && myEvent.numberOfGeneratedLeptons == 0)
   {
   //cout << "leps" << leps << endl;
   //cout << "nrOfGenleps" << myEvent.numberOfGeneratedLeptons << endl;
   //cout << "nrOfGoodleps" << myEvent.ngoodleps << endl;
   }
   if(leps > 1) return;
   if(leps != 0 && idx == 3333) return;
   if(myEvent.ngoodleps != 1) return;


   Double_t dR = 3333;
   
   if(idx != 3333)
   {
   TLorentzVector l;
   TLorentzVector g;
   l.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
   g.SetPtEtaPhiM(myEvent.gen_pt.at(idx), myEvent.gen_eta.at(idx), myEvent.gen_phi.at(idx), myEvent.gen_m.at(idx));
   dR = l.DeltaR(g);
   }
    
   onTheFlyVariables.m_lepDr = dR;

   if(leps == 0)
   {
        onTheFlyVariables.m_fakeLep = 0.5; //fake hadronic
   }
   else if(abs(myEvent.gen_id.at(idx)) != abs(myEvent.lep1_pdgid) || abs(dR) > 0.1)
   {
        onTheFlyVariables.m_fakeLep = 1.5; //fakes
   }
   else
   { 
       onTheFlyVariables.m_fakeLep = 2.5; //reals
   }
}

bool evaluateLeptonAndAk10Overlap()
{

        TLorentzVector l1;
        TLorentzVector j1;

        bool overlap = false;

        for(uint32_t idx = 0; idx < myEvent.ak10pfjets_pruned_mass.size(); idx++)
        {

            if(idx > 0 && overlap == false)
                break;

            l1.SetPtEtaPhiM(myEvent.lep1_pt, myEvent.lep1_eta, myEvent.lep1_phi, myEvent.lep1_mass);
            j1.SetPtEtaPhiM(myEvent.ak10pfjets_pt.at(idx), myEvent.ak10pfjets_eta.at(idx), myEvent.ak10pfjets_phi.at(idx), myEvent.ak10pfjets_pruned_mass.at(idx));

            Double_t dR = l1.DeltaR(j1);

            overlap = dR < 1.0 ? true : false;
        }

        return overlap;

}


#endif
float returnSigCS(float stopm)
{
    if(stopm == 100 ) return 1521.11 ;
    else if(stopm == 125) return 574.981;
    else if(stopm == 150) return 249.409;
    else if(stopm == 175) return 121.416;
    else if(stopm == 200) return 64.5085;
    else if(stopm == 225) return 36.3818;
    else if(stopm == 250) return 21.5949;
    else if(stopm == 275) return 13.3231;
    else if(stopm == 300) return 8.51615;
    else if(stopm == 325) return 5.60471;
    else if(stopm == 350) return 3.78661;
    else if(stopm == 375) return 2.61162;
    else if(stopm == 400) return 1.83537;
    else if(stopm == 425) return 1.31169;
    else if(stopm == 450) return 0.948333;
    else if(stopm == 475) return 0.697075;
    else if(stopm == 500) return 0.51848;
    else if(stopm == 525) return 0.390303;
    else if(stopm == 550) return 0.296128;
    else if(stopm == 600) return 0.174599;
    else if(stopm == 650) return 0.107045;
    else if(stopm == 700) return 0.0670476;
    else if(stopm == 750) return 0.0431418;
    else if(stopm == 800) return 0.0283338;
    else if(stopm == 850) return 0.0189612;
    else if(stopm == 900) return 0.0128895;
    else if(stopm == 950) return 0.00883465;
    else 
    {
        std::cout << "stop mass: " << stopm << std::endl;
        throw std::runtime_error("no cross section for given stop mass"); //@MJ@ TODO make this working!
    }

}



#endif
