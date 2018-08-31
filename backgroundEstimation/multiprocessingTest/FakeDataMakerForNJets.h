#ifndef FakeDataMakerForNJets_h 
#define FakeDataMakerForNJets_h

#include "FakeDataMakerIfc.h"

class FakeDataMakerForNJets: public IFakeDataMaker
{
public:
     virtual TH1* smear(TH1* hist, void* data, uint32_t size)
     {
         Double_t* weights =  static_cast<Double_t*>(data);

         if(size != hist->GetNbinsX())
         {
             cout << "ERROR: incorect number of scale factors" << endl;
             return NULL;
         }
         
         for(uint32_t h = 0; h < size; h++)
         {
             Double_t binCont = hist->GetBinContent(h+1);
             hist->SetBinContent((binCont)*((*weights)[h]));
         }

         return hist;

     }

     virtual float recomputeMET(void* met){cout << "not implemented" << endl;}

     virtual int recomputeNJet(void* nJet){cout << "not implemented" << endl;}
};


#endif // FakeDataMakerForNJets_h

