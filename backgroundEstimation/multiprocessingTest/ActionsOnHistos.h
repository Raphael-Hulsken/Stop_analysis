#ifndef ActionsOnHistos_h 
#define ActionsOnHistos_h

#include <iostream>
#include "TH1.h"

using namespace std;

class ActionsOnHistos
{
public:
     virtual Double_t div32jets(TH1* hist);
     virtual Double_t div42jets(TH1* hist);
     virtual Double_t predictionInBinI(TH1* MCSR, TH1* MCCR, TH1* DATACR, uint32_t bin);
};

#endif //ActionsOnHistos_h

