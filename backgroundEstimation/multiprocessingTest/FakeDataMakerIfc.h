#ifndef FakeDataMakerIfc_h 
#define FakeDataMakerIfc_h

#include <iostream>
#include "TH1.h"

using namespace std;

class IFakeDataMaker
{
public:
     virtual TH1* smear(TH1* hist, void* data, uint32_t size) = 0;
    
     virtual float recomputeMET(void* met) = 0;
   
     virtual int recomputeNJet(void* nJet) = 0;
};

#endif //FakeDataMakerIfc_h
