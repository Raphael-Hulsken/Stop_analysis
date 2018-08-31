#include "ActionsOnHistos.h"

using namespace std;

    Double_t ActionsOnHistos::div32jets(TH1* hist)
    {
        //bins 0,1,2, jets
        Double_t num = hist->GetBinContent(4);
        Double_t denom = hist->GetBinContent(3);
        return num/denom;
    }

    Double_t ActionsOnHistos::div42jets(TH1* hist)
    {
        //bins 0,1,2, jets
        Double_t num = hist->GetBinContent(5);
        Double_t denom = hist->GetBinContent(3);
        return num/denom;
    }

    Double_t ActionsOnHistos::predictionInBinI(TH1* MCSR, TH1* MCCR, TH1* DATACR, uint32_t bin)
    {
        Double_t num = MCSR->GetBinContent(bin);
    }

