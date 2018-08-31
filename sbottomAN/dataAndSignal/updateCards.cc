#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

#include "../../sonicScrewdriver/interface/tables/CombineCardMaker.h"

using namespace std;
using namespace theDoctor;

int main()
{

    theDoctor::CombineCardMaker card;
    card.UpdateCardTable("card.tab", "bkgLostLeptonTF.txt", "bkgLostLepton");
    card.UpdateCardTable("card.tab", "bkgOneLepFromWTF.txt", "bkgOneLepFromW");
    card.ProduceCard("card.tab" ,"../cards");
    //totalSM recompute


    return 0;
}
