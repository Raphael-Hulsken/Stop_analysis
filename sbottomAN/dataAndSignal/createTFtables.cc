#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "../TFFactory.h"

using namespace std;

int main()
{

    vector<string> tfreg = { "_A_250lessMETless350" , "_A_350lessMETless450" , "_A_450lessMETless600" , "_A_600lessMETlessInf" , "_B_250lessMETless450" , "_B_450lessMETless600" , "_B_600lessMETlessInf" , "_C_250lessMETless350" , "_C_350lessMETless450" , "_C_450lessMETless550" , "_C_550lessMETless650" , "_C_650lessMETlessInf" , "_D_250lessMETless350" , "_D_350lessMETless450" , "_D_450lessMETless550" , "_D_550lessMETlessInf" , "_E_250lessMETless350" , "_E_350lessMETless550" , "_E_550lessMETlessInf" , "_F_250lessMETless450" , "_F_450lessMETlessInf" , "_G_250lessMETless350" , "_G_350lessMETless450" , "_G_450lessMETless600" , "_G_600lessMETlessInf" , "_H_250lessMETless450" , "_H_450lessMETlessInf" , "_I_250lessMETless350" , "_I_350lessMETless450" , "_I_450lessMETless550" , "_I_550lessMETlessInf"   };



    vector<string> tfregextr = { "_B_450lessMETless600" , "_E_350lessMETless550" , "_F_250lessMETless450" , "_H_250lessMETless450"  };//only first of the two

    TFProducer prod(tfreg, "bkgLostLepton", "CR2l");
    prod.produceTFTable("sbottomSignal.tab", "yieldsCRs.tab", "bkgLostLeptonTF");
    prod.extrapolateTFTable("bkgLostLeptonTF", "sbottomSignal.tab", "yieldsCRs.tab", tfregextr);
    TFProducer prod2(tfreg, "bkgOneLepFromW", "CR0b");
    prod2.produceTFTable("sbottomSignal.tab", "yieldsCRs.tab" , "bkgOneLepFromWTF");

    return 0;
}
