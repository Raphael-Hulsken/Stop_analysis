#include <exception>
#include <iostream>
#include "stdlib.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFrame.h"
#include "TStyle.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"
#include "interface/Figure.h"

using namespace std;
using namespace theDoctor;

int main()

{

   TF1 *fa2 = new TF1("fa2","TMath::LogNormal(x,0.4)",0,10);
   fa2->Draw();

   TF1 *fa1 = new TF1("fa2","TMath::Poisson(x,5)",0,10);
   fa1->Draw();

}
