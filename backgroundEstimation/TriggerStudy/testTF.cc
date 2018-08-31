#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../../common/TFFactory.h"

#include <iostream>

using namespace std;

int main(){

        //vector<string> label;
        vector<string> tfreg;
        tfreg.push_back("standard");

        TFProducer prod(tfreg);
        prod.produceTFTable("yieldNew.tab", "testOuti2");

}
