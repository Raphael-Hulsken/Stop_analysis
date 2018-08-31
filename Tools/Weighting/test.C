#include "BTagCalibrationStandalone.cc"


int main(){

// setup calibration + reader
//BTagCalibration calib("CSVv2", "files/CSVv2_ichep.csv");
BTagCalibration* calib = new BTagCalibration("CSV", "files/CSV_13TEV_Combined_14_7_2016.csv");
BTagCalibrationReader reader(BTagEntry::OP_LOOSE,  // operating point
                             //"central",             // central sys type
                             "central");             // central sys type
     //                        {"up", "down"});      // other sys types

reader.load(*calib,                // calibration instance
            BTagEntry::FLAV_B,    // btag flavour
            //"comb");               // measurement type
            "fastsim");               // measurement type
}
