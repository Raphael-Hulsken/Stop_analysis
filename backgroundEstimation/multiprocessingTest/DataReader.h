#ifndef DataReader_h 
#define DataReader_h

#include <iostream>
#include <vector>
#include "TString.h"

using namespace std;

class DataReader
{
public:
     void readData(TString input);
 
     void fillInfo();

public:
     vector<TString> inputRootFile;
     vector<TString> variable;
     vector<TString> inputCanvasName;
     vector<TString> inputDir;
     vector<TString> inputLegend;

private:
     vector<TString> in_stream;
};

#endif //DataReader_h

