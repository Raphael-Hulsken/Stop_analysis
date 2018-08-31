#include <fstream>
#include "DataReader.h"

void DataReader::readData(TString input)
{
    ifstream is;
    is.open(input);
    if( is.is_open() )
    {
        char * sline;
        TString tsline("");
        Int_t line = 0;
        string strline;
        while(getline(is, strline)) {
                tsline = strline;
                if(tsline.BeginsWith(" ")) continue;
                if(tsline.BeginsWith("%")) continue;
                if(tsline.Length()==0)     continue;
                in_stream.push_back(tsline);
                line++;
        }
        is.close();
    } 
}

void DataReader::fillInfo()
{
    uint32_t nrOfVars = in_stream.at(0).Atof();
    for(uint32_t f = 0; f < nrOfVars; f++)
    {
        variable.push_back(in_stream.at(1+f));
        inputRootFile.push_back(in_stream.at(1+nrOfVars+f));
        inputCanvasName.push_back(in_stream.at(1+2*nrOfVars+f));
        inputDir.push_back(in_stream.at(1+3*nrOfVars+f));
        inputLegend.push_back(in_stream.at(1+4*nrOfVars+f));
    }
}
