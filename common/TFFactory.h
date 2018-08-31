#ifndef TFFactory_h
#define TFFactory_h
// C/C++ headerrs
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <exception>

//SS headers
#include "../sonicScrewdriver/interface/Table.h"
#include "../sonicScrewdriver/interface/Figure.h"

using namespace std;
using namespace theDoctor;

struct TFnames
{
    string CR;
    string SR;
    string data;
    string MC;
};

class TFFactory
{
public:

TFnames TFmem;
string TFRegion;

private:

Table tab;
Figure dataCR;
Figure MCSR;
Figure MCCR;
//TODO
double additionalError;

public:

TFFactory(string inputTable, string region = "standard", double error = 0)
{
    tab = Table(inputTable);
    dataCR = Figure(0,0);
    MCSR = Figure(0,0);
    MCCR = Figure(0,0);
    additionalError = error;
    TFmem.CR = "CR2l";
    TFmem.SR = "SR1l";
    TFmem.data = "data";
    TFmem.MC = "totalSM";
    TFRegion = region;
}

void defineFactors(string CRtype)
{
    if(TFRegion == "standard")
        return;
    else
    {
        TFmem.CR = CRtype + TFRegion;
        TFmem.SR = "SR1l" + TFRegion;
    }

}

void readFactors(string CRtype)
{

    defineFactors(CRtype);

    dataCR = tab.Get(TFmem.CR,TFmem.data);
    MCSR = tab.Get(TFmem.SR,TFmem.MC);
    MCCR = tab.Get(TFmem.CR,TFmem.MC);

}

Figure returnTF()
{
    if(MCCR.value() != 0)
        return MCSR/MCCR;
    else 
        return Figure(-1,0);
}

Figure returnEstimate()
{
    if(MCCR.value() != 0)
        return dataCR*(MCSR/MCCR);
    else 
        return Figure(-1,0);
}

Figure returnCRdata()
{
    return dataCR;
}
//return TF
};

class TFProducer
{
public:

vector<string> TFRegions;
string TFBkg;
string controlReg;

public:

TFProducer(vector<string> regions, string bkgType, string CRtype)
{
    TFRegions.clear();
    TFRegions = regions;
    TFBkg = bkgType;
    controlReg = CRtype;
}

void produceTFTable(string inputName, string outputName)
{
    vector<string> colId;
    colId.push_back("CR data");
    colId.push_back("TF");
    colId.push_back(TFBkg);

    vector<string> rowId;
    vector<Figure> valueCRdata;
    vector<Figure> valueTF;
    vector<Figure> valueEst;

    for(uint32_t i=0; i<TFRegions.size(); i++)
    {
        TFFactory fact(inputName, TFRegions.at(i));
        fact.TFmem.MC = TFBkg;
        fact.readFactors(controlReg);
        Figure resCRdata = fact.returnCRdata();
        Figure resTF = fact.returnTF();
        Figure resEst = fact.returnEstimate();

        valueCRdata.push_back(resCRdata);
        valueTF.push_back(resTF);
        valueEst.push_back(resEst);
        rowId.push_back(TFRegions.at(i));
    }

    Table outTab(colId, rowId);

    for(uint32_t j=0; j<rowId.size(); j++)
    {
        outTab.Set(colId.at(0), rowId.at(j), valueCRdata.at(j));
        outTab.Set(colId.at(1), rowId.at(j), valueTF.at(j));
        outTab.Set(colId.at(2), rowId.at(j), valueEst.at(j));
    }

    outTab.Print(outputName + ".txt");
    outTab.PrintLatex(outputName + ".tex");
}

};

#endif //TFFactory_h
