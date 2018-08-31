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

//@MJ@ TODO can be used for syst as well -  just the regions will be the ones with name of syst

struct TFnames
{
    string CR;
    string SR;
    string data;
    string MC;
    string MCtot;
};

class TFFactory
{
public:

TFnames TFmem;
string TFRegion;

public:

Table tabSR;
Table tabCR;
Figure dataCR;
Figure MCSR;
Figure MCCR;
//TODO
double additionalError;

public:

TFFactory(string inputTableSR, string inputTableCR, string region = "standard", double error = 0)
{
    tabSR = Table(inputTableSR);
    tabCR = Table(inputTableCR);
    dataCR = Figure(0,0);
    MCSR = Figure(0,0);
    MCCR = Figure(0,0);
    additionalError = error;
    TFmem.CR = "CR2l";
    TFmem.SR = "SR1l";
    TFmem.data = "data1";
    TFmem.MC = "totalSM";
    TFmem.MCtot = "totalSM";
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

    dataCR = tabCR.Get(TFmem.CR,TFmem.data);
    MCSR = tabSR.Get(TFmem.SR,TFmem.MC);
    MCCR = tabCR.Get(TFmem.CR,TFmem.MCtot);

}

void addFactors(Figure dataCR2, Figure MCSR2, Figure MCCR2)
{

    dataCR += dataCR2;
    MCSR += MCSR2; //
    MCCR += MCCR2;
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

void produceTFTable(string inputTabNameSR, string inputTabNameCR, string outputName)
{
    vector<string> colId;
    colId.push_back("CRdata");
    colId.push_back("TF");
    colId.push_back(TFBkg);

    vector<string> rowId;
    vector<Figure> valueCRdata;
    vector<Figure> valueTF;
    vector<Figure> valueEst;

    for(uint32_t i=0; i<TFRegions.size(); i++)
    {
        TFFactory fact(inputTabNameSR, inputTabNameCR, TFRegions.at(i));
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


void extrapolateTFTable(string tableName, string inputTabNameSR, string inputTabNameCR, vector<string> toExtrapolate)
{

    Table TFTable(tableName+".txt");
    vector<string> colIdEx;
    colIdEx = TFTable.colTags;
    colIdEx.push_back("TFextra");
    colIdEx.push_back("TFtotal");

    uint32_t TFIndex = 0;
    uint32_t TFTotIndex = 0;
    uint32_t bkgTypeIndex = 0;
    uint32_t TFExtIndex = 0;
    for(uint32_t t=0; t< colIdEx.size(); t++)
    {
        if(colIdEx.at(t) == "TF")
            TFIndex=t;
        if(colIdEx.at(t) == "TFtotal")
            TFTotIndex=t;
        if(colIdEx.at(t) == TFBkg)
            bkgTypeIndex=t;
        if(colIdEx.at(t) == "TFextra")
            TFExtIndex=t;
    }

    //zkopirovat tabulku do nove 

    vector<string> rowIdEx = TFTable.rowTags;

    Table outTabExt(colIdEx, rowIdEx);

    for(uint32_t j=0; j<rowIdEx.size(); j++)
    {
        for(uint32_t k=0; k<TFTable.colTags.size(); k++)
        {
            outTabExt.Set(colIdEx.at(k), rowIdEx.at(j), TFTable.Get(k,j) );
        }
    }

    uint32_t indexToSkip = 11111;

    for(uint32_t i=0; i<TFRegions.size(); i++)
    {
        for(uint32_t j=0; j<toExtrapolate.size(); j++)
        {
            if( TFRegions.at(i).find(toExtrapolate.at(j))!=std::string::npos)
            {

		TFFactory factExt1(inputTabNameSR, inputTabNameCR, TFRegions.at(i));
		factExt1.TFmem.MC = TFBkg; //background type
		factExt1.readFactors(controlReg);
                Figure dataCR1NoExt = factExt1.dataCR;

		TFFactory factExt2(inputTabNameSR, inputTabNameCR, TFRegions.at(i+1));
		factExt2.TFmem.MC = TFBkg; //background type
		factExt2.readFactors(controlReg);

                cout << " region1 " << TFRegions.at(i) << " dataCR, MCSR, MCCR  " << factExt1.dataCR.value() << " " << factExt1.MCSR.value() << " " << factExt1.MCCR.value() << endl;
                cout << " region2 " << TFRegions.at(i+1) << " dataCR, MCSR, MCCR  " << factExt2.dataCR.value() << " " << factExt2.MCSR.value() << " " << factExt2.MCCR.value() << endl;

                factExt1.addFactors(factExt2.dataCR, factExt2.MCSR, factExt2.MCCR);
               

		Figure resTFExt = factExt1.returnTF();
		Figure resEstExt = factExt1.returnEstimate();
                Figure extrapolationReg1 =  dataCR1NoExt/factExt1.dataCR; 
                Figure extrapolationReg2 =   factExt2.dataCR/factExt1.dataCR; 

                cout << " region1 updated " << TFRegions.at(i) << " dataCR, MCSR, MCCR  " << factExt1.dataCR.value() << " " << factExt1.MCSR.value() << " " << factExt1.MCCR.value() << " TF " << resTFExt.value() << endl;

                outTabExt.Set(colIdEx.at(0), rowIdEx.at(i), factExt1.dataCR ); //data always first
                outTabExt.Set(colIdEx.at(0), rowIdEx.at(i+1), factExt1.dataCR );
                outTabExt.Set(colIdEx.at(TFIndex), rowIdEx.at(i), resTFExt );
                outTabExt.Set(colIdEx.at(TFIndex), rowIdEx.at(i+1), resTFExt );
                outTabExt.Set(colIdEx.at(bkgTypeIndex), rowIdEx.at(i), resEstExt*extrapolationReg1 );
                outTabExt.Set(colIdEx.at(bkgTypeIndex), rowIdEx.at(i+1), resEstExt*extrapolationReg2 );
                outTabExt.Set(colIdEx.at(TFExtIndex), rowIdEx.at(i), extrapolationReg1);
                outTabExt.Set(colIdEx.at(TFExtIndex), rowIdEx.at(i+1), extrapolationReg2);

                cout << " column " << colIdEx.at(TFIndex) << " row " << rowIdEx.at(i) << " value " << resTFExt.value() << endl;

                indexToSkip = i;               
              

 
                //nahazet to do tabulky pro i a i+1
                // + vysledny frakce a vsechny ostatni radky

            }
        }
        if( !(i==indexToSkip || i-1==indexToSkip) )
        {
            outTabExt.Set(colIdEx.at(TFExtIndex), rowIdEx.at(i), Figure(1.0,0.0)  );

        }


    }

    for(uint32_t i=0; i<TFRegions.size(); i++)
    {
        outTabExt.Set(colIdEx.at(TFTotIndex), rowIdEx.at(i), outTabExt.Get(colIdEx.at(TFIndex),rowIdEx.at(i)) * outTabExt.Get(colIdEx.at(TFExtIndex),rowIdEx.at(i)) );
    }

    outTabExt.Print(tableName + "Ext.txt");
    outTabExt.PrintLatex(tableName + "Ext.tex");
}



};

#endif //TFFactory_h
