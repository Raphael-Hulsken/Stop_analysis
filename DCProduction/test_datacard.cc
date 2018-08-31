//USAGE: 3 input parameters:
// 1) table with yield
// 2) file with signal, background names, systematics etc..
// 3) file with signal regions names

#include <exception>
#include "producer/DataCardProducer.h"

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 4)
            throw std::runtime_error("Bad number of arguments!");

        string inputTab = argv[1];
        TString inputFile = argv[2];
        TString inputReg = argv[3];

        std::vector<std::vector<string> >     mydata;
        std::ifstream          myfile(inputFile);

        std::string   myline;
        while(std::getline(myfile, myline))
        {
            std::vector<string>   mylineData;
            std::stringstream  mylineStream(myline);

            string myvalue;

            while(mylineStream >> myvalue)
            {   
                mylineData.push_back(myvalue);
            }

            mydata.push_back(mylineData);
        }

        cout << "#processes: " << mydata.size() << endl; 
        if(mydata.size() > 5)
            throw std::runtime_error("more datasets than expected");
        if(mydata.at(0).size() != 4)
            throw std::runtime_error("more parameters than expected");
       	
        DataCardProducer prod;
	prod.ReadTable(inputTab);
    	
	string datasetName = "data"; //@MJ@ TODO this is good for what?!
	vector<string> names;
	vector<string> table_names;
	vector<float> syst_uncert;
	vector<bool> stat_uncert;
        for(uint32_t i=0; i<mydata.size(); i++)
        {
            names.push_back(mydata.at(i).at(0));
            table_names.push_back(mydata.at(i).at(1));
            syst_uncert.push_back(stof(mydata.at(i).at(2)));
            stat_uncert.push_back(static_cast<bool>(stoi(mydata.at(i).at(3))));
        }

	prod.LoadProcesses(datasetName,  names, table_names, syst_uncert, stat_uncert);

	vector<string> bin_names;

	string line;
        ifstream regfile (inputReg);
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                bin_names.push_back(line);
            }
            regfile.close();
        }
	
        prod.LoadBins(bin_names);
        
        string combineStr = "";
      
        for(uint32_t i = 0; i < bin_names.size(); i++)
        {
            string outStr = table_names.at(0) + string(bin_names.at(i));
            prod.WriteBinCard("cards/" +outStr ,i);

            if(i==0)
                combineStr = "python combineCards.py ";

            combineStr += outStr + ".txt ";
        }

        cout << "output string: " << combineStr << " > " <<  mydata.at(0).at(1) << ".log" << endl; 

}
