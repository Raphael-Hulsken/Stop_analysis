#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <sstream>
#include <thread>
#include <iostream>

using namespace std;

#define DPM_DIR "FlatTree/triggerV1/SingleMuon/Run2016B_PromptReco_v2_MINIAOD/160622_212210/0000/"
#define NUMBER_OF_FILES 499
#define DATASET "SingleElectron"
//#define TREE_NAME "FlatTree/tree"
#define TREE_NAME "babyTuple"


// --------------------------------------------------------
// name of the tree: treename
// selection applied: selection
// name of the file to read: ifilename
// --------------------------------------------------------

void SkimOneFile(string ifilename, string treename, string selection, string ofilename){
	cout<<"@@ START @@ " <<endl;
	//cout<<"taking care of file "<<number<<endl;
	//string dpm_path = "root://sbgse1.in2p3.fr//cms/phedex/store/user/echabert/";
	//stringstream fullfilename;
	//fullfilename<<dpm_path<<DPM_DIR<<"output_"<<number<<".root";
	//fullfilename<<base_ifilename<<number<<".root";
	//cout<<"filename: "<<fullfilename.str()<<endl;
	cout<<"filename: "<<ifilename<<endl;
	TFile* fin = TFile::Open(ifilename.c_str());
	TTree* tree = (TTree*) fin->Get(treename.c_str());
	cout<<"Oringal tree have<<"<<tree->GetEntries()<<endl;
	TTree* selTree = tree->CopyTree(selection.c_str());	
	cout<<"Selected tree have "<<selTree->GetEntries()<<endl;
	//stringstream ofilename;
	//ofilename<<"DATASET"<<"_"<<number<<".root";
	TFile* output = new TFile(ofilename.c_str(),"RECREATE");
	selTree->Write();
	output->Write();
	output->Close();
	cout<<"@@ DONE @@"<<endl;
}

int main(int argc, char** argv){
	if(argc<=4){
		cerr<<"Need to provide in the following order"<<endl;
		cerr<<" 1 - full path and name of the input file "<<endl;
		cerr<<"2 - name of the tree to skim "<<endl;
		cerr<<"3 - name of the output file"<<endl;
		cerr<<"4 - selection to be applied [to not use spaces]"<<endl;
		return -1;
	}
	string ifilename(argv[1]);
	string treename(argv[2]);
	string ofilename(argv[3]);
	string selection;
	for(int i=4;i<argc;i++)
		selection+=string(argv[4]);

	// Needed to be able to run on dpm file
	gSystem->Load("/usr/lib64/libdpm.so");
	SkimOneFile(ifilename, treename, selection, ofilename);

	//std::thread th(SkimOneFile,number);
	//th.join();
	/*
	int nofloops = NUMBER_OF_FILES/nWorkers;
	std::vector<thread> workers;
	for(int i=0;i<NUMBER_OF_FILES;i++){
		workers.push_back(std::thread(SkimOneFile,i));
	}
	std::for_each(workers.begin(),workers.end(),[] (std::thread &t){t.join();});
	*/
	/*
	cout<<"@@@ LAST FILES"<<endl;
	// run on the last files [not //]
	// to be improved
	for(int i=nofloops*nWorkers;i<NUMBER_OF_FILES+1;i++){
		std::thread th(SkimOneFile,i);
		
	}
	*/
}
        /*
	int count150 = 0;
	for(unsigned int i=0;i<tree->GetEntries();i++){
		metBranch->GetEntry(i);
		//cout<<met<<endl;
		if(met>150){
			count50++;
			tree->GetEntry();
			
		}
		if(i%1000 == 0) cout<<count50*1./i<<endl;
	}
	*/



