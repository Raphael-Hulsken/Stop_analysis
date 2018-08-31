
// ########################################################################
// # Credits goes to Wouter Verkerke for original code                    #
// #                                                                      #
// # 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #501      #
// # Using simultaneous p.d.f.s to describe simultaneous fits to multiple #
// # datasets                                                             #
// #                                                                      #
// # 07/2008                                                              #
// ########################################################################

#include "MTTailCorrectionClosureTests.h"
//#include "MTTailCorrectionProducer.h"

#define VERBOSE 0
#define DO_NORM false

std::pair<double,double>  GetSF(RooFitResult* res, string param)
{

    RooRealVar* par_init = (RooRealVar*) res->floatParsInit().find(param.c_str());
    RooRealVar* par_final = (RooRealVar*) res->floatParsFinal().find(param.c_str());
    double vinit = par_init->getVal();
    double vfinal = par_final->getVal();
    double SF = 1.;
    if(par_init!=0) SF = vfinal/vinit;

    double SFerror = 0.;
    if(vfinal!=0) SFerror = SF*par_final->getError()/vfinal;
    if(VERBOSE>0){
   	 cout << "vfinal, parfinalGetErro, SF = " << vfinal << " ; " << par_final->getError() << " ; " << SF << endl;
   	 cout << "SF error = " << SFerror << endl;
    }
    std::pair<double,double> pSF(SF,SFerror);
    return pSF;
}

//----------  Retrieve pdf (norm distrib) -----------------//

TH1F* GetHisto(TFile* fin, string region, string process, string varname, float& norm, bool do_norm, float input_norm)
{
    string cname = CHANNEL_NAME+string("/")+region+"/"+varname;
    TCanvas* c = (TCanvas*) fin->Get(cname.c_str());
    string hname = "v:"+varname+"|p:"+process+"|r:"+region+string("|c:")+CHANNEL_NAME+string("|t:1DEntries");
    TH1F* h = 0;
    if(VERBOSE>0){
 	 cerr<<"cname :"<<cname<<endl;
   	 cerr<<"histo name: "<<hname<<endl;
   	 cerr<<"pointer: "<<c<<endl;
    } 
    TList* l = c->GetListOfPrimitives();
    TPad* pad = (TPad*) l->At(0);
    THStack* stack = (THStack*) pad->GetPrimitive("");
    h = (TH1F*) stack->GetHists()->FindObject(hname.c_str());
    if(do_norm) h->Scale(input_norm/h->Integral());
    norm = h->Integral();
    return (TH1F*) h->Clone();
}

RooHistPdf* GetRooHistPdf(TFile* fin, string region, string process, string varname, RooRealVar* var, float& norm, TH1F*& histo, bool do_mcstat, bool do_norm, float input_norm){
    TH1F* h  =  GetHisto(fin,region,process,varname, norm, do_norm, input_norm);
    histo = h;
    if(do_mcstat)
    {
        //randomisation of the histo ...
        for(int i=1;i<=h->GetNbinsX();i++)
        {
            h->SetBinContent(i,randomnessGenerator->Gaus(h->GetBinContent(i),h->GetBinError(i)));
        }
    }
    string rdhname = "rdh_"+region+"_"+process;
    RooDataHist *rdh  = new RooDataHist(rdhname.c_str(),rdhname.c_str(),RooArgList(*var),Import(*h));
    string pdfname = "pdf_"+region+"_"+process;
    RooHistPdf *pdf  = new  RooHistPdf(pdfname.c_str(),pdfname.c_str(),RooArgSet(*var),*rdh);
    return pdf;
}
//---------------------------------------------------------//

//----------  Retrieve data histo  -----------------
TH1F* GetData(TFile* fin, string region, string varname)
{
    string cname = CHANNEL_NAME+string("/")+region+"/"+varname;
    TCanvas* c = (TCanvas*) fin->Get(cname.c_str());
    TList* l = c->GetListOfPrimitives();
    TPad* pad = (TPad*) l->At(0);
    string hname = "v:"+varname+"|r:"+region+string("|c:")+CHANNEL_NAME+string("|t:1DSumData");
    TH1F* h = (TH1F*) pad->GetPrimitive(hname.c_str());
    return (TH1F*) h->Clone();
}

TH1F* GetPseudoData(TH1F* h_tt1l, TH1F* h_Wjets, TH1F* h_tt2l){
    TH1F* h = (TH1F*) h_tt1l->Clone("data");
    h->Add(h_Wjets);
    h->Add(h_tt2l);
    //randomisation of the histo ...
    for(int i=1;i<=h->GetNbinsX();i++){
           
	    h->SetBinContent(i,randomnessGenerator->Poisson(h->GetBinContent(i)));
    }
    return h;
}

RooDataHist* GetRooData(TH1F* h_tt1l, TH1F* h_Wjets, TH1F* h_tt2l, RooRealVar* var)
{
    //TH1F* h = GetData(fin, region, varname);
    TH1F* h = GetPseudoData(h_tt1l,h_Wjets, h_tt2l);
    //cout<<"nof events in data = "<<h->Integral()<<endl;
    string dname = "data";
    RooDataHist *datah = new RooDataHist(dname.c_str(),dname.c_str(), RooArgList(*var), Import(*h));
    return datah;
}

struct FitSetup
{
    //--name of the root file --//
    string filename;
    //--- variable used to the fit --//
    string varname;
    float  varMin;
    float  varMax;
    //--- region used to the fit ---//
    string region;
    //--- Xsection uncertainties --//
    float  xs_sysfactor;
    bool   do_xs_tt2l_sys;
    bool   do_xs_rare_sys;
    //--- MC stat uncertainties --//
    bool do_mcstat;
    //-- algorithm used to the fit --//
    string type;
    string algo;
    //-- Init value uncert --//
    bool do_init_uncert;
    float init_1ltop;
    float init_Wjets;


    //-- Normalization --//
    float Ndata;
    float rel_norm_1ltop;
    float rel_norm_Wjets;
    float rel_norm_tt2l;

    void Reset()
    {
        filename=string(INPUT_FOLDER)+"/1DDataMCComparison.root";
        varname=OBSERVABLE_FOR_FIT;
        varMin = 0;
        varMax = 600;
        region = "0btag_MTtail";

        //--- Xsection uncertainties --// (or at least, the one who store them later since these variables don't say anything at all anyway...)
        xs_sysfactor=1.;
        do_xs_tt2l_sys = false;
        do_xs_rare_sys = false;

        //--- MC stat uncertainties --//
        do_mcstat = false;

        //-- algorithm used to the fit --//
        type = "Minuit2";
        algo = "MIGRAD";

        //-- Init value uncert --//
        do_init_uncert = false;
        init_1ltop=1.;
        init_Wjets=1.;
    
    	//-- Normalization --//
    	Ndata = 1000;
	rel_norm_1ltop = 0.2;
	rel_norm_Wjets = 0.7;
	rel_norm_tt2l = 0.1;
    }
};

struct FitResult
{
    string conditions;
    pair<float,float> SF_1ltop;
    pair<float,float> SF_Wjets;
    float norm_1ltop;
    float correlation;
    float edm;

    void Reset()
    {
        norm_1ltop = 0;
        edm = 0;
        SF_1ltop = pair<float,float>(0,0);
        SF_Wjets = pair<float,float>(0,0);
        correlation = 0;
    }
    void Print()
    {
        cout<<"#################"<<endl;
        cout<<"# "<<conditions<<"\t SF_1ltop = "<<SF_1ltop.first<<"+/-"<<SF_1ltop.second<<"\t SF_Wjets = "<<SF_Wjets.first<<"+/-"<<SF_Wjets.second<<endl;
        cout<<"#\t edm = "<<edm<<" correlation = "<<correlation<<endl;
        cout<<"#\t init 1ltop: "<<norm_1ltop<<" \t fitted 1ltop: "<<norm_1ltop*SF_1ltop.first<<endl;
        cout<<"#################"<<endl;
    }
};


FitResult  doFit(const FitSetup& setup, string conditions, string fname=string(""))
{
	
    //cerr<<"DO FIT"<<endl;
    string varname = setup.varname;
    RooRealVar var(varname.c_str(),varname.c_str(),setup.varMin, setup.varMax);

    //string region="0btag_MTtail";
    string region= setup.region;

    //should it be an argument ?
    TFile* fin = 0;
    if(fname=="") fin = TFile::Open(setup.filename.c_str());
    //else fin = TFile::Open(fname.c_str());
    else fin = TFile::Open(fname.c_str());

    //-- normalisation in the MC --//
    float mc_norm_1ltop = 0;
    float mc_norm_tt2l = 0;
    float mc_norm_Wjets = 0;
    //float mc_norm_rare = 0;

    // C r e a t e   m o d e l   f o r  CR1_peak_lowM3b
    // -------------------------------------------------------------
    // Construct pdfs for 1ltop, tt2l, Wjets and rare
    TH1F* histo_1ltop = 0;
    TH1F* histo_tt2l = 0;
    TH1F* histo_Wjets = 0;
    
    RooHistPdf *pdf_1ltop  = GetRooHistPdf(fin,region,PROCESS_NAME_TT_1L,varname,&var,mc_norm_1ltop, histo_1ltop, setup.do_mcstat, DO_NORM, setup.Ndata*setup.rel_norm_1ltop);
    RooHistPdf *pdf_tt2l   = GetRooHistPdf(fin,region,PROCESS_NAME_TT_2L,varname,&var,mc_norm_tt2l, histo_tt2l, setup.do_mcstat, DO_NORM, setup.Ndata*setup.rel_norm_tt2l);
    RooHistPdf *pdf_Wjets  = GetRooHistPdf(fin,region,PROCESS_NAME_WJETS,varname,&var,mc_norm_Wjets, histo_Wjets, setup.do_mcstat, DO_NORM, setup.Ndata*setup.rel_norm_Wjets);
    //RooHistPdf *pdf_rare   = GetRooHistPdf(fin,region,PROCESS_NAME_RARE,varname,&var,mc_norm_rare, setup.do_mcstat);

    //cerr<<"TT_1L: "<<mc_norm_1ltop<<endl;
    //cerr<<"TT_2l: "<<mc_norm_tt2l<<endl;
    //cerr<<"WJets: "<<mc_norm_Wjets<<endl;
    //cerr<<"OTHER: "<<mc_norm_rare<<endl;

    // normalization factors (RooRealVar)
    float val_1ltop = mc_norm_1ltop;
    float val_Wjets = mc_norm_Wjets;
    if(setup.do_init_uncert)
    {
        val_1ltop = setup.init_1ltop*mc_norm_1ltop;
        val_Wjets = setup.init_Wjets*mc_norm_Wjets;
    }
    RooRealVar norm_1ltop("norm_1ltop","norm_1ltop",val_1ltop,0.25*mc_norm_1ltop,10.*mc_norm_1ltop);
    RooRealVar norm_Wjets("norm_Wjets","norm_Wjets",val_Wjets,0.25*mc_norm_Wjets,10.*mc_norm_Wjets);
    RooRealVar norm_tt2l("norm_tt2l","norm_tt2l",mc_norm_tt2l,0.25*mc_norm_tt2l,2*mc_norm_tt2l);
    //RooRealVar norm_rare("norm_rare","norm_rare",mc_norm_rare,0.25*mc_norm_rare,2*mc_norm_rare);
    // possibility to study a systematic on it
    if(setup.do_xs_tt2l_sys) mc_norm_tt2l*=setup.xs_sysfactor;
    //if(setup.do_xs_rare_sys) mc_norm_rare*=setup.xs_sysfactor;
    //RooConstVar norm_rare("norm_rare","norm_rare",mc_norm_rare);

    /*
    RooAddPdf model("model","model",
            RooArgList(*pdf_1ltop,*pdf_tt2l,*pdf_Wjets,*pdf_rare),
            RooArgList(norm_1ltop,norm_tt2l,norm_Wjets,norm_rare)) ;
    */
    RooAddPdf model("model","model",
            RooArgList(*pdf_1ltop,*pdf_tt2l,*pdf_Wjets),
            RooArgList(norm_1ltop,norm_tt2l,norm_Wjets)) ;


    //RooDataHist *data_CR1_peak_lowM3b = GetRooData(fin,region,varname,&var);
    RooDataHist *data_CR1_peak_lowM3b = GetRooData(histo_1ltop,histo_Wjets, histo_tt2l,&var);

    fin->Close();


    //--  Constraints on single top and rare --//
    float RelUncert = 0.2;
    // Construct another Gaussian constraint p.d.f on "rare" bkg
    //RooGaussian constr_rare("constr_rare","constr_rare",norm_rare,RooConst(mc_norm_rare),RooConst(RelUncert*mc_norm_rare)) ;
    // Construct another Gaussian constraint p.d.f on "tt2l" bkg
    RooGaussian constr_tt2l("constr_tt2l","constr_tt2l",norm_tt2l,RooConst(mc_norm_tt2l),RooConst(RelUncert*mc_norm_tt2l)) ;

    // P e r f o r m   t em p l a t e   f i t
    // ---------------------------------------------------

    //Minimizer(type,algo) -- Choose minimization package and algorithm to use. Default is MINUIT/MIGRAD through the RooMinimizer
    //                       interface, but rare can be specified (through RooMinimizer interface). Select OldMinuit to use
    //                       MINUIT through the old RooMinuit interface
    //
    //     Type         Algorithm
    //     ------       ---------
    //     OldMinuit    migrad, simplex, minimize (=migrad+simplex), migradimproved (=migrad+improve)
    //     Minuit       migrad, simplex, minimize (=migrad+simplex), migradimproved (=migrad+improve)
    //     Minuit2      migrad, simplex, minimize, scan
    //     GSLMultiMin  conjugatefr, conjugatepr, bfgs, bfgs2, steepestdescent
    //     GSLSimAn     -


    // ---  Perform simultaneous fit of model to data and model_ctl to data_ctl --//
    //RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,Save());
    //RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,ExternalConstraints(constr_rare),ExternalConstraints(constr_tt2l),PrintLevel(-1),Save(),
    RooFitResult* res = model.fitTo(*data_CR1_peak_lowM3b,ExternalConstraints(constr_tt2l),PrintLevel(-1),Save(),
            Minimizer(setup.type.c_str(),setup.algo.c_str()),Verbose(0));

    //--- Writing the results ---///
    FitResult fitRes;
    fitRes.Reset();
    fitRes.norm_1ltop = mc_norm_1ltop;
    fitRes.SF_1ltop = GetSF(res,"norm_1ltop");
    fitRes.SF_Wjets = GetSF(res,"norm_Wjets");
    fitRes.edm = res->edm();
    fitRes.correlation = res->correlationMatrix()[0][1];
    fitRes.conditions = conditions;

    return fitRes;

}

int main()
{
    randomnessGenerator = new TRandom();

    system((string("mkdir -p ")+OUTPUT_FOLDER).c_str());

    // ###########################
    // # Prepare final SFR table #
    // ###########################

    vector<string> columns = { "SFR-1ltop", "SFR-Wjets" };

    vector<string> listAllSignalRegion = listCutAndCounts;

    vector<string> listRawRegion = listIndividualCuts;

    Table tableSFRToBeUsed(columns,listAllSignalRegion);
    Table tableRawSFR(columns,listRawRegion);
    
    columns.clear();
    columns = {"NormPeak_1ltop", "NormTail_1ltop", "SFR-1ltop", "NormPeak_Wjets", "NormTail_Wjets", "SFR-Wjets"};
    Table tableRawNormValues(columns,listRawRegion);

    // Create observables
    RooRealVar var(OBSERVABLE_FOR_FIT,OBSERVABLE_FOR_FIT,0,600) ;
    string varname(OBSERVABLE_FOR_FIT);

    FitSetup setup;
    FitResult res;
    string conditions;


    // ########################
    // #  _____        _____  #
    // # /  __ \ ___  /  __ \ #
    // # | /  \/( _ ) | /  \/ #
    // # | |    / _ \/\ |     #
    // # | \__/\ (_>  < \__/\ #
    // #  \____/\___/\/\____/ #
    // #                      #
    // ########################

    std::map<string,Figure> SFR_CC_1ltop_map;
    std::map<string,Figure> SFR_CC_Wjets_map;

    //Create histos
    TH1F h_SF_MTpeak_CC_1ltop ("h_SF_MTpeak_CC_1ltop",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SF_MTtail_CC_1ltop ("h_SF_MTtail_CC_1ltop",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SFR_CC_1ltop       ("h_SFR_CC_1ltop",        "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SF_MTpeak_CC_Wjets ("h_SF_MTpeak_CC_Wjets",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SF_MTtail_CC_Wjets ("h_SF_MTtail_CC_Wjets",  "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());
    TH1F h_SFR_CC_Wjets       ("h_SFR_CC_Wjets",        "", listIndividualCuts_MTtail.size(), 0, listIndividualCuts_MTtail.size());

    for(unsigned int i=0;i<listIndividualCuts.size();i++)
    {
        cout<<"%%%%%%%%%%%%%%%%%% "<<listIndividualCuts_MTtail[i]<<endl;

        string label = listIndividualCuts[i];
        Figure SFR_1ltop;
        Figure SFR_Wjets;

        //MT tail
        setup.Reset();
        conditions="sigRegions_CC_tail";
        setup.region=listIndividualCuts_MTtail[i];
        setup.varname=varname;
        setup.varMin=0;
        setup.varMax=600;


	//test on Ndata
	
	/*
	vector<int> Ndata = {100,200,500,1000,5000,10000};
	//vector<int> Ndata = {100,200};//,500,1000,5000,10000};
	*/
	TFile* file = new TFile("closure.root","RECREATE");
	//for(unsigned int c=0;c<Ndata.size();c++){
	/*
	TString dirname;
	dirname+=Ndata[c];
	file->mkdir(dirname.Data());
	*/
	TH1F* h_tt1l_Pull = new TH1F("h_tt1l_Pull", "Pull distrib. SF_1ltop", 30, -3, 3);
	TH1F* h_Wjets_Pull = new TH1F("h_Wjets_Pull", "Pull distrib. SF_1ltop", 30, -3, 3);
	TH1F* h_tt1l_uncert = new TH1F("h_tt1l_uncert", "Uncertainty on  SF_1ltop", 100, 0, 5);
	TH1F* h_Wjets_uncert = new TH1F("h_Wjets_uncert", "Uncertainty on SF_1ltop", 100, 0, 1);
        
	//setup.Ndata = Ndata[c];

	for(int np = 0; np<5000; np++){
		res = doFit(setup,conditions);
        	Figure NormTail_1ltop(res.SF_1ltop.first,res.SF_1ltop.second);
      	        Figure NormTail_Wjets=Figure(res.SF_Wjets.first,res.SF_Wjets.second);
		cout<<NormTail_1ltop.Print()<<endl;
		cout<<NormTail_Wjets.Print()<<endl;
		h_tt1l_Pull->Fill((res.SF_1ltop.first-1)/res.SF_1ltop.second);
		h_Wjets_Pull->Fill((res.SF_Wjets.first-1)/res.SF_Wjets.second);
		h_tt1l_uncert->Fill(res.SF_1ltop.second/res.SF_1ltop.first);
		h_Wjets_uncert->Fill(res.SF_Wjets.second/res.SF_Wjets.first);
		//cout<<<<" "<<NormTail_Wjets<<endl;
		//tableRawNormValues.Set("NormTail_1ltop",listIndividualCuts[i],NormTail_1ltop);
		//tableRawNormValues.Set("NormTail_Wjets",listIndividualCuts[i],NormTail_Wjets);
        /*
	h_SF_MTtail_CC_1ltop.SetBinContent(i+1,res.SF_1ltop.first);
        h_SF_MTtail_CC_1ltop.SetBinError(i+1,res.SF_1ltop.second);
        h_SF_MTtail_CC_1ltop.GetXaxis()->SetBinLabel(i+1,label.c_str());
        h_SF_MTtail_CC_Wjets.SetBinContent(i+1,res.SF_Wjets.first);
        h_SF_MTtail_CC_Wjets.SetBinError(i+1,res.SF_Wjets.second);
        h_SF_MTtail_CC_Wjets.GetXaxis()->SetBinLabel(i+1,label.c_str());
	*/
	}
		//file->cd(dirname.Data());
		file->cd();
		h_tt1l_Pull->Write();
		h_Wjets_Pull->Write();
		h_tt1l_uncert->Write();
		h_Wjets_uncert->Write();
		file->cd();
	//}
    	file->Close();
    }

    //Save plots in roofile
    TFile fCR1_CC((string(OUTPUT_FOLDER)+"/results_CC.root").c_str(),"RECREATE");
    h_SF_MTpeak_CC_1ltop.Write();
    h_SF_MTpeak_CC_Wjets.Write();
    h_SF_MTtail_CC_1ltop.Write();
    h_SF_MTtail_CC_Wjets.Write();
    h_SFR_CC_1ltop.Write();
    h_SFR_CC_Wjets.Write();

    //---------------------------------------------
    // Results for C&C
    //---------------------------------------------
    initCutAndCountCuts();
    vector<string> cuts;

    for(unsigned int r=0;r<listCutAndCounts.size();r++)
    {
        cuts = listCutAndCounts_cuts[listCutAndCounts[r]];
        Figure SFR_CC_1ltop;
        Figure SFR_CC_Wjets;
        for(unsigned i=0;i<cuts.size();i++)
        {
            cout<<cuts[i]<<" "<<SFR_CC_1ltop_map[cuts[i]].Print()<<endl;
            if(i == 0) SFR_CC_1ltop = SFR_CC_1ltop_map[cuts[i]];
            else SFR_CC_1ltop = Figure(SFR_CC_1ltop.value(),
                                       sqrt(pow(SFR_CC_1ltop.error(),2)
                                           +pow(SFR_CC_1ltop.value()
                                                -SFR_CC_1ltop_map[cuts[i].c_str()].value(),2)
                                           ));

            if(i == 0) SFR_CC_Wjets = SFR_CC_Wjets_map[cuts[i]];
            else SFR_CC_Wjets = Figure(SFR_CC_Wjets.value(),
                                       sqrt(pow(SFR_CC_Wjets.error(),2)
                                           +pow(SFR_CC_Wjets.value()
                                               -SFR_CC_Wjets_map[cuts[i].c_str()].value(),2)
                                           ));
        }

        // Add 20% of uncertainty for fit itself (JES, MC stat. ...)
        SFR_CC_1ltop *= Figure(1.0,TEMPLATE_FIT_METHOD_UNCERTAINTY);
        SFR_CC_Wjets *= Figure(1.0,TEMPLATE_FIT_METHOD_UNCERTAINTY);

        tableSFRToBeUsed.Set("SFR-1ltop",listCutAndCounts[r],SFR_CC_1ltop);
        tableSFRToBeUsed.Set("SFR-Wjets",listCutAndCounts[r],SFR_CC_Wjets);
    }

    tableRawSFR.Print(string(OUTPUT_FOLDER)+"/rawSFR.tab"   ,4);
    tableRawNormValues.Print(string(OUTPUT_FOLDER)+"/rawNormValues.tab"   ,4);
    tableSFRToBeUsed.Print(string(OUTPUT_FOLDER)+"/SF_MTtail.tab",4);
    
    tableRawSFR.PrintLatex(string(OUTPUT_FOLDER)+"/SF_MTtail.tex",4);
}
