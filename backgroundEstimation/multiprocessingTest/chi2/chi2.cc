#include "chi2.h"

bool debug = false;


//--------------------------------------------------------------------
double fc2 (double c1, double m12, double m22, double m02, bool verbose)
{
  if (verbose) {
    printf("c1: %4.2f\n", c1);
    printf("m12: %4.2f\n", m12);
    printf("m22: %4.2f\n", m22);
    printf("m02: %4.2f\n", m02);
  }

  double a = m22;
  double b = (m02 - m12 - m22) * c1;
  double c = m12 * c1 * c1 - PDG_W_MASS * PDG_W_MASS;

  if (verbose) {
    printf("a: %4.2f\n", a);
    printf("b: %4.2f\n", b);
    printf("c: %4.2f\n", c);
  }

  double num = -1. * b + sqrt(b * b - 4 * a * c);
  double den = 2 * a;

  if (verbose) {
    printf("num: %4.2f\n", num);
    printf("den: %4.2f\n", den);
    printf("num/den: %4.2f\n", num/den);
  }

  return (num/den);
}

//--------------------------------------------------------------------
double fchi2 (void *num, void* denom, uint32_t size){

  double numerator[size];
  double denominator[size];
  //@MJ@ TODO fix this, this can be dangerous with other types than double
  memcpy(numerator, num, size*sizeof(double));
  memcpy(denominator, denom, size*sizeof(double));
  vector<double> members;
  
  for(uint32_t i = 0; i < size; i++)
  {
       //std::cout << "fchi num "  <<numerator[i] << " denom: " << denominator[i] << std::endl;
       double rat = numerator[i]/denominator[i];
       members.push_back(rat*rat);
  }

  double result = 0;

  for(std::vector<double>::iterator it = members.begin(); it != members.end(); ++it) 
  {
    result = result + (*it);
  }
  return result;
}


//--------------------------------------------------------------------
/*void minuitFunction(int&, double* , double &result, double par[], int){
  result=fchi2(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7]);
}
*/


vector<double> calculateChi2ForWMass(vector<TLorentzVector>& jets, vector<float>& sigma_jets, int njets_max = 6, bool new_formula = false)
{

  if(jets.size() != sigma_jets.size())
      throw std::runtime_error("lenght of jets vector is not equal to lenght of its sigmas");

	//check at most first 6 jets
	int n_jets = jets.size();
	if (n_jets>njets_max) n_jets = njets_max;
	//consider at least 3 jets
	if (n_jets<3)
        {
            std::vector<double> nothing;
            return nothing;
        }

	vector<double> v_ij;
	for ( int i=0; i<n_jets; ++i )
        {
		for ( int j=i+1; j<n_jets; ++j )
                {

			//
			//  W
			//
			TLorentzVector hadW = jets[i] + jets[j];
			double pt_w1 = jets[i].Pt();
			double pt_w2 = jets[j].Pt();
                       
                        double sigma_w = 0;
			if(!new_formula)
				sigma_w = sqrt(pow(pt_w1*sigma_jets[i], 2)
				+ pow(pt_w2*sigma_jets[j], 2));
			else		
				//the formula above was wrong
				sigma_w = sqrt( pow(hadW.M()/(2*pt_w1)*sigma_jets[i],2) + pow(hadW.M()/(2*pt_w2)*sigma_jets[j],2)); 

                        double diff = hadW.M() - PDG_W_MASS;

		        //std::cout << "w mass from jets: " << hadW.M() << "PDG mass: " << PDG_W_MASS << " sigma: " << sigma_w << std::endl;
			v_ij.push_back(diff/sigma_w);
                        if(debug) std::cout << "calculateChi2ForWMass chi: " << diff/sigma_w << " = " << diff << " / " << sigma_w << std::endl; 
		}
        }
        return v_ij;
}


vector<double> calculateChi2ForWdR(vector<TLorentzVector>& jets, vector<float>& sigma_jets, int njets_max = 6)
{

  if(jets.size() != sigma_jets.size())
      throw std::runtime_error("lenght of jets vector is not equal to lenght of its sigmas");

	//check at most first 6 jets
	int n_jets = jets.size();
	if (n_jets>njets_max) n_jets = njets_max;
	//consider at least 3 jets
	if (n_jets<3)
        {
            std::vector<double> nothing;
            return nothing;
        }

	vector<double> v_ij;
	for ( int i=0; i<n_jets; ++i )
        {
		for ( int j=i+1; j<n_jets; ++j )
                {

			//
			//  W
			//
			TLorentzVector hadW = jets[i] + jets[j];
                        double hadWPt = hadW.Pt();
                        double hadWdR = jets[i].DeltaR(jets[j]);
                        
			double pt_w1 = jets[i].Pt();
			double pt_w2 = jets[j].Pt();
                        
			/*
                        double sigma_w = sqrt(pow(pt_w1*sigma_jets[i], 2)
				+ pow(pt_w2*sigma_jets[j], 2));
			*/
			// the formula above was wrong

			float sigma_w = hadWdR* sqrt( pow(pt_w1/hadWPt*sigma_jets[i],2) + pow(pt_w2/hadWPt*sigma_jets[j],2) );

                        double diff = hadWdR*hadWPt - (2*PDG_W_MASS);
			
			v_ij.push_back(diff/sigma_w);
			if(debug) std::cout << "calculateChi2FordR_W chi: " << diff/sigma_w << " = " << diff << " / " << sigma_w << std::endl; 
		}
        }
        return v_ij;
}


vector< vector<double> > calculateChi2ForTopMass(vector<TLorentzVector>& jets, vector<float>& sigma_jets, int njets_max = 6, bool new_formula = false)
{

  if(jets.size() != sigma_jets.size())
      throw std::runtime_error("lenght of jets vector is not equal to lenght of its sigmas");

	//check at most first 6 jets
	int n_jets = jets.size();
	if (n_jets>njets_max) n_jets = njets_max;
	//consider at least 3 jets
	if (n_jets<3)
        {
            vector< vector<double> > nothing;
            return nothing;
        }

	vector< vector<double> > v_ij;
	for ( int i=0; i<n_jets; ++i ) //@MJ@ TODO not sure about this!
        {
		for ( int j=i+1; j<n_jets; ++j )
                {
	            vector<double> v_ijk;
		    for ( int k=0; k<n_jets; ++k )
                    {    
                        if(k==i || k==j)
                            continue;

			//
			//  W
			//
			TLorentzVector hadW = jets[i] + jets[j];
			TLorentzVector hadTop = hadW + jets[k];
			double pt_w1 = jets[i].Pt();
			double pt_w2 = jets[j].Pt();
			double pt_w3 = jets[k].Pt();
                        
                        double sigma_top = 0;
			
			if(!new_formula){
			   sigma_top = sqrt(pow(pt_w1*sigma_jets[i], 2)
				+ pow(pt_w2*sigma_jets[j], 2)
				+ pow(pt_w3*sigma_jets[k], 2));
			}
			else{ // the formula above was wrong
			  sigma_top = pow(sigma_jets[i]*(jets[j].Pt()*(1-cos(jets[i].Theta()-jets[j].Theta()))/(sin(jets[i].Theta())*sin(jets[j].Theta()))
				+jets[k].Pt()*(1-cos(jets[i].Theta()-jets[k].Theta())))/(sin(jets[i].Theta())*sin(jets[k].Theta())),2);
			  sigma_top += pow(sigma_jets[j]*(jets[i].Pt()*(1-cos(jets[i].Theta()-jets[j].Theta()))/(sin(jets[i].Theta())*sin(jets[j].Theta()))
				+jets[k].Pt()*(1-cos(jets[j].Theta()-jets[k].Theta())))/(sin(jets[j].Theta())*sin(jets[k].Theta())),2);
			  sigma_top += pow(sigma_jets[k]*(jets[i].Pt()*(1-cos(jets[i].Theta()-jets[k].Theta()))/(sin(jets[i].Theta())*sin(jets[k].Theta()))
				+jets[j].Pt()*(1-cos(jets[j].Theta()-jets[k].Theta())))/(sin(jets[j].Theta())*sin(jets[k].Theta())),2);
			  sigma_top = sqrt(sigma_top/pow(hadTop.M(),2));
			}
			double diff = hadTop.M() - PDG_TOP_MASS;
			
			v_ijk.push_back(diff/sigma_top);
			if(debug) std::cout << "calculateChi2ForMass_top chi: " << diff/sigma_top << " = " << diff << " / " << sigma_top << std::endl; 
                    }
                    v_ij.push_back(v_ijk);
		}
        }
        return v_ij;
}

vector< vector<double> > calculateChi2ForTopdR(vector<TLorentzVector>& jets, vector<float>& sigma_jets, int njets_max = 6)
{

  if(jets.size() != sigma_jets.size())
      throw std::runtime_error("lenght of jets vector is not equal to lenght of its sigmas");

	//check at most first 6 jets
	int n_jets = jets.size();
	if (n_jets>njets_max) n_jets = njets_max;
	//consider at least 3 jets
	if (n_jets<3)
        {
            vector< vector<double> > nothing;
            return nothing;
        }

	vector< vector<double> > v_ij;
	for ( int i=0; i<n_jets; ++i ) //@MJ@ TODO not sure about this!
        {
		for ( int j=i+1; j<n_jets; ++j )
                {
	            vector<double> v_ijk;
		    for ( int k=0; k<n_jets; ++k )
                    {    
                        if(k==i || k==j)
                            continue;

			//
			//  W
			//
			TLorentzVector hadW = jets[i] + jets[j];
			TLorentzVector hadTop = hadW + jets[k];
                        double hadTopPt = hadTop.Pt();
			double hadTopdR = hadW.DeltaR(jets[k]);
			double pt_w1 = jets[i].Pt();
			double pt_w2 = jets[j].Pt();
			double pt_w3 = jets[k].Pt();
                        
                        /*
			double sigma_top = sqrt(pow(pt_w1*sigma_jets[i], 2)
				+ pow(pt_w2*sigma_jets[j], 2)
				+ pow(pt_w3*sigma_jets[k], 2));
			*/
			//the formula above was wrong
			double sigma_top = pow(sigma_jets[i]*(jets[i].Pt()+jets[j].Pt()*(1-cos(jets[i].Theta()-jets[j].Theta()))+jets[k].Pt()*(1-cos(jets[i].Theta()-jets[k].Theta()))),2);
			sigma_top += pow(sigma_jets[j]*(jets[j].Pt()+jets[i].Pt()*(1-cos(jets[i].Theta()-jets[j].Theta()))+jets[k].Pt()*(1-cos(jets[j].Theta()-jets[k].Theta()))),2);
			sigma_top += pow(sigma_jets[k]*(jets[k].Pt()+jets[i].Pt()*(1-cos(jets[k].Theta()-jets[i].Theta()))+jets[j].Pt()*(1-cos(jets[j].Theta()-jets[k].Theta()))),2);

			sigma_top = hadTopdR*sqrt(sigma_top/pow(hadTopPt,2));

                        double diff = hadTopdR*hadTopPt - (2*PDG_TOP_MASS)/hadTopPt;
                        
			if(debug) std::cout << "calculateChi2FordR_top chi: " << diff/sigma_top << " = " << diff << " / " << sigma_top << std::endl; 
			
			v_ijk.push_back(diff/sigma_top);
                        }
                    v_ij.push_back(v_ijk);
		}
        }
        return v_ij;
}

void minuitFunction(int&, double* , double &result, double par[], int){
}

// This function calculates the hadronic chi2 - SNT version
double calculateChi2(vector<LorentzVector>& jets, vector<float>& sigma_jets, vector<bool>& btag){

	assert(jets.size() == sigma_jets.size());
	assert(jets.size() == btag.size());

	//check at most first 6 jets
	int n_jets = jets.size();
	if (n_jets>6) n_jets = 6;
	//consider at least 3 jets
	if (n_jets<3) return 999999.;

	vector<int> v_i, v_j;
	vector<double> v_k1, v_k2;
	for ( int i=0; i<n_jets; ++i )
		for ( int j=i+1; j<n_jets; ++j ){

			//
			//  W
			//
			LorentzVector hadW = jets[i] + jets[j];

			//
			//  W Mass Constraint.
			//
			TFitter *minimizer = new TFitter();
			double p1 = -1;

			minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);
			minimizer->SetFCN(minuitFunction);
			minimizer->SetParameter(0 , "c1"     , 1.1             , 1 , 0 , 0);
			minimizer->SetParameter(1 , "pt1"    , 1.0             , 1 , 0 , 0);
			minimizer->SetParameter(2 , "sigma1" , sigma_jets[i]   , 1 , 0 , 0);
			minimizer->SetParameter(3 , "pt2"    , 1.0             , 1 , 0 , 0);
			minimizer->SetParameter(4 , "sigma2" , sigma_jets[j]   , 1 , 0 , 0);
			minimizer->SetParameter(5 , "m12"    , jets[i].mass2() , 1 , 0 , 0);
			minimizer->SetParameter(6 , "m22"    , jets[j].mass2() , 1 , 0 , 0);
			minimizer->SetParameter(7 , "m02"    , hadW.mass2()    , 1 , 0 , 0);

			for (unsigned int k = 1; k < 8; k++)
				minimizer->FixParameter(k);

			minimizer->ExecuteCommand("SIMPLEX", 0, 0);
			minimizer->ExecuteCommand("MIGRAD", 0, 0);

			double c1 = minimizer->GetParameter(0);
			if (c1!=c1) {
				//cout<<"[PartonCombinatorics::recoHadronicTop] ERROR: c1 parameter is NAN! Skipping this parton combination"
					//<<endl;
				continue;
			}
			double c2 = fc2(c1, jets[i].mass2(), jets[j].mass2(), hadW.mass2());

			delete minimizer;


			//     * W Mass check :)
			//     *  Never trust a computer you can't throw out a window.
			//      *  - Steve Wozniak

			// cout << "c1 = " <<  c1 << "  c1 = " << c2 << "   M_jj = "
			// 	   << ((jets[i] * c1) + (jets[j] * c2)).mass() << endl;

			v_i.push_back(i);
			v_j.push_back(j);
			v_k1.push_back(c1);
			v_k2.push_back(c2);
		}

	//Apply b-consistency requirement
	int n_btag = 0;
	for( int i = 0 ; i < n_jets ; i++ )
		if( btag.at(i) ) n_btag++;

	double chi2min = 99999.;

	//consider b-jet in leading 3 jets
	for ( int b=0; b<n_jets; ++b ) {    

		//if not tagged, consider only 3 leading jets
		if( (!btag.at(b)) && b>2 ) continue;

		//require b-tagging if have more than 1 b-tag
		if( n_btag>1 && (!btag.at(b)) ) continue;
		double pt_b = jets[b].Pt();

		for (unsigned int w = 0; w < v_i.size() ; ++w ) {
			int i = v_i[w];
			int j = v_j[w];
			if ( i==b || j==b ) continue;
			//count number of b-tagged Ws
			int nwb = 0;
			if (btag.at(i)) nwb++;
			if (btag.at(j)) nwb++;
			//no btagged jets in W if have few btags
			if ( n_btag<3  && nwb>0 ) continue;
			//In 3 b-tag case, allow for 1 W jet to be tagged
			// If have more b-tags then btagging information not useful
			if ( n_btag==3 && nwb>1 ) continue;

			double pt_w1 = jets[i].Pt();
			double pt_w2 = jets[j].Pt();

			///
			//  W Mass.
			///
			LorentzVector hadW = jets[i] + jets[j];
			double massW = hadW.mass();

			double c1 = v_k1[w];
			double c2 = v_k2[w];

			///
			// Top Mass.
			///
			LorentzVector hadT = (jets[i] * c1) + (jets[j] * c2) + jets[b];
			double massT = hadT.mass();

			double pt_w = hadW.Pt();
			double sigma_w2 = pow(pt_w1*sigma_jets[i], 2)
				+ pow(pt_w2*sigma_jets[j], 2);
			double smw2 = (1. + 2.*pow(pt_w,2)/pow(massW,2))*sigma_w2;
			double pt_t = hadT.Pt();
			double sigma_t2 = pow(c1*pt_w1*sigma_jets[i],2)
				+ pow(c2*pt_w2*sigma_jets[j],2)
				+ pow(pt_b*sigma_jets[b],2);
			double smtop2 = (1. + 2.*pow(pt_t,2)/pow(massT,2))*sigma_t2;

			double c_chi2 = pow(massT-PDG_TOP_MASS, 2)/smtop2
				+ pow(massW-PDG_W_MASS, 2)/smw2;
			if (c_chi2<chi2min) chi2min = c_chi2;

		}
	}

	return chi2min;
}

// This function is interfacing python code and  hadronic chi2 - SNT version
//double Chi2Interface(int njets, double* jets_pt, double* jets_eta, double* jets_phi, double* jets_E, double* sigma, bool* jets_btag){
/*extern "C" double Chi2Interface(int njets, double* jets_pt, double* jets_eta, double* jets_phi, double* jets_E, bool* jets_btag){
	vector<LorentzVector> jets;
	vector<float> sigmas;
	vector<bool> btag;
	for (int i=0;i< njets;i++){
		TLorentzVector j;
		j.SetPtEtaPhiE(jets_pt[i], jets_eta[i], jets_phi[i], jets_E[i]);
		LorentzVector lv(j.Px(), j.Py(), j.Pz(), j.E());
		//fill jets
		jets.push_back(lv);
		//fill sigmas
		sigmas.push_back(0.1);
		//fill btag
		btag.push_back(jets_btag[i]);
	}
	return calculateChi2(jets, sigmas, btag);
}
*/

// method:
// 1: Wmass+top-mass
// 2: Wmass + DRW
// 3: DRW + DRTop
// 4: Wmass + top-mass + DRW
// 5: Wmass + top-mass + DRW + DRTop

double Chi2(vector<float> jets_pt, vector<float> jets_eta, vector<float> jets_phi, vector<float> jets_m, int method = 1, int njets_max = 6, bool new_formula = false)
{
	vector<TLorentzVector> jets;
	vector<float> sigmas;
	vector<bool> btag;
	for (uint32_t i=0;i< jets_pt.size();i++){
		TLorentzVector j;
		j.SetPtEtaPhiM(static_cast<Double_t>(jets_pt.at(i)), static_cast<Double_t>(jets_eta.at(i)), static_cast<Double_t>(jets_phi.at(i)), static_cast<Double_t>(jets_m.at(i)));
		//fill jets
		jets.push_back(j);
		//fill sigmas
		sigmas.push_back(0.1 * jets_pt.at(i));
	}

        vector<double> m1;
        vector<double> m2;
        vector< vector<double> > m3;
        vector< vector<double> > m4;
	m1 =  calculateChi2ForWMass(jets, sigmas, njets_max, new_formula);
	m2 =  calculateChi2ForWdR(jets, sigmas, njets_max);
	m3 =  calculateChi2ForTopMass(jets, sigmas, njets_max, new_formula);
	m4 =  calculateChi2ForTopdR(jets, sigmas, njets_max);

        if(m1.size() != m2.size() || m3.size() != m4.size() || m1.size() != m3.size())
            throw std::runtime_error("1: sizes of chi 2 memebers are different");

        vector<double> chi2s;
        for(uint32_t ch = 0; ch < m1.size(); ch++)
        {
	    ///*
            vector<double> mm3 = m3.at(ch);
            vector<double> mm4 = m4.at(ch);
            if(mm3.size() != mm4.size())
                throw std::runtime_error("2: sizes of chi 2 memebers are different");
            for(uint32_t h = 0; h < mm3.size(); h++)
            {
                /*
		double arr[] = {m1.at(ch), m2.at(ch), mm3.at(h), mm4.at(h)};
                //std::cout << "m1.at(ch) " << m1.at(ch) << " m2.at(ch) " << m2.at(ch) << "  mm3.at(h) " << mm3.at(h) << " mm4.at(h) " << mm4.at(h)<< std::endl;
                double arrSig[] = {1,1,1,1};
                chi2s.push_back(fchi2(arr, arrSig, sizeof(arr)/sizeof(double)));
	        */
		switch(method){
	      	  case 1:{
			double arr1[] = {m1.at(ch), mm3.at(h)};
                	double arrSig1[] = {1,1};
                	chi2s.push_back(fchi2(arr1, arrSig1, sizeof(arr1)/sizeof(double)));
			}break;
		  case 2:{  
			double arr2[] = {m1.at(ch), m2.at(ch)};
                	double arrSig2[] = {1,1};
                	chi2s.push_back(fchi2(arr2, arrSig2, sizeof(arr2)/sizeof(double)));
			}break;
		  case 3:{
			double arr3[] = {m2.at(ch), mm4.at(h)};
                	double arrSig3[] = {1,1};
                	chi2s.push_back(fchi2(arr3, arrSig3, sizeof(arr3)/sizeof(double)));
			}break;
		  case 4:{
			double arr4[] = {m1.at(ch), m2.at(ch), mm3.at(h)};
                	double arrSig4[] = {1,1,1};
                	chi2s.push_back(fchi2(arr4, arrSig4, sizeof(arr4)/sizeof(double)));
			}break;
		  case 5:{
			double arr5[] = {m1.at(ch), m2.at(ch), mm3.at(h),mm4.at(h)};
                	double arrSig5[] = {1,1,1,1};
                	chi2s.push_back(fchi2(arr5, arrSig5, sizeof(arr5)/sizeof(double)));
			}break;
	       }
            }
	    //*/
        }
        
    
        double res = chi2s.size() !=0 ? *min_element(chi2s.begin(), chi2s.end()) : -1; 
        return res;
}
