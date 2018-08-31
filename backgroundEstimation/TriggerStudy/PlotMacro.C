

int main(){

  TFile* fmuon = new TFile("","READ");
  TFile* felectron = new TFile("","READ");
  TFile* fjetht = new TFile("","READ");

  //color

  TFile* fou = new TFile("eff_plots.root","RECREATE");
  TCanvas* c;
  TGraphAsymmErrors* graph;

  // -----------------
  // Muon efficiency
  // -----------------
  string 
  TCanvas c_muonEff("c_muonEff");
  // -- from muon sample
  c = (TCanvas*) fmuon->Get(cname.c_str());
  graph = (TGraphAsymmErrors*) c->GetPrimitive("eff_graph"); 
  c_muonEff.cd();
  graph->SetLineColor(kBlue);
  graph->Draw("AP");
  // -- from electron sample
  c = (TCanvas*) felectron->Get(cname.c_str());
  graph = (TGraphAsymmErrors*) c->GetPrimitive("eff_graph"); 
  graph->SetLineColor(kGreen);
  c_muonEff.cd();
  graph->Draw("same");
  // -- from JetHT sample
  c = (TCanvas*) fjetht->Get(cname.c_str());
  graph = (TGraphAsymmErrors*) c->GetPrimitive("eff_graph"); 
  graph->SetLineColor(kRef);
  c_muonEff.cd();
  graph->Draw("same");
  
  c_muonEff.Write();
  
  
  // -----------------
  // Electron efficiency
  // -----------------
  

  // -----------------
  // MET or MHT efficiency
  // -----------------
  

}
