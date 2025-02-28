void pull_plot(string file_name, string title){
  auto file_ = new TFile(file_name.c_str(), "READ");
  auto c_ = new TCanvas("c","c", 500, 500);
  auto tree = (TTree*)file_->Get("tree_fit_sb");
  c_->cd();
  double r, rHiErr, rLoErr;
  int status;
  tree->SetBranchAddress("r", &r);
  tree->SetBranchAddress("rHiErr", &rHiErr);
  tree->SetBranchAddress("rLoErr", &rLoErr);
  tree->SetBranchAddress("fit_status", &status);

  auto h = new TH1F("h", "h", 80, -4, 4);
  //tree->Draw("r/(0.5*(rHiErr+rLoErr))>>h");
  for (int i(0); i<tree->GetEntries(); i++){
    tree->GetEntry(i);
    if (status>0) h->Fill(r/(0.5*(rHiErr+rLoErr)));
  }
  h->Fit("gaus");
  h->SetTitle(title.c_str());
  h->GetXaxis()->SetTitle("Pull");
  gStyle->SetOptFit(1);
  c_->Update();
}


  
