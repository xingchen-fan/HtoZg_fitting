using namespace RooFit;
using namespace std;

double br_left  = 120;
double br_right = 130;


std::vector<string> filename_vector(string prefix, std::vector<string> &filenames){
  std::vector<string> ret_vector = {};
  for(unsigned int idx = 0; idx < filenames.size(); idx++){
    ret_vector.push_back(prefix + filenames[idx]);
  }
  return ret_vector;
}


TH1D* get_hist(string infile, string hist_name){
  TFile *file = new TFile(infile.c_str());
  TCanvas* canvas_hist = static_cast<TCanvas*>(file->Get("canvas"));
  TH1D *hist_return = nullptr;
  //Finds histogram in file
  for(int idx_h=0; idx_h<50000; idx_h++){
    if(hist_return == nullptr){
      hist_return = static_cast<TH1D*>( canvas_hist->FindObject( (hist_name+std::to_string(idx_h)).c_str() ) );
    }

    if(hist_return != nullptr){break;}
  }

  return hist_return;
}

TH1D* get_bkg_hist(string infile){ return get_hist(infile,"bkg_Background_");}
TH1D* get_sig_hist(string infile){ return get_hist(infile,"sig_Signal_");}
TH1D* get_mismatch_hist(string infile){ return get_hist(infile,"bkg_H#rightarrowZ#gamma_");}


RooGenericPdf* fit_cthsq_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, RooFitResult **fit_result){ 
  RooRealVar* alpha = new RooRealVar("alpha","alpha", 0.99, 0., 5.); alpha -> setError(0.1);
  RooRealVar* beta  = new RooRealVar("beta", "beta", 1.02, 0., 5.);  beta  -> setError(0.1);
  RooGenericPdf *one_plus_costhetasq = new RooGenericPdf(name, "1 + cos^2#theta", "@1 + @2*@0*@0", RooArgList(x,*alpha,*beta));

  *fit_result = (one_plus_costhetasq -> fitTo(hist_to_fit, SumW2Error(kTRUE), Save(kTRUE)));
  return one_plus_costhetasq;
}

//This function will plot just one function on the histogram with a residual plot
void fit_and_plot_cth(string infile, string outfile, string func){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins = 72;
  double ctheta_low = 0.0;
  double ctheta_high = 0.9;

  //Debug 1
  cout << "Top of function." << endl;

  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar ctheta("ctheta" ,"ctheta", ctheta_low, ctheta_high, "[Bin Width = 0.01]");

  //Debug 1.25
  cout << "..." << endl;

  //RooDataHist hist_total_rdh("hist_total_rdh", "ctheta_dh", ctheta,Import(*(get_bkg_hist(infile))));
  RooDataHist hist_total_rdh("hist_total_rdh", "ctheta_dh", ctheta, Import(*(get_mismatch_hist(infile))));

  //Debug 1.5
  cout << "Really? Here?" << endl;

  //RooCrystalBall* signal_fit = fit_sig_model(ctheta,hist_signal_rdh,"dscb");
  //RooGaussian *signal_fit = fit_sig_gauss(ctheta, hist_signal_rdh,"dscb");


  //Calls whichever function specified by the input variable "func"
  RooFitResult  *fit_result = nullptr;
  RooGenericPdf* fit_onepcthetasq = fit_cthsq_model(ctheta, hist_total_rdh, func.c_str(), &fit_result);   
  cout << fit_result << endl;
  //Debug 2
  cout << "model returned." << endl;

  //Define different plots for the data/fit plot and the residual plot
  //Both will be placed on the same canvas  
  RooPlot* data_and_fit  = ctheta.frame(ctheta_low, ctheta_high);
  hist_total_rdh.plotOn( data_and_fit,Binning(nbins), RooFit::Name("hist_total_rdh"));
  //splusb_fit -> plotOn( data_and_fit, RooFit::Components("splusb_fit"), LineStyle(kSolid),  LineColor(TColor::GetColor("#0000ff")), RooFit::Name("splusb_fit_fr"));
  //splusb_fit -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#0000ff")), RooFit::Name((func+"_fr").c_str()));
  //splusb_fit -> plotOn( data_and_fit, RooFit::Components("dscb"),       LineStyle(kDashed), LineColor(TColor::GetColor("#ff0000")), RooFit::Name("dscb_fr"));
  fit_onepcthetasq -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#df9f1f")), RooFit::Name(func.c_str()));

  cout << fit_result << endl;

  cout << "everything before the ploton" <<endl;
  fit_onepcthetasq -> plotOn( data_and_fit, RooFit::VisualizeError(*fit_result,1,false), RooFit::FillColor(TColor::GetColor("#df9f1f")), RooFit::Components(func.c_str()));

 //RooFit::Normalization(1.0,RooAbsReal::RelativeExpected),
  data_and_fit -> getAttFill() -> SetFillColorAlpha(TColor::GetColor("#df9f1f"),0.3);
  //Debug 3
  cout << "data_and_fit plot made." << endl;


  //Residual plot created from the data_and_fit frame
  RooHist *residual_rdh = data_and_fit -> residHist("hist_total_rdh",func.c_str()); 

  //Debug 4
  cout << "Residual plot made." << endl;


  TString pdf= outfile + ".pdf";
  TCanvas *canvas_plot     = new TCanvas("canvas_plot", "outplot", 1920, 1440);
  TLegend *legend          = new TLegend(0.6,0.8,0.9,0.9);
  TLine   *resid_plot_line = new TLine(ctheta_low,0,ctheta_high,0);

  //Add only entry to legend, make transparent?
  legend->AddEntry("hist_total_rdh"  , "Bkg. MC");//,"L");
  //legend->AddEntry("splusb_fit_fr"      , "S + B Fit", "L");//,"L");
  legend->AddEntry((func).c_str() , "Bkg. Fit", "L");//,"L");
  //legend->AddEntry("dscb_fr"            , "Sig. Fit", "L");//,"L");
  legend->SetFillStyle(0);
  legend->SetLineColorAlpha(kWhite,0);

  //Divide canvas into two pieces
  canvas_plot -> Divide(1,2,0.0,0.0);

  //Set minimum number of events to 0 (quality of life change)
  data_and_fit -> SetMinimum(0.0);

  //A couple of options that should make plotting look nice
  canvas_plot -> Update();
  gStyle      -> SetOptFit();
  
  cout << "plotting on canvas here" << endl;
  double margin_lr  = 0.1;
  double margin_top = 0.05;
  double margin_bot = 0.1;

  double split_point = 0.35;
  double min_x = 0.05;
  double max_x = 1.00;
  double min_y = 0.15;
  
  
  //plot upper plot in upper 7/10ths of plot
  canvas_plot -> cd(1);
  gPad        -> SetPad(min_x, split_point, max_x, 1.0);
  gPad        -> SetMargin(margin_lr, margin_lr, 0.0, margin_top);

  data_and_fit -> GetXaxis() -> SetLabelSize(0);
  data_and_fit -> SetTitle("");
  data_and_fit -> Draw();
  legend       -> Draw();


  //plot lower plot in lower 3/10ths of plot
  canvas_plot -> cd(2);
  gPad        -> SetPad(min_x, min_y,  1, split_point);
  gPad        -> SetMargin(margin_lr, margin_lr, margin_bot, 0);

  //residual_rdh.Draw();
  //residual_rdh -> SetTitle("");
  //residual_rdh -> GetXaxis() -> SetTitle("m_{ll#gamma}");
  residual_rdh -> SetTitle(";m_{ll#gamma};data - fit");

  data_and_fit -> GetXaxis() -> SetLabelSize(0);
  data_and_fit -> SetTitle("");
  data_and_fit -> Draw();
  legend       -> Draw();


  //plot lower plot in lower 3/10ths of plot
  canvas_plot -> cd(2);
  gPad        -> SetPad(min_x, min_y,  1, split_point);
  gPad        -> SetMargin(margin_lr, margin_lr, margin_bot, 0);

  //residual_rdh.Draw();
  //residual_rdh -> SetTitle("");
  //residual_rdh -> GetXaxis() -> SetTitle("m_{ll#gamma}");
  residual_rdh -> SetTitle(";m_{ll#gamma};data - fit");
  residual_rdh -> GetXaxis() -> SetTitleSize(0.15);
  residual_rdh -> GetXaxis() -> SetLabelSize(0.1);
  residual_rdh -> GetXaxis() -> SetLimits(ctheta_low, ctheta_high);
  residual_rdh -> GetXaxis() -> SetTitleOffset(1);

  //residual_rdh -> GetYaxis() -> SetTitle("data - fit");
  residual_rdh -> GetYaxis() -> SetTitleSize(0.15);
  residual_rdh -> GetYaxis() -> SetLabelSize(0.10);
  residual_rdh -> GetYaxis() -> SetTickLength(0.01);
  residual_rdh -> GetYaxis() -> SetTitleOffset(1);
  residual_rdh -> Draw("AP");

  resid_plot_line -> SetLineStyle(kDashed);
  resid_plot_line -> SetLineColorAlpha(kBlue,0.75);
  resid_plot_line -> Draw("SAME");


  //Debug 4
 //Debug 4
  cout << "File plotted, will save and return." << endl;

  canvas_plot -> Print(pdf);

  //delete fit_result;

  return;
}


std::vector<std::string> plot_list(std::string name_segment_one, std::string name_segment_two, bool all=true){
  //Loops to cover
  std::vector<std::string> flavor   = {"_ll"};//"_ee","_mumu"};
  std::vector<std::string> topology = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};

  if(!all){
    flavor   = {"_ee","_mumu"};
    topology = {"WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  }

  std::vector<std::string> list_o_hists = {};
  for(std::string flav : flavor){
    for(std::string top : topology){
      list_o_hists.push_back(name_segment_one + top + flav + name_segment_two);
      cout << list_o_hists[list_o_hists.size() -1] << endl;
    }
  }
  return list_o_hists;
}


void costheta_1d(){
  //String with paths and plots of interest
  //std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AfterMichaelAugustPush/draw_pico/plots/run3_cat_mlly/";
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AddCategoriesAndKinFit/draw_pico/plots/test_plots_2d/";
  //TestPlots_cat_ttH_had_ll_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt__lumi_nonorm_lin.pdf 

  std::string name_seg_one = "TestPlots_cat_";
  std::string name_seg_two = "_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin";
  //test fitting code with one plot
  std::string test_file = "TestPlots_cat_ttH_had_ll_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin.root";

  //std::string test_file = "cat_ttH_lep_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  //fit_and_plot_cth(file_path + test_file, "output/test_plot_VBF_ee.pdf", "don");
  //plot_one_func(file_path + test_file, "output/test_plot_VBF_ee.pdf", "pow3");
  //plot_one_func(file_path + test_file, "output/test_plot", "exmg");
  //fit_and_plot_cth(file_path + test_file, "output/test_plot", "dont_matter");



  //Create list of all the plots
  std::vector<std::string> list_of_plots = plot_list(name_seg_one, name_seg_two);//,false);
  //for(string plot_name : list_of_plots){ fit_and_plot(file_path + plot_name + ".root", "output/" + plot_name); }
  for(string plot_name : list_of_plots){ fit_and_plot_cth(file_path + plot_name + ".root", "output/costheta_1d_fit_" + plot_name, "dont_matter"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/pow3_fit_" + plot_name, "pow3"); }

  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/modg_fit_" + plot_name, "modg"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/exmg_fit_" + plot_name, "exmg"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/gam_fit_" + plot_name, "gam"); }
  //for(string plot_name : list_of_plots){ plot_one_func(file_path + plot_name + ".root", "output/agg_fit_" + plot_name, "agg"); }

}
