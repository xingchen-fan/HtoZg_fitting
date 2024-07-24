#include "../APFitting/APF_utilities.cxx"
//#include "../Utilities/HZGRooPdfs.cxx"
#include "../Utilities/RooGaussStepBernstein.h"

using namespace APFUtilities;

double br_left  = 120;
double br_right = 130;

TH1D* get_mismatch_hist(string infile){ return get_hist(infile,"bkg_H#rightarrowZ#gamma_");}
TH1D* get_fake_hist(string infile){ return get_hist(infile,"bkg_Z+FakePhoton_");}

TH2D* get_mismatch_hist_2D(string infile){ return get_hist_2D(infile,"bkg_H#rightarrowZ#gamma_");}
TH2D* get_noname_hist_2D(string infile){ 
  TFile *file = new TFile(infile.c_str());
  TCanvas* canvas_hist = static_cast<TCanvas*>(file->Get("canvas"));
  TH2D *hist_return = static_cast<TH2D*>( canvas_hist->FindObject("") );
  return hist_return;
}

  /*
  //Bernstein functions
  //Bern 5
  RooRealVar p0("p0", "p0", 15.);
  RooRealVar b5p1("b5p1", "b5p1", 0.3, -200., 200.);
  RooRealVar b5p2("b5p2", "b5p2", 0.3,-200.,200.);
  RooRealVar b5p3("b5p3", "b5p3", 0.3,-200.,200.);
  RooRealVar b5p4("b5p4", "b5p4", 0.3,-200.,200.);
  RooRealVar b5p5("b5p5", "b5p5", 0.3,-200.,200.);

  RooRealVar sigma_bern5("sigma_bern5","sigma_bern5"       ,4.7647,  1., 20.);
  RooRealVar stepval_bern5("stepval_bern5", "stepval_bern5", 103.46, 100., 120.);
  RooGaussStepBernstein bern5_model("bern5_model", "Bernstein 5(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2,b5p3, b5p4, b5p5));
  RooFitResult* ber5_fit = ber5_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  //Bern 2
  RooGaussStepBernstein bern2_model("bern2_model", "Bernstein 2(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2));
  RooFitResult* ber2_fit = ber2_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  //Bern 3
  RooGaussStepBernstein bern3_model("bern3_model", "Bernstein 3(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2, b5p3));
  RooFitResult* ber3_fit = ber3_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

  // Bern 4 
  RooGaussStepBernstein bern4_model("bern4_model", "Bernstein 4(X) gauss", m_lly, meanx, sigma_bern5, stepval_bern5, RooArgList(p0,b5p1,b5p2,b5p3, b5p4));
  RooFitResult* ber4_fit = ber4_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));
  //End of Bernstein functions
  */

/*
RooGaussStepBernstein *fit_bern4_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, RooFitResult *fit_result=nullptr){
  RooRealVar *mux  = new RooRealVar("mux", "mean of gaussian x", 0);

  RooRealVar* p0   = new RooRealVar("p0", "p0", 15.);
  RooRealVar* b5p1 = new RooRealVar("b5p1", "b5p1", 0.3, -200., 200.);
  RooRealVar* b5p2 = new RooRealVar("b5p2", "b5p2", 0.3,-200.,200.);
  RooRealVar* b5p3 = new RooRealVar("b5p3", "b5p3", 0.3,-200.,200.);
  RooRealVar* b5p4 = new RooRealVar("b5p4", "b5p4", 0.3,-200.,200.);

  RooRealVar* sigma_bern5   = new RooRealVar("sigma_bern5",   "sigma_bern5",   5.0,  0.1, 15.0);
  RooRealVar* stepval_bern5 = new RooRealVar("stepval_bern5", "stepval_bern5", 108.0,  95.0, 115.0);

  RooGaussStepBernstein *bern4_model = new RooGaussStepBernstein("bern4_model", "Bernstein 4(X) gauss", x, *mux, *sigma_bern5, *stepval_bern5, RooArgList(*p0, *b5p1, *b5p2, *b5p3, *b5p4));
  fit_result = bern4_model -> fitTo(hist_to_fit,SumW2Error(kTRUE) ,Save(kTRUE));

  return bern4_model;

  //RooFitResult* ber4_fit = ber4_model.fitTo(hist_total_rdh,SumW2Error(kTRUE) ,Save(kTRUE));

}
*/


RooProdPdf* fit_2D_signal(RooRealVar &x, RooRealVar &y, RooDataHist &hist_to_fit, RooDataHist &hist_x, RooDataHist &hist_y, const char *name, RooFitResult *result=nullptr){
  cout << "------------------------------------ S FIT --------------------------------" << endl;
  //RooFitResult fit_result_x = nullptr; RooFitResult fit_result_y = nullptr; RooFitResult fit_result_2d = nullptr;
  bool fix = true;

  //RooGaussBernstein *fit_x_axis = fit_bern4_model(x,hist_x,"xdim_fit",result);
  RooCrystalBall *sig_fit_x_axis = APFUtilities::fit_sig_model(x,hist_x,"sig_fit_x");//,fit_result_x);
  RooGenericPdf  *sig_fit_y_axis = fit_cthsq_model(y, hist_y, "sig_fit_y", fix); 
  RooProdPdf *sig_fit_twodim = new RooProdPdf("sig_fit_2D","sig_fit_2D", RooArgList(*sig_fit_x_axis,*sig_fit_y_axis));
  //result = sig_fit_twodim -> fitTo(hist_to_fit, SumW2Error(kTRUE) ,Save(kTRUE));
 
  return sig_fit_twodim;
}

RooProdPdf* fit_2D_bkg(RooRealVar &x, RooRealVar &y, RooDataHist &hist_to_fit, RooDataHist &hist_x, RooDataHist &hist_y, const char *name, RooFitResult *result=nullptr){
  cout << "------------------------------------B FIT --------------------------------" << endl;
  //RooFitResult *fit_result_x = nullptr; RooFitResult *fit_result_y = nullptr; RooFitResult *fit_result_2d = nullptr;
  bool fix = false;

  //RooGaussStepBernstein *bkg_fit_x_axis = fit_bern4_model(x,hist_x,name);//,result);
  //RooFFTConvPdf *bkg_fit_x_axis = fit_pow3_model(x, hist_to_fit, name);
  EXModGaus     *bkg_fit_x_axis = fit_exmg_model(x, hist_to_fit, name);
  RooGenericPdf *bkg_fit_y_axis = fit_cthsq_model(y, hist_y, "bkg_fit_y", false); 
  RooProdPdf *sig_fit_twodim = new RooProdPdf("bkg_fit_2D","bkg_fit_2D", RooArgList(*bkg_fit_x_axis,*bkg_fit_y_axis));
  //result = sig_fit_twodim -> fitTo(hist_to_fit, SumW2Error(kTRUE) ,Save(kTRUE));
 
  return sig_fit_twodim;
}


RooAddPdf* fit_2D_splusb(RooProdPdf *sig_2d_fit, RooProdPdf *bkg_2d_fit, RooDataHist &hist_to_fit, const char *name, double nbkg_init = 0, double nsig_init=0, RooFitResult *result=nullptr){
  cout << "------------------------------------ S + B FIT --------------------------------" << endl;
  RooRealVar *Nbkg = new RooRealVar("nbkg", "Nbkg", nbkg_init, 0, 100000000000);
  RooRealVar *Nsig = new RooRealVar("nsig", "Nsig", nsig_init, -100000, 100000);
  cout << Nbkg -> getValV() << endl;
  cout << Nsig -> getValV() << endl;

  cout << "fit def" << endl;
  RooAddPdf *spb_fit_twodim = new RooAddPdf("spb_fit_twodim", "spb_fit_twodim", RooArgList(*bkg_2d_fit,*sig_2d_fit), RooArgList(*Nbkg,*Nsig));

  cout << "fitting func" << endl;
  result = spb_fit_twodim -> fitTo(hist_to_fit, SumW2Error(kTRUE) ,Save(kTRUE));
 
  return spb_fit_twodim;
}

//This function will plot just one function on the histogram with a residual plot
void fit_and_plot_mllycth(vector<string> infile, string outfile, string func){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins_ctheta = 72;
  int nbins_mlly = 40;

  double ctheta_low = 0.0;
  double ctheta_high = 0.9;
  double m_lly_low = 100;
  double m_lly_high = 180;


  //Debug 1
  cout << "Top of function." << endl;

  //Creates fit variables and converts histogram to RooDataHist
  RooRealVar ctheta("ctheta" ,"ctheta", ctheta_low, ctheta_high, "[Bin Width = 0.01]");
  RooRealVar m_lly("m_lly" ,"m_lly", m_lly_low, m_lly_high, "GeV");

  //Debug 1.25
  cout << "..." << endl;

  //RooDataHist hist_total_rdh("hist_total_rdh", "ctheta_dh", ctheta,Import(*(get_bkg_hist(infile))));
  //infile idces: 0 - mlly_costheta sig, 1 - mlly_costheta bkg, 2 - mlly sig, 3 - costheta sig, 4 - mlly bkg, 5 - costheta bkg

  //TH2D *h2d_sig = get_mismatch_hist_2D(infile[0]);
  TH2D *h2d_sig = get_noname_hist_2D(infile[0]);

  RooDataHist hist_sig_2d("hist_sig_2d", "sig_twodim_dh", RooArgList(m_lly,ctheta), Import(*(h2d_sig)));

  RooDataHist hist_sig_mlly(  "hist_sig_mlly",   "sig_mlly_dh",   m_lly,  Import(*(get_mismatch_hist(infile[2]))));
  RooDataHist hist_sig_ctheta("hist_sig_ctheta", "sig_ctheta_dh", ctheta, Import(*(get_mismatch_hist(infile[3]))));

  //TH2D *h2d_bkg = get_bkg_hist_2D(infile[1]);
  TH2D *h2d_bkg = get_noname_hist_2D(infile[1]);
  RooDataHist hist_bkg_2d("hist_bkg_2d", "bkg_twodim_dh", RooArgList(m_lly,ctheta), Import(*(h2d_bkg)));

  RooDataHist hist_bkg_mlly(  "hist_bkg_mlly",   "bkg_mlly_dh",   m_lly,  Import(*(get_fake_hist(infile[4]))));
  RooDataHist hist_bkg_ctheta("hist_bkg_ctheta", "bkg_ctheta_dh", ctheta, Import(*(get_fake_hist(infile[5]))));

  //Debug 1.5
  cout << "Really? Here?" << endl;

  //Calls whichever function specified by the input variable "func"  
  RooFitResult *signal_result; RooFitResult *bkg_result; RooFitResult *splusb_result;
  RooProdPdf *signal_fit_2d = fit_2D_signal(m_lly, ctheta, hist_sig_2d, hist_sig_mlly, hist_sig_ctheta, "signal_2D", signal_result);
  RooProdPdf *bkg_fit_2d    = fit_2D_bkg(m_lly, ctheta, hist_bkg_2d, hist_bkg_mlly, hist_bkg_ctheta, "bkg_2D", bkg_result);

  //Initial guesses for splusb fit
  double nsig_init = hist_sig_mlly.sumEntries();
  double nbkg_init = hist_bkg_mlly.sumEntries();
  
  cout << "N_B (from 2D plot): " << hist_bkg_2d.sumEntries() << endl; //N_{S+B}
  cout << "N_S,true (2D plot): " << hist_sig_2d.sumEntries() << endl; //N_{S,true}

  //Add the signal histogram to the background
  hist_bkg_2d.add(hist_sig_2d);

  //Perform S+B fit
  RooAddPdf  *spb_fit_2D  = fit_2D_splusb(signal_fit_2d, bkg_fit_2d, hist_bkg_2d, "splusb_fit", nbkg_init, nsig_init, splusb_result);
  //Debug 2
  cout << "model returned." << endl;

  //Define different plots for the data/fit plot and the residual plot
  //Both will be placed on the same canvas  
  //RooPlot* data_and_fit  = ctheta.frame(ctheta_low, ctheta_high);
  //hist_total_rdh.plotOn( data_and_fit,Binning(nbins), RooFit::Name("hist_total_rdh"));
  //splusb_fit -> plotOn( data_and_fit, RooFit::Components("splusb_fit"), LineStyle(kSolid),  LineColor(TColor::GetColor("#0000ff")), RooFit::Name("splusb_fit_fr"));
  //splusb_fit -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#0000ff")), RooFit::Name((func+"_fr").c_str()));
  //splusb_fit -> plotOn( data_and_fit, RooFit::Components("dscb"),       LineStyle(kDashed), LineColor(TColor::GetColor("#ff0000")), RooFit::Name("dscb_fr"));
  //fit_onepcthetasq -> plotOn( data_and_fit, RooFit::Components(func.c_str()), LineStyle(kDashed), LineColor(TColor::GetColor("#df9f1f")), RooFit::Name(func.c_str()));

  RooRealVar* nsig_ext = static_cast<RooRealVar*>((spb_fit_2D -> coefList()).find("nsig") );
  cout << "N_B (from mlly plot): "    << hist_bkg_mlly.sumEntries()   << endl; //N_{S+B}
  cout << "N_B (from ctheta plot): "  << hist_bkg_ctheta.sumEntries() << endl; //N_{S+B}
  cout << "N_S,true (1D plot): "      << nsig_init                    << endl; //N_{S,true}
  cout << "N_S,true (2D plot): "      << hist_sig_2d.sumEntries()     << endl; //N_{S,true}
  cout << "---------------------------------------------------------" << endl;
  cout << "N_{S + B} (from 2D plot): "<< hist_bkg_2d.sumEntries() << endl; //N_{S+B}
  cout << "N_S,fit: "                 << nsig_ext -> getValV()    << endl; //N_{S,fit}
  cout << "sigma N_S,fit: "           << nsig_ext -> getError()   << endl; //sigma_{N_{s,fit}}
  cout << "sigma upper N_S,fit: "     << nsig_ext -> getErrorHi() << endl; //sigma_{N_{s,fit}}
  cout << "sigma lower N_S,fit: "     << nsig_ext -> getErrorLo() << endl; //sigma_{N_{s,fit}}

  cout << "--------------------------------------------" << endl;
  
  delete nsig_ext;
  delete spb_fit_2D;
  delete signal_result; delete bkg_result; delete splusb_result;
  delete signal_fit_2d;
  delete bkg_fit_2d;
  delete h2d_sig;
  delete h2d_bkg;

  return;
}


std::vector<std::vector<std::string>> twodim_plot_list(std::string name_segment_one, std::vector<std::string> topologies, std::vector<std::string> name_segment_two){
  std::vector<std::vector<std::string>> list_o_hists = {};
  for(std::string top : topologies){
    std::vector<std::string> temp_vec = {};
    for(std::string name_seg_two : name_segment_two){
      temp_vec.push_back(name_segment_one + top + name_seg_two);
    }
    list_o_hists.push_back(temp_vec);
  }
  return list_o_hists;
}



void costheta_2d(){
  //String with paths and plots of interest
  //std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AfterMichaelAugustPush/draw_pico/plots/run3_cat_mlly/";
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AddCategoriesAndKinFit/draw_pico/plots/test_plots_2d/";
  //TestPlots_cat_ttH_had_ll_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt__lumi_nonorm_lin.pdf 

  std::string name_seg_one = "TestPlots_cat_";
  std::string name_seg_two = "_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin";
  //test fitting code with one plot
  std::string top = "ggF";

  std::string test_file_mlly_bkg = file_path + "TestPlots_cat_" + top + "_ll_TIGHT_BASELINE_bkg_llphoton_m__llphoton_m0__wgt_nodiff__lumi_nonorm_lin.root";
  std::string test_file_cth_bkg  = file_path + "TestPlots_cat_" + top + "_ll_TIGHT_BASELINE_bkg_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin.root";
  std::string test_file_2d_bkg   = file_path + "TestPlots_cat_" + top + "_ll_TIGHT_BASELINE_bkg_costheta_mlly__llphoton_m0__wgt_nodiff__lumi_lin.root";

  std::string test_file_mlly_sig = file_path + "TestPlots_cat_" + top + "_ll_TIGHT_BASELINE_sig_llphoton_m__llphoton_m0__wgt_nodiff__lumi_nonorm_lin.root";
  std::string test_file_cth_sig  = file_path + "TestPlots_cat_" + top + "_ll_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin.root";
  std::string test_file_2d_sig   = file_path + "TestPlots_cat_" + top + "_ll_TIGHT_BASELINE_sig_costheta_mlly__llphoton_m0__wgt_nodiff__lumi_lin.root";
  std::vector<std::string> test_file_names = {test_file_2d_sig, test_file_2d_bkg, test_file_mlly_sig, test_file_cth_sig, test_file_mlly_bkg, test_file_cth_bkg};
  fit_and_plot_mllycth(test_file_names, "output/TEST_twodim_fit.txt", "dont_matter");

  //infile idces: 0 - mlly_costheta sig, 1 - mlly_costheta bkg, 2 - mlly sig, 3 - costheta sig, 4 - mlly bkg, 5 - costheta bkg  
  vector<string> plot_endings = {"_ll_TIGHT_BASELINE_bkg_llphoton_m__llphoton_m0__wgt_nodiff__lumi_nonorm_lin.root","_ll_TIGHT_BASELINE_bkg_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin.root",
                                 "_ll_TIGHT_BASELINE_bkg_costheta_mlly__llphoton_m0__wgt_nodiff__lumi_lin.root","_ll_TIGHT_BASELINE_sig_llphoton_m__llphoton_m0__wgt_nodiff__lumi_nonorm_lin.root",
                                 "_ll_TIGHT_BASELINE_sig_costheta__abs_costheta__wgt_nodiff__lumi_nonorm_lin.root","_ll_TIGHT_BASELINE_sig_costheta_mlly__llphoton_m0__wgt_nodiff__lumi_lin.root"};
 
  //2D fit for each topology
  std::vector<std::string> topologies = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  std::vector<std::vector<std::string>> plot_names = twodim_plot_list(file_path + "TestPlots_cat_", topologies, plot_endings);
  for(size_t idx = 0; idx<topologies.size(); idx++){
    fit_and_plot_mllycth(plot_names[idx], "output/twodimfit_ctheta_mlly_" + top + ".txt", "doesntmatter");
  }

  //std::string test_file = "cat_ttH_lep_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  //fit_and_plot_cth(file_path + test_file, "output/test_plot_VBF_ee.pdf", "don");
  //plot_one_func(file_path + test_file, "output/test_plot_VBF_ee.pdf", "pow3");
  //plot_one_func(file_path + test_file, "output/test_plot", "exmg");
  //fit_and_plot_cth(file_path + test_file, "output/test_plot", "dont_matter");

  //Create list of all the plots
  //std::vector<std::vector<std::string>> list_of_plots = plot_list(name_seg_one, name_seg_two);//,false);
  //for(string plot_name : list_of_plots){ fit_and_plot_mllycth(file_path + plot_name + ".root", "output/twodim_fit_" + plot_name, "dont_matter"); }
}


/*

  cout << test_file_mlly_bkg << endl;
  cout << test_file_cth_bkg << endl;
  cout << test_file_2d_bkg << endl;

  cout << test_file_mlly_sig << endl;
  cout << test_file_cth_sig << endl;
  cout << test_file_2d_sig << endl;


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



*/
