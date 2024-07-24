#include "../APFitting/APF_utilities.cxx"

using namespace APFUtilities;

AsymGenGaussian* fit_agg_model_sideband(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180, bool perform_fit = true, RooFitResult **result=nullptr){
  //x.setRange("Range1",100,122);
  //x.setRange("Range2",128,180);

  //Naming the vars used in the fit
  RooRealVar* kappa = new RooRealVar(("kappa" + static_cast<string>(name)).c_str(), "kappa", -0.5,  -100,  100);
  RooRealVar* alpha = new RooRealVar(("alpha" + static_cast<string>(name)).c_str(), "alpha",20.01,   0, 1000);
  RooRealVar* xsi   = new RooRealVar(("xsi" + static_cast<string>(name)).c_str(),   "xsi",    110,  0,  180);

  AsymGenGaussian*  agg_model = new AsymGenGaussian(name,"Asymm. Gen. Gaussian", x, *kappa, *alpha, *xsi, low_d, high_d);
  if(!perform_fit){return agg_model;}
  
/*  if(*result==nullptr){
    agg_model -> fitTo(hist_to_fit,SumW2Error(kTRUE), Save(kTRUE), Range("Range1,Range2"));
  } else {
    *result = agg_model -> fitTo(hist_to_fit,SumW2Error(kTRUE), Save(kTRUE), Range("Range1,Range2"));
  }
*/
  return agg_model;
}

EXModGaus* fit_exmg_model_sideband(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180){
  //x.setRange("Range1",100,122);
  //x.setRange("Range2",128,180);
 
  //Naming the vars used in the fit
  RooRealVar* mu     = new RooRealVar("exmg_mu",    "peak g",  110, 100, 180);//   0, -1.0, 1.0); mu -> setError(1);
  RooRealVar* sig    = new RooRealVar("exmg_sigma", "width g",  10,   0,  80); // 1,   0.,   5); sig -> setError(1);
  RooRealVar* lambda = new RooRealVar("exmg_lambda","lambda",  0.5,  -5,  10); //0.5,   0.,   5); lambda -> setError(1);
  EXModGaus* exmg_model = new EXModGaus(name,"Analytic Exp. (N=1)", x, *mu, *sig, *lambda, low_d, high_d);
  //exmg_model -> fitTo(hist_to_fit,SumW2Error(kTRUE), Save(kTRUE), Range("Range1,Range2"));
  return exmg_model;
}

RooAbsPdf* bkg_fit_sideband(RooRealVar &x, RooDataHist bkg_rdh, std::string func, int x_low, int x_high){
  RooAbsPdf* bkg_fit_model = nullptr;
  if (func=="exmg"){
    bkg_fit_model = static_cast<RooAbsPdf*>(fit_exmg_model_sideband(x,bkg_rdh,func.c_str(),x_low,x_high));
  } else if (func=="agg"){
    bkg_fit_model = static_cast<RooAbsPdf*>(fit_agg_model_sideband(x,bkg_rdh,func.c_str(),x_low,x_high));
  } else {
    std::cout << "Please use one of the allowed set of functions" << std::endl;
  }

  return bkg_fit_model;
}


//This function will plot just one function on the histogram with a residual plot
void plot_sideband_fit(string infile, string outfile, string func){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins = 80; //160
  double m_lly_low = 100;
  double m_lly_high = 180;

  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar m_lly("m_lly" ,"m_lly", m_lly_low, m_lly_high, "GeV");
  RooRealVar m_lly_test("m_lly_test" ,"m_lly_test", m_lly_low, m_lly_high, "GeV");

  //Gets histograms and scales the signal by signal boost (defaulted to 1!)
  TH1D *data_hist = APFUtilities::get_data_hist(infile);

  //Complete fit giving the two histograms as inputs. Saves the result to the pointer result
  RooFitResult *result = new RooFitResult("fit_result","fit_result");

  RooDataHist hist_sideband_rdh("hist_signal_rdh", "m_lly_dh", m_lly, Import(*(data_hist)));
  double nsideband_d = data_hist -> Integral();
  RooRealVar *nsideband = new RooRealVar("nsig", "nsig", nsideband_d, 0.0, 1000000.0);

  RooAbsPdf* sideband_bkg = bkg_fit_sideband(m_lly, hist_sideband_rdh, func, m_lly_low, m_lly_high);
  RooAddPdf* sideband_fit = new RooAddPdf("signal_fit", "signal_fit", RooArgList(*sideband_bkg),   RooArgList(*nsideband));
  m_lly.setRange("Range1",100,122);
  m_lly.setRange("Range2",128,180);
  sideband_fit -> fitTo(hist_sideband_rdh,SumW2Error(kTRUE), Save(kTRUE), Range("Range1,Range2"));
  m_lly.removeRange("Range1");
  m_lly.removeRange("Range2");

  //Make plot for every flat integer
  make_sideband_plot_with_residual(m_lly, hist_sideband_rdh, sideband_fit, func.c_str(), outfile, m_lly_low, m_lly_high, nbins);

  /*
  delete sideband_fit;
  delete sideband_bkg;
  delete nsideband;
  delete result;
  delete data_hist;
  */
  return;
}

int sideband_data_fit(){
  //Right now all histograms are located in a common area. Will try to move them to be more central
  //std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AfterMichaelAugustPush/draw_pico/plots/run3_cat_mlly/";
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AddCategoriesAndKinFit/draw_pico/plots/mlly_dist_for_fitting/";
  std::string name_seg_one = "cat_";
  std::string name_seg_two = "_ph_all_lly_m_data__mlly__wgt__lumi_nonorm_lin";
  std::string name_seg_three = "_ph_all_refit_lly_m_data__llphoton_refit_m__wgt__lumi_nonorm_lin";



  //Various test files
  //std::string test_file = "cat_VBF_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_WH_3l_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ZH_MET_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  std::string test_file = "cat_ttH_had_ll_ph_all_lly_m_data__mlly__wgt__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  std::vector<std::string> topology      = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  std::vector<std::string> func_list     = {"exmg", "agg"};
  std::vector<std::string> list_of_plots = plot_list(name_seg_one, name_seg_two, true, true);//,false);
  std::vector<std::string> list_of_plots_refit = plot_list(name_seg_one, name_seg_three, true, true);//,false);

  cout << file_path + test_file;
  string func_test = "agg";
  plot_sideband_fit(file_path + test_file , "./output/TEST_ttH_had_datasideband", func_test);
  
  gErrorIgnoreLevel = kWarning;
  RooMsgService::instance().setStreamStatus(2,false);
  list_of_plots.erase(list_of_plots.begin());
  list_of_plots_refit.erase(list_of_plots.begin());
  topology.erase(topology.begin());
 
  
  for(unsigned int idx_func = 0; idx_func < func_list.size(); idx_func++){
    for(unsigned int idx_top = 0; idx_top < list_of_plots.size(); idx_top++){
      string top = topology[idx_top];
      string plot = list_of_plots[idx_top];
      string func = func_list[idx_func];

      cout << "input file: " << file_path + list_of_plots[idx_top] << endl;
      cout << "output file: " << "output/sideband_fit_plot_" + top + "_" + func << endl; 

      plot_sideband_fit(file_path + plot + ".root" , "output/sideband_fit_plot_" + top + "_" + func, func); 
    }

    for(unsigned int idx_top = 0; idx_top < list_of_plots_refit.size(); idx_top++){
      string top = topology[idx_top];
      string plot = list_of_plots_refit[idx_top];
      string func = func_list[idx_func];

      cout << "input file: " << file_path + list_of_plots[idx_top] << endl;
      cout << "output file: " << "output/sideband_fit_plot_" + top + "_" + func << endl; 

      plot_sideband_fit(file_path + plot + ".root" , "output/sideband_fit_plot_" + top + "_" + func + "_refit", func); 
    }
  }
  
  return 0;
}

