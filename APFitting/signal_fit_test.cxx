#include "../APFitting/APF_utilities.cxx"

using namespace APFUtilities;

//This function will plot just one function on the histogram with a residual plot
void plot_signal_test(string infile, string outfile, string func, double signal_boost=1.0){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins = 60; //160
  double m_lly_low = 110;
  double m_lly_high = 140;

  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar m_lly("m_lly" ,"m_lly", m_lly_low, m_lly_high, "GeV");

  //Here signal is assumed to come in with normalization of 10x which gives the 0.1 scale factor. 
  TH1D *sig_hist = APFUtilities::get_sig_hist(infile);
  double nsig_injected = (sig_hist -> Integral())*0.1*signal_boost;
  
  RooDataHist hist_signal_rdh("hist_signal_rdh", "m_lly_sdh", m_lly, Import(*(sig_hist)));
  double nsigtrue = sig_hist -> Integral();
  RooRealVar *nsig = new RooRealVar("nsig", "nsig", nsigtrue, -1000000.0, 1000000.0);

  //Complete fit giving the two histograms as inputs. Saves the result to the pointer result
  RooFitResult *result = new RooFitResult("fit_result","fit_result");
  //SIGNAL FIT HERE

  RooCrystalBall* sig_dscb_fit = nullptr;
  RooGaussian* sig_gauss_fit = nullptr;
  RooAddPdf* signal_fit = nullptr;
  if(func=="dscb"){
   sig_dscb_fit = fit_sig_model(m_lly, hist_signal_rdh, func.c_str()); 
   signal_fit = new RooAddPdf("signal_fit", "signal_fit", RooArgList(*sig_dscb_fit),   RooArgList(*nsig));
  } else {
   sig_gauss_fit = fit_sig_gauss(m_lly, hist_signal_rdh, func.c_str());
   signal_fit = new RooAddPdf("signal_fit", "signal_fit", RooArgList(*sig_gauss_fit),   RooArgList(*nsig));
  }


  //hist_total_rdh.add(hist_signal_rdh); //Have to do this twice. Once in function and once outside. Maybe rethink how this is being done?

  //Make plot for every flat integer
  //if(signal_boost < 10 || static_cast<int>(signal_boost)%10==0){
  //if( (static_cast<int>(signal_boost)*10)%10==0 ){
    //make_plot_with_residual(m_lly, hist_total_rdh, splusb_fit, func.c_str(),outfile,m_lly_low,m_lly_high,nbins,&hist_signal_rdh,result);
    make_signal_plot_with_residual(m_lly, hist_signal_rdh, signal_fit, func.c_str(), outfile, m_lly_low, m_lly_high, nbins);
  //}
  
  RooRealVar* nsig_ext = static_cast<RooRealVar*>((signal_fit -> coefList()).find("nsig") ); 

  //Outputs # of signal events, fit values for nsig
  //Need to get uncertainties!
  /*out_file << signal_boost                 << "," << flush; //signal boost factor
  out_file << hist_total_rdh.sumEntries()  << "," << flush; //N_{S+B}
  out_file << nsig_injected                << "," << flush; //N_{S,true}
  out_file << nsig_ext -> getValV()        << "," << flush; //N_{S,fit}
  out_file << nsig_ext -> getError()       << "," << flush; //sigma_{N_{s,fit}}
  out_file << nsig_ext -> getErrorHi()     << "," << flush; //sigma_{N_{s,fit}} High
  out_file << nsig_ext -> getErrorLo()     << "," << flush; //sigma_{N_{s,fit}} Low
  out_file << result -> minNll()           << flush; //-log(L)
  out_file << endl;                                         //End of saved information

  cout << "signal injection factor: " << signal_boost << endl; //signal boost factor
  cout << "S + B: "                   << hist_total_rdh.sumEntries()  << endl; //N_{S+B}
  cout << "N_S,true: "                << nsig_injected << endl; //N_{S,true}
  cout << "N_S,true: "                << nsigtrue << endl; //N_{S,true}
  cout << "N_S,fit: "                 << nsig_ext -> getValV()        << endl; //N_{S,fit}
  cout << "sigma N_S fit: "           << nsig_ext -> getError()       << endl; //sigma_{N_{s,fit}}
  cout << "sigma N_S fit: "           << nsig_ext -> getErrorHi()       << endl; //sigma_{N_{s,fit}}
  cout << "sigma N_S fit: "           << nsig_ext -> getErrorLo()       << endl; //sigma_{N_{s,fit}}

  cout << "-log(L): "                 << result -> minNll()           << endl; //-log(L)
  cout << "------------------------------------" << endl;                                         //End of saved information
  */

  delete nsig_ext;
  delete signal_fit;
  delete sig_dscb_fit;
  delete sig_gauss_fit;
  delete result;

  return;

}

int signal_fit_test(){
  //Right now all histograms are located in a common area. Will try to move them to be more central
  //std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AfterMichaelAugustPush/draw_pico/plots/run3_cat_mlly/";
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AddCategoriesAndKinFit/draw_pico/plots/mlly_dist_for_fitting/";
  std::string name_seg_one = "cat_";
  std::string name_seg_two = "_ph_all_lly_m__mlly__wgt__lumi_nonorm_lin";
  std::string name_seg_three = "_ph_all_refit_lly_m__llphoton_refit_m__wgt__lumi_nonorm_lin";

  //test fitting code with one plot
  //std::string test_file = "cat_WH_3l_mumu_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";


  //Various test files
  //std::string test_file = "cat_VBF_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_WH_3l_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ZH_MET_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  std::string test_file = "cat_ttH_had_ll_ph_all_lly_m__mlly__wgt__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  std::vector<std::string> topology      = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  std::vector<std::string> func_list     = {"dscb", "gauss"};
  std::vector<std::string> list_of_plots = plot_list(name_seg_one, name_seg_two, true, true);//,false);
  std::vector<std::string> list_of_refit_plots = plot_list(name_seg_one, name_seg_three, true, true);//,false);

  cout << file_path + test_file;
  string func_test = "gauss";
  ofstream out_file; out_file.open("./sig_fit_test_ttH_" + func_test + ".txt");
  plot_signal_test(file_path + test_file , "./output/TEST_ttH_1_sig_fit", func_test, 1.0);
  out_file.close();

  
  vector<double> signal_yields = {};
  for (double idx=-5; idx <= 10; idx += 0.1){ signal_yields.push_back(idx); }
  
  gErrorIgnoreLevel = kWarning;
  RooMsgService::instance().setStreamStatus(2,false);
  //list_of_plots.erase (list_of_plots.begin());
  //topology.erase (topology.begin());
 
  for(unsigned int idx_func = 0; idx_func < func_list.size(); idx_func++){
    for(unsigned int idx_top = 0; idx_top < list_of_plots.size(); idx_top++){
      string top = topology[idx_top];
      string plot = list_of_plots[idx_top];
      cout << file_path + list_of_plots[idx_top] << endl; 
      string func = func_list[idx_func]; 
      plot_signal_test(file_path + plot + ".root" , "output/sig_fit_test_" + top + "_" + func, func, 1.0); 
    }

    for(unsigned int idx_top = 0; idx_top < list_of_refit_plots.size(); idx_top++){
      string top = topology[idx_top];
      string plot = list_of_refit_plots[idx_top];
      cout << file_path + list_of_refit_plots[idx_top] << endl; 
      string func = func_list[idx_func];
      plot_signal_test(file_path + plot + ".root" , "output/sig_fit_test_" + top + "_" + func + "_refit", func, 1.0); 
    }
  }
  
  return 0;

}

