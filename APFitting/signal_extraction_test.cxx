#include "../APFitting/APF_utilities.cxx"

using namespace APFUtilities;

//This function will plot just one function on the histogram with a residual plot
void plot_signal_extraction(string infile, string outfile, string func, ofstream &out_file, double signal_boost=1.0){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins = 80; //160
  double m_lly_low = 100;
  double m_lly_high = 180;

  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar m_lly("m_lly" ,"m_lly", m_lly_low, m_lly_high, "GeV");

  //Gets histograms and scales the signal by signal boost (defaulted to 1!)
  TH1D *bkg_hist = APFUtilities::get_bkg_hist(infile);
  double min= bkg_hist -> GetMinimum();
  while(min < 0){
    bkg_hist -> SetBinContent( bkg_hist->GetMinimumBin(), 0 );
    min = bkg_hist -> GetMinimum();
  }
  //Here signal is assumed to come in with normalization of 10x which gives the 0.1 scale factor. 
  TH1D *sig_hist = APFUtilities::get_sig_hist(infile);
  double nsig_injected = (sig_hist -> Integral())*0.1*signal_boost;

  //Complete fit giving the two histograms as inputs. Saves the result to the pointer result
  RooFitResult *result = new RooFitResult("fit_result","fit_result");
  RooAddPdf *splusb_fit = fit_splusb_signalinjection(m_lly, bkg_hist, sig_hist, func.c_str(), m_lly_low, m_lly_high, signal_boost, false, &result);

  RooDataHist hist_total_rdh("hist_total_rdh", "m_lly_dh", m_lly, Import(*(bkg_hist)));
  double nbkg_init = hist_total_rdh.sumEntries();

  //sig_hist -> Scale(0.1*signal_boost);
  RooDataHist hist_signal_rdh("hist_signal_rdh", "m_lly_sdh", m_lly, Import(*(sig_hist)));
  double nsigtrue = sig_hist -> Integral();
  hist_total_rdh.add(hist_signal_rdh); //Have to do this twice. Once in function and once outside. Maybe rethink how this is being done?

  //Make plot for every flat integer
  //if(signal_boost < 10 || static_cast<int>(signal_boost)%10==0){
  if(static_cast<int>((signal_boost)*10)%10==0){
    make_plot_with_residual(m_lly, hist_total_rdh, splusb_fit, func.c_str(),outfile,m_lly_low,m_lly_high,nbins,&hist_signal_rdh,result);
  }
  
  RooRealVar* nsig_ext = static_cast<RooRealVar*>((splusb_fit -> coefList()).find("nsig") ); 

  //Outputs # of signal events, fit values for nsig
  //Need to get uncertainties!
  out_file << signal_boost                 << "," << flush; //signal boost factor
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


  delete nsig_ext;
  delete splusb_fit;
  delete result;

  return;

}

int signal_extraction_test(){
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AfterMichaelAugustPush/draw_pico/plots/run3_cat_mlly/";
  std::string name_seg_one = "cat_";
  std::string name_seg_two = "_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin";
  
  //test fitting code with one plot
  //std::string test_file = "cat_WH_3l_mumu_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";


  //Various test files
  //std::string test_file = "cat_VBF_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_WH_3l_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ZH_MET_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  std::string test_file = "cat_ttH_had_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  std::vector<std::string> topology      = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  std::vector<std::string> list_of_plots = plot_list(name_seg_one, name_seg_two, true, true);//,false);

  string func_test = "agg";
  ofstream out_file; out_file.open("./sig_test_ttH_" + func_test + ".txt");
  plot_signal_extraction(file_path + test_file , "./output/TEST_ttH_1", func_test, out_file, 5.5);
  out_file.close();


  vector<double> signal_yields = {};
  for (double idx=-5; idx <= 10; idx += 0.1){ signal_yields.push_back(idx); }
  
  gErrorIgnoreLevel = kWarning;
  RooMsgService::instance().setStreamStatus(2,false);
  list_of_plots.erase (list_of_plots.begin());
  topology.erase (topology.begin());
  
  for(unsigned int idx_top = 0; idx_top < list_of_plots.size(); idx_top++){

  //int idx_top=1;
    string top = topology[idx_top];
    string plot = list_of_plots[idx_top];
    cout << file_path + list_of_plots[idx_top] << endl; 
    string func = "agg";
    ofstream out_file; out_file.open("./sig_ext_" + top + "_" + func + ".txt");
    out_file << "signal_factor,spb_yield,Ns_true,Ns_fit,Ns_fit_unc,Ns_fit_unc_high,Ns_fit_unc_low,NLL" <<endl;
    for(double yield_factor : signal_yields){ plot_signal_extraction(file_path + plot + ".root" , "output/sig_inj_plot_" + top + "_" + func + "_" + to_string(yield_factor), func, out_file, yield_factor); }
    out_file.close();
  }
  

  return 0;

}

