#include "../APFitting/APF_utilities.cxx"

using namespace APFUtilities;



RooWorkspace* fit_workspace = new RooWorkspace("generated_blinded_windows","generated_blinded_windows");
RooCategory* categories = new RooCategory("VHttH_categories", "VHttH_categories", {{"ggF",0},{"VBF",1},{"WH_3l",2},{"ZH_MET",3},{"ttH_had",4},{"ttH_lep",5}});

void add_to_workspace(RooDataHist data_hist, RooAbsPdf* sig_pdf, RooAbsPdf* bkg_pdf){
  fit_workspace -> import(*bkg_pdf);
  fit_workspace -> import(*sig_pdf);
  fit_workspace -> import(data_hist);
}

//Is not needed? Produce workspace says to create these classes in 
void add_custom_class(RooWorkspace* workspace, string function ){
  //gROOT -> ProcessLine((".x ./Utilities/" + function + ".cxx").c_str());

  //TClass* import_class = TClass::GetClass(function.c_str());
  TClass* import_class = new TClass(function.c_str());

  import_class -> SetDeclFile(("./Utilities/" + function + ".h").c_str(), 0);
  import_class -> SetImplFileName(("./Utilities/" + function + ".cxx").c_str());

  workspace -> addClassDeclImportDir("./");
  workspace -> addClassImplImportDir("./");
  workspace -> importClassCode(import_class);

}


//This function will plot just one function on the histogram with a residual plot
void plot_generated_blinded_region(string infile_one, string infile_two, string outfile, string func, string category, double signal_boost=1.0){
  //Sets a couple of constants used. Maybe shift to input variables with default values
  int nbins = 80; //160
  double m_lly_low = 100;
  double m_lly_high = 180;

  //Creates fit x-axis variable and converts histogram to RooDataHist
  RooRealVar m_lly("m_lly" ,"m_lly", m_lly_low, m_lly_high, "GeV");

  //Gets histograms and scales the signal by signal boost (defaulted to 1!)
  TH1D *sideband_data = APFUtilities::get_data_hist(infile_one);
  RooDataHist hist_sideband_rdh( (category + "_hist_dsideband_genblinded").c_str()  , "hist_sideband_rdh", m_lly, Import(*(sideband_data)));
  double nev_sideband = sideband_data -> Integral();

  //Here signal is assumed to come in with normalization of 10x which gives the 0.1 scale factor. 
  TH1D *sig_hist = APFUtilities::get_sig_hist(infile_two);
  RooDataHist hist_signal_rdh("hist_signal_rdh", "hist_signal_rdh", m_lly, Import(*(sig_hist)));
  double nsig_injected = (sig_hist -> Integral())*signal_boost*0.1;

  cout << nsig_injected << endl;

  //Complete fit giving the two histograms as inputs. Saves the result to the pointer result
  //RooFitResult *result = new RooFitResult("fit_result","fit_result");

  //Sideband fit
  RooAddPdf* sideband_fit = fit_sideband(m_lly, hist_sideband_rdh, nev_sideband, func, m_lly_low, m_lly_high, nbins, category + "_sideband_fit");

  //Setting the global random seed to be a new value to ensure each time a different set of psuedo-random data is produced
  gRandom = new TRandom3(0);

  //Scale factor for number of events. 1 = same number as expected in data. 
  int scale_factor = 1;
  bool isAsimov = false;

  //Generated dataset from sideband data
  
  RooDataHist* generated_bkg_dataset = generate_bkg_data_from_fit(m_lly, sideband_fit, nev_sideband, scale_factor, isAsimov);

  //Setting the global random seed to be a new value to ensure each time a different set of psuedo-random data is produced
  gRandom = new TRandom3(0);

  //Creating the random signal dataset. 
  RooAbsPdf* signal_fit_model = signal_fit(m_lly, hist_signal_rdh, "dscb", false, category + "_sig_fit_");
  RooAddPdf* signal_fit_model_add = new RooAddPdf("signal_fit_model_add", "signal_fit_model_add", RooArgList(*signal_fit_model),   RooArgList(nsig_injected));

  //Generate signal datasets
  RooDataHist* generated_signal_dataset = generate_signal_data_from_fit(m_lly, signal_fit_model, nsig_injected, scale_factor, isAsimov);

  //Generated bkg + generated signal
  hist_sideband_rdh.add(*generated_bkg_dataset, "m_lly > 122 && m_lly < 128");
  hist_sideband_rdh.add(*generated_signal_dataset, "m_lly > 122 && m_lly < 128");
  //generated_bkg_dataset -> add(*generated_signal_dataset);

  //Expected yields for backgrounds 
  RooRealVar *nbkg = new RooRealVar("nbkg", "nbkg", nev_sideband,  0         , 1000000.0); //bkg_fraction_from_sideband(m_lly, sideband_fit)*nev_sideband, 0, 1000000.0);
  RooRealVar *nsig = new RooRealVar("nsig", "nsig", nsig_injected, -1000000.0, 1000000.0);

  //S + B Fit, but this is not really used for anything right now
  RooAddPdf *splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", RooArgList(*sideband_fit,*signal_fit_model), RooArgList(*nbkg,*nsig));
 
  //make plot with the generated signal plus background
  //make_sideband_plot_with_residual
  //make_fit_plot_with_residual(m_lly, hist_sideband_rdh, splusb_fit, func.c_str(), outfile, m_lly_low, m_lly_high, nbins);

  //Add datasets, signal model, and background model
  add_to_workspace(hist_sideband_rdh, signal_fit_model_add, sideband_fit);

  //Test code for fitting 
  sideband_fit -> fitTo(hist_sideband_rdh, SumW2Error(kTRUE), Save(kTRUE) );
  splusb_fit   -> fitTo(hist_sideband_rdh, SumW2Error(kTRUE), Save(kTRUE) );

  RooAbsReal* nll_bkg_only = sideband_fit -> createNLL(hist_sideband_rdh);
  RooAbsReal* nll_sig_bkg  = splusb_fit -> createNLL(hist_sideband_rdh);

  make_fit_plot_with_residual(m_lly, hist_sideband_rdh, sideband_fit, splusb_fit, func.c_str(), outfile, m_lly_low, m_lly_high, nbins);


  cout << "NLL bkg only: "  << nll_bkg_only -> getValV() << endl;
  cout << "NLL S+B only: "  <<  nll_sig_bkg -> getValV() << endl;
  //-------------------------------------------------------------------------------//



  //RooAddPdf *splusb_fit = fit_splusb_signalinjection(m_lly, bkg_hist, sig_hist, func.c_str(), m_lly_low, m_lly_high, signal_boost, false, &result);
  //RooRealVar* nsig_ext = static_cast<RooRealVar*>((splusb_fit -> coefList()).find("nsig") ); 
  //delete nsig_ext; 
  
  delete splusb_fit;
  delete nsig; delete nbkg; delete generated_signal_dataset;
  delete signal_fit_model_add; delete signal_fit_model;
  delete generated_bkg_dataset; delete sideband_fit;
  //delete result;

  return;

}

int fit_with_generated_blinded_region(){
  //Add to the workspace
  gSystem -> Load("./Utilities/AsymGenGaussian.cxx");

  //gSystem->AddLinkedLibs("-L./Utilities/AsymGenGaussian.cxx -l*");
  //fit_workspace -> import(categories);

  //Right now all histograms are located in a common area. Will try to move them to be more central
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AddCategoriesAndKinFit/draw_pico/plots/mlly_dist_for_fitting/";
  std::string name_seg_one = "cat_";
  std::string name_seg_two = "_ph_all_lly_m_data__mlly__wgt__lumi_nonorm_lin";
  std::string name_seg_three = "_ph_all_lly_m__mlly__wgt__lumi_nonorm_lin";
  std::string name_seg_refit_two = "_ph_all_refit_lly_m_data__llphoton_refit_m__wgt__lumi_nonorm_lin";
  std::string name_seg_refit_three = "_ph_all_refit_lly_m__llphoton_refit_m__wgt__lumi_nonorm_lin";

 
  //Various test files
  std::string test_file_1 = "cat_ttH_had_ll_ph_all_lly_m_data__mlly__wgt__lumi_nonorm_lin.root";
  std::string test_file_2 = "cat_ttH_had_ll_ph_all_lly_m__mlly__wgt__lumi_nonorm_lin.root";

  //std::string test_file = "cat_ttH_lep_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";

  std::vector<std::string> topology      = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  std::vector<std::string> func_list     = {"agg"};
                                 //options: "agg","swoods","exmg"};

  //Add the classes to the rooworkspace - not needed for now?
  //add_custom_class(fit_workspace,"AsymGenGaussian");
  //add_custom_class(fit_workspace,"ExModGaus");


  //fit_workspace -> importClassCode("AsymGenGaussian");
  //fit_workspace -> importClassCode("EXModGaus");

  std::vector<std::string> list_of_plots_1 = plot_list(name_seg_one, name_seg_two, true, true);//,false);
  std::vector<std::string> list_of_plots_2 = plot_list(name_seg_one, name_seg_three, true, true);//,false);

  
  cout << file_path + test_file_1 << endl;
  cout << file_path + test_file_2 << endl;
  string func_test = "agg";//"agg";
  plot_generated_blinded_region(file_path + test_file_1, file_path + test_file_2 , "./output/test_generate_blinded_ttH_had", func_test, "ttH_had", 1.0);
 
  gErrorIgnoreLevel = kWarning;
  RooMsgService::instance().setStreamStatus(2,false);

  /*   
  for(unsigned int idx_func = 0; idx_func < func_list.size(); idx_func++){
    string func = func_list[idx_func];

    for(unsigned int idx_top = 0; idx_top < list_of_plots_1.size(); idx_top++){
      string top = topology[idx_top];
      string plot_1 = file_path + list_of_plots_1[idx_top] + ".root";
      string plot_2 = file_path + list_of_plots_2[idx_top] + ".root";
      string out_plot = "output/generated_blinded_region_" + top + "_" + func + to_string(1);

      cout << "Plot (data): " << plot_1 << endl;
      cout << "Plot (MC): " << plot_2 << endl;
          
      plot_generated_blinded_region(plot_1, plot_2, out_plot, func, top, 1);
    }

    TFile* workspace_out = new TFile(("workspaces/vhtth_workspace_" + func + ".root").c_str(), "RECREATE");
    fit_workspace ->Print();
    fit_workspace -> Write();
  

    //delete workspace_out;
  }
  
   // delete workspace_out;
  */

  return 0;

}


  //test fitting code with one plot
  //std::string test_file = "cat_WH_3l_mumu_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ttH_lep_ee_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_VBF_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_WH_3l_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";
  //std::string test_file = "cat_ZH_MET_ll_ph_all_lly_m__mlly__wgt_1invfb__lumi_nonorm_lin.root";


  //list_of_plots.erase (list_of_plots.begin());
  //topology.erase (topology.begin());

  //Outputs # of signal events, fit values for nsig
  //Need to get uncertainties!
/*  out_file << signal_boost                 << "," << flush; //signal boost factor
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
      //out_file.close();

     //ofstream out_file; out_file.open("dontcarerightnow.txt");//./sig_ext_" + top + "_" + func + ".txt");
      //out_file << "signal_factor,spb_yield,Ns_true,Ns_fit,Ns_fit_unc,Ns_fit_unc_high,Ns_fit_unc_low,NLL" <<endl;
      //for(double yield_factor : signal_yields){  }
  
