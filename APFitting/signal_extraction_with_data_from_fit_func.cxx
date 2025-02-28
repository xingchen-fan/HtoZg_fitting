#include "../APFitting/APF_utilities.cxx"

#include <sstream>
#include <iostream>
using namespace APFUtilities;

RooWorkspace* fit_workspace = new RooWorkspace("generated_blinded_windows","generated_blinded_windows");
RooCategory* categories = new RooCategory("VHttH_categories", "VHttH_categories", {{"ggF",0},{"VBF",1},{"WH_3l",2},{"ZH_MET",3},{"ttH_had",4},{"ttH_lep",5}});

//Combines the background and signal generation to make it easier to loop over!
//Right now the full ranges are generated because I got that to work but it is inefficient
RooDataHist mlly_plot_generated_whole_fit_region(RooRealVar &m_lly, RooAbsPdf *signal_fit, RooAddPdf *background_fit, int nbackground, double nsignal, ostream &out_file, int scale_factor=1, bool isAsimov=false){
  //Setting the global random seed to be a new value to ensure each time a different set of psuedo-random data is produced
  gRandom = new TRandom3(0);
    
  //Generated background
  RooDataHist generated_bkg_dataset = *generate_bkg_data_from_fit(m_lly, background_fit, nbackground, scale_factor, isAsimov);

  //Generate signal datasets
  RooDataHist generated_signal_dataset = *generate_signal_data_from_fit(m_lly, signal_fit, nsignal, scale_factor, isAsimov);

  if(nsignal > 0.01){
    generated_bkg_dataset.add(generated_signal_dataset);
    out_file << generated_signal_dataset.sumEntries() << "," << flush;
  } else {
    out_file << "0," << flush;
  }

  out_file << generated_bkg_dataset.sumEntries() << "," << flush;

  //At least in one spot we can stop some memory leaks.
  //delete generated_bkg_dataset;
  //delete generated_signal_dataset;
  
  //return the copy of the sideband data with the blinded region filled in
  return generated_bkg_dataset;
  
}


//This function will plot just one function on the histogram with a residual plot
void signal_extraction_fullfitreg_generate(string infile_one, string infile_two, string outfile, string func, string category, double signal_boost=1.0, int ntrials=100){
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
  double nsig_expected = (sig_hist -> Integral())*0.1;

  cout << nsig_injected << endl;

  //Fits to the histogram
  //RooAddPdf* sideband_fit = fit_sideband(m_lly, hist_sideband_rdh, nev_sideband, func, m_lly_low, m_lly_high, nbins, category + "_sideband_fit");
  RooAddPdf sideband_fit = *fit_sideband(m_lly, hist_sideband_rdh, nev_sideband, func, m_lly_low, m_lly_high, nbins, category + "_sideband_fit");

  RooAbsPdf* signal_fit_model = signal_fit(m_lly, hist_signal_rdh, "dscb", false, category + "_sig_fit_");

  
  //Expected yields for backgrounds and signal
  double nbkg_fullrange = bkg_fraction_from_sideband(m_lly, &sideband_fit)*nev_sideband;
  double nsig_trial = nsig_expected*signal_boost;

  //Open the outfile that will hold the results of these toys
  ofstream out_file; out_file.open("./siginjection_csv/csv_ffrgen" + category + "_sig_strength_" + to_string(static_cast<int>(signal_boost)) + "_extended.csv");
  out_file << "Nsiginj,Nbkginj,Nsigfit,Nsig_err,Nbkgfit,Nbkgfit_blinded,Nbkgfit_err,Nbonly_fit,NBonly_err,Nll_splusb,Nll_bkgonly" << endl;
  for(int idx = 1; idx <= ntrials; idx++){
    //Here I grab the RooAbsPdf out of the RooAddPdf. Probably could change fit_sideband to return a RooAbsPdf? I do want the events in the sideband though... Idk leave for now
    //RooArgSet *sideband_fit_comps = (sideband_fit -> getComponents());
    //RooAbsPdf *abspdf_sideband_fit = static_cast<RooAbsPdf*>(sideband_fit_comps -> FindObject((func).c_str()));
    //RooArgSet *sideband_fit_comps = (sideband_fit.getComponents());
    RooAddPdf sideband_fit = *fit_sideband(m_lly, hist_sideband_rdh, nev_sideband, func, m_lly_low, m_lly_high, nbins, category + "_sideband_fit");

    RooAddPdf* abspdf_sideband_fit = (static_cast<RooAddPdf*>((sideband_fit).Clone()));
    //RooAddPdf* abspdf_sideband_fit_ptr; sideband_fit-> Copy(abspdf_sideband_fit_ptr); 
    //RooAddPdf  abspdf_sideband_fit_nptr; sideband_fit.Copy(abspdf_sideband_fit_nptr);// =  *abspdf_sideband_fit_ptr;
    //RooAddPdf* abspdf_sideband_fit = &abspdf_sideband_fit_nptr; 

    //Should create copies of the objects without disrupting the originals (it does! yay!)
    RooRealVar nbkg_toy = RooRealVar("nbkg_toy", "nbkg", nbkg_fullrange, 0.75*nbkg_fullrange, 1.25*nbkg_fullrange);
    RooRealVar nsig_toy = RooRealVar("nsig_toy", "nsig", nsig_expected, -100*nsig_expected, 100.0*nsig_expected);
    RooRealVar nbonly_toy = RooRealVar("nbkg_toy", "nbkg", nbkg_fullrange, 0.5*nbkg_fullrange, 1.5*nbkg_fullrange);

    RooAddPdf* splusb_fit_two = new RooAddPdf("splusb_fit_two", "splusb_fit", RooArgList(*abspdf_sideband_fit,*signal_fit_model), RooArgList(nbkg_toy, nsig_toy));
    RooAddPdf* bonly_fit = new RooAddPdf("bonly_fit", "bonly_fit", RooArgList(*abspdf_sideband_fit), RooArgList(nbonly_toy) );
    RooDataHist full_range_with_bdata_two = (mlly_plot_generated_whole_fit_region(m_lly, signal_fit_model, abspdf_sideband_fit, nev_sideband, nsig_trial, out_file));

    //Fit and plot
    RooFitResult *result_spb = splusb_fit_two -> fitTo(full_range_with_bdata_two, SumW2Error(kFALSE), Save(kTRUE), Extended());
    RooFitResult *result_bo  = bonly_fit -> fitTo(full_range_with_bdata_two, SumW2Error(kFALSE), Save(kTRUE), Extended());
    
    //Add the information to the .csv output file
    out_file << nsig_toy.getValV()     << "," << flush; //N_{S,fit}
    out_file << nsig_toy.getError()    << "," << flush; //sigma_{N_{s,fit}}
    out_file << nbkg_toy.getValV()     << "," << flush;
    out_file << (nbkg_toy.getValV())*bkg_fraction_in_blinded_region(m_lly, splusb_fit_two, 122, 128) << "," << flush;
    out_file << nbkg_toy.getError()    << "," << flush;
    out_file << nbonly_toy.getValV()     << "," << flush;
    out_file << nbonly_toy.getError()    << "," << flush;
    out_file << result_spb -> minNll() << "," << flush; //-log(L)
    out_file << result_bo -> minNll()  << flush; //-log(L)
    out_file << endl;                                     //End of saved information

    //Every 100th fit plot the result 
    if(idx%100==0){
      make_fit_plot_with_residual(m_lly, full_range_with_bdata_two, splusb_fit_two, &sideband_fit, func, outfile + "_"  + to_string(idx),m_lly_low, m_lly_high, nbins);
    }


    //out_file << nbkg_toy.getErrorHi()  << "," << flush;
    //out_file << nbkg_toy.getErrorLo()  << "," << flush;
    //out_file << nsig_toy.getErrorHi()  << "," << flush; //sigma_{N_{s,fit}} High
    //out_file << nsig_toy.getErrorLo()  << "," << flush; //sigma_{N_{s,fit}} Low


    //Try to not leak all of the memory!
    delete result_spb;
    delete result_bo;
    delete splusb_fit_two;
    delete bonly_fit;
    //delete abspdf_sideband_fit_ptr;
    //delete sideband_fit_comps;
    
  }

  //Close the output file
  out_file.close();

  //the trials are done!
  return;

}
int signal_extraction_with_data_from_fit_func(string category, double signal_strength){

  //Right now all histograms are located in a common area. Will try to move them to be more central
  std::string file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/AddCategoriesAndKinFit/draw_pico/plots/mlly_dist_for_fitting/";
  std::string name_seg_one = "cat_";
  std::string name_seg_two = "_ph_all_lly_m_data__mlly__wgt__lumi_nonorm_lin";
  std::string name_seg_three = "_ph_all_lly_m__mlly__wgt__lumi_nonorm_lin";
  std::string name_seg_refit_two = "_ph_all_refit_lly_m_data__llphoton_refit_m__wgt__lumi_nonorm_lin";
  std::string name_seg_refit_three = "_ph_all_refit_lly_m__llphoton_refit_m__wgt__lumi_nonorm_lin";


  std::vector<std::string> topology      = {"ggF", "VBF", "WH_3l", "ZH_MET", "ttH_had", "ttH_lep"};
  std::vector<std::string> func_list     = {"agg"};
                                 //options: "agg","swoods","exmg"};
  
  string plot_one = name_seg_one + category + "_ll" + name_seg_two + ".root";
  string plot_two = name_seg_one + category + "_ll" + name_seg_three + ".root";
 
  string plot_refit_one = name_seg_one + category + "_ll" + name_seg_refit_two + ".root";
  string plot_refit_two = name_seg_one + category + "_ll" + name_seg_refit_three + ".root";

  cout << "file path for the plots: " << file_path << endl;
  cout << "plot with sideband data: " << plot_refit_one << endl;
  cout << "plot with signal MC: "     << plot_refit_two << endl;

  int toys = 1000;
  for(string func : func_list){
    string out_plot_name = "./output/test_sig_ext_ffrgen_" + category + "_" + func + "_sigstr_" + to_string(static_cast<int>(signal_strength*100)) +"percent"; 
    signal_extraction_fullfitreg_generate(file_path + plot_one, file_path + plot_two, out_plot_name, func, category + "_" + func, signal_strength, toys);
  }

  //gErrorIgnoreLevel = kWarning;
  //RooMsgService::instance().setStreamStatus(2,false);

  return 0;

}

  //Add to the workspace
  //gSystem -> Load("./Utilities/AsymGenGaussian.cxx");

  //gSystem->AddLinkedLibs("-L./Utilities/AsymGenGaussian.cxx -l*");
  //fit_workspace -> import(categories);


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
          
      signal_extraction_blinded_region(plot_1, plot_2, out_plot, func, top, 1);
    }

    TFile* workspace_out = new TFile(("workspaces/vhtth_workspace_" + func + ".root").c_str(), "RECREATE");
    fit_workspace ->Print();
    fit_workspace -> Write();
  

    //delete workspace_out;
  }
  */
  
   // delete workspace_out;
 
  //RooArgSet *abspdf_components = abspdf_sideband_fit -> getComponents();
  //vector<string> pdf_components = get_component_vector(abspdf_components -> contentsString());

  /*RooRealVar *nbkg = new RooRealVar("nbkg", "nbkg", nbkg_fullrange, 0.5*nbkg_fullrange, 1.5*nbkg_fullrange);
  RooRealVar *nsig = new RooRealVar("nsig", "nsig", nsig_expected, -1*nsig_expected, 10.0*nsig_expected);

  //S + B Fit, but this is not really used for anything right now
  RooAddPdf *splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", RooArgList(*abspdf_sideband_fit,*signal_fit_model), RooArgList(*nbkg,*nsig));


  RooDataHist* full_range_with_bdata = mlly_plot_generated_blinded_window(m_lly, hist_sideband_rdh, signal_fit_model, sideband_fit, nev_sideband, 10);
  make_sideband_plot_with_residual(m_lly, *full_range_with_bdata, splusb_fit, func.c_str(), outfile + "_before", m_lly_low, m_lly_high, nbins);


  cout << "--------------------TEST BEFORE ADD TO WS -------------------------------" <<endl;
  (splusb_fit -> getVariables()) -> Print();
  cout << "-------------------------------------------------------------------------" <<endl;
  */


  //cout << *((signal_fit_model -> getComponents()) -> << endl;
  /*
  RooArgSet *dscb_argset     = signal_fit_model -> getComponents();
  vector<string>dscb_components = get_component_vector(dscb_argset -> contentsString());
  for(string component : dscb_components){
    RooRealVar *curr_parameter = static_cast<RooRealVar*>(dscb_argset -> find(component.c_str()));
    double curr_val = curr_parameter -> getValV();
    curr_parameter -> setRange("",0.9*curr_val,1.1*curr_val);
  }
  */

  //RooArgSet *sideband_fit_comps = sideband_fit -> getObservables(*(sideband_fit -> getComponents()));
  //cout << sideband_fit_comps -> FindObject((category + "_sideband_fit").c_str()) << endl;

  /*
    RooWorkspace *toy_workspace = new RooWorkspace("toy_workspace","toy_workspace");
    RooAddPdf *splusb_fit = new RooAddPdf("splusb_fit", "splusb_fit", RooArgList(*abspdf_sideband_fit,*signal_fit_model), RooArgList(*nbkg,*nsig));

    toy_workspace -> import(splusb_fit);
    toy_workspace -> import(full_range_with_bdata);
  */

  /*
  cout << "---------------------TEST AFTER ADD TO WS -------------------------------" <<endl;
  (splusb_fit -> getVariables()) -> Print();
  cout << "-------------------------------------------------------------------------" <<endl;


  make_sideband_plot_with_residual(m_lly, *full_range_with_bdata, splusb_fit, func.c_str(), outfile + "_after", m_lly_low, m_lly_high, nbins);
  

  //splusb_fit -> fitTo(*full_range_with_bdata, SumW2Error(kTRUE), Save(kTRUE));
  //cout << nsig_injected     << endl;
  //cout << nsig -> getValV() << endl;


  //make plot with the generated signal plus background
  make_sideband_plot_with_residual(m_lly, *full_range_with_bdata, splusb_fit, func.c_str(), outfile, m_lly_low, m_lly_high, nbins);
  */
  /*
  //Add datasets, signal model, and background model
  add_to_workspace(hist_sideband_rdh, signal_fit_model_add, sideband_fit);
  */
  
/*
  //Test code for fitting 
  sideband_fit -> fitTo(hist_sideband_rdh, SumW2Error(kTRUE), Save(kTRUE) );
  splusb_fit   -> fitTo(hist_sideband_rdh, SumW2Error(kTRUE), Save(kTRUE) );

  RooAbsReal* nll_bkg_only = sideband_fit -> createNLL(hist_sideband_rdh);
  RooAbsReal* nll_sig_bkg  = splusb_fit -> createNLL(hist_sideband_rdh);
  
  cout << "NLL bkg only: "  << nll_bkg_only -> getValV() << endl;
  cout << "NLL S+B only: "  <<  nll_sig_bkg -> getValV() << endl;
  //-------------------------------------------------------------------------------//
*/


  //RooAddPdf *splusb_fit = fit_splusb_signalinjection(m_lly, bkg_hist, sig_hist, func.c_str(), m_lly_low, m_lly_high, signal_boost, false, &result);
  //RooRealVar* nsig_ext = static_cast<RooRealVar*>((splusb_fit -> coefList()).find("nsig") ); 
  //delete nsig_ext; 
  
  //delete splusb_fit;
  //delete nsig; delete nbkg; delete generated_signal_dataset;
  //delete signal_fit_model_add; delete signal_fit_model;
  //delete generated_bkg_dataset; delete sideband_fit;
  //delete result;
 

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
  
