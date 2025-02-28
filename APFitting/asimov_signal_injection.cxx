#include "../APFitting/APF_utilities.cxx"

#include <sstream>
#include <iostream>
using namespace APFUtilities;

RooWorkspace* fit_workspace = new RooWorkspace("generated_blinded_windows","generated_blinded_windows");
RooCategory* categories = new RooCategory("VHttH_categories", "VHttH_categories", {{"ggF",0},{"VBF",1},{"WH_3l",2},{"ZH_MET",3},{"ttH_had",4},{"ttH_lep",5}});

//Combines the background and signal generation to make it easier to loop over!
//Right now the full ranges are generated because I got that to work but it is inefficient
RooDataHist* mlly_plot_generated_blinded_window(RooRealVar &m_lly, RooDataHist hist_sideband_rdh, RooAbsPdf *signal_fit, RooAddPdf *background_fit, int nbackground, double nsignal, ofstream &out_file, int scale_factor=1, bool isAsimov=false){
  //Setting the global random seed to be a new value to ensure each time a different set of psuedo-random data is produced
  gRandom = new TRandom3(0);
    
  //Generated background 
  RooDataHist* generated_bkg_dataset = generate_bkg_data_from_fit(m_lly, background_fit, nbackground, scale_factor, isAsimov);

  //Generate signal datasets
  RooDataHist* generated_signal_dataset = generate_signal_data_from_fit(m_lly, signal_fit, nsignal, scale_factor, isAsimov);

  out_file << generated_signal_dataset -> sumEntries("m_lly > 122 && m_lly < 128") << "," << flush;
  out_file << generated_bkg_dataset -> sumEntries("m_lly > 122 && m_lly < 128") << "," << flush;

  //Generated bkg + generated signal
  hist_sideband_rdh.add(*generated_bkg_dataset, "m_lly > 122 && m_lly < 128");
  
  if(nsignal > 0.01){
    hist_sideband_rdh.add(*generated_signal_dataset, "m_lly > 122 && m_lly < 128");
  }

 
  if(isAsimov){
    generated_bkg_dataset -> add(*generated_signal_dataset);
    return generated_bkg_dataset;
  }

  //In one spot at least we can stop some memory leaks.
  delete generated_bkg_dataset;
  delete generated_signal_dataset;


  //return the copy of the sideband data with the blinded region filled in
  return (new RooDataHist(hist_sideband_rdh));
}


//This function will plot just one function on the histogram with a residual plot
void asimov_signal_extraction_blinded_region(string infile_one, string infile_two, string outfile, string func, string category, ofstream &out_file, double signal_boost=1.0){
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
  RooAddPdf* sideband_fit = fit_sideband(m_lly, hist_sideband_rdh, nev_sideband, func, m_lly_low, m_lly_high, nbins, category + "_sideband_fit");
  RooAbsPdf* signal_fit_model = signal_fit(m_lly, hist_signal_rdh, "dscb", false, category + "_sig_fit_");

  //Here I grab the RooAbsPdf out of the RooAddPdf. Probably could change fit_sideband to return a RooAbsPdf? I do want the events in the sideband though... Idk leave for now
  RooArgSet *sideband_fit_comps = (sideband_fit -> getComponents());
  RooAbsPdf *abspdf_sideband_fit = static_cast<RooAbsPdf*>(sideband_fit_comps -> FindObject((func).c_str()));
  
  //Expected yields for backgrounds and signal
  double nbkg_fullrange = bkg_fraction_from_sideband(m_lly, sideband_fit)*nev_sideband;
  double nsig_trial = nsig_expected*signal_boost;

  RooRealVar nbkg_toy = RooRealVar("nbkg_toy", "nbkg", nbkg_fullrange, 0.9*nbkg_fullrange, 1.1*nbkg_fullrange);
  RooRealVar nsig_toy = RooRealVar("nsig_toy", "nsig", nsig_expected, -100*nsig_expected, 100.0*nsig_expected);
  RooRealVar nbonly_toy = RooRealVar("nbkg_toy", "nbkg", nbkg_fullrange, 0.9*nbkg_fullrange, 1.1*nbkg_fullrange);

  RooAddPdf* splusb_fit_two = new RooAddPdf("splusb_fit_two", "splusb_fit", RooArgList(*abspdf_sideband_fit,*signal_fit_model), RooArgList(nbkg_toy, nsig_toy));
  RooAddPdf* bonly_fit = new RooAddPdf("bonly_fit", "bonly_fit", RooArgList(*abspdf_sideband_fit), RooArgList(nbonly_toy) );

  RooDataHist full_range_with_bdata_two = *(mlly_plot_generated_blinded_window(m_lly, hist_sideband_rdh, signal_fit_model, sideband_fit, nev_sideband, nsig_trial, out_file, 1, true));

  //Fit and plot
  RooFitResult *result_spb = splusb_fit_two -> fitTo(full_range_with_bdata_two, SumW2Error(kFALSE), Save(kTRUE), Extended());
  RooFitResult *result_bo  = bonly_fit -> fitTo(full_range_with_bdata_two, SumW2Error(kFALSE), Save(kTRUE), Extended());
 
  //make_sideband_plot_with_residual(m_lly, full_range_with_bdata_two, splusb_fit_two, func.c_str(), outfile + "_"  + to_string(static_cast<int>(signal_boost*100)) + "percent", m_lly_low, m_lly_high, nbins); 
  make_fit_plot_with_residual(m_lly, full_range_with_bdata_two, splusb_fit_two, sideband_fit, func, outfile,m_lly_low, m_lly_high, nbins);

  out_file << nsig_toy.getValV()   << "," << flush;
  out_file << nsig_toy.getError()  << "," << flush; //sigma_{N_{s,fit}}
  out_file << nbkg_toy.getValV()   << "," << flush;
  out_file << nbkg_toy.getError()  << "," << flush;
  out_file << (nbkg_toy.getValV())*bkg_fraction_in_blinded_region(m_lly, splusb_fit_two, 122, 128) << "," << flush;
  out_file << nbonly_toy.getValV() << "," << flush;
  out_file << nbonly_toy.getError()  << "," << flush;
  out_file << result_spb -> minNll() << "," << flush; //-log(L)
  out_file << result_bo -> minNll()  << flush; //-log(L)
  out_file << endl;                                     //End of saved information


  delete splusb_fit_two; delete bonly_fit; delete result_spb; delete result_bo;

  //the trials are done!
  return;

}
int asimov_signal_injection(){

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
 
  for(string category : topology){
    for(string func : func_list){

      string plot_one = name_seg_one + category + "_ll" + name_seg_two + ".root";
      string plot_two = name_seg_one + category + "_ll" + name_seg_three + ".root";
     
      string plot_refit_one = name_seg_one + category + "_ll" + name_seg_refit_two + ".root";
      string plot_refit_two = name_seg_one + category + "_ll" + name_seg_refit_three + ".root";

      cout << "file path for the plots: " << file_path << endl;
      cout << "plot with sideband data: " << plot_refit_one << endl;
      cout << "plot with signal MC: "     << plot_refit_two << endl;

      ofstream out_file; out_file.open("./siginjection_csv/asimov_siginj_" + category + "_" + func + ".csv");
      out_file << "Nsiginj,Nbkginj,Nsigfit,Nsig_err,Nbkgfit,Nbkg_err,Nbonly,Nbonly_err,Nll_splusb,Nll_bkgonly" << endl;
      for(double signal_strength = 1; signal_strength < 10; signal_strength++){
        string out_plot_name = "./output/asimov_sig_ext_blinded_" + category + "_" + func + "_sigstr_" + to_string(static_cast<int>(signal_strength*100)) +"percent"; 
        asimov_signal_extraction_blinded_region(file_path + plot_one, file_path + plot_two, out_plot_name, func, category + "_" + func, out_file, signal_strength);
      }
      out_file.close();

    }
  }
  //gErrorIgnoreLevel = kWarning;
  //RooMsgService::instance().setStreamStatus(2,false);

  return 0;

}

