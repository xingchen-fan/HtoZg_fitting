namespace APFUtilities{ 

//
std::vector<string> filename_vector(string prefix, std::vector<string> &filenames);
TH1D* get_hist(string infile, string hist_name);
TH1D* get_bkg_hist(string infile);
TH1D* get_sig_hist(string infile);

//Signal fit functions
RooCrystalBall* fit_sig_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);
RooGaussian* fit_sig_gauss(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);

//Bkg. fit functions
ModGaus* fit_modgauss_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);
EXModGaus* fit_exmg_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180);
AsymGenGaussian* fit_agg_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180);
RooGamma* fit_gam_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name, double low_d=100, double high_d=180);
RooFFTConvPdf* fit_pow3_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);
RooFFTConvPdf* fit_pow1_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);
RooFFTConvPdf* fit_exp1_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);
RooFFTConvPdf* fit_lau1_model(RooRealVar &x, RooDataHist &hist_to_fit, const char *name);

//Functions are ones to actually call to make plots. First one only plots one function with a residual. Second plots all relevant background functions
//Maybe expand to be able to choose which functions you want?
void plot_one_func(string infile, string outfile, string func);
void fit_and_plot(string infile,string outfile);

//This function will create a list of strings with each category's plot to make looping easier.
std::vector<std::string> plot_list(std::string name_segment_one, std::string name_segment_two, bool all=true);

RooAddPdf* fit_splusb_siginjection(RooRealVar &x, TH1D* bkg_rdh, TH1D* signal_rdh, string func, int x_low, int x_high, double nbkg_init=10000, double nsig_init=10, double scale_factor = 1);


}
