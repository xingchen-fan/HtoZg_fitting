import ROOT
ROOT.gInterpreter.AddIncludePath('../Utilities/AsymGenGaussian.h')
ROOT.gSystem.Load('../Utilities/AsymGenGaussian_cxx.so')

ROOT.gInterpreter.AddIncludePath('../Utilities/EXModGaus.h')
ROOT.gSystem.Load('../Utilities/EXModGaus_cxx.so')


VHttH_file_path = "/net/cms27/cms27r0/abarzdukas/drawPico/gitChangesDirectory/PushControlRegions/draw_pico/plots/an_int_fitting/"

def get_hist(infile, hist_name):
  file = ROOT.TFile(infile)
  canvas_hist = file.Get("canvas")
  hist_return = None
  
  #Finds histogram in file  
  for idx_h in range(0,50000):
    if(hist_return is None):
      hist_return = canvas.FindObject(hist_name+str(idx_h))
    else:
      break
  return hist_return  

def get_hist_data(infile):
  file = ROOT.TFile(infile)
  canvas_hist = file.Get("canvas")
  hist_return = canvas_hist.FindObject("dat_Data")

  return hist_return;

def get_bkg_hist(infile):
  return get_hist(infile,"bkg_Background_")

def get_sig_canvas(infile):
  file = ROOT.TFile(infile)
  return file.Get("canvas")

def get_hist_from_canvas(canvas, hist_name):
  hist_return = None
  
  #Finds histogram in file  
  for idx_h in range(0,50000):
    if(hist_return is None):
      hist_return = canvas.FindObject(hist_name+str(idx_h))
    else:
      break
  return hist_return  


def get_sig_hist(infile):
  return get_hist(infile,"sig_Signal_")

def get_data_hist(infile):
  return get_hist_data(infile)

def get_data_rdh(infile, name, x):
  histogram = get_data_hist(infile)
  return ROOT.RooDataHist(name, name, x, histogram)

def define_mlly():
  mlly_low = 100.0
  mlly_high = 180.0
  mlly_sideband_low = 120.0
  mlly_sideband_high = 130.0

  x = ROOT.RooRealVar("x", "mllg", mlly_low, mlly_high)
  x.setRange('left', mlly_low, mlly_sideband_low)
  x.setRange('right', mlly_sideband_high, mlly_high)
  x.setRange('full', mlly_low, mlly_high)
  
  return x


