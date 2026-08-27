import ROOT
import json
import argparse
parser = argparse.ArgumentParser(description = "Make workspace")
parser.add_argument('-con', '--config', help = 'Configuration', default='../Config/chi2_config_xgboost_nodrop.json')
parser.add_argument('-c', '--cat', help="category")
args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
lowx = setting["Range"]

#----------------------------------------------------------------------------------------
# This script reads a toy from a combine output root file and saves the toy histogram into a wkspace
#----------------------------------------------------------------------------------------

ROOT.gInterpreter.Declare("""
RooDataSet readToy(TString filename, string num){
auto file_ = new TFile(filename, "READ");
TDirectoryFile *dir = (TDirectoryFile *)file_->Get("toys");
RooDataSet *oneToy  = (RooDataSet *)dir->Get(("toy_" + num).c_str());
return *oneToy;
}
""")

x = ROOT.RooRealVar("CMS_hzg_mass_"+CAT, "CMS_hzg_mass_"+CAT, lowx, lowx + 65.)
x.setBins(260)
file_ = 'higgsCombine.'+CAT+'.oneToy.Significance.mH125.123456.root'
hist_toy = ROOT.readToy(file_, '2')
#hist_toy.SetNameTitle('hist_'+CAT, 'hist_'+CAT)
hist_save = ROOT.RooDataHist('hist_'+CAT, 'hist_'+CAT, x)
for i in range(262):
    hist_toy.get(i-1)
    hist_save.get(i-1)
    hist_save.set(hist_toy.weight())
    #print('bin ', i-2, ' = ', hist_toy.weight())

if hist_save.isNonPoissonWeighted(): print('weighted')
f_out2 = ROOT.TFile('workspaces/workspace_toy_'+CAT+'.root', "RECREATE")
w_bkg = ROOT.RooWorkspace("workspace_bkg","workspace_bkg")
getattr(w_bkg, "import")(hist_save)
w_bkg.Print()
w_bkg.Write()
f_out2.Close()
