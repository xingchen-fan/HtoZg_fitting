import ROOT
import os
import sys
import json
import argparse
sys.path.append(os.path.abspath("../Utilities/"))
from bkg_functions_fit import *
from bkg_functions_class import *
from sig_functions_class import *
from Xc_Minimizer import *
from plot_utility import *
from sample_reader import *
from profile_class import *
from sig_functions_class import *

x = ROOT.RooRealVar("x", "mllg", 100, 165)
MH = ROOT.RooRealVar("MH","MH"       ,124.7, 120., 130.)

for CAT in ['ggf1','ggf2','ggf3','ggf4']:
    combine_model = combineSignal(x, MH, CAT, '../Config/config_DSCB.json')
    print(CAT, ' = ', combine_model.ntot)
