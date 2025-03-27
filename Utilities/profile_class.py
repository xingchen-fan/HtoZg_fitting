import ROOT
import json
from bkg_functions_class import *

class profileClass:
    def __init__(self, x, mu_gauss, cat='', config=''):
        jfile_ = open(config, 'r')
        configs_ = json.load(jfile_)
        setting = configs_[cat]
        self.setting_ = setting
        self.bern2_model = Bern2Class(x, mu_gauss, cat, setting["bern2"]['p0'], setting["bern2"]['p1_init'], setting["bern2"]['p2_init'], setting["bern2"]['bond'],setting["bern2"]['sigma_init'],setting["bern2"]['step_init'], setting["bern2"]['fix_sigma'])
        self.bern3_model = Bern3Class(x, mu_gauss, cat, setting["bern3"]['p0'], setting["bern3"]['p1_init'], setting["bern3"]['p2_init'], setting["bern3"]['p3_init'], setting["bern3"]['bond'],setting["bern3"]['sigma_init'],setting["bern3"]['step_init'], setting["bern3"]['fix_sigma'])
        self.bern4_model = Bern4Class(x, mu_gauss, cat, setting["bern4"]['p0'], setting["bern4"]['p1_init'], setting["bern4"]['p2_init'], setting["bern4"]['p3_init'], setting["bern4"]['p4_init'], setting["bern4"]['bond'],setting["bern4"]['sigma_init'],setting["bern4"]['step_init'], setting["bern4"]['fix_sigma'])
        self.bern5_model = Bern5Class(x, mu_gauss, cat, setting["bern5"]['p0'], setting["bern5"]['p1_init'], setting["bern5"]['p2_init'], setting["bern5"]['p3_init'], setting["bern5"]['p4_init'], setting["bern5"]['p5_init'], setting["bern5"]['bond'],setting["bern5"]['sigma_init'],setting["bern5"]['step_init'], setting["bern5"]['fix_sigma'])
        
        self.pow1_model = Pow1Class(x, mu_gauss, cat, setting["pow1"]['sigma_init'], setting["pow1"]['step_init'], setting["pow1"]['p_init'], setting["pow1"]['p_low'], setting["pow1"]['p_high'], setting["pow1"]['di_gauss'], setting["pow1"]['fix_sigma'])
        self.pow2_model = Pow2Class(x, mu_gauss, cat, setting["pow2"]['sigma_init'], setting["pow2"]['step_init'], setting["pow2"]['p1_init'], setting["pow2"]['p1_low'], setting["pow2"]['p1_high'], setting["pow2"]['p2_init'], setting["pow2"]['p2_low'], setting["pow2"]['p2_high'], setting["pow2"]['f1_init'], setting["pow2"]['f2_init'], setting["pow2"]['xmax'], setting["pow2"]['const_f1'], setting["pow2"]['di_gauss'], setting["pow2"]['fix_sigma'])
        self.pow3_model = Pow3Class(x, mu_gauss, cat,setting["pow3"]['sigma_init'], setting["pow3"]['step_init'], setting["pow3"]['p1_init'], setting["pow3"]['p1_low'], setting["pow3"]['p1_high'], setting["pow3"]['p2_init'], setting["pow3"]['p2_low'], setting["pow3"]['p2_high'], setting["pow3"]['p3_init'], setting["pow3"]['p3_low'], setting["pow3"]['p3_high'], setting["pow3"]['f1_init'], setting["pow3"]['f2_init'], setting["pow3"]['f3_init'], setting["pow3"]['xmax'], setting["pow3"]['const_f1'], setting["pow3"]['di_gauss'], setting["pow3"]['fix_sigma'])
        
        self.exp1_model = Exp1Class(x, mu_gauss, cat, setting["exp1"]['sigma_init'], setting["exp1"]['step_init'], setting["exp1"]['p_init'], setting["exp1"]['p_low'], setting["exp1"]['p_high'], setting["exp1"]['di_gauss'], setting["exp1"]['fix_sigma'])
        self.exp2_model = Exp2Class(x, mu_gauss, cat, setting["exp2"]['sigma_init'], setting["exp2"]['step_init'], setting["exp2"]['p1_init'], setting["exp2"]['p1_low'], setting["exp2"]['p1_high'], setting["exp2"]['p2_init'], setting["exp2"]['p2_low'], setting["exp2"]['p2_high'], setting["exp2"]['f1_init'], setting["exp2"]['f2_init'], setting["exp2"]['xmax'], setting["exp2"]['const_f1'], setting["exp2"]['di_gauss'], setting["exp2"]['fix_sigma'])
        self.exp3_model = Exp3Class(x, mu_gauss, cat, setting["exp3"]['sigma_init'], setting["exp3"]['step_init'], setting["exp3"]['p1_init'], setting["exp3"]['p1_low'], setting["exp3"]['p1_high'], setting["exp3"]['p2_init'], setting["exp3"]['p2_low'], setting["exp3"]['p2_high'], setting["exp3"]['p3_init'], setting["exp3"]['p3_low'], setting["exp3"]['p3_high'], setting["exp3"]['f1_init'], setting["exp3"]['f2_init'], setting["exp3"]['f3_init'], setting["exp3"]['xmax'], setting["exp3"]['const_f1'], setting["exp3"]['di_gauss'], setting["exp3"]['fix_sigma'])
        
        self.lau2_model = Lau2Class(x, mu_gauss, cat, setting["lau2"]['sigma_init'], setting["lau2"]['step_init'], setting["lau2"]['p1'], setting["lau2"]['p2'], setting["lau2"]['f_init'], setting["lau2"]['xmax'], setting["lau2"]['const_f1'], setting["lau2"]['di_gauss'], setting["lau2"]['fix_sigma'])
        self.lau3_model = Lau3Class(x, mu_gauss, cat, setting["lau3"]['sigma_init'], setting["lau3"]['step_init'], setting["lau3"]['p1'], setting["lau3"]['p2'], setting["lau3"]['p3'], setting["lau3"]['f_init'], setting["lau3"]['xmax'], setting["lau3"]['const_f1'], setting["lau3"]['di_gauss'], setting["lau3"]['fix_sigma'])
        self.lau4_model = Lau4Class(x, mu_gauss, cat, setting["lau4"]['sigma_init'], setting["lau4"]['step_init'], setting["lau4"]['p1'], setting["lau4"]['p2'], setting["lau4"]['p3'], setting["lau4"]['p4'], setting["lau4"]['f_init'], setting["lau4"]['xmax'], setting["lau4"]['const_f1'], setting["lau4"]['di_gauss'], setting["lau4"]['fix_sigma'])

        self.modg_model = ModGausClass(x, cat, x.getMin(), x.getMax(), setting["modg"]['m0'], setting["modg"]['sl'], setting["modg"]['sh'], setting["modg"]['vl'], setting["modg"]['vr'])

    def testSelection(self, test=""):
        profile_ = []
        if "bern2" in self.setting_[test]: profile_.append(self.bern2_model)
        if "bern3" in self.setting_[test]: profile_.append(self.bern3_model)
        if "bern4" in self.setting_[test]: profile_.append(self.bern4_model)
        if "bern5" in self.setting_[test]: profile_.append(self.bern5_model)
        if "pow1" in self.setting_[test]: profile_.append(self.pow1_model)
        if "pow2" in self.setting_[test]: profile_.append(self.pow2_model)
        if "pow3" in self.setting_[test]: profile_.append(self.pow3_model)
        if "exp1" in self.setting_[test]: profile_.append(self.exp1_model)
        if "exp2" in self.setting_[test]: profile_.append(self.exp2_model)
        if "exp3" in self.setting_[test]: profile_.append(self.exp3_model)
        if "lau2" in self.setting_[test]: profile_.append(self.lau2_model)
        if "lau3" in self.setting_[test]: profile_.append(self.lau3_model)
        if "lau4" in self.setting_[test]: profile_.append(self.lau4_model)
        if "modg" in self.setting_[test]: profile_.append(self.modg_model)
        return profile_
