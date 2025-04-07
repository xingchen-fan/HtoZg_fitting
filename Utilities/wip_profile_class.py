import ROOT
import json
from bkg_functions_class import *

class profileClass:
    def __init__(self, x, mu_gauss, cat='', config=''):
        jfile_ = open(config, 'r')
        configs_ = json.load(jfile_)
        setting = configs_[cat]
        self.cat = cat
        self.config_file = config
        self.config_dict = configs_
        self.setting_ = setting
        self.bern2_model = Bern2Class(x, mu_gauss, cat, setting["bern2"]['p0'], setting["bern2"]['p1'], setting["bern2"]['p2'], setting["bern2"]['bond'],setting["bern2"]['sigma'],setting["bern2"]['step'], setting["bern2"]['fix_sigma'])
        self.bern3_model = Bern3Class(x, mu_gauss, cat, setting["bern3"]['p0'], setting["bern3"]['p1'], setting["bern3"]['p2'], setting["bern3"]['p3'], setting["bern3"]['bond'],setting["bern3"]['sigma'],setting["bern3"]['step'], setting["bern3"]['fix_sigma'])
        self.bern4_model = Bern4Class(x, mu_gauss, cat, setting["bern4"]['p0'], setting["bern4"]['p1'], setting["bern4"]['p2'], setting["bern4"]['p3'], setting["bern4"]['p4'], setting["bern4"]['bond'],setting["bern4"]['sigma'],setting["bern4"]['step'], setting["bern4"]['fix_sigma'])
        self.bern5_model = Bern5Class(x, mu_gauss, cat, setting["bern5"]['p0'], setting["bern5"]['p1'], setting["bern5"]['p2'], setting["bern5"]['p3'], setting["bern5"]['p4'], setting["bern5"]['p5'], setting["bern5"]['bond'],setting["bern5"]['sigma'],setting["bern5"]['step'], setting["bern5"]['fix_sigma'])
        
        self.pow1_model = Pow1Class(x, mu_gauss, cat, setting["pow1"]['sigma'], setting["pow1"]['step'], setting["pow1"]['p'], setting["pow1"]['p_low'], setting["pow1"]['p_high'], setting["pow1"]['di_gauss'], setting["pow1"]['fix_sigma'])
        self.pow2_model = Pow2Class(x, mu_gauss, cat, setting["pow2"]['sigma'], setting["pow2"]['step'], setting["pow2"]['p1'], setting["pow2"]['p1_low'], setting["pow2"]['p1_high'], setting["pow2"]['p2'], setting["pow2"]['p2_low'], setting["pow2"]['p2_high'], setting["pow2"]['f1'], setting["pow2"]['f2'], setting["pow2"]['xmax'], setting["pow2"]['const_f1'], setting["pow2"]['di_gauss'], setting["pow2"]['fix_sigma'])
        self.pow3_model = Pow3Class(x, mu_gauss, cat,setting["pow3"]['sigma'], setting["pow3"]['step'], setting["pow3"]['p1'], setting["pow3"]['p1_low'], setting["pow3"]['p1_high'], setting["pow3"]['p2'], setting["pow3"]['p2_low'], setting["pow3"]['p2_high'], setting["pow3"]['p3'], setting["pow3"]['p3_low'], setting["pow3"]['p3_high'], setting["pow3"]['f1'], setting["pow3"]['f2'], setting["pow3"]['f3'], setting["pow3"]['xmax'], setting["pow3"]['const_f1'], setting["pow3"]['di_gauss'], setting["pow3"]['fix_sigma'])
        
        self.exp1_model = Exp1Class(x, mu_gauss, cat, setting["exp1"]['sigma'], setting["exp1"]['step'], setting["exp1"]['p'], setting["exp1"]['p_low'], setting["exp1"]['p_high'], setting["exp1"]['di_gauss'], setting["exp1"]['fix_sigma'])
        self.exp2_model = Exp2Class(x, mu_gauss, cat, setting["exp2"]['sigma'], setting["exp2"]['step'], setting["exp2"]['p1'], setting["exp2"]['p1_low'], setting["exp2"]['p1_high'], setting["exp2"]['p2'], setting["exp2"]['p2_low'], setting["exp2"]['p2_high'], setting["exp2"]['f1'], setting["exp2"]['f2'], setting["exp2"]['xmax'], setting["exp2"]['const_f1'], setting["exp2"]['di_gauss'], setting["exp2"]['fix_sigma'])
        self.exp3_model = Exp3Class(x, mu_gauss, cat, setting["exp3"]['sigma'], setting["exp3"]['step'], setting["exp3"]['p1'], setting["exp3"]['p1_low'], setting["exp3"]['p1_high'], setting["exp3"]['p2'], setting["exp3"]['p2_low'], setting["exp3"]['p2_high'], setting["exp3"]['p3'], setting["exp3"]['p3_low'], setting["exp3"]['p3_high'], setting["exp3"]['f1'], setting["exp3"]['f2'], setting["exp3"]['f3'], setting["exp3"]['xmax'], setting["exp3"]['const_f1'], setting["exp3"]['di_gauss'], setting["exp3"]['fix_sigma'])
        
        self.lau2_model = Lau2Class(x, mu_gauss, cat, setting["lau2"]['sigma'], setting["lau2"]['step'], setting["lau2"]['p1'], setting["lau2"]['p2'], setting["lau2"]['f'], setting["lau2"]['xmax'], setting["lau2"]['const_f1'], setting["lau2"]['di_gauss'], setting["lau2"]['fix_sigma'])
        self.lau3_model = Lau3Class(x, mu_gauss, cat, setting["lau3"]['sigma'], setting["lau3"]['step'], setting["lau3"]['p1'], setting["lau3"]['p2'], setting["lau3"]['p3'], setting["lau3"]['f'], setting["lau3"]['xmax'], setting["lau3"]['const_f1'], setting["lau3"]['di_gauss'], setting["lau3"]['fix_sigma'])
        self.lau4_model = Lau4Class(x, mu_gauss, cat, setting["lau4"]['sigma'], setting["lau4"]['step'], setting["lau4"]['p1'], setting["lau4"]['p2'], setting["lau4"]['p3'], setting["lau4"]['p4'], setting["lau4"]['f'], setting["lau4"]['xmax'], setting["lau4"]['const_f1'], setting["lau4"]['di_gauss'], setting["lau4"]['fix_sigma'])
        self.modg_model = ModGausClass(x, cat, x.getMin(), x.getMax(), setting["modg"]['m0'], setting["modg"]['sl'], setting["modg"]['sh'], setting["modg"]['vl'], setting["modg"]['vr'])
        self.exmg_model = EXMGClass(x, cat, setting["exmg"]['mu'], setting["exmg"]['sig'], setting["exmg"]['xsi'])
        self.agg_model = AGGClass(x, cat, setting["agg"]['kappa'], setting["agg"]['alpha'], setting["agg"]['zeta'])

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
        if "exmg" in self.setting_[test]: profile_.append(self.exmg_model)
        if "agg"  in self.setting_[test]: profile_.append(self.agg_model)
        return profile_

    #Code to update the JSON file
    def write_config_file(self, histogram, test, config_out=None):
      #Update setting dictionary with updated

      test_models = self.testSelection(test)
      
      slimmed_setting = self.setting_
      
      for model in test_models:
        parameters = model.pdf.getParameters(histogram)
        model_name = model.name.replace('_' + self.cat,'')
        #slimmed_setting[model_name] = {}

        for param in parameters:
          #Format name
          param_name = (param.GetName()).replace('_' + self.cat,'') #Check GetName command
          param_name = param_name.replace(model_name + "_",'')
          if param_name=="t":
            param_name = "step"

          #Get value
          param_val = param.getValV()
          #self.config_dict[self.cat][model_name][param_name] = param.getValV()
          #self.setting_[model_name][param_name] = param.getValV()

          #set correct dictionary
          slimmed_setting[model_name][param_name] = param.getValV()

      #Correctly identify which file to write to
      config_to_update = self.config_file
      if config_out:
        config_to_update = config_out

      config_to_dump = {}
      #config_to_dump[self.cat] = self.setting_
      config_to_dump[self.cat] = slimmed_setting

      #Dump the file to a JSON to update CONFIG
      with open(config_to_update, 'w') as output_file:
        json.dump(config_to_dump, output_file, indent=4, separators=(',',': '))


