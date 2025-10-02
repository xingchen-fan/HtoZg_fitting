"""Simple utility to validate signal fits are reasonable and output some summary info
"""
#!/usr/bin/env python3
from argparse import ArgumentParser
import json
from os.path import isfile, abspath
import sys
sys.path.append(abspath("../Utilities/"))
from constants import *

def get_args():
    """Parses arguments

    Returns:
        Namespace with arguments
    """
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', help = 'Config file', default='')
    return parser.parse_args()

if __name__=="__main__":
    args = get_args()
    config = {}
    if isfile(args.config):
        with open(args.config, "r") as config_file:
            config = json.load(config_file)
    else:
        raise RuntimeError("Invalid json")
    procs = ["Htozg_el","Htozg_mu","Htomm"]
    warnings = []
    for proc in procs:
        for cat in CATEGORIES:
            dmh = config[f"{proc}_cat_{cat}_nominal"]["dMH"]
            sigmal = config[f"{proc}_cat_{cat}_nominal"]["sigmaL"]
            sigmar = config[f"{proc}_cat_{cat}_nominal"]["sigmaR"]
            dmh_eup = config[f"{proc}_cat_{cat}_CMS_scale_eUp"]["dMH"]-dmh
            dmh_edn = config[f"{proc}_cat_{cat}_CMS_scale_eDown"]["dMH"]-dmh
            dmh_gup = config[f"{proc}_cat_{cat}_CMS_scale_gUp"]["dMH"]-dmh
            dmh_gdn = config[f"{proc}_cat_{cat}_CMS_scale_gDown"]["dMH"]-dmh
            sigmal_eup = (config[f"{proc}_cat_{cat}_CMS_res_eUp"]["sigmaL"]
                          -sigmal)
            sigmar_eup = (config[f"{proc}_cat_{cat}_CMS_res_eUp"]["sigmaR"]
                          -sigmar)
            sigmal_edn = (config[f"{proc}_cat_{cat}_CMS_res_eDown"]["sigmaL"]
                          -sigmal)
            sigmar_edn = (config[f"{proc}_cat_{cat}_CMS_res_eDown"]["sigmaR"]
                          -sigmar)
            sigmal_gup = (config[f"{proc}_cat_{cat}_CMS_res_gUp"]["sigmaL"]
                          -sigmal)
            sigmar_gup = (config[f"{proc}_cat_{cat}_CMS_res_gUp"]["sigmaR"]
                          -sigmar)
            sigmal_gdn = (config[f"{proc}_cat_{cat}_CMS_res_gDown"]["sigmaL"]
                          -sigmal)
            sigmar_gdn = (config[f"{proc}_cat_{cat}_CMS_res_gDown"]["sigmaR"]
                          -sigmar)
            print(f"process {proc} in category {cat}:")
            print(f"dmh = {dmh:.2f}+{dmh_eup:.2f}(eup)+{dmh_edn:.2f}(edn)"
                  +f"+{dmh_gup:.2f}(gup)+{dmh_gdn:.2f}(gdn)")
            print(f"sigmal = {sigmal:.2f}+{sigmal_eup:.2f}(eup)"
                  +f"+{sigmal_edn:.2f}(edn)"
                  +f"+{sigmal_gup:.2f}(gup)+{sigmal_gdn:.2f}(gdn)")
            print(f"sigmar = {sigmar:.2f}+{sigmar_eup:.2f}(eup)"
                  +f"+{sigmar_edn:.2f}(edn)"
                  +f"+{sigmar_gup:.2f}(gup)+{sigmar_gdn:.2f}(gdn)")
            if (dmh < -1.0 or dmh > 1.0):
                warnings.append(f"MH<124 or >126 for {proc} in {cat}")
            if ((sigmal < 0.5 or sigmal > 3.0) and not proc=="Htomm"):
                warnings.append(f"sigmal<0.5 or >3.0 for {proc} in {cat}")
            if ((sigmar < 0.5 or sigmar > 3.0) and not proc=="Htomm"):
                warnings.append(f"sigmar<0.5 or >3.0 for {proc} in {cat}")
            if (abs(dmh_eup)>0.4 or abs(dmh_edn)>0.4 or abs(dmh_gup)>0.4 
                or abs(dmh_gdn)>0.4):
                warnings.append(f"dmh unc>0.4 for {proc} in {cat}")
            if (abs(sigmal_eup)/sigmal>0.5 or abs(sigmal_edn)/sigmal>0.5 
                or abs(sigmal_gup)/sigmal>0.5 or abs(sigmal_gdn)/sigmal>0.5):
                warnings.append(f"sigmal rel unc>0.5 for {proc} in {cat}")
            if (abs(sigmar_eup)/sigmar>0.5 or abs(sigmar_edn)/sigmar>0.5 
                or abs(sigmar_gup)/sigmar>0.5 or abs(sigmar_gdn)/sigmar>0.5):
                warnings.append(f"sigmar rel unc>0.5 for {proc} in {cat}")
                
    if len(warnings)==0:
        print("All fits within specifications.")
    else:
        print("WARNINGS (may be okay, but user should verify):")
        for warning in warnings:
            print(warning)
        
