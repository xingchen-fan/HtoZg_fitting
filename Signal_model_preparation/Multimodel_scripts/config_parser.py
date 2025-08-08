#!/usr/bin/env python3
import json
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-p', '--prod', help = 'Production mode')
parser.add_argument('-y', '--year', help = 'Year')
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-log', '--log', help = 'Fitting log')
args = parser.parse_args()
CAT = args.cat
PROD = args.prod
jfile = open(args.config, 'r')
configs = json.load(jfile)
model_name_el = "Htozg_el_"+CAT+"_"+args.year+"_"+PROD
model_name_mu = "Htozg_mu_"+CAT+"_"+args.year+"_"+PROD
setting_el = configs[model_name_el]
setting_mu = configs[model_name_mu]
log = open(args.log)
disigma_el = False
disigma_mu = False

for line in log:
    if "sigmaL" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting_el["sigmaL"] = float(r[1])
    elif "sigmaR" in line and "_el" in line:
        disigma_el = True
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting_el["sigmaR"] = float(r[1])
    elif "sigma_eff_err" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting_el["sigma_eff_err " + PROD] = float(r[0])
    elif "alphaL" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting_el["alphaL"] = float(r[0])
        elif len(r) >1:
            setting_el["alphaL"] = float(r[1])
    elif "alphaR" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting_el["alphaR"] = float(r[0])
        elif len(r) >1:
            setting_el["alphaR"] = float(r[1])
    elif "nL" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting_el["nL"] = float(r[0])
        elif len(r) >1:
            setting_el["nL"] = float(r[1])
    elif "nR" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting_el["nR"] = float(r[0])
        elif len(r) >1:
            setting_el["nR"] = float(r[1])
    elif "nexp" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting_el["nexp"] = float(r[0])
    elif "dMH" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting_el["dMH"] = float(r[1])
        
    elif "sigmaL" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting_mu["sigmaL"] = float(r[1])
    elif "sigmaR" in line and "_mu" in line:
        disigma_mu = True
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting_mu["sigmaR"] = float(r[1])
    elif "sigma_eff_err" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting_mu["sigma_eff_err " + PROD] = float(r[0])
    elif "alphaL" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting_mu["alphaL"] = float(r[0])
        elif len(r) >1:
            setting_mu["alphaL"] = float(r[1])
    elif "alphaR" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting_mu["alphaR"] = float(r[0])
        elif len(r) >1:
            setting_mu["alphaR"] = float(r[1])
    elif "nL" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting_mu["nL"] = float(r[0])
        elif len(r) >1:
            setting_mu["nL"] = float(r[1])
    elif "nR" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting_mu["nR"] = float(r[0])
        elif len(r) >1:
            setting_mu["nR"] = float(r[1])
    elif "nexp" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting_mu["nexp"] = float(r[0])
    elif "dMH" in line and "_mu" in line: 
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting_mu["dMH"] = float(r[1])
if disigma_el:
    setting_el["disigma"] = 1
else:
    setting_el["disigma"] = 0

if disigma_mu:
    setting_mu["disigma"] = 1
else:
    setting_mu["disigma"] = 0
jfile_mod =  open(args.config, 'w')
modify = json.dump(configs, jfile_mod, indent=4)
#f = open('chi2.log')
