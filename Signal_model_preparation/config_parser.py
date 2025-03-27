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
jfile = open('../Config/'+ args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT][args.year]
log = open(args.log)
disigma_el = False
disigma_mu = False
for line in log:
    if "sigmaL" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting["el"]["sigmaL "+PROD] = float(r[1])
    elif "sigmaR" in line and "_el" in line:
        disigma_el = True
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting["el"]["sigmaR "+PROD] = float(r[1])
    elif "alphaL" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting["el"]["alphaL "+PROD] = float(r[0])
        elif len(r) >1:
            setting["el"]["alphaL "+PROD] = float(r[1])
    elif "alphaR" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting["el"]["alphaR "+PROD] = float(r[0])
        elif len(r) >1:
            setting["el"]["alphaR "+PROD] = float(r[1])
    elif "nL" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting["el"]["nL "+PROD] = float(r[0])
        elif len(r) >1:
            setting["el"]["nL "+PROD] = float(r[1])
    elif "nR" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        if len(r) == 1:
            setting["el"]["nR "+PROD] = float(r[0])
        elif len(r) >1:
            setting["el"]["nR "+PROD] = float(r[1])
    elif "nexp" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting["el"]["nexp "+PROD] = float(r[0])
    elif "MH" in line and "_el" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
        setting["el"]["MH "+PROD] = float(r[1])
        
    elif "sigmaL" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting["mu"]["sigmaL "+PROD] = float(r[1])
    elif "sigmaR" in line and "_mu" in line:
        disigma_mu = True
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting["mu"]["sigmaR "+PROD] = float(r[1])
    elif "alphaL" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting["mu"]["alphaL "+PROD] = float(r[0])
        elif len(r) >1:
            setting["mu"]["alphaL "+PROD] = float(r[1])
    elif "alphaR" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting["mu"]["alphaR "+PROD] = float(r[0])
        elif len(r) >1:
            setting["mu"]["alphaR "+PROD] = float(r[1])
    elif "nL" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting["mu"]["nL "+PROD] = float(r[0])
        elif len(r) >1:
            setting["mu"]["nL "+PROD] = float(r[1])
    elif "nR" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        if len(r) == 1:
            setting["mu"]["nR "+PROD] = float(r[0])
        elif len(r) >1:
            setting["mu"]["nR "+PROD] = float(r[1])
    elif "nexp" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting["mu"]["nexp "+PROD] = float(r[0])
    elif "MH" in line and "_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_mu")+3:])
        setting["mu"]["MH "+PROD] = float(r[1])
if disigma_el:
    setting["el"]["disigma "+PROD] = 1
else:
    setting["el"]["disigma "+PROD] = 0

if disigma_mu:
    setting["mu"]["disigma "+PROD] = 1
else:
    setting["mu"]["disigma "+PROD] = 0
jfile_mod =  open('../Config/'+args.config, 'w')
modify = json.dump(configs, jfile_mod, indent=4)
#f = open('chi2.log')
