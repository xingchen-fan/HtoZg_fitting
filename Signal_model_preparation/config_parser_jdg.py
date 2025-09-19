#!/usr/bin/env python3
import json
import argparse
import re

def configParser(configFile, logFile, CAT, FLAV, YEAR, PROD):
    datString = '_'+CAT+'_'+FLAV+'_'+YEAR+'_'+PROD
    lenDatString = len(datString)
    newLabel = 'Htozg_' + FLAV + '_' + CAT + '_' + YEAR + '_' + PROD + '_' +'nominal'
    jfile = open(configFile, 'r')
    configs = json.load(jfile)

    disigma_el = False
    disigma_mu = False
    log = open(logFile)

    try:
        test = configs[newLabel]
    except KeyError:
        configs[newLabel]={}

    for line in log:
        if "sigmaL" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
        #check if sigmaL in line. Then we will also have el, mu, or neither (comb)
            configs[newLabel]["sigmaL"] = float(r[1])
        elif "sigmaR" in line:
            disigma_el = True
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            configs[newLabel]["sigmaR"] = float(r[1])
        elif "alphaL" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            if len(r) == 1:
                configs[newLabel]["alphaL"] = float(r[0])
            elif len(r) >1:
                configs[newLabel]["alphaL"] = float(r[1])
        elif "alphaR" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            if len(r) == 1:
                configs[newLabel]["alphaR"] = float(r[0])
            elif len(r) >1:
                configs[newLabel]["alphaR"] = float(r[1])
        elif "nL" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            if len(r) == 1:
                configs[newLabel]["nL"] = float(r[0])
            elif len(r) >1:
                configs[newLabel]["nL"] = float(r[1])
        elif "nR" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            if len(r) == 1:
                configs[newLabel]["nR"] = float(r[0])
            elif len(r) >1:
                configs[newLabel]["nR"] = float(r[1])
        elif "nexp" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
            configs[newLabel]["nexp"] = float(r[0])
        elif "MH" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find("_el")+3:])
            configs[newLabel]["MH"] = float(r[1])

    if disigma_el:
        configs[newLabel]["disigma"] = 1
    else:
        configs[newLabel]["disigma"] = 0

    if disigma_mu:
        configs[newLabel]["disigma"] = 1
    else:
        configs[newLabel]["disigma"] = 0
    jfile_mod =  open(configFile, 'w')
    modify = json.dump(configs, jfile_mod, indent=4)

    log.close()
    jfile_mod.close()

