#!/usr/bin/env python3
import json
import argparse
import re

def configParser(configFile, logFile, CAT, FLAV, YEAR, PROD):
    datString = '_'+CAT+'_'+FLAV+'_'+YEAR+'_'+PROD
    lenDatString = len(datString)
    newLabel = 'Htozg_' + CAT + '_' + FLAV + '_' + YEAR + '_' + PROD + '_' +'nominal'
    jfile = open(configFile, 'r')
    configs = json.load(jfile)
    log = open(logFile)

    disigma = False

    try:
        test = configs[newLabel]
    except KeyError:
        configs[newLabel]={}

    for line in log:
        if "sigmaL" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
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
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            configs[newLabel]["nexp"] = float(r[0])
        elif "dMH" in line:
            r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(datString)+lenDatString:])
            configs[newLabel]["dMH"] = float(r[1])

    if disigma:
        configs[newLabel]["disigma"] = 1
    else:
        configs[newLabel]["disigma"] = 0

    jfile_mod =  open(configFile, 'w')
    modify = json.dump(configs, jfile_mod, indent=4)

    log.close()
    jfile_mod.close()

