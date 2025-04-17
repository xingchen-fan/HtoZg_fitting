#!/usr/bin/env python3
import json
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', help = 'Category')
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-log', '--log', help = 'Fitting log')
args = parser.parse_args()
CAT = args.cat
jfile = open('../Config/'+ args.config, 'r')
configs = json.load(jfile)
setting = configs[CAT]
log = open(args.log)
for line in log:
    if "bern2_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern2"]['p1'] = float(r[1])
    elif "bern2_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern2"]['p2'] = float(r[1])
    elif "bern2_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern2"]["sigma"] = float(r[0])
            setting["bern2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern2"]["sigma"] = float(r[1])
            setting["bern2"]["fix_sigma"] = 0
    elif "bern2_step" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern2"]["step"] = float(r[1])

    if "bern3_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]['p1'] = float(r[1])
    elif "bern3_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]['p2'] = float(r[1])
    elif "bern3_p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]['p3'] = float(r[1])
    elif "bern3_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern3"]["sigma"] = float(r[0])
            setting["bern3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern3"]["sigma"] = float(r[1])
            setting["bern3"]["fix_sigma"] = 0
    elif "bern3_step" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]["step"] = float(r[1])

    if "bern4_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p1'] = float(r[1])
    elif "bern4_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p2'] = float(r[1])
    elif "bern4_p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p3'] = float(r[1])
    elif "bern4_p4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p4'] = float(r[1])
    elif "bern4_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern4"]["sigma"] = float(r[0])
            setting["bern4"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern4"]["sigma"] = float(r[1])
            setting["bern4"]["fix_sigma"] = 0
    elif "bern4_step" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]["step"] = float(r[1])

    if "bern5_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p1'] = float(r[1])
    elif "bern5_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p2'] = float(r[1])
    elif "bern5_p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p3'] = float(r[1])
    elif "bern5_p4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p4'] = float(r[1])
    elif "bern5_p5" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p5'] = float(r[1])
    elif "bern5_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern5"]["sigma"] = float(r[0])
            setting["bern5"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern5"]["sigma"] = float(r[1])
            setting["bern5"]["fix_sigma"] = 0
    elif "bern5_step" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]["step"] = float(r[1])
        
    if "pow1_t_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["step"]= float(r[1])
    elif "pow1_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["pow1"]["sigma"] =  float(r[0])
            setting["pow1"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["pow1"]["sigma"] =  float(r[1])
            setting["pow1"]["fix_sigma"] = 0
    elif "pow1_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["sigma2"] =  float(r[1])
        setting["pow1"]["di_gauss"] = 1
    elif "pow1_p_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["p"] = float(r[1])
    elif "pow1_gc" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["gc"] = float(r[1])
    
    if "pow2_t_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["step"]= float(r[1])
    elif "pow2_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["pow2"]["sigma"] =  float(r[0])
            setting["pow2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["pow2"]["sigma"] =  float(r[1])
            setting["pow2"]["fix_sigma"] = 0
    elif "pow2_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["sigma2"] =  float(r[1])
        setting["pow2"]["di_gauss"] = 1
    elif "pow2_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["p1"] = float(r[1])
    elif "pow2_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["p2"] = float(r[1])
    elif "pow2_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["f1"] = float(r[0])
    elif "pow2_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["f2"] = float(r[1])
    elif "pow2_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["gc"] = float(r[1])

    if "pow3_t_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["step"]= float(r[1])
    elif "pow3_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["pow3"]["sigma"] =  float(r[0])
            setting["pow3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["pow3"]["sigma"] =  float(r[1])
            setting["pow3"]["fix_sigma"] = 0
    elif "pow3_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["sigma2"] =  float(r[1])
        setting["pow3"]["di_gauss"] = 1
    elif "pow3_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["p1"] = float(r[1])
    elif "pow3_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["p2"] = float(r[1])
    elif "pow3_p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["p3"] = float(r[1])
    elif "pow3_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["f1"] = float(r[0])
    elif "pow3_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["f2"] = float(r[1])
    elif "pow3_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["f3"] = float(r[1])
    elif "pow3_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["gc"] = float(r[1])


    if "exp1_t_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["step"]= float(r[1])
    elif "exp1_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["exp1"]["sigma"] =  float(r[0])
            setting["exp1"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["exp1"]["sigma"] =  float(r[1])
            setting["exp1"]["fix_sigma"] = 0
    elif "exp1_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["sigma2"] =  float(r[1])
        setting["exp1"]["di_gauss"] = 1
    elif "exp1_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["p"] = float(r[1])
    elif "exp1_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["gc"] = float(r[1])
        
    if "exp2_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["step"]= float(r[1])
    elif "exp2_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["exp2"]["sigma"] =  float(r[0])
            setting["exp2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["exp2"]["sigma"] =  float(r[1])
            setting["exp2"]["fix_sigma"] = 0
    elif "exp2_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["sigma2"] =  float(r[1])
        setting["exp2"]["di_gauss"] = 1
    elif "exp2_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["p1"] = float(r[1])
    elif "exp2_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["p2"] = float(r[1])
    elif "exp2_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["f1"] = float(r[0])
    elif "exp2_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["f2"] = float(r[1])
    elif "exp2_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["gc"] = float(r[1])

    if "exp3_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["step"]= float(r[1])
    elif "exp3_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["exp3"]["sigma"] =  float(r[0])
            setting["exp3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["exp3"]["sigma"] =  float(r[1])
            setting["exp3"]["fix_sigma"] = 0
    elif "exp3_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["sigma2"] =  float(r[1])
        setting["exp3"]["di_gauss"] = 1
    elif "exp3_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["p1"] = float(r[1])
    elif "exp3_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["p2"] = float(r[1])
    elif "exp3_p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["p3"] = float(r[1])
    elif "exp3_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["f1"] = float(r[0])
    elif "exp3_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["f2"] = float(r[1])
    elif "exp3_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["f3"] = float(r[1])
    elif "exp3_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["gc"] = float(r[1])
        
    if "lau2_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["step"]= float(r[1])
    elif "lau2_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["f1"] = float(r[0])
    elif "lau2_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["f2"] = float(r[1])
    elif "lau2_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["lau2"]["sigma"] =  float(r[0])
            setting["lau2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["lau2"]["sigma"] =  float(r[1])
            setting["lau2"]["fix_sigma"] = 0
    elif "lau2_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["sigma2"] =  float(r[1])
        setting["lau2"]["di_gauss"] = 1
    elif "lau2_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["gc"] = float(r[1])
        
    if "lau3_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["step"]= float(r[1])
    elif "lau3_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["f1"] = float(r[0])
    elif "lau3_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["f2"] = float(r[1])
    elif "lau3_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["f3"] = float(r[1])
    elif "lau3_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["lau3"]["sigma"] =  float(r[0])
            setting["lau3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["lau3"]["sigma"] =  float(r[1])
            setting["lau3"]["fix_sigma"] = 0
    elif "lau3_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["sigma2"] =  float(r[1])
        setting["lau3"]["di_gauss"] = 1
    elif "lau3_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["gc"] = float(r[1])
        
    if "lau4_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["step"]= float(r[1])
    elif "lau4_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f1"] = float(r[0])
    elif "lau4_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f2"] = float(r[1])
    elif "lau4_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f3"] = float(r[1])
    elif "lau4_f4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f4"] = float(r[1])
    elif "lau4_sigma_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["lau4"]["sigma"] =  float(r[0])
            setting["lau4"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["lau4"]["sigma"] =  float(r[1])
            setting["lau4"]["fix_sigma"] = 0
    elif "lau4_sigma2_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["sigma2"] =  float(r[1])
        setting["lau4"]["di_gauss"] = 1
    elif "lau4_gc_" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["gc"] = float(r[1])

    if "modg_m0" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["m0"] = float(r[1])
    elif "modg_nuL" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["vl"] = float(r[1])
    elif "modg_nuRange" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["vr"] = float(r[1])
    elif "modg_sigmaL" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["sl"] = float(r[1])
    elif "modg_sigmaH" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["sh"] = float(r[1])


     if "agg_kappa" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["agg"]["kappa"] = float(r[1])
    elif "agg_alpha" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["agg"]["alpha"] = float(r[1])
    elif "agg_zeta" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["agg"]["zeta"] = float(r[1])


    if "exmg_mu" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exmg"]["mu"] = float(r[1])
    elif "exmg_sig" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exmg"]["sig"] = float(r[1])
    elif "exmg_xsi" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exmg"]["xsi"] = float(r[1])

jfile_mod =  open('../Config/'+args.config, 'w')
modify = json.dump(configs, jfile_mod, indent=4)
#f = open('chi2.log')
