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
    if "b2p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern2"]['p1_init'] = float(r[1])
    elif "b2p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern2"]['p2_init'] = float(r[1])
    elif "sigma_bern2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern2"]["sigma_init"] = float(r[0])
            setting["bern2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern2"]["sigma_init"] = float(r[1])
            setting["bern2"]["fix_sigma"] = 0
    elif "stepval_bern2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern2"]["step_init"] = float(r[1])

    if "b3p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]['p1_init'] = float(r[1])
    elif "b3p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]['p2_init'] = float(r[1])
    elif "b3p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]['p3_init'] = float(r[1])
    elif "sigma_bern3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern3"]["sigma_init"] = float(r[0])
            setting["bern3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern3"]["sigma_init"] = float(r[1])
            setting["bern3"]["fix_sigma"] = 0
    elif "stepval_bern3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern3"]["step_init"] = float(r[1])

    if "b4p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p1_init'] = float(r[1])
    elif "b4p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p2_init'] = float(r[1])
    elif "b4p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p3_init'] = float(r[1])
    elif "b4p4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]['p4_init'] = float(r[1])
    elif "sigma_bern4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern4"]["sigma_init"] = float(r[0])
            setting["bern4"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern4"]["sigma_init"] = float(r[1])
            setting["bern4"]["fix_sigma"] = 0
    elif "stepval_bern4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern4"]["step_init"] = float(r[1])

    if "b5p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p1_init'] = float(r[1])
    elif "b5p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p2_init'] = float(r[1])
    elif "b5p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p3_init'] = float(r[1])
    elif "b5p4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p4_init'] = float(r[1])
    elif "b5p5" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]['p5_init'] = float(r[1])
    elif "sigma_bern5" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["bern5"]["sigma_init"] = float(r[0])
            setting["bern5"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["bern5"]["sigma_init"] = float(r[1])
            setting["bern5"]["fix_sigma"] = 0
    elif "stepval_bern5" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["bern5"]["step_init"] = float(r[1])
        
    if "pow1t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["step_init"]= float(r[1])
    elif "sigma_pow1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["pow1"]["sigma_init"] =  float(r[0])
            setting["pow1"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["pow1"]["sigma_init"] =  float(r[1])
            setting["pow1"]["fix_sigma"] = 0
    elif "sigma2_pow1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["sigma2_init"] =  float(r[1])
        setting["pow1"]["di_gauss"] = 1
    elif "pow1p" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["p_init"] = float(r[1])
    elif "gc_pow1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow1"]["gc_init"] = float(r[1])
    
    if "pow2t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["step_init"]= float(r[1])
    elif "sigma_pow2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["pow2"]["sigma_init"] =  float(r[0])
            setting["pow2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["pow2"]["sigma_init"] =  float(r[1])
            setting["pow2"]["fix_sigma"] = 0
    elif "sigma2_pow2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["sigma2_init"] =  float(r[1])
        setting["pow2"]["di_gauss"] = 1
    elif "pow2p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["p1_init"] = float(r[1])
    elif "pow2p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["p2_init"] = float(r[1])
    elif "pow2f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["f1_init"] = float(r[0])
    elif "pow2f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["f2_init"] = float(r[1])
    elif "gc_pow2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow2"]["gc_init"] = float(r[1])

    if "pow3t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["step_init"]= float(r[1])
    elif "sigma_pow3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["pow3"]["sigma_init"] =  float(r[0])
            setting["pow3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["pow3"]["sigma_init"] =  float(r[1])
            setting["pow3"]["fix_sigma"] = 0
    elif "sigma2_pow3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["sigma2_init"] =  float(r[1])
        setting["pow3"]["di_gauss"] = 1
    elif "pow3p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["p1_init"] = float(r[1])
    elif "pow3p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["p2_init"] = float(r[1])
    elif "pow3p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["p3_init"] = float(r[1])
    elif "pow3f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["f1_init"] = float(r[0])
    elif "pow3f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["f2_init"] = float(r[1])
    elif "pow3f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["f3_init"] = float(r[1])
    elif "gc_pow3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["pow3"]["gc_init"] = float(r[1])


    if "exp_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["step_init"]= float(r[1])
    elif "sigma_exp1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["exp1"]["sigma_init"] =  float(r[0])
            setting["exp1"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["exp1"]["sigma_init"] =  float(r[1])
            setting["exp1"]["fix_sigma"] = 0
    elif "sigma2_exp1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["sigma2_init"] =  float(r[1])
        setting["exp1"]["di_gauss"] = 1
    elif "exp_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["p_init"] = float(r[1])
    elif "gc_exp1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp1"]["gc_init"] = float(r[1])
        
    if "exp2_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["step_init"]= float(r[1])
    elif "sigma_exp2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["exp2"]["sigma_init"] =  float(r[0])
            setting["exp2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["exp2"]["sigma_init"] =  float(r[1])
            setting["exp2"]["fix_sigma"] = 0
    elif "sigma2_exp2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["sigma2_init"] =  float(r[1])
        setting["exp2"]["di_gauss"] = 1
    elif "exp2_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["p1_init"] = float(r[1])
    elif "exp2_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["p2_init"] = float(r[1])
    elif "exp2_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["f1_init"] = float(r[0])
    elif "exp2_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["f2_init"] = float(r[1])
    elif "gc_exp2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp2"]["gc_init"] = float(r[1])

    if "exp3_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["step_init"]= float(r[1])
    elif "sigma_exp3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["exp3"]["sigma_init"] =  float(r[0])
            setting["exp3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["exp3"]["sigma_init"] =  float(r[1])
            setting["exp3"]["fix_sigma"] = 0
    elif "sigma2_exp3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["sigma2_init"] =  float(r[1])
        setting["exp3"]["di_gauss"] = 1
    elif "exp3_p1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["p1_init"] = float(r[1])
    elif "exp3_p2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["p2_init"] = float(r[1])
    elif "exp3_p3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["p3_init"] = float(r[1])
    elif "exp3_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["f1_init"] = float(r[0])
    elif "exp3_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["f2_init"] = float(r[1])
    elif "exp3_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["f3_init"] = float(r[1])
    elif "gc_exp3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["exp3"]["gc_init"] = float(r[1])
        
    if "lau2_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["step_init"]= float(r[1])
    elif "lau2_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["f1_init"] = float(r[0])
    elif "lau2_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["f2_init"] = float(r[1])
    elif "sigma_lau2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["lau2"]["sigma_init"] =  float(r[0])
            setting["lau2"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["lau2"]["sigma_init"] =  float(r[1])
            setting["lau2"]["fix_sigma"] = 0
    elif "sigma2_lau2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["sigma2_init"] =  float(r[1])
        setting["lau2"]["di_gauss"] = 1
    elif "gc_lau2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau2"]["gc_init"] = float(r[1])
        
    if "lau3_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["step_init"]= float(r[1])
    elif "lau3_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["f1_init"] = float(r[0])
    elif "lau3_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["f2_init"] = float(r[1])
    elif "lau3_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["f3_init"] = float(r[1])
    elif "sigma_lau3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["lau3"]["sigma_init"] =  float(r[0])
            setting["lau3"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["lau3"]["sigma_init"] =  float(r[1])
            setting["lau3"]["fix_sigma"] = 0
    elif "sigma2_lau3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["sigma2_init"] =  float(r[1])
        setting["lau3"]["di_gauss"] = 1
    elif "gc_lau3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau3"]["gc_init"] = float(r[1])
        
    if "lau4_t" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["step_init"]= float(r[1])
    elif "lau4_f1" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f1_init"] = float(r[0])
    elif "lau4_f2" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f2_init"] = float(r[1])
    elif "lau4_f3" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f3_init"] = float(r[1])
    elif "lau4_f4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["f4_init"] = float(r[1])
    elif "sigma_lau4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        if len(r) == 1:
            setting["lau4"]["sigma_init"] =  float(r[0])
            setting["lau4"]["fix_sigma"] = 1
        elif len(r) >1:
            setting["lau4"]["sigma_init"] =  float(r[1])
            setting["lau4"]["fix_sigma"] = 0
    elif "sigma2_lau4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["sigma2_init"] =  float(r[1])
        setting["lau4"]["di_gauss"] = 1
    elif "gc_lau4" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["lau4"]["gc_init"] = float(r[1])

    if "m0" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["m0"] = float(r[1])
    elif "nuL" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["vl"] = float(r[1])
    elif "nuRange" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["vr"] = float(r[1])
    elif "sigmaL" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["sl"] = float(r[1])
    elif "sigmaH" in line:
        r = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",line[line.find(CAT)+4:])
        setting["modg"]["sh"] = float(r[1])



jfile_mod =  open('../Config/'+args.config, 'w')
modify = json.dump(configs, jfile_mod, indent=4)
#f = open('chi2.log')
