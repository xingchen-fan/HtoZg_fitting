#!/usr/bin/env python3
import numpy as np 
import matplotlib.pyplot as plt
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description = "Bias test heat map")
parser.add_argument('-c', '--cat', help="category")
parser.add_argument('-con', '--config', help = 'Configuration')
parser.add_argument('-s', '--sig', help = 'Signal injected', default='0')
parser.add_argument('-b', '--bad', help = 'Bad toy heat map?', type=bool, default=False)

args = parser.parse_args()
jfile = open(args.config, 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
func_list = setting["CMSBias"]
N_func = len(func_list)
if args.bad:
    keyword = "bad funct"
else:
    keyword = "good func"

def fill_best(best_list_, func, file_):
    found = False
    for i, entry in enumerate(func_list):
        if func == entry + "_" + CAT +"_model":
            best_list_[N_func - 1 - i]+=1
            found = True
    if not found:
        print('Not found', func, file_)
        

def read_dir(func, N_func):
    best_list = [0]* N_func
    for i in range(200):
        #print ('file ',func, ' ',i)
        file_path = Path("condor/"+ func + "_" +CAT +"_"+args.sig+"sig/output"+str(i)+".txt")
        if not file_path.is_file():
            continue
        f =  open("condor/"+ func + "_" +CAT +"_"+args.sig+"sig/output"+str(i)+".txt")
        for line in f:
            best = ""
            if line[:9] == keyword:
                #print("the first 9 = ", line[:9])
                brac = line.find("[")
                comma = line.find(",")
                if comma == -1:
                    best = line[brac+2:line.find("]")-1]
                    fill_best(best_list, best, i)
                    #print (best)
                    continue
                #print("brac = ", brac)
                #print("comma = ", comma)
                best = line[brac+2:comma-1]
                #print (best)
                fill_best(best_list, best, i)
                theEnd = False
                while not theEnd:
                    left = comma
                    comma = line.find(",", left+1)
                    if comma == -1:
                        best = line[left+3: line.find("]")-1]
                        #print (best)
                        fill_best(best_list, best, i)
                        theEnd = True
                    else:
                        best = line[left+3: comma-1]
                        #print (best)
                        fill_best(best_list, best, i)

        f.close()
    return best_list


inverse_list = [0]*N_func
for i in range(N_func):
    inverse_list[i] = func_list[N_func-1-i]
arr = [[0]*N_func]*N_func
Ntoys = [0]*N_func
for j in range(len(func_list)):
    arr[j] = read_dir(func_list[j], N_func)
    Ntoys[j] = sum(arr[j])
print(arr)
print(Ntoys)
for i in range(N_func):
  for j in range(N_func):
      arr[i][j] /= Ntoys[i]

plt.imshow(arr)
plt.title(CAT + " Bias Test "+ args.sig + " Signal")
plt.colorbar()
for i in range(N_func): 
    for j in range(N_func): 
        plt.annotate(format(arr[i][j], '.2f'), xy=(j+0.05, i+0.05), ha='center', va='center', color='white')

plt.xlabel("Fit") 
plt.ylabel("Truth") 
plt.xticks(range(N_func), inverse_list, rotation=90)
plt.yticks(range(N_func), func_list)
if args.sig == '0' and not args.bad:
    plt.savefig("plots/heatmap_" + CAT + ".pdf", format="pdf", bbox_inches="tight")
elif args.sig == '0' and args.bad:
    plt.savefig("plots/bad_heatmap_" + CAT + ".pdf", format="pdf", bbox_inches="tight")
elif args.sig != '0' and not args.bad:
    plt.savefig("plots/heatmap_" + CAT + "_" + args.sig+"sig.pdf", format="pdf", bbox_inches="tight")
elif args.sig != '0' and args.bad:
    plt.savefig("plots/bad_heatmap_" + CAT + "_" + args.sig+"sig.pdf", format="pdf", bbox_inches="tight")
plt.show()

