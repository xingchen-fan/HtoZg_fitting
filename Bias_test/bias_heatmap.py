import numpy as np 
import matplotlib.pyplot as plt
import argparse
import json

parser = argparse.ArgumentParser(description = "Bias test heat map")
parser.add_argument('-c', '--cat', help="category")
args = parser.parse_args()
jfile = open('../Config/config.json', 'r')
configs = json.load(jfile)
CAT = args.cat
setting = configs[CAT]
func_list = setting["FT"]
N_func = len(func_list)

def fill_best(best_list_, func):
    for i, entry in enumerate(func_list):
        if func == entry + "_" + CAT +"_model": best_list_[N_func - 1 - i]+=1 

def read_dir(func, N_func):
    best_list = [0]* N_func
    for i in range(200):
        f = open("condor/"+ func + "_" +CAT +"/output"+str(i)+".txt")
        for line in f:
            best = ""
            if line[:9] == "best func":
                #print("the first 9 = ", line[:9])
                brac = line.find("[")
                comma = line.find(",")
                #print("brac = ", brac)
                #print("comma = ", comma)
                best = line[brac+2:comma-1]
                fill_best(best_list, best)
                theEnd = False
                while not theEnd:
                    left = comma
                    comma = line.find(",", left+1)
                    best = line[left+3: comma-1]
                    fill_best(best_list, best)
                    if line.find(",", comma+1) == -1:
                        theEnd = True
                best = line[comma+3: line.find("]")-1]
                fill_best(best_list, best)
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
plt.title(CAT + " Bias Test")
plt.colorbar()
for i in range(N_func): 
    for j in range(N_func): 
        plt.annotate(format(arr[i][j], '.2f'), xy=(j+0.05, i+0.05), ha='center', va='center', color='white')

plt.xlabel("Fit") 
plt.ylabel("Truth") 
plt.xticks(range(N_func), inverse_list, rotation=90)
plt.yticks(range(N_func), func_list)
plt.savefig("plots/heatmap_" + CAT + ".pdf", format="pdf", bbox_inches="tight")
plt.show()
