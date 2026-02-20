#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import matplotlib.patches as mpatches
import math

x_el = [-1.33, -1.00, 3.97, -0.13, 0.99, -4.97]
x_mu = [-0.64, -1.66, 2.3, 1.15, 1.31, -1.59]
y = [13, 11, 9, 6, 4, 2]
y = [3*entry for entry in y]

xerrud_el = [[1.90, 2.28, 3.26, 0.97, 2.58, 6.48], [2.02, 2.23, 3.29, 1.12, 2.66, 6.54]]
xerrud_mu = [[1.37, 1.44, 2.06, 1.01, 1.83, 5.17], [1.47, 1.50, 3.18, 1.15, 2.48, 5.24]]
x = [x_el[i] - x_mu[i] for i in range(len(x_el))]

xerr1_el = [(xerrud_el[0][i] + xerrud_el[1][i])/2 for i in range(len(x_el))]
xerr1_mu = [(xerrud_mu[0][i] + xerrud_mu[1][i])/2 for i in range(len(x_mu))]

xerr1_ = [math.sqrt(xerr1_el[i]*xerr1_el[i] + xerr1_mu[i]*xerr1_mu[i]) for i in range(len(x_el))]
xerr2_ = [2*x for x in xerr1_]

xlow = -10
xhigh = 10
ylow = 0
yhigh = 45
red_patch = mpatches.Patch(color='red', label='$\pm$1 $\sigma$')
blue_patch = mpatches.Patch(color='blue', label='$\pm$2 $\sigma$')

plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)

plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.cms.text("Toy Simulation with Systematics", fontsize=12)

plt.errorbar(x, y,xerr=xerr2_, fmt='sk',ecolor='b',elinewidth=2,markersize = 0.3)
plt.errorbar(x, y,xerr=xerr1_, fmt='sk',ecolor='r',elinewidth=4,markersize = 4)
plt.legend(handles=[red_patch, blue_patch], loc='upper left', prop={'size': 12})

for i in range(len(x)):
  plt.annotate('%.2f'%x[i]+'$\pm$'+ '%.2f'%xerr1_[i], xy=(x[i], y[i]), xytext=(6, y[i]+1), size = 12)

plt.yticks(y, ['ggF1&2', 'ggF3', 'ggF4', 'VBF1&2', 'VBF3', 'VBF4'], size=12)
plt.tick_params(axis="x",direction="in")
plt.tick_params(axis="y",direction="in")
plt.xticks(np.arange(xlow, xhigh, step=1))

plt.vlines(0, ylow, yhigh, colors='k').set_linewidth(0.5)
#plt.vlines(x[0], ylow, yhigh, colors='k', linestyles ="dashed", dashes=(0, (5, 10))).set_linewidth(0.5)
#plt.annotate('SM', xy=(0, 0.2),size=12)

#plt.title('CMS', loc='left', fontname="sans serif", fontstyle='italic', fontsize=18)
#plt.title('Toy Simulation', x=0.22,y=1, fontname="sans serif", fontstyle='italic',fontsize=11)
#plt.title('137.61 fb$^{-1}$ (13 TeV) + 62.32 fb$^{-1}$ (13.6 TeV)', loc = 'right', fontstyle='normal', fontname="sans serif", fontsize=10)
plt.xlabel('$\Delta\mu$ between electron and muon channels')
plt.savefig('toy_lep_comp.png')
