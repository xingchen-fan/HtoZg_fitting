#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import math

x_el = [4.545, 1.114, 3.070, -0.237, 2.099, -1.590, 4.181]
x_mu = [0.082, 2.488, 1.707, -0.585, 0.510, 1.157, 7.949]
y = [1, 2, 3, 4, 6, 7, 8]
xerrud_el = [[2.659, 2.009, 2.439, 5.546, 1.504, 2.158, 5.802], [2.877, 2.037, 2.463, 5.545, 1.713, 2.291, 5.918]]
xerrud_mu = [[2.002, 1.604, 1.975, 3.825, 0.953, 1.903, 6.424], [2.178, 1.604, 1.992, 9.239, 1.385, 2.022, 4,616]]
x = [x_el[i] - x_mu[i] for i in range(len(x_el))]

xerr1_el = [(xerrud_el[0][i] + xerrud_el[1][i])/2 for i in range(len(x_el))]
xerr1_mu = [(xerrud_mu[0][i] + xerrud_mu[1][i])/2 for i in range(len(x_mu))]

xerr1_ = [math.sqrt(xerr1_el[i]*xerr1_el[i] + xerr1_mu[i]*xerr1_mu[i]) for i in range(len(x_el))]
xerr2_ = [2*x for x in xerr1_]

xlow = -10
xhigh = 10
ylow = 0
yhigh = 10
red_patch = mpatches.Patch(color='red', label='$\pm$1 $\sigma$')
blue_patch = mpatches.Patch(color='blue', label='$\pm$2 $\sigma$')

plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)

plt.errorbar(x, y,xerr=xerr2_, fmt='sk',ecolor='b',elinewidth=2,markersize = 0.3)
plt.errorbar(x, y,xerr=xerr1_, fmt='sk',ecolor='r',elinewidth=4,markersize = 4)
plt.legend(handles=[red_patch, blue_patch], loc='upper left', prop={'size': 12})

for i in range(len(x)):
  plt.annotate('%.3f'%x[i]+'$\pm$'+ '%.3f'%xerr1_[i], xy=(x[i], y[i]), xytext=(6, y[i]), size = 10)

plt.yticks(y, ['ggF1', 'ggF2', 'ggF3', 'ggF4', 'VBF1&2', 'VBF3', 'VBF4'], size=12)
plt.tick_params(axis="x",direction="in")
plt.tick_params(axis="y",direction="in")
plt.xticks(np.arange(xlow, xhigh, step=1))

plt.vlines(0, ylow, yhigh, colors='k').set_linewidth(0.5)
#plt.vlines(x[0], ylow, yhigh, colors='k', linestyles ="dashed", dashes=(0, (5, 10))).set_linewidth(0.5)
#plt.annotate('SM', xy=(0, 0.2),size=12)

plt.title('CMS', loc='left', fontname="sans serif", fontstyle='italic', fontsize=18)
plt.title('Toy Simulation', x=0.22,y=1, fontname="sans serif", fontstyle='italic',fontsize=11)
plt.title('137.61 fb$^{-1}$ (13 TeV) + 62.32 fb$^{-1}$ (13.6 TeV)', loc = 'right', fontstyle='normal', fontname="sans serif", fontsize=10)
plt.xlabel('$\Delta\mu$ between electron and muon channels')
plt.savefig('example.png')
