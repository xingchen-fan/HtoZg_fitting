#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

x = [1.557, 1.53, 1.427, 2.266, 2.502, 2.159, 1.122, -0.528, 5.25]
y = [1, 3, 4, 5, 6, 8, 9, 10, 11]
xerr1_ = [[0.606, 1.717, 1.296, 1.535, 3.484, 1.335, 0.993, 1.432, 3.489], [0.547, 2.213, 1.237, 1.563, 3.414, 1.628, 1.065, 1.497, 3.489]]
xerr2_ = [[2*x for x in xerr1_[0]], [2*x for x in xerr1_[1]]]

xlow = -8
xhigh = 12
ylow = 0
yhigh = 12
red_patch = mpatches.Patch(color='red', label='$\pm$1 $\sigma$')
blue_patch = mpatches.Patch(color='blue', label='$\pm$2 $\sigma$')

plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)

plt.errorbar(x, y,xerr=xerr2_, fmt='sk',ecolor='b',elinewidth=2,markersize = 0.3)
plt.errorbar(x, y,xerr=xerr1_, fmt='sk',ecolor='r',elinewidth=4,markersize = 4)
plt.legend(handles=[red_patch, blue_patch], loc='upper left', prop={'size': 12})

for i in range(len(x)):
  #plt.annotate(str(x[i])+'$\pm$'+str(xerr1_[i]), xy=(x[i], y[i]), xytext=(8, y[i]), size = 12)
  plt.annotate('%.3f'%x[i], xy=(x[i], y[i]), xytext=(8, y[i]), size = 12)
  plt.annotate('-' + str(xerr1_[0][i]), xy=(x[i], y[i]), xytext=(10.2, y[i]-0.1), size = 6)
  plt.annotate('+' + str(xerr1_[1][i]), xy=(x[i], y[i]), xytext=(10.2, y[i]+0.3), size = 6)

plt.yticks(y, ['Combine', 'ggF1', 'ggF2', 'ggF3', 'ggF4', 'VBF1', 'VBF2', 'VBF3', 'VBF4'], size=12)
plt.tick_params(axis="x",direction="in")
plt.tick_params(axis="y",direction="in")
plt.xticks(np.arange(xlow, xhigh, step=1))

plt.vlines(1, ylow, yhigh, colors='k').set_linewidth(0.5)
plt.vlines(x[0], ylow, yhigh, colors='k', linestyles ="dashed", dashes=(0, (5, 10))).set_linewidth(0.5)
plt.annotate('SM', xy=(0, 0.2),size=12)

plt.title('CMS', loc='left', fontname="sans serif", fontstyle='italic', fontsize=18)
plt.title('Toy Simulation', x=0.22,y=1, fontname="sans serif", fontstyle='italic',fontsize=11)
plt.title('137.61 fb$^{-1}$ (13 TeV) + 62.32 fb$^{-1}$ (13.6 TeV)', loc = 'right', fontstyle='normal', fontname="sans serif", fontsize=10)
plt.xlabel('$\mu$')
plt.savefig('toy_signal.png')
