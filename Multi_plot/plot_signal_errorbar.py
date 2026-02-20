#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import mplhep as hep

x = [1.26, -1.24, -1.51, 2.89, 2.16, 0.87, 1.22, -2.56, 5.84, -16.59, 2.22, 6.78, -9.06, 0.73]
y = [30, 28, 26, 24, 21, 19, 17, 15, 12, 10, 8, 6, 4, 1]
y = [3*entry for entry in y]
xerr1_ = [[1.86, 1.37, 1.28, 1.74, 1.12, 0.91, 1.54, 6.04, 5.89, 9.10, 3.52, 7.97, 6.05, 0.43], [2.27, 1.40, 1.26, 1.76, 1.51, 1.03, 1.61, 4.25, 7.07, 8.88, 4.67, 9.52, 5.74, 0.46]]
xerr2_ = [[2*x for x in xerr1_[0]], [2*x for x in xerr1_[1]]]

xlow = -20
xhigh = 20
ylow = 0
yhigh = 105
red_patch = mpatches.Patch(color='red', label='$\pm$1 $\sigma$')
blue_patch = mpatches.Patch(color='blue', label='$\pm$2 $\sigma$')

plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)

plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.cms.text("Toy Simulation with Systematics", fontsize=12)
#hep.cms.label(loc=0)

plt.errorbar(x, y,xerr=xerr2_, fmt='sk',ecolor='b',elinewidth=2,markersize = 0.3)
plt.errorbar(x, y,xerr=xerr1_, fmt='sk',ecolor='r',elinewidth=4,markersize = 4)
plt.legend(handles=[red_patch, blue_patch], loc='upper left', prop={'size': 12})

for i in range(len(x)):
  #plt.annotate(str(x[i])+'$\pm$'+str(xerr1_[i]), xy=(x[i], y[i]), xytext=(8, y[i]), size = 12)
  if i == 13:
    plt.annotate('%.2f'%x[i], xy=(x[i], y[i]), xytext=(13, y[i]+1), size = 12, color='magenta')
  else:
    plt.annotate('%.2f'%x[i], xy=(x[i], y[i]), xytext=(13, y[i]+1), size = 12)
  plt.annotate('- ' + str(xerr1_[0][i]), xy=(x[i], y[i]), xytext=(17, y[i]+0.4), size = 8)
  plt.annotate('+' + str(xerr1_[1][i]), xy=(x[i], y[i]), xytext=(17, y[i]+3.4), size = 8)

plt.yticks(y, ['ggF1', 'ggF2', 'ggF3', 'ggF4', 'VBF1', 'VBF2', 'VBF3', 'VBF4', 'VH3l', 'VHMET', 'ttH lep', 'ttH had', 'Untag', 'Combine'], size=12)
plt.tick_params(axis="x",direction="in")
plt.tick_params(axis="y",direction="in")
plt.xticks(np.arange(xlow, xhigh, step=1))

plt.vlines(1, ylow, yhigh, colors='k', linestyles ='dashed', dashes=(0, (5, 10))).set_linewidth(0.5)
plt.vlines(x[13], ylow, yhigh, colors='magenta').set_linewidth(1)
plt.annotate('SM', xy=(2, 3),size=12)

#plt.title('CMS', loc='left', fontname="sans serif", fontstyle='italic', fontsize=18)
#plt.title('Toy Simulation', x=0.22,y=1, fontname="sans serif", fontstyle='italic',fontsize=11)
#plt.title('137.61 fb$^{-1}$ (13 TeV) + 62.32 fb$^{-1}$ (13.6 TeV)', loc = 'right', fontstyle='normal', fontname="sans serif", fontsize=10)
plt.xlabel('$\mu$')
plt.savefig('toy_signal.png')
