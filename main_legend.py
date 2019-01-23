import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm       # color map
from matplotlib.ticker import FuncFormatter



densities = np.array([
    0.10059
    ,.11101
    ,.12143
    ,.13184
    ,.14226
    ,.14962
    ,.15699
    ,.16436
    ,.17172
    ,.17909
    ,.18645
    ,.19382
    ,.20118
    ,.22202
    ,.24285
    ,.26368
    ,.28452
    ])

gradient = np.linspace(.10, .285, 256)
gradient = np.vstack((gradient, gradient))


cmap = cm.nipy_spectral

fig = plt.figure()
ax = plt.gca()
#fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
#ax.set_title('colormaps', fontsize=14)

ax.set_ylim(.10, .285)
ax.set_yticks(densities)
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.2%}'.format(y)))

ax.imshow(gradient.T, aspect='auto', extent=[.1,.285,.1,.285], cmap=cmap, origin='lower')
pos = list(ax.get_position().bounds)

ax.set_xticks([])
# Turn off *all* ticks & spines, not just the ones with colormaps.
#ax.set_axis_off()

plt.show()