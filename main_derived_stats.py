
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')

import numpy as np

from matplotlib import cm       # color map
from matplotlib import pyplot   # plotting
import matplotlib.ticker        # defining ticks on the axes

test_case_data_filename = '/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments/analysis/all_statistics.csv'

with open(test_case_data_filename) as test_case_data_file:
    data = np.genfromtxt(test_case_data_filename, delimiter=',', skip_header=1)



pyplot.clf()

x_range = [np.min(data[:,2]), np.max(data[:,2])]

figsize = (4, 3)

symbol = {3:'E', 4:'\\sigma'}
name = {3: 'Youngs modulus', 4: 'Plateau stress'}
file = {3: 'youngs_modulus', 4: 'plateau_stress'}

def plotStatistic(stat_row):
    for dithering in [True, False]:
        top_color = 'blue'  # if dithering else 'darkcyan'
        side_color = 'orange'  # if dithering else 'red'

        dither_str = ('_no' if not dithering else '') + '_dither'

        fig = pyplot.figure(dithering * 2 + stat_row, figsize=figsize)
        fig.tight_layout()
        ax = pyplot.gca()
        # ax.autoscale()
        for density in [10.0592231765546, 14.2258898432212, 20.1184463531091, 28.4517796864425]:
            ax.axvline(density, 0, 1, linestyle='dashed', color='grey')  # plot all whole densities
        ax.set_xlabel('Density (%)')
        ax.set_ylabel(name[stat_row] + ' ($MPa$)')
        if stat_row == 3:
            ax.set_ylim([0.0, 6.0])
        else:
            ax.set_ylim([0.0, 0.62])

        for top in [True, False]:
            color = top_color if top else side_color
            label_stress = symbol[stat_row] + ('_z' if top else '_{xy}')
            # label = label + dither_str
            data_here = data[np.logical_and(data[:, 0] == dithering, data[:, 1] == top), :]
            densities = data_here[:, 2]

            stats = data_here[:, stat_row]
            ax.scatter(densities, stats, color=color, label='$' + label_stress + '$')
            if dithering:
                fit = np.polyfit(densities, stats, 1)
                fit_fn = np.poly1d(fit)
                fit_label = '\\hat' + label_stress
                pyplot.plot(x_range, fit_fn(x_range), color=color, label='$' + fit_label + '$')
            else:
                stats_uniq = np.unique(densities)
                avg_stats = []
                for density in stats_uniq:
                    stats_ = data_here[densities == density, stat_row]
                    avg_stats.append(np.mean(stats_))

                avg_label = '\\bar' + label_stress
                pyplot.plot(stats_uniq, avg_stats, color=color, label='$' + avg_label + '$', drawstyle='steps-mid')

        ax.legend(loc=2)
        fig.savefig('analysis/results/' + file[stat_row] + '_top_and_side' + dither_str + '.pdf', bbox_inches='tight')
        pyplot.clf()
        pyplot.close(fig)

plotStatistic(3)
plotStatistic(4)

pyplot.close()