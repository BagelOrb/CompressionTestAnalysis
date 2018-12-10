
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

for dithering in [True, False]:
    top_color = 'blue' # if dithering else 'darkcyan'
    side_color = 'orange' # if dithering else 'red'

    dither_str = ('_no' if not dithering else '') + '_dither'

    fig = pyplot.figure(int(dithering))
    fig.tight_layout()
    ax = pyplot.gca()
    ax.autoscale()

    for top in [True, False]:
        color = top_color if top else side_color
        label_stress = 'E_z' if top else 'E_{xy}'
        # label = label + dither_str
        data_here = data[np.logical_and(data[:,0] == dithering, data[:,1] == top),:]
        densities = data_here[:, 2]


        pyplot.figure(int(dithering))
        ax = pyplot.gca()
        plateau_stresses = data_here[:, 4]
        ax.scatter(densities, plateau_stresses, color=color, label='$'+label_stress+'$')
        if dithering:
            fit = np.polyfit(densities, plateau_stresses, 1)
            fit_fn = np.poly1d(fit)
            fit_label = '\\hat' + label_stress
            pyplot.plot(x_range, fit_fn(x_range), color=color, label='$'+fit_label+'$')
        else:
            densities_uniq = np.unique(densities)
            avg_stresses = []
            for density in densities_uniq:
                stresses = data_here[densities == density, 4]
                avg_stresses.append(np.mean(stresses))

            avg_label = '\\bar' + label_stress
            pyplot.plot(densities_uniq, avg_stresses, color=color, label='$'+avg_label+'$', drawstyle ='steps-mid')

    pyplot.figure(int(dithering))
    ax = pyplot.gca()

    for density in [10.0592231765546, 14.2258898432212, 20.1184463531091, 28.4517796864425]:
        ax.axvline(density, 0, 1, linestyle='dashed', color='grey')  # plot all whole densities
    ax.set_xlabel('Density (%)')
    ax.set_ylabel('Plateau stress ($MPa$)')
    ax.set_ylim([0.0, 0.62])
    pyplot.legend()
    fig.savefig('analysis/results/plateau_height_top_and_side' + dither_str +'.pdf', bbox_inches='tight')

# pyplot.show()



pyplot.clf()

x_range = [np.min(data[:,2]), np.max(data[:,2])]

for dithering in [True, False]:
    top_color = 'blue' # if dithering else 'darkcyan'
    side_color = 'orange' # if dithering else 'red'

    dither_str = ('_no' if not dithering else '') + '_dither'

    fig = pyplot.figure(10 + int(dithering))
    fig.tight_layout()
    ax = pyplot.gca()
    ax.autoscale()

    for top in [True, False]:
        color = top_color if top else side_color
        label_modulus = '\\sigma_z' if top else '\\sigma_{xy}'
        # label = label + dither_str
        data_here = data[np.logical_and(data[:,0] == dithering, data[:,1] == top),:]
        densities = data_here[:, 2]

        pyplot.figure(10 + int(dithering))
        ax = pyplot.gca()
        youngs_moduli = data_here[:,3]
        ax.scatter(densities, youngs_moduli, color=color, label='$'+label_modulus+'$')

        if dithering:
            fit = np.polyfit(densities, youngs_moduli, 1)
            fit_fn = np.poly1d(fit)
            fit_label = '\\hat' + label_modulus
            pyplot.plot(x_range, fit_fn(x_range), color=color, label='$'+fit_label+'$')
        else:
            densities_uniq = np.unique(densities)
            avg_moduli = []
            for density in densities_uniq:
                moduli = data_here[densities == density, 3]
                avg_moduli.append(np.mean(moduli))

            avg_label = '\\bar' + label_modulus
            pyplot.plot(densities_uniq, avg_moduli, color=color, label='$'+avg_label+'$', drawstyle ='steps-mid')


    pyplot.figure(10 + int(dithering))
    ax = pyplot.gca()
    for density in [10.0592231765546, 14.2258898432212, 20.1184463531091, 28.4517796864425]:
        ax.axvline(density, 0, 1, linestyle='dashed', color='grey')  # plot all whole densities
    ax.set_xlabel('Density (%)')
    ax.set_ylabel('Youngs modulus ($MPa$)')
    ax.set_ylim([0.0, 6.0])
    pyplot.legend()
    fig.savefig('analysis/results/youngs_modulus_top_and_side' + dither_str +'.pdf', bbox_inches='tight')

# pyplot.show()