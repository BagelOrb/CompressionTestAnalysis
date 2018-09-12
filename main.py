
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')

import TestCases
from CompSlowDecompTest import CompSlowDecompTest
from TangentModulus import TangentModulus
import PolyFitting
import PlottingUtil
import MathUtils
from Plateau import Plateau

from typing import List
from typing import Dict

import numpy as np
from scipy import signal    # savgol_filter

import functools

from matplotlib import cm       # color map
from matplotlib import pyplot   # plotting
import matplotlib.ticker        # defining ticks on the axes
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.axes3d import Axes3D

import pickle



cmap = cm.nipy_spectral

def getColor(test):
    if density_colors:
        return cmap((test.limit_density / 100 - .1) / (.285 - .1))
    else:
        return cmap(test.printer_number / 6)

def getLabel(test):
    if density_colors:
        return format(test.limit_density, '2.2f')
    else:
        return test.printer_name

#



def gatherStatistics():
    for test in tests:
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]
        strain_decomp = test.strain[test.decompression_range]
        stress_decomp = test.stress[test.decompression_range]

        compression_energy = abs(MathUtils.integral(strain_comp, stress_comp)[-1])
        decompression_energy = abs(MathUtils.integral(strain_decomp, stress_decomp)[-1])
        densities.append(test.limit_density)
        compression_energies.append(compression_energy)
        energy_diffs.append(compression_energy - decompression_energy)
        energy_ratios.append((compression_energy - decompression_energy) / compression_energy)
        stress_t = stress_comp[100:]
        strain_t = strain_comp[100:] - strain_comp[0]
        min_secant_moduli.append(min(stress_t / strain_t))
        max_tangent_moduli.append(TangentModulus.getTangentModulus(test).val)
        end_tangent_moduli.append(TangentModulus.getEndTangentModulus(test).val)


def plotCompressions():
    if plot_tangent:
        for test in tests:
            color = PlottingUtil.lighten_color(getColor(test), .25)
            tangent_modulus: TangentModulus = TangentModulus.getTangentModulus(test)
            (xs, ys) = TangentModulus.getTangentLine(tangent_modulus)
            pyplot.figure(0)
            # pyplot.scatter(xs[0:1], ys[0:1])
            pyplot.plot(xs, ys, color=color)
            if enable_3D_plot:
                ax.plot(xs, [test.limit_density, test.limit_density], ys, color=color)

    for test in tests:
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]

        color = getColor(test)
        if enable_3D_plot:
            ax.plot(strain_comp, np.ones(len(test.compression_range)) * test.limit_density, stress_comp, color = color, label = getLabel(test))

        pyplot.figure(0)
        pyplot.plot(strain_comp, stress_comp, color = color, label = getLabel(test))

        # fit = PolyFitting.fit(strain_comp, stress_comp, 25)
        # fit_map = np.vectorize(fit)
        # pyplot.plot(strain_comp, fit_map(strain_comp), color = getColor(test), label = getLabel(test))


def plotCompressionDerivativesFitting():
    for test in tests:
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]

        der = PolyFitting.getDerivative(strain_comp, stress_comp, 1, 25)
        der_map = np.vectorize(der)
        pyplot.plot(strain_comp, der_map(strain_comp), color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))


def plotCompressionDerivativesSmoothing():
    for test in tests:
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]

        (strain_comp, stress_comp) = signal.savgol_filter((strain_comp, stress_comp), window_length = 51, polyorder = 3)

        pyplot.plot(strain_comp[:-1], MathUtils.derivative(strain_comp, stress_comp),
                    color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))


def plotPlateaus():
    for i in range(len(tests)):
        test = tests[i]
        strain_comp = test.strain[test.compression_range]
        # stress_comp = test.stress[test.compression_range]
        plateau = Plateau(test)
        line_x = np.array([strain_comp[plateau.start_idx], strain_comp[plateau.end_idx]])
        line_y = np.array([plateau.stress, plateau.stress])
        # plateau_xs[i] = strain_comp[plateau.]
        pyplot.plot(line_x, line_y, color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))

        if enable_3D_plot:
            line_height = np.array([test.limit_density, test.limit_density])
            ax.plot(line_x, line_height, line_y, color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))



def plotPlateaus3D():
    plateaus_per_density: Dict[float, List[Plateau]] = {}
    for test in tests:
        plateaus_per_density[test.limit_density] = []

    min_stress = 99999
    max_stress = 0
    for test in tests:
        plateau = Plateau(test)
        min_stress = min(min_stress, plateau.stress)
        max_stress = max(max_stress, plateau.stress)
        plateaus_per_density[test.limit_density].append(plateau)

    plateau_vertex_strains: List[List[float]] = []
    plateau_vertex_stresses: List[List[float]] = []
    plateau_vertex_densities: List[List[float]] = []
    plateau_strains_start: List[float] = []
    plateau_strains_end: List[float] = []
    plateau_stresses: List[float] = []
    plateau_densities: List[float] = []
    for density, plateaus in sorted(plateaus_per_density.items()):
        strain_start = functools.reduce(lambda tot, plat: tot + plat.strain_start, plateaus, 0.0) / len(plateaus)
        strain_end = functools.reduce(lambda tot, plat: tot + plat.strain_end, plateaus, 0.0) / len(plateaus)
        stress = functools.reduce(lambda tot, plat: tot + plat.stress, plateaus, 0.0) / len(plateaus)

        plateau_strains_start.append(strain_start)
        plateau_strains_end.append(strain_end)
        plateau_stresses.append(stress)
        plateau_densities.append(density)

        strains = [strain_start, strain_end]
        stresses = [stress, stress]
        densities = [density, density]
        color = PlottingUtil.darken_color(cmap((stress - min_stress) / (max_stress - min_stress)), .5)
        ax.plot(np.array(strains), np.array(densities), np.array(stresses), color=color)

        plateau_vertex_strains.append(strains)
        plateau_vertex_densities.append(densities)
        plateau_vertex_stresses.append(stresses)

    xs = np.array(plateau_vertex_strains)
    ys = np.array(plateau_vertex_densities)
    zs = np.array(plateau_vertex_stresses)

    ax.plot_surface(xs, ys, zs, cmap=cmap, alpha=0.5)
    ax.plot(plateau_strains_start, plateau_densities, plateau_stresses, color='black', alpha=0.15)
    ax.plot(plateau_strains_end, plateau_densities, plateau_stresses, color='black', alpha=0.15)

    ax.plot_surface(xs, ys, zs * 0, cmap=cmap, alpha=0.15)
    ax.plot(plateau_strains_start, plateau_densities, np.zeros(len(plateau_densities)), color='black', alpha=0.15)
    ax.plot(plateau_strains_end, plateau_densities, np.zeros(len(plateau_densities)), color='black', alpha=0.15)

    max_strain = 0.8 # max(tests, key = lambda test: np.maximum(test.strain[test.compression_range]))
    ax.plot(np.ones(len(plateau_densities)) * max_strain, plateau_densities, plateau_stresses, color='black', alpha=0.25)

    # plot projection lines for start and end
    for i in [0, -1]:
        ax.plot(
            [plateau_strains_end[i], plateau_strains_end[i], max_strain],
            [plateau_densities[i], plateau_densities[i], plateau_densities[i]],
            [0, plateau_stresses[i], plateau_stresses[i]],
            color='black', alpha=0.15, linestyle='--')
        ax.plot(
            [plateau_strains_start[i], plateau_strains_start[i]],
            [plateau_densities[i], plateau_densities[i]],
            [0, plateau_stresses[i]],
            color='black', alpha=0.15, linestyle='--')


#

# test_top = False
enable_stress_strain_plot = True
enable_3D_plot = True
plot_plateaus = True
plot_tangent = False
density_colors = True
export = True


all_densities: List[List[float]] = []
all_max_tangent_moduli: List[List[float]] = []

for test_top in {True, False}:



    test_file_names = TestCases.top_file_names if test_top else TestCases.side_file_names

    tests = []

    densities: List[float] = []
    compression_energies: List[float] = []
    energy_diffs: List[float] = []
    energy_ratios: List[float] = []
    min_secant_moduli: List[float] = []
    max_tangent_moduli: List[float] = []
    end_tangent_moduli: List[float] = []

    for test_file_name in test_file_names:
        print('gathering info for: ' + test_file_name)
        test = CompSlowDecompTest("test_results/Top/" if test_top else "test_results/Side/", test_file_name)
        tests.append(test)

    tests = sorted(tests, key = lambda test: - test.limit_density)

    ticks_x = matplotlib.ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x * 100)) # axes tick formatter for percentages

    if enable_3D_plot:
        fig, ax = pyplot.subplots(subplot_kw = {'projection': '3d'})
        ax.set_xlabel('Strain ε (%)')
        ax.set_ylabel('Density (%)')
        ax.set_zlabel('Stress σ (MPa)')

    gatherStatistics()

    all_densities.append(densities)
    all_max_tangent_moduli.append(max_tangent_moduli)

    if enable_3D_plot and plot_plateaus:
        plotPlateaus3D()

    if enable_stress_strain_plot:
        fig2D = pyplot.figure(0)
        fig2D.tight_layout()
        plotCompressions()
        pyplot.legend()
        pyplot.xlabel('Strain ε (%)')
        pyplot.ylabel('Stress σ (MPa)')
        pyplot.subplot().xaxis.set_major_formatter(ticks_x)

    if enable_3D_plot:
        fig.legend()
        fig.tight_layout()
        ax.view_init(elev=22, azim=162)
        ax.xaxis.set_major_formatter(ticks_x)
        ax.invert_yaxis()

    # pyplot.figure(2)
    # plotPlateaus()

    # pyplot.figure(2)
    # plotCompressionDerivativesFitting()

    # pyplot.figure(2)
    # plotCompressionDerivativesSmoothing()




    # pyplot.figure(2)
    # pyplot.plot(densities, energy_diffs)
    # pyplot.xlabel('Structure density (%)')
    # pyplot.ylabel('Energy consumption (J?)')

    # pyplot.figure(3)
    # pyplot.autoscale()
    # pyplot.scatter(densities, compression_energies)
    # pyplot.xlabel('Structure density (%)')
    # pyplot.ylabel('Compressive energy (J?)')
    
    # pyplot.figure(4)
    # pyplot.plot(densities, energy_ratios)
    # pyplot.xlabel('Structure density (%)')
    # pyplot.ylabel('energy absorption ratio (J?)')
    
    # pyplot.figure(5)
    # pyplot.autoscale()
    # pyplot.scatter(densities, min_secant_moduli)
    # pyplot.xlabel('density')
    # pyplot.ylabel('Minimal secant modulus')
    
    pyplot.figure(6)
    pyplot.autoscale()
    colors = list(map(lambda test: getColor(test), tests))
    pyplot.scatter(densities, max_tangent_moduli, color=colors)
    pyplot.xlabel('Density (%)')
    pyplot.ylabel('Youngs modulus (MPa)')
    
    # pyplot.figure(7)
    # pyplot.autoscale()
    # pyplot.scatter(densities, end_tangent_moduli)
    # pyplot.xlabel('density')
    # pyplot.ylabel('tangent modulus at 2 kN')




    if export:
        config = 'top' if test_top else 'side'

        pyplot.figure(6).savefig('analysis/results/youngs_modulus_' + config + '.pdf', bbox_inches='tight')

        # pickle.dump(pyplot.figure(0), open('analysis/results/stress_strain_' + config + '.pickle', 'wb'))
        pyplot.figure(0).savefig('analysis/results/stress_strain_' + config + '.pdf', bbox_inches='tight')
        if enable_3D_plot:
            # pickle.dump(pyplot.figure(1), open('analysis/results/stress_strain_' + config + '_3D.pickle', 'wb'))
            pyplot.figure(1).savefig('analysis/results/stress_strain_' + config + '_3D.pdf', bbox_inches='tight')

    # pyplot.show()



pyplot.figure(16)
pyplot.autoscale()
pyplot.scatter(all_densities[0], all_max_tangent_moduli[0], color='blue', label='Top')
pyplot.scatter(all_densities[1], all_max_tangent_moduli[1], color='orange', label='Side')
pyplot.xlabel('Density (%)')
pyplot.ylabel('Youngs modulus (MPa)')
pyplot.legend()
pyplot.figure(16).savefig('analysis/results/youngs_modulus_top_and_side.pdf', bbox_inches='tight')


pyplot.show()

