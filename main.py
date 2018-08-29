
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
from typing import NamedTuple

import numpy as np
from scipy import signal    # savgol_filter

from matplotlib import cm       # color map
from matplotlib import pyplot   # plotting
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.axes3d import Axes3D



def getColor(test):
    if density_colors:
        return cm.rainbow((test.limit_density / 100 - .1) / .3)
    else:
        return cm.rainbow(test.printer_number  / 6)

def getLabel(test):
    if density_colors:
        return str(int(test.limit_density))
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
        line_y = np.array([plateau.height, plateau.height])
        # plateau_xs[i] = strain_comp[plateau.]
        pyplot.plot(line_x, line_y, color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))

        if enable_3D_plot:
            line_height = np.array([test.limit_density, test.limit_density])
            ax.plot(line_x, line_height, line_y, color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))



def plotPlateaus3D():
    plateau_densities: List[List[float]] = []
    plateau_strains: List[List[float]] = []
    plateau_stresses: List[List[float]] = []
    for test in tests:
        strain_comp = test.strain[test.compression_range]
        plateau = Plateau(test)

        line_x = np.array([strain_comp[plateau.start_idx], strain_comp[plateau.end_idx]])
        plateau_strains.append(line_x)

        plateau_densities.append([test.limit_density, test.limit_density])

        plateau_stresses.append([plateau.height, plateau.height])

    xs = np.array(plateau_strains)
    ys = np.array(plateau_densities)
    zs = np.array(plateau_stresses)

    ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm, alpha=0.5)



#

test_top = False
enable_3D_plot = True
plot_tangent = False
density_colors = True



test_file_names = TestCases.top_file_names if test_top else TestCases.side_file_names

tests = []

densities: List[float] = []
compression_energies: List[float] = []
energy_diffs: List[float] = []
energy_ratios: List[float] = []
min_secant_moduli: List[float] = []
max_tangent_moduli: List[float] = []
end_tangent_moduli: List[float] = []

if enable_3D_plot:
    fig, ax = pyplot.subplots(subplot_kw = {'projection': '3d'})
    ax.set_xlabel('Strain ε (%)')
    ax.set_ylabel('Density (%)')
    ax.set_zlabel('Stress σ (MPa)')

for test_file_name in test_file_names:
    print('gathering info for: ' + test_file_name)
    test = CompSlowDecompTest("test_results/Top/" if test_top else "test_results/Side/", test_file_name)
    tests.append(test)





pyplot.figure(0)
plotCompressions()
pyplot.legend()
pyplot.xlabel('Strain ε (%)')
pyplot.ylabel('Stress σ (MPa)')
if enable_3D_plot:
    fig.legend()
    plotPlateaus3D()

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
'''
pyplot.figure(3)
pyplot.autoscale()
pyplot.scatter(densities, compression_energies)
pyplot.xlabel('Structure density (%)')
pyplot.ylabel('Compressive energy (J?)')

# pyplot.figure(4)
# pyplot.plot(densities, energy_ratios)
# pyplot.xlabel('Structure density (%)')
# pyplot.ylabel('energy absorption ratio (J?)')

pyplot.figure(5)
pyplot.autoscale()
pyplot.scatter(densities, min_secant_moduli)
pyplot.xlabel('density')
pyplot.ylabel('Minimal secant modulus')

pyplot.figure(6)
pyplot.autoscale()
pyplot.scatter(densities, max_tangent_moduli)
pyplot.xlabel('density')
pyplot.ylabel('Max tangent modulus')

pyplot.figure(7)
pyplot.autoscale()
pyplot.scatter(densities, end_tangent_moduli)
pyplot.xlabel('density')
pyplot.ylabel('tangent modulus at 2 kN')
'''



pyplot.show()











