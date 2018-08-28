
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')

import TestCases
from CompSlowDecompTest import CompSlowDecompTest
from TangentModulus import TangentModulus
import PolyFitting
import PlottingUtil

from typing import List
from typing import NamedTuple

import numpy as np

from matplotlib import cm       # color map
from matplotlib import pyplot   # plotting
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d.axes3d import Axes3D


def integral(disp: np.array, force: np.array) -> np.array:
    ret = np.zeros(disp.size)

    total_volume: float = 0
    for i in range(disp.size - 1):
        volume = (disp[i+1] - disp[i]) * 0.5 * (force[i+1] + force[i])
        total_volume += volume
        ret[i+1] = total_volume

    return ret





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

        compression_energy = abs(integral(strain_comp, stress_comp)[-1])
        decompression_energy = abs(integral(strain_decomp, stress_decomp)[-1])
        densities.append(test.limit_density)
        compression_energies.append(compression_energy)
        energy_diffs.append(compression_energy - decompression_energy)
        energy_ratios.append((compression_energy - decompression_energy) / compression_energy)
        stress_t = stress_comp[100:-1]
        strain_t = strain_comp[100:-1] - strain_comp[0]
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


def plotCompressionDerivatives():
    for test in tests:
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]

        pyplot.figure(2)
        der = PolyFitting.getDerivative(strain_comp, stress_comp, 1, 25)
        der_map = np.vectorize(der)
        pyplot.plot(strain_comp, der_map(strain_comp), color = PlottingUtil.lighten_color(getColor(test), .7), label = getLabel(test))


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



plotCompressions()
plotCompressionDerivatives()




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

pyplot.figure(0)
pyplot.legend()
pyplot.xlabel('Strain ε (%)')
pyplot.ylabel('Stress σ (MPa)')

if enable_3D_plot:
    fig.legend()

pyplot.show()











