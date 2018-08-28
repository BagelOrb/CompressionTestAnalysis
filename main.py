
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')

import TestCases
from CompSlowDecompTest import CompSlowDecompTest
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


Tangent = NamedTuple('Tangent', [('source_x', float), ('source_y', float), ('val', float)])


def getTangentModulus(test: CompSlowDecompTest, polynomial_order: int = 10, overall_polynomial_order: int = 25) -> Tangent:
    strain_comp = test.strain[test.compression_range]
    stress_comp = test.stress[test.compression_range]
    derivative_optima_indices = PolyFitting.getRoots(strain_comp, stress_comp, 2, overall_polynomial_order) # where the derivative is maximal or minimal

    max_data_index = strain_comp.size - 1

    # determine the data on which to fit a polynomial
    der_3 = PolyFitting.getDerivative(strain_comp, stress_comp, 3, overall_polynomial_order)
    derivative_optima = strain_comp[derivative_optima_indices]
    mapping_der_3 = np.vectorize(der_3) # so that we can map the function over an np.array
    derivative_optimum_direction = mapping_der_3(derivative_optima) # positive for minima, negative for maxima
    der_minima_indices = np.where(np.logical_and(derivative_optimum_direction > 0, derivative_optima > 0.05))[0]
    if len(der_minima_indices) > 0:
        max_data_index = derivative_optima_indices[der_minima_indices[0]]

    # fit curve to first part of the graph and get the max derivative from that
    strain_comp_sub = strain_comp[0:max_data_index]
    stress_comp_sub = stress_comp[0:max_data_index]
    derivative_optima_indices_sub = PolyFitting.getRoots(strain_comp_sub, stress_comp_sub, 2, polynomial_order) # where the derivative is maximal or minimal
    der_sub = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 1, polynomial_order)
    der_3_sub = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 3, polynomial_order)
    first_optimum_is_max = der_3_sub(strain_comp_sub[derivative_optima_indices_sub[0]]) < 0 # to determine whether it's a top or a bottom peek in the derivative
    tangent: float = 0
    source_x: float = 0
    if first_optimum_is_max:
        strain_at_derivative_max = strain_comp_sub[derivative_optima_indices_sub[0]-1]
        assert(strain_at_derivative_max > 0)
        tangent = der_sub(strain_at_derivative_max)
        source_x = strain_at_derivative_max
    else:
        tangent = der_sub(0)
        source_x = 0
    assert(tangent > 0)
    polynomial = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 0, polynomial_order)
    source_y = polynomial(source_x)

    #
    # pyplot.clf()
    # pyplot.plot(strain_comp, stress_comp)
    # pyplot.plot(strain_comp_sub, stress_comp_sub)
    # pyplot.plot(strain_comp_sub, polynomial(strain_comp_sub))
    # pyplot.show()
    # exit(1)

    return Tangent(source_x, source_y, tangent)


# get the tangent modulus at the end of the graph
def getEndTangentModulus(test: CompSlowDecompTest, polynomial_order: int = 10, overall_polynomial_order: int = 25) -> Tangent:
    strain_comp = test.strain[test.compression_range]
    stress_comp = test.stress[test.compression_range]
    derivative_optima_indices = PolyFitting.getRoots(strain_comp, stress_comp, 2, overall_polynomial_order) # where the derivative is maximal or minimal

    min_data_index = 0

    # determine the data on which to fit a polynomial
    der_3 = PolyFitting.getDerivative(strain_comp, stress_comp, 3, overall_polynomial_order)
    mapping_der_3 = np.vectorize(der_3)
    derivative_optima = strain_comp[derivative_optima_indices]
    derivative_optimum_direction = mapping_der_3(derivative_optima) # positive for minima, negative for maxima
    der_minima_indices = np.where(derivative_optimum_direction > 0)[0]
    if len(der_minima_indices) > 0:
        min_data_index = derivative_optima_indices[der_minima_indices[-1]]

    # fit curve to first part of the graph and get the max derivative from that
    strain_comp_sub = strain_comp[min_data_index:-1]
    stress_comp_sub = stress_comp[min_data_index:-1]
    polynomial = PolyFitting.fit(strain_comp_sub, stress_comp_sub, polynomial_order)
    derivative = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 1, polynomial_order)
    source_x = strain_comp_sub[-1]
    source_y = polynomial(source_x)
    tangent = derivative(source_x)

    # pyplot.clf()
    # pyplot.plot(strain_comp, stress_comp)
    # pyplot.plot(strain_comp_sub, stress_comp_sub)
    # pyplot.plot(strain_comp_sub, polynomial(strain_comp_sub))
    # pyplot.show()

    return Tangent(source_x, source_y, tangent)


def getTangentLine(tangent: Tangent) -> (np.array, np.array):
    tangent_max = np.array([.4, .4 * tangent.val])
    if tangent_max[1] > 0.45:
        tangent_max = tangent_max * 0.45 / (.4 * tangent.val)
    xs = np.array([0, tangent_max[0]]) + tangent.source_x
    ys = np.array([0, tangent_max[1]]) + tangent.source_y
    return (xs, ys)


def getColor(test):
    if density_colors:
        return cm.rainbow((test.limit_density / 100 - .1) / .3)
    else:
        return cm.rainbow(test.printer_number  / 6)

def getLabel(test):
    if density_colors:
        return str(int(test.limit_density * 100))
    else:
        return test.printer_name

#


test_top = False
enable_3D_plot = False
plot_tangent = True

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

if plot_tangent:
    for test in tests:
        color = PlottingUtil.lighten_color(getColor(test), .25)
        tangent_modulus: Tangent = getTangentModulus(test)
        (xs, ys) = getTangentLine(tangent_modulus)
        pyplot.figure(0)
        # pyplot.scatter(xs[0:1], ys[0:1])
        pyplot.plot(xs, ys, color = color)
        if enable_3D_plot:
            ax.plot(xs, [test.limit_density, test.limit_density], ys, color = color)

if True:
    for test in tests:

        color = getColor(test)

        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]
        strain_decomp = test.strain[test.decompression_range]
        stress_decomp = test.stress[test.decompression_range]

        if enable_3D_plot:
            ax.plot(strain_comp, np.ones(len(test.compression_range)) * test.limit_density, stress_comp, color = color, label = getLabel(test))

        pyplot.figure(0)
        pyplot.plot(strain_comp, stress_comp, color = color, label = getLabel(test))

        compression_energy = abs(integral(strain_comp, stress_comp)[-1])
        decompression_energy = abs(integral(strain_decomp, stress_decomp)[-1])
        densities.append(test.limit_density)
        compression_energies.append(compression_energy)
        energy_diffs.append(compression_energy - decompression_energy)
        energy_ratios.append((compression_energy - decompression_energy) / compression_energy)
        stress_t = stress_comp[100:-1]
        strain_t = strain_comp[100:-1] - strain_comp[0]
        min_secant_moduli.append(min(stress_t / strain_t))
        max_tangent_moduli.append(getTangentModulus(test).val)
        end_tangent_moduli.append(getEndTangentModulus(test).val)


if False:
    i = 0
    for test in tests:
        pyplot.figure(20 + i)
        i += 1
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]
        PolyFitting.fitPolynomial(strain_comp, stress_comp)
    pyplot.show()

'''
pyplot.figure(2)
pyplot.plot(densities, energy_diffs)
pyplot.xlabel('Structure density (%)')
pyplot.ylabel('Energy consumption (J?)')

pyplot.figure(3)
pyplot.plot(densities, compression_energies)
pyplot.xlabel('Structure density (%)')
pyplot.ylabel('Compressive energy (J?)')

pyplot.figure(4)
pyplot.plot(densities, energy_ratios)
pyplot.xlabel('Structure density (%)')
pyplot.ylabel('energy absorption ratio (J?)')

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











