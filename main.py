
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')

from DataCutting import DataCutting
from ElastoPlasticDeformationCutter import ElastoPlasticDeformationCutter
from CompSlowDecompTest import CompSlowDecompTest
import PolyFitting
import PlottingUtil

from typing import List
from typing import NamedTuple

import numpy as np

from matplotlib import cm       # color map
from matplotlib import pyplot   # plotting
from mpl_toolkits.mplot3d.axes3d import Axes3D


def computeDeformationRecovery(data: np.array, cycle_count: int = 5) -> np.array:
    ret = np.zeros(cycle_count - 1, 2)

    disp = data[:,0]
    force = data[:,1]

    cutter: DataCutting = DataCutting.separateData(disp, cycle_count)

    disp_here = disp[cutter.compression_ranges[0]]
    force_here = force[cutter.compression_ranges[0]]
    start_disp_cutoff_index = ElastoPlasticDeformationCutter.getNongrippedDisplacementIndex(disp_here,
                                                                                 force_here, True)
    start_disp = disp_here[start_disp_cutoff_index]
    start_time = data[start_disp_cutoff_index:,2]
    for i in range(cycle_count - 1):
        decompression_range = cutter.decompression_ranges[i]
        disp_here = disp[decompression_range]
        force_here = force[decompression_range]
        cutoff_index = ElastoPlasticDeformationCutter.getNongrippedDisplacementIndex(disp_here, force_here, False)
        ... # not done yet

    return ret


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

    return Tangent(source_x, source_y, tangent)


def getTangentLine(tangent: Tangent) -> (np.array, np.array):
    tangent_max = np.array([.4, .4 * tangent.val])
    if tangent_max[1] > 0.45:
        tangent_max = tangent_max * 0.45 / (.4 * tangent.val)
    xs = np.array([0, tangent_max[0]]) + tangent.source_x
    ys = np.array([0, tangent_max[1]]) + tangent.source_y
    return (xs, ys)



top_file_names: List[str] = [
    "3.04_100_9_PW3",
    "3.04_100_11_D",
    "0.76_11.10_6_PW1",
    "0.76_11.10_2_T",
    "0.76_12.14_11_T",
    "0.76_12.14_8_PW2",
    "0.76_13.18_2_PW2",
    "0.76_13.18_1_PW2",
    "2.15_100_14_PW2",
    "2.15_100_14_PW4",
    "2.15_100_13_T",
    "2.15_100_3_T",
    "2.15_100_11_PW3",
    "0.76_15.7_6_PW4",
    "0.76_15.70_2_PW3",
    "0.76_17.17_7_T",
    "0.76_17.17_7_PW3",
    "0.76_18.65_4_D",
    "0.76_18.65_2_D",
    "1.52_100_4_PW3",
    "1.52_100_4_PW2",
    "0.76_22.20_13_PW1",
    "0.76_22.20_7_PW2",
    "0.76_24.29_7_PW3",
    "0.76_24.29_10_PW2",
    "0.76_26.37_5_d",
    "0.76_26.37_3_D",
    "1.075_100_3_PW3",
    "1.075_100_15_D",
    "1.075_100_3_PW2",
    "1.075_100_15_PW4",
    "1.075_100_15_T",
    "1.075_100_15_PW3"
    ]

# "0.76_22.20_7_PW1",# --> this test seems to have been performed wrongly!!!

side_file_names: List[str] = [
    "3.04_100_11_PW4",
    "3.04_100_12_PW4",
    "0.76_11.10_8_PW3",
    "0.76_11.10_11_PW1",
    "0.76_12.14_14_D",
    "0.76_12.14_1_PW1",
    "0.76_13.18_8_PW4",
    "0.76_13.18_13_PW4",
    "2.15_100_8_T",
    "2.15_100_14_PW3",
    "2.15_100_15_PW1",
    "0.76_15.7_1_PW3",
    "0.76_15.70_8_D",
    "0.76_17.17_13_PW3",
    "0.76_17.17_2_PW4",
    "0.76_18.65_1_D",
    "0.76_18.65_10_PW4",
    "1.52_100_10_D",
    "1.52_100_13_T",
    "0.76_22.20_15_PW2",
    "0.76_22.20_10_PW3",
    "0.76_24.29_7_T",
    "0.76_24.29_12_PW2",
    "0.76_26.37_13_D",
    "0.76_26.37_10_T",
    "1.075_100_7_D",
    "1.075_100_7_PW4",
    ]


test_top = True
enable_3D_plot = True
plot_tangent = True

test_file_names = top_file_names if test_top else side_file_names

tests = []

densities: List[float] = []
compression_energies: List[float] = []
energy_diffs: List[float] = []
energy_ratios: List[float] = []
min_secant_moduli: List[float] = []
max_tangent_moduli: List[float] = []

pyplot.figure(0)
pyplot.xlabel('Strain ε (%)')
pyplot.ylabel('Stress σ (MPa)')


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
        color = cm.rainbow((test.limit_density / 100 - .1) / .3)
        tangent_modulus: Tangent = getTangentModulus(test)
        (xs, ys) = getTangentLine(tangent_modulus)
        pyplot.figure(0)
        pyplot.plot(xs, ys, color = PlottingUtil.lighten_color(color, .25))
        if enable_3D_plot:
            ax.plot(xs, [test.limit_density, test.limit_density], ys, color = PlottingUtil.lighten_color(color, .25))

if True:
    for test in tests:

        color = cm.rainbow((test.limit_density / 100 - .1) / .3)

        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]
        strain_decomp = test.strain[test.decompression_range]
        stress_decomp = test.stress[test.decompression_range]

        if enable_3D_plot:
            ax.plot(strain_comp * 100, np.ones(len(test.compression_range)) * test.limit_density, stress_comp * 1000, color = color)

        pyplot.figure(0)
        pyplot.plot(strain_comp * 100, stress_comp * 1000, color = color)

        tangent_modulus: Tangent = getTangentModulus(test)

        max_tangent_moduli.append(tangent_modulus.val)

        compression_energy = abs(integral(strain_comp, stress_comp)[-1])
        decompression_energy = abs(integral(strain_decomp, stress_decomp)[-1])
        densities.append(test.limit_density)
        compression_energies.append(compression_energy)
        energy_diffs.append(compression_energy - decompression_energy)
        energy_ratios.append((compression_energy - decompression_energy) / compression_energy)
        stress_t = stress_comp[100:-1] - stress_comp[0]
        strain_t = strain_comp[100:-1] - strain_comp[0]
        min_secant_moduli.append(min(stress_t / strain_t))


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
pyplot.plot(densities, min_secant_moduli)
pyplot.xlabel('density')
pyplot.ylabel('Minimal secant modulus')

pyplot.figure(6)
pyplot.autoscale(tight=True)
pyplot.scatter(densities, max_tangent_moduli)
pyplot.xlabel('density')
pyplot.ylabel('Max tangent modulus')
'''


pyplot.show()











