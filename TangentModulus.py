
import PolyFitting
import CompSlowDecompTest

# from __future__ import annotations # forward declaration of type >> only python 3.7+

import numpy as np

class TangentModulus:
    def __init__(self, source_x: float, source_y: float, val: float):
        self.source_x = source_x
        self.source_y = source_y
        self.val = val

    def getTangentModulus(test: CompSlowDecompTest, polynomial_order: int = 10, overall_polynomial_order: int = 25) -> 'TangentModulus':
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]
        derivative_optima_indices = PolyFitting.getRoots(strain_comp, stress_comp, 2, overall_polynomial_order)  # where the derivative is maximal or minimal

        max_data_index = strain_comp.size - 1

        # determine the data on which to fit a polynomial
        der_3 = PolyFitting.getDerivative(strain_comp, stress_comp, 3, overall_polynomial_order)
        derivative_optima = strain_comp[derivative_optima_indices]
        mapping_der_3 = np.vectorize(der_3)  # so that we can map the function over an np.array
        derivative_optimum_direction = mapping_der_3(derivative_optima)  # positive for minima, negative for maxima
        der_minima_indices = np.where(np.logical_and(derivative_optimum_direction > 0, derivative_optima > 0.05))[0]
        if len(der_minima_indices) > 0:
            max_data_index = derivative_optima_indices[der_minima_indices[0]]

        # fit curve to first part of the graph and get the max derivative from that
        strain_comp_sub = strain_comp[0:max_data_index]
        stress_comp_sub = stress_comp[0:max_data_index]
        derivative_optima_indices_sub = PolyFitting.getRoots(strain_comp_sub, stress_comp_sub, 2, polynomial_order)  # where the derivative is maximal or minimal
        der_sub = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 1, polynomial_order)
        der_3_sub = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 3, polynomial_order)
        first_optimum_is_max = der_3_sub(strain_comp_sub[derivative_optima_indices_sub[0]]) < 0  # to determine whether it's a top or a bottom peek in the derivative
        tangent: float = 0
        source_x: float = 0
        if first_optimum_is_max:
            strain_at_derivative_max = strain_comp_sub[derivative_optima_indices_sub[0] - 1]
            assert (strain_at_derivative_max > 0)
            tangent = der_sub(strain_at_derivative_max)
            source_x = strain_at_derivative_max
        else:
            tangent = der_sub(0)
            source_x = 0
        assert (tangent > 0)
        polynomial = PolyFitting.getDerivative(strain_comp_sub, stress_comp_sub, 0, polynomial_order)
        source_y = polynomial(source_x)

        #
        # pyplot.clf()
        # pyplot.plot(strain_comp, stress_comp)
        # pyplot.plot(strain_comp_sub, stress_comp_sub)
        # pyplot.plot(strain_comp_sub, polynomial(strain_comp_sub))
        # pyplot.show()
        # exit(1)

        return TangentModulus(source_x, source_y, tangent)

    # get the tangent modulus at the end of the graph
    def getEndTangentModulus(test: CompSlowDecompTest, polynomial_order: int = 10, overall_polynomial_order: int = 25) -> 'TangentModulus':
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]
        derivative_optima_indices = PolyFitting.getRoots(strain_comp, stress_comp, 2, overall_polynomial_order)  # where the derivative is maximal or minimal

        min_data_index = 0

        # determine the data on which to fit a polynomial
        der_3 = PolyFitting.getDerivative(strain_comp, stress_comp, 3, overall_polynomial_order)
        mapping_der_3 = np.vectorize(der_3)
        derivative_optima = strain_comp[derivative_optima_indices]
        derivative_optimum_direction = mapping_der_3(derivative_optima)  # positive for minima, negative for maxima
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

        return TangentModulus(source_x, source_y, tangent)

    def getTangentLine(tangent: 'TangentModulus') -> (np.array, np.array):
        tangent_max = np.array([.4, .4 * tangent.val])
        if tangent_max[1] > 0.45:
            tangent_max = tangent_max * 0.45 / (.4 * tangent.val)
        xs = np.array([0, tangent_max[0]]) + tangent.source_x
        ys = np.array([0, tangent_max[1]]) + tangent.source_y
        return (xs, ys)
