
import numpy as np

from typing import List

from matplotlib import pyplot

class ElastoPlasticDeformationCutter(object):

    ## get the displacement when the object is just/still fully gripped
    # when the 2nd order derivative is zero
    @staticmethod
    def getNongrippedDisplacementIndex(disp: np.array, force: np.array, compression_vs_decomp: bool, max_nongripped_displacement: float = 10) -> int:
        out_of_disp_range_indices = np.where(disp > max_nongripped_displacement)
        max_non_gripped_disp_idx = out_of_disp_range_indices[0][0] if compression_vs_decomp else out_of_disp_range_indices[0][-1]
        window = range(0, max_non_gripped_disp_idx) if compression_vs_decomp else range(max_non_gripped_disp_idx, disp.size)  # the window within which we are looking for the start/end of gripping

        fitted_polynomial_coefficients = np.polyfit(disp, force, 21)
        der_coeff = np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], 2)[::-1]  # polynomial lib handles coefficients in reverse order)
        roots = np.polynomial.polynomial.polyroots(der_coeff[::-1])

        pyplot.plot(disp[window], force[window])
        disp_range = np.arange(0,disp[window[-1]], .1)
        pyplot.plot(disp_range, np.poly1d(fitted_polynomial_coefficients)(disp_range))
        pyplot.show()

        # get where the second order derivative is zero: it's roots
        real_roots: List[float] = []  # roots which only have a rel part and no imaginary part
        for root in roots:
            if (root.imag == 0 and root.real > 0):  # check if root is a proper complex number
                real_roots.append(root.real)

        # get the data point closest to the root
        non_gripped_displacement = min(real_roots)
        best_disp_dist = float('inf')
        best_index = -1
        for i in range(disp.size):
            dist = abs(disp[i] - disp[best_index])
            if dist < best_disp_dist:
                best_disp_dist = dist
                best_index = i

        return best_index
