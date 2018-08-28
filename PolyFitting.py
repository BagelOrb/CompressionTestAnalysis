
import numpy as np

from typing import List

from matplotlib import pyplot


def fitPolynomial(disp: np.array, force: np.array) -> None:
    max_disp: float = max(disp)
    fitted_polynomial_coefficients = np.polyfit(disp, force, 21)
    print(disp.size)
    fitted_polynomial = np.poly1d(fitted_polynomial_coefficients)
    #print(fitted_polynomial)
    #print(fitted_polynomial_coefficients)

    order = 2 # the number of time to differentiate
    der_coeff = np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], order)[::-1]  # polynomial lib handles coefficients in reverse order
    der = np.poly1d(der_coeff)

    der_1 = np.poly1d(np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], order - 1)[::-1])

    x = np.arange(0, max_disp, max_disp / 1000)
    pyplot.plot(x, der_1(x) / max(der_1(x)) * max(force))
    pyplot.plot(disp, force)
    pyplot.plot(x, fitted_polynomial(x))

    roots: np.ndarray = np.polynomial.polynomial.polyroots(der_coeff[::-1])
    for root in roots:
        if root.imag == 0 and root.real > min(disp) and root.real < max(disp):  # check if root is a proper complex number and within the range of the data
            pyplot.axvline(root.real)
            print('root: ' + str(root.real))

    pyplot.plot(disp, force)

def getDerivative(x_data: np.array, y_data: np.array, derivate_nr: int, polynomial_order: int) -> np.lib.polynomial.poly1d:
    fitted_polynomial_coefficients = np.polyfit(x_data, y_data, polynomial_order)
    der_coeff = np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], derivate_nr)[::-1]  # polynomial lib handles coefficients in reverse order
    der = np.poly1d(der_coeff)
    return der


##
# get the positions where a given order derivative has a maximum or minimum
def getRoots(x_data: np.array, y_data: np.array, derivate_nr: int, polynomial_order: int) -> List[int]:
    fitted_polynomial_coefficients = np.polyfit(x_data, y_data, polynomial_order)
    der_coeff = np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], derivate_nr)[::-1]  # polynomial lib handles coefficients in reverse order

    real_roots: List[int] = []
    roots: np.ndarray = np.polynomial.polynomial.polyroots(der_coeff[::-1])
    for root in roots:
        if root.imag == 0 and root.real > min(x_data) and root.real < max(x_data):  # check if root is a proper complex number and within the range of the data
            disp_indices = np.where(x_data > root.real)[0]
            if len(disp_indices) > 0:
                real_roots.append(disp_indices[0])

    return real_roots

##
# Get the index of the datum of the first peak, or None if there is no first peak
def findFirstTop(disp: np.array, force: np.array, min_force: float = 0.02) -> int:
    max_force: float = 0.9 * np.max(force)
    fitted_polynomial_coefficients = np.polyfit(disp, force, 21)
    fitted_polynomial = np.poly1d(fitted_polynomial_coefficients)
    der_coeff = np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], 1)[::-1]  # polynomial lib handles coefficients in reverse order
    roots = np.polynomial.polynomial.polyroots(der_coeff[::-1])

    for root in roots:
        if root.imag == 0 and root.real > 0 and fitted_polynomial(root.real) > min_force and fitted_polynomial(root.real) < max_force:  # check if root is a proper complex number
            indices_beyond_top = np.where(disp > root.real)
            return indices_beyond_top[0][0]

    return None
