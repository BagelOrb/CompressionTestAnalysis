
import numpy as np


from matplotlib import pyplot


def fitPolynomial(disp: np.array, force: np.array) -> None:
    max_disp: float = max(disp)
    fitted_polynomial_coefficients = np.polyfit(disp, force, 21)
    fitted_polynomial = np.poly1d(fitted_polynomial_coefficients)
    #print(fitted_polynomial)
    #print(fitted_polynomial_coefficients)

    der_coeff = np.polynomial.polynomial.polyder(fitted_polynomial_coefficients[::-1], 2)[::-1]  # polynomial lib handles coefficients in reverse order
    der = np.poly1d(der_coeff)

    x = np.arange(0, max_disp, max_disp / 1000)
    pyplot.plot(x, der(x))
    pyplot.plot(x, fitted_polynomial(x))

    roots = np.polynomial.polynomial.polyroots(der_coeff[::-1])
    for root in roots:
        if (root.imag == 0 and root.real > 0):  # check if root is a proper complex number
            pyplot.axvline(root.real)
            print(root.real)

    pyplot.plot(disp, force)
    pyplot.show()
