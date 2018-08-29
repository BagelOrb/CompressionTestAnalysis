
import numpy as np
from scipy import signal    # savgol_filter


from CompSlowDecompTest import CompSlowDecompTest
import MathUtils

class Plateau:
    def __init__(self, test: CompSlowDecompTest, cutoff_modulus: float = 0.4):
        strain_comp = test.strain[test.compression_range]
        stress_comp = test.stress[test.compression_range]

        (strain, stress) = signal.savgol_filter((strain_comp, stress_comp), window_length = 51, polyorder = 3)

        ders = MathUtils.derivative(strain, stress)

        plateaud_indices = np.where(np.logical_and(ders < cutoff_modulus, strain_comp[:-1] > 0.05))[0]

        if plateaud_indices.size > 0:
            self.start_idx = plateaud_indices[0]
            self.end_idx = plateaud_indices[-1]
            self.stress = np.mean(stress[self.start_idx:self.end_idx])
        else:
            # there is no plateau, so we return the point with lowest derivative
            strain_start_idx = np.where(strain > 0.05)[0][0] # assume that there are location beyond that strain
            min_derivative = np.min(ders[strain_start_idx:-1])
            min_der_idx = np.where(ders == min_derivative)[0][0] # assume we can find the minimum
            self.start_idx = min_der_idx
            self.end_idx = min_der_idx
            self.stress = stress_comp[min_der_idx]

        self.strain_start = strain_comp[self.start_idx]
        self.strain_end = strain_comp[self.end_idx]
