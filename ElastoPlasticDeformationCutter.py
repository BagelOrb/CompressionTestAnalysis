
import numpy as np

from typing import List

from matplotlib import pyplot

class ElastoPlasticDeformationCutter(object):

    ## get the displacement when the object is just/still fully gripped
    # when the 2nd order derivative is zero
    @staticmethod
    def getNongrippedDisplacementIndex(disp: np.array, force: np.array, compression_vs_decomp: bool, force_cutoff: float = .01) -> int:
        out_of_force_range_indices = np.where(force > force_cutoff)
        return out_of_force_range_indices[0][0] if compression_vs_decomp else out_of_force_range_indices[0][-1]
