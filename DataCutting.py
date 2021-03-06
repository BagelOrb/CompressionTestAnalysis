
from typing import List

import numpy as np


class DataCutting(object):

    def _separateData(self, disp: np.array, cycle_count: int) -> None:
        max_indices: List[int] = [0]  # indices where disp is maximal
        min_indices: List[int] = [0]  # indices where disp is minimal

        approx_cycle_size = int(disp.shape[0] / cycle_count)
        for compression_nr in range(0, cycle_count):
            start_idx = min_indices[-1]
            end_idx = min(disp.shape[0] - 1, start_idx + approx_cycle_size)
            max_idx = start_idx + np.argmax(disp[start_idx: end_idx])
            max_indices.append(max_idx)
            self.compression_ranges.append(range(start_idx, max_idx))

            start_idx = max_indices[-1]
            end_idx = min(disp.shape[0] - 1, start_idx + approx_cycle_size)
            min_idx = start_idx + np.argmin(disp[start_idx: end_idx])
            min_indices.append(min_idx)
            self.decompression_ranges.append(range(start_idx, min_idx))
            #print('max: ' + str(max_indices[-1]))
            #print('min: ' + str(min_indices[-1]))


    def __init__(self, disp: np.array, cycle_count: int):
        self.compression_ranges: List[range] = []
        self.decompression_ranges: List[range] = []
        self._disp = disp
        self._separateData(disp, cycle_count)