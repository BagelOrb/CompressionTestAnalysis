

from DataCutting import DataCutting
from ElastoPlasticDeformationCutter import ElastoPlasticDeformationCutter

import math         # sqrt
import typing       # special types such as List
import numpy as np
import csv

import typing
from typing import Dict

class CompSlowDecompTest:
    compression_count: int = 1

    gcode_dimensions: typing.List[float] = [48.67, 48.48, 48.48]
    gcode_area: float = 48.48 * 48.48

    def __init__(self, folder_name: str, file_name_base: str):
        test_case_data_filename: str = folder_name + file_name_base + '_1.csv'
        # test_case_data = np.genfromtxt(test_case_data_file, delimiter=',', skip_header=1, skip_footer=3)

        test_variables: Dict[str, str] = {}
        with open(test_case_data_filename) as test_case_data_file:
            test_case_data_reader = csv.reader(test_case_data_file, delimiter=',')
            vars = next(test_case_data_reader)
            vals = next(test_case_data_reader)

            for var, val in list(zip(vars, vals)):
                test_variables[var] = val

            self.print_cycle_nr = int(test_variables['Sample number inputs : Print cycle number'])
            self.dimensions_before = [
                float(test_variables['Sample number inputs : Height (mm)']),
                float(test_variables['Sample number inputs : Side Height (mm)']),
                float(test_variables['Sample number inputs : Width (mm)'])]
            self.dithering_density_setting = float(test_variables['Sample number inputs : dithering density'])
            self.line_dist_setting = float(test_variables['Sample number inputs : line dist (mm)'])

            self.limit_density = self.dithering_density_setting if self.dithering_density_setting != 100 else 0.38 / (self.line_dist_setting / 100) * (1/3 + 1/3*math.sqrt(2))

        test_data_file: str = folder_name + file_name_base + '.is_ccyclic_Exports/' + file_name_base + '_1.csv'
        data = np.genfromtxt(test_data_file, delimiter=',', skip_header=2)
        self.disp = data[:, 0]
        self.force = data[:, 1]
        if np.isnan(np.min(self.disp)):
            raise SyntaxError('Couldn\'t read CSV data: encountered nan.\n Check whether data contains quotation marks etc.')
        cutter: DataCutting = DataCutting(self.disp, CompSlowDecompTest.compression_count)
        assert(len(cutter.compression_ranges) == 1)

        start_disp_cutoff_index = ElastoPlasticDeformationCutter.getNongrippedDisplacementIndex(self.disp[cutter.compression_ranges[0]],
                                                                                                self.force[cutter.compression_ranges[0]],
                                                                                                compression_vs_decomp = True,
                                                                                                force_cutoff = 0.005)
        self.compression_range = cutter.compression_ranges[0][start_disp_cutoff_index:-1]

        disp_at_start_cutoff = self.disp[self.compression_range[0]]
        disp_indices_below_cutoff: np.array[int] = np.where(self.disp[cutter.decompression_ranges[0][0]:-1] < disp_at_start_cutoff)
        end_cutoff_index: int = disp_indices_below_cutoff[0][0]
        self.decompression_range = cutter.decompression_ranges[0][0:end_cutoff_index]

        self.disp -= disp_at_start_cutoff

        self.strain = self.disp / self.gcode_dimensions[0]
        self.stress = self.force / self.gcode_area * 1000