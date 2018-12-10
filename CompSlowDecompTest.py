

from DataCutting import DataCutting
from ElastoPlasticDeformationCutter import ElastoPlasticDeformationCutter

import math                 # sqrt
from typing import List     # special types such as List
import numpy as np
import csv

import traceback
import logging

import typing
from typing import Dict

class CompSlowDecompTest:
    my_printer_names: List[str] = [
            'Tims barfer',
            'PW001',
            'PW002',
            'PW003',
            'PW004',
            'Davids printer',
            'TUD'
        ]

    gcode_dimensions: typing.List[float] = [48.67, 48.48, 48.48]
    gcode_area: float = 48.48 * 48.48

    def __init__(self, folder_name: str, file_name_base: str, compression_count: int = 1, use_gcode_dimensions: bool = False):

        try:
            test_case_data_filename: str = folder_name + file_name_base + '_1.csv'
            # test_case_data = np.genfromtxt(test_case_data_file, delimiter=',', skip_header=1, skip_footer=3)
            with open(test_case_data_filename) as test_case_data_file:
                test_variables: Dict[str, str] = {}
                test_case_data_reader = csv.reader(test_case_data_file, delimiter=',')
                vars = next(test_case_data_reader)
                vals = next(test_case_data_reader)

                for var, val in list(zip(vars, vals)):
                    test_variables[var] = val

                try:
                    self.printer_name = test_variables['Sample choice input : Printer origin']
                except KeyError:
                    if 'T' in file_name_base:
                        self.printer_name = CompSlowDecompTest.my_printer_names[0]
                    elif 'PW1' in file_name_base:
                        self.printer_name = CompSlowDecompTest.my_printer_names[1]
                    elif 'PW2' in file_name_base:
                        self.printer_name = CompSlowDecompTest.my_printer_names[2]
                    elif 'PW3' in file_name_base:
                        self.printer_name = CompSlowDecompTest.my_printer_names[3]
                    elif 'PW4' in file_name_base:
                        self.printer_name = CompSlowDecompTest.my_printer_names[4]
                    elif 'D' in file_name_base:
                        self.printer_name = CompSlowDecompTest.my_printer_names[5]
                    else:
                        self.printer_name = CompSlowDecompTest.my_printer_names[6]
                        logging.warning('Printer origin of \'' + file_name_base + '\' defaulted to ' + self.printer_name)
                self.printer_number = CompSlowDecompTest.my_printer_names.index(self.printer_name)
                self.print_cycle_nr = int(test_variables['Sample number inputs : Print cycle number'])
                self.dimensions_before = [
                    float(test_variables['Sample number inputs : Height (mm)']),
                    float(test_variables['Sample number inputs : Side Height (mm)']),
                    float(test_variables['Sample number inputs : Width (mm)'])]
                self.dithering_density_setting = float(test_variables['Sample number inputs : dithering density'])
                self.line_dist_setting = float(test_variables['Sample number inputs : line dist (mm)'])

                self.limit_density = self.dithering_density_setting if self.dithering_density_setting != 100 else 0.38 / (self.line_dist_setting / 100) * (1/3 + 1/3*math.sqrt(2))
        except Exception as e:
            logging.error(traceback.format_exc())
            logging.warning('For file \'' + test_case_data_filename + '\'')
            exit(1)

        try:

            test_data_file: str = folder_name + file_name_base + '.is_ccyclic_Exports/' + file_name_base + '_1.csv'
            data = np.genfromtxt(test_data_file, delimiter=',', skip_header=2)
            self.disp = data[:, 0] # in mm
            self.force = data[:, 1] # in kN
            if np.isnan(np.min(self.disp)):
                raise SyntaxError('Couldn\'t read CSV data: encountered nan.\n Check whether data contains quotation marks etc.')
            cutter: DataCutting = DataCutting(self.disp, compression_count)
            assert(len(cutter.compression_ranges) == 1)

            start_disp_cutoff_index = ElastoPlasticDeformationCutter.getNongrippedDisplacementIndex(self.disp[cutter.compression_ranges[0]],
                                                                                                    self.force[cutter.compression_ranges[0]],
                                                                                                    compression_vs_decomp = True,
                                                                                                    force_cutoff = 0.005)
            self.compression_range = cutter.compression_ranges[0][start_disp_cutoff_index:]

            disp_at_start_cutoff = self.disp[self.compression_range[0]]
            if len(cutter.decompression_ranges) == 0 or len(cutter.decompression_ranges[0]) == 0:
                self.decompression_range = range(self.compression_range[1], self.compression_range[-1])
            else:
                disp_indices_below_cutoff: np.array[int] = np.where(self.disp[cutter.decompression_ranges[0][0]:] < disp_at_start_cutoff)
                end_cutoff_index: int = disp_indices_below_cutoff[0][0]
                self.decompression_range = cutter.decompression_ranges[0][0:end_cutoff_index]

            self.disp -= disp_at_start_cutoff

            mm_to_m = 0.001
            if use_gcode_dimensions:
                self.strain = self.disp / self.gcode_dimensions[0] # in ratio (mm/mm)
                self.stress = self.force / (self.gcode_dimensions[1] * self.gcode_dimensions[2]) / (mm_to_m * mm_to_m * 1000) # in MPa (10^6 * N/m^2 = 10^3 * kN/m^2)
            else:
                self.strain = self.disp / self.dimensions_before[0] # in ratio (mm/mm)
                self.stress = self.force / (self.dimensions_before[1] * self.dimensions_before[2]) / (mm_to_m * mm_to_m * 1000) # in MPa (10^6 * N/m^2 = 10^3 * kN/m^2)
        except Exception as e:
            logging.error(traceback.format_exc())
            logging.warning('For file \'' + test_data_file + '\'')
            exit(1)