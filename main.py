
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')

import PolyFitting
from DataCutting import DataCutting

from ElastoPlasticDeformationCutter import ElastoPlasticDeformationCutter

import numpy as np


from matplotlib import pyplot


compression_count = 5



def addFileToPlot(filename: str, compression_vs_decomp: bool, colorspec: str) -> None:
    data = np.genfromtxt(filename, delimiter=',', skip_header=2)
    disp = data[:,0]
    force = data[:,1]
    if np.isnan(np.min(disp)):
        raise SyntaxError('Couldn\'t read CSV data: encountered nan.\n Check whether data contains quotation marks etc.')
    cutter: DataCutting = DataCutting.separateData(disp, compression_count)
    for i in {1}: #range(compression_count):
        data_range = cutter.compression_ranges[i] if compression_vs_decomp else cutter.decompression_ranges[i]
        disp_here = disp[data_range]
        force_here = force[data_range]
        pyplot.axvline(disp[ElastoPlasticDeformationCutter.getNongrippedDisplacementIndex(disp_here, force_here, compression_vs_decomp)])
        pyplot.plot(disp_here, force_here, colorspec)


def compare8top(compression_vs_decomp: bool = True):
    addFileToPlot('preliminaries/results/8.4_top.is_ccyclic_Exports/8.4_top_1.csv', compression_vs_decomp, 'g')
    addFileToPlot('preliminaries/results/8.3_top.is_ccyclic_Exports/8.3_top_1.csv', compression_vs_decomp, 'b')
    #addFileToPlot('preliminaries/results/5_top.is_ccyclic_Exports/5_top_1.csv', compression_vs_decomp, 'r')
    pyplot.show()


def compare8side(compression_vs_decomp: bool = True):
    addFileToPlot('preliminaries/results/8.4_side.is_ccyclic_Exports/8.4_side_1.csv', compression_vs_decomp, 'g')
    addFileToPlot('preliminaries/results/8.3_side.is_ccyclic_Exports/8.3_side_1.csv', compression_vs_decomp, 'b')
    #addFileToPlot('preliminaries/results/5_side.is_ccyclic_Exports/5_side_1.csv', compression_vs_decomp, 'r')
    pyplot.show()

def compare9side(compression_vs_decomp: bool = True):
    addFileToPlot('preliminaries/results/9_side.is_ccyclic_Exports/9_side_1.csv', compression_vs_decomp, 'g')
    addFileToPlot('preliminaries/results/7_side.is_ccyclic_Exports/7_side_1.csv', compression_vs_decomp, 'b')
    pyplot.show()


#compare8top()
#compare8side(False)


addFileToPlot('preliminaries/results/8.4_top.is_ccyclic_Exports/8.4_top_1.csv', True, 'g')
#addFileToPlot('preliminaries/results/8.4_top.is_ccyclic_Exports/8.4_top_1.csv', False, 'b')
pyplot.show()

'''
data = np.genfromtxt('preliminaries/results/8.4_top.is_ccyclic_Exports/8.4_top_1.csv', delimiter=',', skip_header=2)
disp = data[:,0]
force = data[:,1]
if np.isnan(np.min(disp)):
    raise SyntaxError('Couldn\'t read CSV data: encountered nan.\n Check whether data contains quotation marks etc.')
cutter: DataCutting = DataCutting.separateData(disp, compression_count)
PolyFitting.fitPolynomial(disp[cutter.compression_ranges[0]], force[cutter.compression_ranges[0]])
'''















