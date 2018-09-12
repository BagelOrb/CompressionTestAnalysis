
import os
os.chdir('/home/t.kuipers/Documents/PhD/Fractal Dithering project/experiments')


import pickle

from matplotlib import pyplot

fig_0 = pickle.load(open('analysis/results/top_3D.pickle', 'rb'))
fig_1 = pickle.load(open('analysis/results/side_3D.pickle', 'rb'))

pyplot.show()