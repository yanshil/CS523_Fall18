import pandas    
df = pandas.read_csv('test', header=None)
df = df.values
df.size

import numpy as np

X = range(16)
Y = range(16)
X, Y = np.meshgrid(X, Y, sparse=False)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
Z= df
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)

import pylab
pylab.savefig('foo.png')