
# coding: utf-8

# In[113]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import cm


# In[114]:


import csv
results = []
with open('test', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        results.append(row)


# In[117]:


import numpy as np
X = [i for i in np.linspace(1,32,32)]
Y = [i for i in np.linspace(1,32,32)]
X, Y = np.meshgrid(X, Y)

Z = np.array(results)
Z = Z.astype(np.float)


# In[118]:


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap = cm.coolwarm, linewidth = 0)
plt.show()

