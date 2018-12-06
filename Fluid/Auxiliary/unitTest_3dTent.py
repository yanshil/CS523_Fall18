
# coding: utf-8

# In[9]:


dim1 = 128
dim2 = 128
dim3 = 128


# In[10]:


import pandas    
df = pandas.read_csv('test', header=None)
df = df.values
df.size  # 128 × 128 = 16384


# In[11]:


import numpy as np
X = range(dim1)
Y = range(dim2)
Z = range(dim3)
X, Y, Z = np.meshgrid(X, Y, Z, sparse=False)
P = df
P = P.reshape(dim1, dim2, dim3)


# In[12]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface((X[:,:,1]), (Y[:,:,1]), (P[:,:,59]))


# In[ ]:


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, P)
