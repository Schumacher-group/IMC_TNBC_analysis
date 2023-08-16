#!/usr/bin/env python
# coding: utf-8

# In[1]:


import napari
import numpy as np
from skimage.io import imread
import pandas as pd
'''This script reads the file img data, and shows in napari with the protein tag'''

# In[2]:
#change folder to location of Leap 

folder = '../IMC_SegmentationResults/Leap001/'
img = imread(folder+'img/Leap001_008.tiff')


# In[4]:

#assuming panel file is in Leap***/panel.csv
panel = pd.read_csv(folder+'panel.csv')
panel = panel[panel.keep==1]#select only the keep files
panel.name = panel.name.str.split('-',n=1).str.get(1)
panel.reset_index(drop = True,inplace=True)


# In[5]:


viewer = napari.Viewer()


# In[6]:


layer = viewer.add_image(img, name=panel.name,channel_axis=0)
napari.run() 








