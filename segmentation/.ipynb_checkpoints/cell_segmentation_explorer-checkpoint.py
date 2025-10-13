#!/usr/bin/env python
# coding: utf-8

# In[1]:

import tifffile

#import napari
import numpy as np
from skimage.io import imread
import pandas as pd
import os
import glob
'''This script reads the file img data, and  creates the split_channels'''

# In[2]:
# Search for one.tiff files
#change folder to location of Leap 
base_dir = '../../combined_tiff'
file_pattern = '*.tiff'  # Change to '*.tif' if your files have the '.tif' extension

# Create a pattern to search for subdirectories with names starting with 'Leap'
sub_dir_pattern = os.path.join(base_dir, 'Leap*')

# Use glob to find all TIFF files in the subdirectories
tiff_files = glob.glob(os.path.join(sub_dir_pattern, file_pattern), recursive=True)
panel_folders = ['/'.join(el.split('/')[:-1]) for el in tiff_files]# 1 level up tiff files




# In[4]:
for tiff_file,panel_folder in zip(tiff_files,panel_folders):
    #assuming panel file is in Leap***/panel.csv
    panel = pd.read_csv(panel_folder+'/panel.csv',index_col=0)
    panel = panel[panel.keep==1]#select only the keep files
    panel.reset_index(drop = True,inplace=True)
    img = imread(tiff_file)

    # In[5]:
    Leap_name = panel_folder.split('/')[-1]
    img_folder = '/'.join(tiff_file.split('/')[:-1]) # ed '../combined_tiff/Leap006/
    tiff_file_no_path = tiff_file.split('/')[-1].rstrip('.tiff')#Leap006_4
    output_img_folder = '../split_channels/'+tiff_file_no_path+'/'
    if not os.path.exists(output_img_folder):
        #if /img folder does not exist, create it
        os.makedirs(output_img_folder)
    for channel,marker in list(zip(img,panel['name'].values)):
        tifffile.imwrite( output_img_folder+'/'+marker+'.tiff',
                    data=channel,
                    )
'''
viewer = napari.Viewer()


# In[6]:


layer = viewer.add_image(img, name=panel.name,channel_axis=0)
napari.run() 
'''







