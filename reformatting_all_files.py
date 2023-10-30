import numpy as np
import os
from os import PathLike
from pathlib import Path
from readimc import MCDFile, TXTFile
from collections import defaultdict
import pandas as pd
from typing import Any, Dict, Generator, List, Optional, Sequence, Tuple, Union
from tqdm import tqdm
from readimc.data.acquisition import Acquisition, AcquisitionBase
from readimc.data.panorama import Panorama
from readimc.data.slide import Slide
import shutil
from scipy.ndimage import maximum_filter
import tifffile
#import imcsegpipe
#from imcsegpipe.utils import sort_channels_by_mass
import  imcsegpipe._imcsegpipe as imcsegpipe

def preprocess_image(img: np.ndarray, hpf: Optional[float] = None) -> np.ndarray:
    img = img.astype(np.float32)
    if hpf is not None:
        img = filter_hot_pixels(img, hpf)
    return img    
    return io._to_dtype(img, io.img_dtype)
def filter_hot_pixels(img: np.ndarray, thres: float) -> np.ndarray:
    kernel = np.ones((1, 3, 3), dtype=bool)
    kernel[0, 1, 1] = False
    max_neighbor_img = maximum_filter(img, footprint=kernel, mode="mirror")
    return np.where(img - max_neighbor_img > thres, max_neighbor_img, img)
def create_panel_from_acquisition(acquisition):
    '''Take as input the acquisition from the mcd file'''
    panel = pd.DataFrame(
                    data={
                        "channel": pd.Series(
                            data=acquisition.channel_names,
                            dtype=pd.StringDtype(),
                        ),
                        "name": pd.Series(
                            data=acquisition.channel_labels,
                            dtype=pd.StringDtype(),
                        ),
                    },
                )
    name = panel.name.str.split('-',n=1).str.get(1)
    panel['keep'] = ~name.isna()
    panel['name'] = name.where(panel.keep,'-')    
    panel.loc[panel['channel'] == 'Pt195', ['name', 'keep']] = ['Carboplatin', True]#add carboplatin to panel

    return panel

def save_tiff_from_mcd(mcd_file,Leap_folder):
    '''Save in ../combined_tiff and ../split_channels'''
    output_path = '../combined_tiff/'+Leap_folder+'/'#
    
    if not os.path.exists(output_path):
        #if /img folder does not exist, create it
        os.makedirs(output_path)

    with MCDFile(mcd_file) as f:
        for slide in f.slides:
            for acquisition in slide.acquisitions:
                if int(acquisition.metadata['DataEndOffset'])==0:
                    #corrupted acquisition
                    continue
                try:
                    img = f.read_acquisition(acquisition)
                except OSError:
                    continue # if mcd corrupt, move to next acquisition
                panel = create_panel_from_acquisition(acquisition)
                panel.to_csv(output_path.rstrip('/img')+'/panel.csv')
                img = img[panel.keep.values] #select only channels associated to antibody
                img = preprocess_image(img,hpf=50)
                filename = Leap_folder+'_'+acquisition.metadata['ID']
                if np.any(np.array(img.shape[1:])<64):
                    #if image is not at least of 64 pixel per side, ignore it
                    continue
                tifffile.imwrite( 
                output_path+filename+'.tiff',
                data=img[np.newaxis, np.newaxis, :, :, :, np.newaxis],
                imagej=img.dtype in (np.uint8, np.uint16, np.float32),
                )
                #save in  the split_channels
                output_path_split = '../split_channels/'+filename+'/'# eg. ../split_channels/Leap001_10'               
                if not os.path.exists(output_path_split):
                    #if folder does not exist, create it
                    os.makedirs(output_path_split)
                for channel,marker in list(zip(img,panel[panel.keep.values]['name'].values)):
                    tifffile.imwrite( 
                    output_path_split+marker+'.tiff',
                    data=channel
                    )

def main():
    path0 = '../IMC_data'
    mcd_files_list = [str(el) for el in sorted(Path(path0).rglob("[!.]*.mcd"))]
    data_folder = ['/'.join(mcd_file.split('/')[:-1])for mcd_file in mcd_files_list]#take the path up to a level of the mcd file, it should contain the txt files

    tiff_files = [str(el) for el in sorted(Path('../split_channels').rglob("[!.]*.tiff"))]
    Leap_existing_files =        pd.Series(tiff_files).str.lstrip('../split_channels').str.split('/').str[0].str.split('_').str[0].unique()
    acquisition_metadatas = []    
    for mcd_file,main_folder in tqdm(list(zip(mcd_files_list,data_folder))):
        #main_folder is where the mcd files sit
        Leap_folder = mcd_file.replace(path0,'.').split('/',maxsplit = 2)[1] 
        mcd_folder = '/'.join(mcd_file.split('/')[:-1])
        if mcd_file.lstrip(path0).split('/')[0] in Leap_existing_files:
            print('skipping'+mcd_file)
            continue
        #load metadata        
        acquisition_metadata = imcsegpipe.extract_mcd_file(mcd_file,acquisition_dir= mcd_folder)
        #extract_mcd_file
        acquisition_metadatas.append(acquisition_metadata)
        save_tiff_from_mcd(mcd_file,Leap_folder)
    acquisition_metadata = pd.concat(acquisition_metadatas, copy=False)
    acquisition_metadata.to_csv(path0 +"/acquisition_metadata.csv",mode = 'a')
if __name__ == "__main__":
    main()