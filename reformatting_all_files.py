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
import re
from scipy.ndimage import maximum_filter
import tifffile
import logging
#import imcsegpipe
#from imcsegpipe.utils import sort_channels_by_mass
import  imcsegpipe._imcsegpipe as imcsegpipe

def preprocess_image(img: np.ndarray, hpf: Optional[float] = None) -> np.ndarray:
    img = img.astype(np.float32)
    if hpf is not None:
        img = filter_hot_pixels(img, hpf)
    return img    

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
'''
def save_tiff_from_mcd(mcd_file,Leap_folder):
    Save in ../combined_tiff and ../split_channels
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
                #img = preprocess_image(img,hpf=50)
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
                output_path_split = '../split_channels_nohpf/'+filename+'/'# eg. ../split_channels/Leap001_10'               
                if not os.path.exists(output_path_split):
                    #if folder does not exist, create it
                    os.makedirs(output_path_split)
                for channel,marker in list(zip(img,panel[panel.keep.values]['name'].values)):
                    tifffile.imwrite( 
                    output_path_split+marker+'.tiff',
                    data=channel
                    )
    '''
def mcd_2_ome_tiff():
    '''This function extracts anything there is to extract from the mcd files and generate a metadata csv file.'''
    reprocess_mcd = True
    path0 = '../IMC_data'
    mcd_files_list = [str(el) for el in sorted(Path(path0).rglob("[!.]*.mcd"))]
    mcd_files_list = [str(el) for el in sorted(Path(path0).rglob("[!.]*.mcd"))]
    mcd_files_list = pd.Series(mcd_files_list)
    mcd_files_list = mcd_files_list[~mcd_files_list.str.contains(r"Large Pano MCD files|PanoramasCRUCKCI")]# these are all empty without acquisitions


    biobank = pd.read_excel(os.path.join(path0,'ExtraDocs','biobank list.xlsx'))
    biobank.rename({'BIOBANK ID':'BIOBANK_ID'},inplace=True,axis = 1)
    biobank.dropna(axis = 0,inplace=True)
    biobank['code'] = biobank.BIOBANK_ID.str.split('-').str[0]
    code_2_Leap = biobank[['LEAP ID','code']].drop_duplicates().set_index('code')

    data_folder = ['/'.join(mcd_file.split('/')[:-1])for mcd_file in mcd_files_list]#take the path up to a level of the mcd file, it should contain the txt files

    tiff_files = [str(el) for el in sorted(Path('../split_channels').rglob("[!.]*.tiff"))]
    if reprocess_mcd:
        Leap_existing_files = []

    else:
        Leap_existing_files =        pd.Series(tiff_files).str.lstrip('../split_channels').str.split('/').str[0].str.split('_').str[0].unique()

    acquisition_metadatas = []
    file_saved = []
    for mcd_file,main_folder in tqdm(list(zip(mcd_files_list,data_folder))):
        #main_folder is where the mcd files sit
        Leap_folder = mcd_file.replace(path0,'.').split('/',maxsplit = 2)[1] 
        mcd_folder = '/'.join(mcd_file.split('/')[:-1])
        if mcd_file.lstrip(path0).split('/')[0] in Leap_existing_files:
            print('skipping'+mcd_file)
            continue
        #load metadata        
        try:
            acquisition_metadata = imcsegpipe.extract_mcd_file(mcd_file,acquisition_dir= os.path.join(path0,Leap_folder))#acquisition_dir is the path to the folder: Leap_ID
            #extract_mcd_file
            acquisition_metadatas.append(acquisition_metadata)
            
            file_saved+=(list(Leap_folder +'_'+acquisition_metadata.id.astype(str).values))
        except OSError:
            continue
    acquisition_metadata = pd.concat(acquisition_metadatas, copy=False)
    acquisition_metadata.to_csv(path0 +"/acquisition_metadata.csv",mode = 'w')

def ome_tiff_2_tiff():
    '''Extracts and saves the tiff files. It also correct the names whenever appropriate'''
    tqdm.pandas()
    path0 = '../IMC_data'
    biobank = pd.read_excel(os.path.join(path0,'ExtraDocs','biobank list.xlsx'))
    biobank.rename({'BIOBANK ID':'BIOBANK_ID'},inplace=True,axis = 1)
    biobank.dropna(axis = 0,inplace=True)
    biobank['code'] = biobank.BIOBANK_ID.str.split('-').str[0]
    code_2_Leap = biobank[['LEAP ID','code']].drop_duplicates().set_index('code')

    ome_tiff_paths = list(Path(path0).rglob("acquisition/[!.]*.ome.tiff"))
    path_tb = pd.DataFrame(ome_tiff_paths,columns = ['path'])
    path_tb['filename'] = path_tb['path'].apply(lambda x:x.name).str.rstrip('_ac.ome.tiff')
    path_tb['root'] = path_tb.path.astype(str).str.lstrip('../IMC_data/').str.split('/acquisition/').str[0]
    new_names = []
    metadata = pd.read_csv(path0+'/acquisition_metadata.csv')
    metadata ['root'] = metadata.source_path.str.lstrip('/home/giuseppe/devices/Delta_Tissue/IMC/IMC_analysis../IMC_data').str.split('/').str[0]
    metadata['description'] = metadata['description'].str.replace('170469269','17046926')#correct code for Leap54
    panel_file = path_tb.iloc[0].path
    panel_file = panel_file.with_name(panel_file.name[:-9]+'.csv')
    panel = pd.read_csv(panel_file)
    panel['marker'] = panel[panel.channel_label.str.contains('-')].channel_label.str.split('-',n = 1).str[1]
    panel.loc[panel['channel_name'] == 'Pt195','marker']='Carboplatin'
    for _,file_row in tqdm(path_tb.iterrows()):
        #old_name = row['AcSession']+'_'+str(row['id'])
        # Pattern to match the specific format
        pattern = r"(.*)_s(\d)_a(\d{1,2})"
        
        # Find all occurrences of the pattern
        Leap_ID,_,acquisition_id = np.squeeze(re.findall(pattern, file_row.filename))
        row = metadata[(metadata.root==file_row.root)&(metadata.id==int(acquisition_id))&(metadata.AcSession == Leap_ID)]
        if len(row)>1:
            warnings.warn('Multiple matches')
        row = row.iloc[0]
        description = row.description.lstrip('ROI_')
        if any( x in description.upper() for x in ['LASER','TEST']):
            # this is a test file, ignore acquisition 
            continue

        if '_' in description:
            code = description.split('_')[1]
            if any([ el in description for el in ['Leap067', 'Leap068']]):
                code = description.split('_')[0]  
            
            if Leap_ID == 'Leap009_010_011':
                leap_9_10_11_mapper = {'19005858':'LEAP009','19005859':'LEAP011','19005860':'LEAP010'}#the sample id in the description is wrong
                Leap_ID = leap_9_10_11_mapper[code].capitalize()
            else:
                if re.match("^\d{8}$",code):
                    #it is a 8 digits, looks like the code we want to use from biobank
                    try:
                        Leap_ID = code_2_Leap.loc[code]['LEAP ID'].capitalize()
                    except KeyError:
                        warning.warn('Potential problem, skipping')
                        pass
        if Leap_ID == 'Leap015_016':
            if 'TOP' in description:
                Leap_ID = 'Leap015'
            else:
                Leap_ID = 'Leap016'
        if Leap_ID == 'Leap017_018':
            if 'TOP' in description:
                Leap_ID = 'Leap017'
            else:
                Leap_ID = 'Leap018'
        if Leap_ID == 'Slide 42_MK_ROI':
            Leap_ID,row['id'] = description.split('_')
        if Leap_ID =='Leap091':
            #Claudia said that is to be removed
            continue
        if (Leap_ID =='Leap091') and (str(row['id']) in ['8','9','10']):
            Leap_ID ='Leap092'#Claudia mention on the Whatsapp this is to be renamed
        new_name = Leap_ID+'_'+str(row['id'])
        new_names+=[new_name]

        img = tifffile.imread(file_row.path)
        if np.any(np.array(img.shape[1:])<128):
            #if image is not at least of 128 pixel per side, ignore it
            continue
        img= img[~panel.marker.isna()]
        output_path = '../combined_tiff/'+Leap_ID+'/'#
        if not os.path.exists(output_path):
            #if folder does not exist, create it
            os.makedirs(output_path)

        tifffile.imwrite( 
        output_path+new_name+'.tiff',
        data=img[np.newaxis, np.newaxis, :, :, :, np.newaxis],
        imagej=img.dtype in (np.uint8, np.uint16, np.float32),
        )
        #save in  the split_channels
        output_path_split = '../split_channels_nohpf/'+new_name+'/'# eg. ../split_channels/Leap001_10'               
        if not os.path.exists(output_path_split):
            #if folder does not exist, create it
            os.makedirs(output_path_split)
        for channel,marker in list(zip(img,panel.dropna(axis = 0).marker.values)):
            tifffile.imwrite( 
            output_path_split+marker+'.tiff',
            data=channel
            )
        panel.to_csv(path0+'/panel.csv')
        #rename files and correct according to Leor table
def rename_leap_id():
    def leap3_4():
        swap_from_leap4_to_leap3 =[14,15,16]
        swap_from_leap3_to_leap4 =[11,12,13]
        #old_name:new_name
        file_2_rename = {'Leap003_'+str(id):'Leap004_'+str(id) for id in swap_from_leap3_to_leap4}|{'Leap004_'+str(id):'Leap003_'+str(id) for id in swap_from_leap4_to_leap3}
        files_2_delete = []    
        return files_2_delete,file_2_rename

    def leap17_18():
        files_2_delete = ['Leap018_'+str(x)for x in range(14,25)] 
        return files_2_delete,{}
    def leap9_11():
        files_2_delete = ['Leap009_'+x for x in ['1','2','3','4','5']]
        #swap_from_leap11_to_leap10 =['6','8','9','10','11','12','13']
        #swap_from_leap10_to_leap11 =['1','2','3','4','5']
        #file_2_rename = {'Leap011_'+str(id):'Leap010_'+str(id) for id in swap_from_leap11_to_leap10}|{'Leap010_'+str(id):'Leap011_'+str(id) for id in swap_from_leap10_to_leap11}
        return files_2_delete,{}
    def leap24_25():
        files_2_delete = ['Leap024_'+str(x) for x in range(1,7)]+['Leap025_7']
        return files_2_delete,{}    
    def leap40():
        files_2_delete = ['Leap040_'+str(x) for x in [2,3,4,5,7]]
        return files_2_delete,{}
    def leap91_92():
        file_2_rename = {'Leap092_'+str(id):'Leap091_'+str(id) for id in [8,9,10]}
        return [],file_2_rename
    def aggregate_changes():
        functions = [leap3_4,leap17_18,leap9_11,leap24_25,leap40,leap91_92]
        files_2_rename = {}
        files_2_delete = []
        for f in functions:
            a,b = f()
            files_2_delete+=a
            files_2_rename|=b
        return files_2_delete,files_2_rename

    files_2_delete,files_2_rename = aggregate_changes()
    output_path = '../split_channels_nohpf/'
    for file in files_2_delete:
        try:
            shutil.rmtree(output_path+file)
        except FileNotFoundError:
            logging.warn(output_path+file+' not found')
    for old,new in files_2_rename.items():
        try:
            os.rename(output_path+old,output_path+new)
        except:
            logging.warn(output_path+old+' not found')
def main():
    mcd_2_ome_tiff()
    ome_tiff_2_tiff()
    rename_leap_id()
if __name__=='main':
    main()