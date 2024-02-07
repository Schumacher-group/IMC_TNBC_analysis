import numpy as np
import pandas as pd
#import umap
import matplotlib.pyplot as plt
import scanpy as sc
import pickle 
import squidpy as sq

#from muon import prot as pt
from sklearn import preprocessing
import os
import sys
sys.path.append('../phenotyping/')
from collections import defaultdict
from tqdm import tqdm
import seaborn as sns
input_dir = 'pre_processed_files/'
output_fig = '../figures'
def strip_extension(self):
    '''strip the file extension from a pd.Series of file names'''
    return self.apply(lambda x: os.path.splitext(x)[0])

from phenotyping_utils import generate_anndata_from_ark_analysis
# Add the method to the Pandas Series class
pd.Series.strip_extension = strip_extension

adata = generate_anndata_from_ark_analysis(cell_table_path='../../segmentation/cell_table_Denoised/cell_table_size_normalized_cell_labels.csv',biosamples_path='../../IMC_data/ExtraDocs/processed_response.csv')
dic = defaultdict(list)


types_to_inspect = np.array(adata.obs['Pixie'].cat.categories)
def distance_dictionary(key0,adata,types_to_inspect=None,by_response_status = True):
    interval = np.linspace(0,1200,30)
    def make_dic(adata,acquisition_list,key0,types_to_inspect):
        dic =defaultdict(list)
        for ID in  tqdm(acquisition_list):
            sub_adata = adata[adata.obs.acquisition_ID==ID].copy()
            a = np.array(sub_adata.obs['Pixie'].cat.categories)#label list

            if ((key0 in a)& np.isin(types_to_inspect,a).any()): 
                #if there is at least a cell to inspect
                sq.gr.co_occurrence(sub_adata, cluster_key="Pixie",n_splits=1,interval=interval,show_progress_bar=False)
                x = sub_adata.uns['Pixie_co_occurrence']['interval'][1:]
                cooccur = np.squeeze(sub_adata.uns['Pixie_co_occurrence']['occ'][np.arange(len(a))[a==key0]])
                for key in types_to_inspect:
                    if (key in a):
                        y = np.squeeze(cooccur[a ==key])
                        dic[key]+=[y]

        return dic

    if types_to_inspect is None:
        types_to_inspect = adata.obs.Pixie.unique()
   

    if by_response_status == True:
        a = {}
        response_acq_ID = pd.DataFrame([*adata.obs[['Response','acquisition_ID']].value_counts().index],columns = ['Response','acquisition_ID'])
        for response_status,index in response_acq_ID.groupby('Response').groups.items():
            a[response_status] = make_dic(adata,response_acq_ID.loc[index].acquisition_ID,key0,types_to_inspect)
        return a['pCR'],a['Non-Responder'],interval
    else:
        return make_dic(adata,adata.obs.acquisition_ID.unique(),key0,types_to_inspect),interval
core = adata[adata.obs['SAMPLE_TYPE_(CORE/RESECTION)']=='CORE']
resection = adata[adata.obs['SAMPLE_TYPE_(CORE/RESECTION)']=='RESECTION']
dic_pCR = {}
dic_nR = {}
dic_resection={}
for key0 in types_to_inspect:
    dic_pCR[key0],dic_nR[key0],x =  distance_dictionary(key0,core,types_to_inspect,by_response_status = True)
    dic_resection[key0],x =  distance_dictionary(key0,resection,types_to_inspect,by_response_status = False)

with open('distance_distribution.pkl', 'wb') as f:
    pickle.dump({'resection':dic_resection,'pCR':dic_pCR,'non_responder':dic_nR,'x':x}, f)