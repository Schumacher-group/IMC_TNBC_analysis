import scanpy as sc
import pandas as pd
import os
import numpy as np
import logging
logger = logging.getLogger()
def discretise(data,thr):
    '''binarise the counts to 0/1'''
    if isinstance(data, sc.AnnData):
        return data.X>thr
    else:
        return data>thr
def normalise(adata,quantile):
    if isinstance(adata, sc.AnnData):
        data =  adata.X
    else:
        data = adata.copy()
    if np.all(data<1):
        logger.warning('data seem already normalised, skipping normalisation')
        return adata
    q = np.quantile(data,q = quantile,axis = 0)
    data = data/q
    data[data>1] = 1
    if isinstance(adata, sc.AnnData):
        adata.X = data
        return adata
    else:
        return data


def quality_control(intensities,low_gene_active = 0.2,high_gene_active = 0.5):
    if 'pass_qc' in intensities.columns:
        return intensities['pass_qc']
    intensities_protein = intensities.iloc[:,1:intensities.columns.get_loc('label')]#proteins are from the second columns up to the column called label
    #For phenotyping , we  want the cells to express some markers, but not all together at the same time.
    #We are gonna consider the set of markers that we use for phenotyping

    markers_4_phenotyping = ['CD38', 'CD14', 'Tbet', 'CD16', 'CD163',
       'Pan-keratin', 'CD11b', 'CD107a', 'CD45', 'CD44', 'CD366',
       'FOXP3', 'CD4', 'E-Cadherin', 'CD68', 'HLA-DR-DQ-DP', 'CD20',
       'CD8a', 'Beta-Catenin', 'B7-H4', 'Granzyme-B',
       'CD3', 'CD27', 'CD45RO',
       'Alpha-SMA', 'Vimentin', 'CD31' ]
    intensities_protein =    normalise(intensities_protein,0.95)
    A = intensities_protein.loc[:,markers_4_phenotyping]
    cond_few_genes_in_cell  = np.sum(discretise(A,thr = low_gene_active),axis=1)>0
    cond_many_genes_in_cell  = np.sum(discretise(A,thr = high_gene_active),axis=1)<11
    dna_count = intensities_protein[['DNA1', 'DNA2']].sum(axis = 1)
    dna_thr = np.quantile(dna_count,0.05)
    dna_cond = dna_count>dna_thr
    cond = cond_few_genes_in_cell& cond_many_genes_in_cell&dna_cond
    return cond
def generate_anndata_from_cell_table(cell_table_path = None,biosamples_path = None):
    '''
    Here I load the spatial data, which consists of the protein intensity per cell, and the geometry location of the cell. I use the cell type annotation from Pixie. I filter out images with less than 1000 cells 
    Also I remove cells with the lowest 5%
    
    '''
    base_dir = "../../"

    if cell_table_path is None:
        cell_table_path = os.path.join(base_dir, 'segmentation', 'cell_table', 'cell_table_size_normalized_cell_labels.csv')
    if biosamples_path is None:
        biosamples_path = base_dir+'IMC_data/ExtraDocs/processed_response.csv'
    cell_table = pd.read_csv(cell_table_path)
    if 'qc_pass' not in cell_table.columns:
        logger.info('Generating quality control mask')
        qc_pass = quality_control(cell_table)
        cell_table['qc_pass'] = qc_pass
        logger.info('writing the cell table file at '+cell_table_path)
        cell_table.to_csv(cell_table_path,index = False)

    if 'cell_meta_cluster' in cell_table.columns:
        cell_table = cell_table[cell_table['cell_meta_cluster']!='Unassigned']#remove cells that have not been assigned yet
    biosamples =pd.read_csv(biosamples_path)
    intensities_protein = cell_table.iloc[:,1:cell_table.columns.get_loc('label')]#proteins are from the second columns up to the column called label
    intensities_protein['Carboplatin'] = cell_table['Carboplatin_nuclear']

    logger.info('Finished loading, now create the anndata object')
    adata = sc.AnnData(intensities_protein, obsm={"spatial": cell_table[['centroid-0', 'centroid-1']].values})
    try:
        adata.obs['Pixie'] = pd.Categorical(cell_table.cell_meta_cluster.values.astype(str))
    except:
        print('cell type label not present')
    adata.obs['acquisition_ID'] = cell_table.fov.values
    adata.obs['Leap_ID'] = adata.obs.acquisition_ID.str.split('_',n = 1).str[0].str.upper()
    adata.obs = adata.obs.reset_index().merge(biosamples,left_on='Leap_ID',right_on= 'LEAP_ID').drop(['LEAP_ID'],axis = 1).set_index('index')
    adata.obs['qc_pass'] = cell_table['qc_pass'].values
    adata = adata[~((adata.obs.Response == 'Responder')&(adata.obs['SAMPLE_TYPE_(CORE/RESECTION)']=='RESECTION'))]#remove cases of resection of responders

    # get fovs having more than 1000 cells
    fovs = adata.obs.acquisition_ID.value_counts()[adata.obs.acquisition_ID.value_counts()>=1000].index
    adata = adata[adata.obs.acquisition_ID.isin(fovs)]
    adata.raw = adata#raw data are unfiltered and unnormalised
    

    #Normalise each channel independently by quantile
    adata = normalise(adata,quantile=0.95)

    
    return adata
generate_anndata_from_ark_analysis = generate_anndata_from_cell_table