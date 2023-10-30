import scanpy as sc
import pandas as pd
import os
def generate_anndata_from_ark_analysis(cell_table_path = None,biosamples_path = None):
    '''
    Here I load the spatial data, which consists of the protein intensity per cell, and the geometry location of the cell. I use the cell type annotation from Pixie. I filter out images with less than 1000 cells 
    '''
    base_dir = "../../"
    if cell_table_path is None:
        cell_table_path = os.path.join(base_dir, 'segmentation', 'cell_table', 'cell_table_size_normalized_cell_labels.csv')
    if biosamples_path is None:
        biosamples_path = base_dir+'IMC_data/ExtraDocs/processed_response.csv'
    intensities = pd.read_csv(cell_table_path)
    intensities = intensities[intensities['cell_meta_cluster']!='Unassigned']#remove cells that have not been assigned yet
    biosamples =pd.read_csv(biosamples_path)
    intensities_protein = intensities.iloc[:,1:intensities.columns.get_loc('label')]#proteins are from the second columns up to the column called label
    adata = sc.AnnData(intensities_protein, obsm={"spatial": intensities[['centroid-0', 'centroid-1']].values})
    adata.obs['Pixie'] = pd.Categorical(intensities.cell_meta_cluster.values.astype(str))
    adata.obs['acquisition_ID'] = intensities.fov.values
    adata.obs['Leap_ID'] = adata.obs.acquisition_ID.str.split('_',n = 1).str[0].str.upper()
    adata.obs = adata.obs.merge(biosamples,left_on='Leap_ID',right_on= 'LEAP_ID').drop(['LEAP_ID'],axis = 1)
    # get fovs having more than 1000 cells
    fovs = adata.obs.acquisition_ID.value_counts()[adata.obs.acquisition_ID.value_counts()>=1000].index
    adata = adata[adata.obs.acquisition_ID.isin(fovs)].copy()
    return adata