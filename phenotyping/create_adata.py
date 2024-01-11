import scanpy as sc
from phenotyping_itils import generate_anndata_from_cell_table
base_dir = '../../'
cell_table_path=base_dir+'segmentation/cell_table_Denoised/cell_table_size_normalized.csv'
biosamples_path=base_dir+'IMC_data/ExtraDocs/processed_response.csv'
### TO DO: It should  read the adata inside phenograph and pixie, and combine together into a unique adata