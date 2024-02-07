import pandas as pd
cell_table_path='../../segmentation/cell_table_Denoised/cell_table_size_normalized_cell_labels.csv'
biosamples_path='../../IMC_data/ExtraDocs/processed_response.csv'
cell_table = pd.read_csv(cell_table_path,index_col=0)
biosamples =pd.read_csv(biosamples_path)
a = pd.DataFrame()
a['acquisition_ID'] = cell_table.fov.values
a['Leap_ID'] = a.acquisition_ID.str.split('_',n = 1).str[0].str.upper()
a['Leap_ID'] = a.Leap_ID.str[:7]#leap_ID should be Leap123, anything more is stripped
a = a.reset_index().merge(biosamples,left_on='Leap_ID',right_on= 'LEAP_ID').drop(['LEAP_ID'],axis = 1).set_index('index')
cell_table[(a['SAMPLE_TYPE_(CORE/RESECTION)']=='CORE')*(a.Response == 'Non-Responder')].to_csv('../../segmentation/cell_table_Denoised/cell_tb_by_group/nonresponder.csv')

cell_table[(a['SAMPLE_TYPE_(CORE/RESECTION)']=='CORE')*(a.Response == 'pCR')].to_csv('../../segmentation/cell_table_Denoised/cell_tb_by_group/pCR.csv')

cell_table[a['SAMPLE_TYPE_(CORE/RESECTION)']=='RESECTION'].to_csv('../../segmentation/cell_table_Denoised/cell_tb_by_group/resection.csv')
