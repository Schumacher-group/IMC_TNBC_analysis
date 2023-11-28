import pandas as pd
import numpy as np
import warnings
import re
def unravel_response_df(response):     
    '''Transform the response Excel file to machine readable'''      
    def process_row(row):
        def check_keywords(input_text):
            keywords = ["bilateral", "same patient"]
            for keyword in keywords:
                if keyword in input_text:
                    return True
            return False

        leap_id = row.LEAP_ID
        prefix = 'LEAP'
        num_part = leap_id.lstrip(prefix)
        parts = num_part.split('/')
        if len(parts) == 2:
            #print(parts)
            #split the row into two
            new_row_1 = row.copy()
            new_row_1.LEAP_ID = prefix+parts[0]
            new_row_2 = row.copy()
            id1,id2 = list(map(int,parts))
            #if comments distinguis between core/section
            if ('core' in row.COMMENTS) and ('resection' in row.COMMENTS):    
                mapper = re.findall(r'(\d+)\s*=\s*(\S+)', row.COMMENTS)#key is leap number, value is the corresponding comment (core or resection)
                mapper = {int(el[0]):el[1] for el in mapper}#transform mapper to dict and make the key an integer to help indexing few lines below
            elif check_keywords(row.COMMENTS):
                # if two samples for the same patient
                if 'core' in row.COMMENTS:
                    mapper = {id1:'core',id2:'core'}
                elif 'resection' in row.COMMENTS:
                    mapper = {id1:'resection',id2:'resection'}
            else:
                raise ValueError('Failed to interpret the comment: ', row.COMMENTS)
            try:
                new_row_1.COMMENTS = mapper[id1]
                new_row_2.COMMENTS = mapper[id2]
            except KeyError:
                raise KeyError('Failed to interpret the comment: ', row.COMMENTS)

            #print(id1,id2)
            new_row_2.LEAP_ID = prefix+parts[0][:-len(parts[1])] + parts[1]
            return [new_row_1,new_row_2]
        else:
            return [row]

    transformed = []
    for _,el in response.iterrows():

        transformed+=[process_row(el)]
    return pd.DataFrame([b for el in transformed for b in el])
path0 = '../'
response = pd.read_excel(path0+"IMC_data//ExtraDocs/LEAP code response data 020523.xlsx")
response.columns = response.columns.str.replace(' ','_')#use underscore rather than space
response = unravel_response_df(response)

biosamples =pd.read_excel(path0+'IMC_data/ExtraDocs/LEAP biobank samples.xlsx')
biosamples.columns = biosamples.columns.str.replace(' ','_')#use underscore rather than space
biosamples.drop('OUTCOME____(RESPONDER,_NON-RESPONDER)',inplace = True,axis = 1)# response status is in response
response = pd.read_excel(path0+"IMC_data//ExtraDocs/processed_response.xlsx")
response.columns = response.columns.str.replace(' ','_')#use underscore rather than space
metadata = biosamples.merge(response,on='LEAP_ID')
if ~np.all(metadata['SAMPLE_TYPE_(CORE/RESECTION)']==metadata['COMMENTS'].str.upper()):
    warnings.warn(metadata[['LEAP_ID','SAMPLE_TYPE_(CORE/RESECTION)']][metadata['SAMPLE_TYPE_(CORE/RESECTION)']!=metadata['COMMENTS'].str.upper()],+'not matching')
metadata.to_csv(path0+'IMC_data//ExtraDocs/processed_response.csv',index = False)