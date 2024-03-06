import pandas as pd
import numpy as np
import warnings
import re
def split_text(text):
    '''
    Text is match at the beginning and the string is shortened repeatedely
    '''
    regex = r'^\d{1,3}=core|resection'
    text = text.replace(' ','')#remove any space in the string
    data = {}
    while True:
        a = re.search(regex,text,flags=re.IGNORECASE)
        if a is None:
            break
        temp =  text[:a.end()]
        key,value = temp.split('=')
        data[int(key)] = value
        text= text[a.end():]
    return data
def unravel_response_df(response):     
    '''Transform the response Excel file to machine readable'''      
    def process_row(row):
        def check_keywords(input_text):
            keywords = ["bilateral", "same patient"]
            for keyword in keywords:
                if keyword.upper() in input_text.upper():
                    return True
            return False

        leap_id = row.LEAP_ID
        prefix = 'LEAP'
        num_part = leap_id.replace(prefix,'')#remove LEAP
        parts = num_part.split('/')
        if len(parts) > 1:
            #if comments distinguis between core/section
            if ('core'.upper() in row.COMMENTS.upper()) or ('resection'.upper() in row.COMMENTS.upper()):  
                if len(parts) == row.COMMENTS.count('='):
                    #the number of = characters  matches the number of split in the LEAP_ID column, that is the normal cases
                    #key is leap number, value is the corresponding comment (core or resection)
                    mapper = split_text(row.COMMENTS)#transform mapper to dict and make the key an integer to help indexing of the dictionary below
                elif (('core'.upper() in row.COMMENTS.upper()) ^('resection'.upper() in row.COMMENTS.upper())):
                    # we are in the case where we have more Leaps that information about CORE/RESECTION.
                    #If Comments column contains either CORE or Resection, but not both, this will be assigned to all leaps

                    if 'core' in row.COMMENTS:
                        name = 'core'
                    else:
                        name = 'resection'
                    mapper = {}
                    for part in parts:
                        ind = int(part)
                        mapper[ind] = name.upper()
                else:
                    raise ValueError('Multiple Leap ID do not correspond but not the same number of info in COMMENTS column.\n Failed to interpret the comment: ', row.COMMENTS)
                
            else:
                raise ValueError('Failed to interpret the comment: ', row.COMMENTS)
            try:

                new_rows= []
                for part in parts:
                    #split the row into two
                    new_row = row.copy()
                    new_row.LEAP_ID = prefix+"%03d" %int(part)
                    ind = int(part)
                    new_row.COMMENTS = mapper[ind].upper()
                    new_rows +=[new_row]
            except KeyError:
                raise KeyError('Failed to interpret the comment: ', row.COMMENTS, 'Cannot find key',ind)
            return new_rows
        else:
            if 'CORE' in row.COMMENTS.upper():
                name = 'core'
            elif 'RESECTION' in row.COMMENTS.upper():
                name = 'resection'
            else:
                raise ValueError('Failed to interpret the comment: ', row.COMMENTS)
            row.LEAP_ID = prefix+"%03d" %int(num_part)
            row.COMMENTS = name.upper()
            return [row]

    transformed = []
    for _,el in response.iterrows():
        #print(el)
        transformed+=[process_row(el)]
    return pd.DataFrame([b for el in transformed for b in el])
path0 = '../'
response = pd.read_excel(path0+"IMC_data//ExtraDocs/LEAP code response data 05122023.xlsx")
#response = pd.read_excel(path0+"IMC_data//ExtraDocs/FLIM LEAP TNBC records.xlsx")
response.columns = response.columns.str.replace(' ','_')#use underscore rather than space
response.COMMENTS = response.COMMENTS.str.replace(',',' ')#the file mix space and , as seprators, using only whitespace for simplicity

response = unravel_response_df(response[['LEAP_ID', 'Response', 'COMMENTS']])

biosamples =pd.read_excel(path0+'IMC_data/ExtraDocs/LEAP biobank samples.xlsx')
biosamples.columns = biosamples.columns.str.replace(' ','_')#use underscore rather than space
biosamples.drop('OUTCOME____(RESPONDER,_NON-RESPONDER)',inplace = True,axis = 1)# response status is in response
metadata = biosamples.merge(response,on='LEAP_ID')
#check that the core/resection of the 2 dataframes is consistent whenever data are provided
cond1 = ~metadata['SAMPLE_TYPE_(CORE/RESECTION)'].isna()
cond2 = ~metadata['COMMENTS'].isna()

if not (metadata['SAMPLE_TYPE_(CORE/RESECTION)']==metadata['COMMENTS'].str.upper())[cond1*cond2].all():
    warnings.warn('The response data and biobank sample dataframes have different core/resection annotation')
metadata.to_csv(path0+'IMC_data//ExtraDocs/processed_response.csv',index = False)
