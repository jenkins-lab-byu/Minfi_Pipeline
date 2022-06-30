import argparse
from operator import index
import pandas as pd
import re
import os
import glob

def natural_sort(mylist):
    #got code from https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(mylist, key = alphanum_key)

def get_files(path,file_suffix):
    os.chdir(path)

    my_files=glob.glob('*'+file_suffix)

    return(natural_sort(my_files))

def check_path(path):
    if path[-1]!='/':
        path+='/'

    return(path)

if __name__ == "__main__" :
    parser=argparse.ArgumentParser(description='Combined all with a given suffix in a given directory')
    parser.add_argument('--dir', help='Full path to directory holding files to concatenate',type=str, required=True)
    parser.add_argument('--suffix', help='Suffix of files to concatenate',type=str, required=True)
    parser.add_argument('--name', help='Base name of output file (without extension)',type=str, required=True)

    args=parser.parse_args()

    dir=args.dir
    suffix=args.suffix
    name=args.name

    dir=check_path(dir)

    files_to_concat=get_files(dir,suffix)


    combined_data=pd.DataFrame()

    for file_path in files_to_concat:
        curr_file=pd.read_csv(file_path,header=0)

        combined_data=pd.concat([combined_data,curr_file],ignore_index=True)

    combined_data['Specimen_Type']='sperm'

    combined_data.loc[combined_data.Sample_Name.str[:2]=='IN','Specimen_Type']='blood'
    combined_data.loc[combined_data.Sample_Name.str[:6]=='20UHCQ','Specimen_Type']='blood'
    combined_data.loc[combined_data.Sample_Name.str[-6:]=='_blood','Specimen_Type']='blood'

    combined_data.loc[combined_data.Sample_Name.str[:8]=='Positive','Specimen_Type']='positive_control'


    combined_data.to_csv('{}.csv'.format(name),header=True,index=None)
    combined_data.to_excel('{}.xlsx'.format(name),header=True,index=None)
