import argparse
import pandas as pd 
import os
import glob
import re

def natural_sort(mylist):
    #got code from https://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(mylist, key = alphanum_key)


def get_files(path,file_suffix):
    os.chdir(path)
    my_files=glob.glob('*'+file_suffix)
    return(natural_sort(my_files))

def check_path(my_path):
    if my_path[-1]!='/':
        my_path+='/'
    return(my_path)

if __name__ == "__main__" :
    parser=argparse.ArgumentParser(description='Merge methylation analysis data')
    parser.add_argument('--input_dir', help='Full path to input directory',type=str, required=True)
    parser.add_argument('--output_dir', help='Full path to output directory',type=str, required=True)
    parser.add_argument('--name', help='Name of microarray run or batch',type=str, required=True)
    parser.add_argument('--stability_suffix', help='Suffix of stability score files (e.g. "_stability_scores.csv")',type=str, required=True)
    parser.add_argument('--minfi_suffix', help='Suffix of minfi results file (e.g. "_minfi_results.csv")',type=str, required=True)
    parser.add_argument('--sample_sheet', help='Full path to sample sheet',type=str, required=True)

    args=parser.parse_args()

    input_dir=args.input_dir
    output_dir=args.output_dir
    batch_name=args.name
    stability_suffix=args.stability_suffix
    minfi_suffix=args.minfi_suffix
    sample_sheet_path=args.sample_sheet

    input_dir=check_path(input_dir)
    output_dir=check_path(output_dir)

    # get minfi file (everything gets merged into 'minfi_file')
    minfi_file_name=get_files(input_dir,minfi_suffix)[0]
    minfi_file=pd.read_csv(minfi_file_name,header=0)
        
    # get stability score files
    stability_files=get_files(input_dir,stability_suffix)
    for stability_file_name in stability_files:
        stability_file=pd.read_csv(stability_file_name,header=0,index_col=0)
        minfi_file=pd.merge(minfi_file,stability_file,left_index=True,right_index=True)

    # get info file
    sample_sheet=pd.read_csv(sample_sheet_path,header=0)

    ## remove unnamed last column of info file (if necessary)
    if sample_sheet.columns[-1].split(':')[0]=='Unnamed':
        sample_sheet.drop(sample_sheet.columns[len(sample_sheet.columns)-1], axis=1, inplace=True)
    sample_sheet['Sentrix_Info']=sample_sheet.Sentrix_ID.astype(str)+'_'+sample_sheet.Sentrix_Position.astype(str)
    
    merged_data=pd.merge(minfi_file,sample_sheet,left_index=True,right_on='Sentrix_Info')

    first_cols=['Sample_Name','Sentrix_ID','Sentrix_Position','Sentrix_Info']
    middle_cols=[]
    snp_cols=[]

    for col in merged_data.columns:
        if col[:2]=='rs':
            snp_cols.append(col)
        elif col not in first_cols:
            middle_cols.append(col)

    merged_data=merged_data.loc[:,(first_cols+middle_cols+snp_cols)]
    
    merged_data.to_csv('{}{}_methylation_analysis_results.csv'.format(output_dir,batch_name),header=True,index=False,sep=',')