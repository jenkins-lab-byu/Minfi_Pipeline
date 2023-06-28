
### Pathnames to change before running pipeline
# Pathname to the minfi pipeline directory that was downloaded from github. (Make sure to include a '/' at the end')
script_base=/Volumes/Research_Data/Jenkins_Lab_Github/Minfi_Pipeline/EPICv2_Pipeline/

# Pathname to the basedirectory that contains all files relevant to the study
base_dir="/Volumes/Research_Data/Research_Datasets/Neuro_Panel/Jenkins_MethylationEPICv2_20230609"

#Name of the basedirectory that contains all files relevant to the study (Make sure it is just the name of the directory, not the pathname)
name=Jenkins_MethylationEPICv2_20230609

# Name of the Sample Sheet inside your basedirectory that contains all the sample metadata. (Make sure the first three columns are 'Sample_Name','Sentrix_ID','Sentrix_Position')
sample_sheet=sample_sheet.csv

# You dont need to touch anything after this point!!!!

# Pathname to the basedirectory that contains all 
additional_files_base="${script_base}additional_files/"
minfi_script="${script_base}src/array_preprocessing.R"
merge_script="${script_base}src/merge_data.py"


# batch info
input_dir=$base_dir
output_dir=$base_dir
sample_sheet_path=$base_dir$sample_sheet
mvals="${base_dir}${name}_m_values.csv"

###########
## MINFI ##
###########
echo ""
echo "**** RUNNING MINFI SCRIPT ****"
echo ""

Rscript $minfi_script --input_dir $input_dir --output_dir $output_dir --sample_sheet $sample_sheet --name $name --betas --mvals

# Outputs included: beta_values, m_values, p_values, methylated_points, unmethylated_points, QC_report