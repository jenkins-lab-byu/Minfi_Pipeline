### Pathnames to change before running pipeline
# Pathname to the minfi pipeline directory that was downloaded from github. (Make sure to include a '/' at the end')
script_base=/Volumes/Research_Data/My_Github/Minfi_Pipeline/

# Pathname to the basedirectory that contains all files relevant to the study
base_dir="/Volumes/Research_Data/Research_Datasets/Sperm/Blood_Sperm/Jenkins-BYU_MethylationEPIC_20201221_bad_samples_removed"

#Name of the basedirectory that contains all files relevant to the study (Make sure it is just the name of the directory, not the pathname)
name=Jenkins-BYU_MethylationEPIC_20201221_bad_samples_removed

# Name of the Sample Sheet inside your basedirectory that contains all the sample metadata. (Make sure the first three columns are 'Sample_Name','Sentrix_ID','Sentrix_Position')
sample_sheet=sample_sheet.csv



# Pathname to the basedirectory that contains all 
additional_files_base="${script_base}additional_files/"
minfi_script="${script_base}src/array_preprocessing.R"
sperm_age_script="${script_base}src/calc_sperm_age.R"
# stability_script="${script_base}src/calculate_variability.py"
merge_script="${script_base}src/merge_data.py"
sperm_age_model="${additional_files_base}GLA_Model.rds"


# batch info
input_dir=$base_dir
output_dir=$base_dir
#sample_sheet=Inherent_Infinium_MethyEPIC_QC_Report_${name}.csv
sample_sheet_path=$base_dir$sample_sheet
mvals="${base_dir}${name}_m_values.csv"

###########
## MINFI ##
###########
echo ""
echo "**** RUNNING MINFI SCRIPT ****"
echo ""

Rscript $minfi_script --model $sperm_age_model --sperm_age_script $sperm_age_script --input_dir $input_dir --output_dir $output_dir --sample_sheet $sample_sheet --name $name --betas --mvals --analyze

# Outputs included: beta_values, m_values, p_values, QC_report, minfi_results
# Outputs included in minfi_results: Sperm_Age, DLK1 Mean, Methylation Signal Intensities, SNP rs values