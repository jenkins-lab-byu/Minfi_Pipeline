# batch info
name=Batch51_06042022
sample_sheet=Inherent_Infinium_MethyEPIC_QC_Report_Batch_51_06042022.csv

all_analysis_files_dir=/data/EPIC_Summary_Data/Analyses/
base_dir="/data/${name}/" # use '/' at end
#base_dir=/data/
zip_dir="/data/"
chip_dir="/data/${name}/Raw_Data/"

# define file names
sample_sheet=${base_dir}${sample_sheet}
mvals="${base_dir}${name}_m_values.csv"
betas="${base_dir}${name}_beta_values.csv"
intensity_plot="${base_dir}${name}_intensity_plot.pdf"
qc_report="${base_dir}${name}_qc_report.pdf"
analysis="${base_dir}${name}_methylation_analysis_results.csv"
analysis_with_qc="${base_dir}${name}_methylation_analysis_results-with_qc.csv"
zip="${zip_dir}${name}.zip"

# transfer to s3
aws s3 cp $sample_sheet s3://xfer-inherentbio/Microarray_Data/Sample_Sheets/
aws s3 cp $mvals s3://xfer-inherentbio/Microarray_Data/M_Values/
aws s3 cp $betas s3://xfer-inherentbio/Microarray_Data/Beta_Values/
aws s3 cp $intensity_plot s3://xfer-inherentbio/Microarray_Data/QC/
aws s3 cp $qc_report s3://xfer-inherentbio/Microarray_Data/QC/
aws s3 cp $analysis s3://xfer-inherentbio/Microarray_Data/Analysis/
aws s3 cp $analysis_with_qc s3://xfer-inherentbio/Microarray_Data/Analysis/
aws s3 cp $zip s3://xfer-inherentbio/Microarray_Data/Original_Data/

cp $analysis_with_qc $all_analysis_files_dir

#cp $intensity_plot ${base_dir}qc/
#cp $qc_report ${base_dir}qc/
#cp $analysis ${base_dir}analyses/

#aws s3 cp $zip s3://xfer-inherentbio/Microarray_Data/Original_Data/

cd $chip_dir
# transfer chip data (.idats)
#for chip in 203219670028 203219670145 203219670153 203219670191 203219730092 203219750033 203220070080 203220070153 203225140050 203225140082 203225140170 
for chip in 206129770168
do
        cd $chip_dir$chip

        for idat in *.idat
        do

                aws s3 cp $idat s3://xfer-inherentbio/Microarray_Data/Chip_Data/${chip}/  
        done
        #break

done
