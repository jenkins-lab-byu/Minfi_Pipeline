
cat("\n**LOADING PACKAGES**\n\n")
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(argparse))
#source("calc_sperm_age.R")

## set up parser
parser = ArgumentParser(description='script to normalize microarray data with minfi and calculate sperm age')
parser$add_argument('--model', type='character', help='full path to sperm age model RDS file')
parser$add_argument('--sperm_age_script', type='character', help='full path to sperm age script')
parser$add_argument('--input_dir', type='character', help='full path input directory',required=TRUE)
parser$add_argument('--output_dir', type='character', help='full path output directory',required=TRUE)
parser$add_argument('--sample_sheet', type='character', help='name of sample sheet',required=TRUE)
parser$add_argument('--name', type='character', help='name you want to give to this run',required=TRUE)
# parser$add_argument('--array', type='character', help='use "450k" or "epic"',default='epic')
parser$add_argument("--betas", action='store_true', help='use this flag if you want to write beta values to .csv file', default=FALSE)
parser$add_argument("--mvals", action='store_true', help='use this flag if you want to write m-values to .csv file', default=FALSE)
#parser$add_argument("--sperm_age", action='store_true', help='use this flag if you want to calculate sperm age and write to file', default=FALSE)
#parser$add_argument("--dlk1", action='store_true', help='use this flag if you want to calculate and write DLK1 means to file', default=FALSE)
parser$add_argument("--analyze", action='store_true', help='use this flag if want to calculate 1) sperm age 2) DLK1 mean and 3) mean signal intensities and write all to file', default=FALSE)

args = parser$parse_args()

## get arguments
model_path=args$model
sperm_age_script=args$sperm_age_script
input_dir=args$input_dir
output_dir=args$output_dir
sample_sheet=args$sample_sheet
run_name=args$name
write_betas=args$betas
write_mvals=args$mvals
analyze=args$analyze

source(sperm_age_script)

## this fxn comes from: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
## fxn to grab last N characters of string
substr_right = function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## make sure last char of path is '/' to help with naming output file
clean_path = function(my_path){
  if(substr_right(my_path,1)!='/'){
    my_path=paste(my_path,'/',sep='')
  }
  return(my_path)
}

get_dlk1=function(df,cgs){
  dlk1_means=data.frame(rowMeans(df[,cgs],na.rm=T))
  colnames(dlk1_means)=c("DLK1_Mean")
  return(dlk1_means)
}

dlk1_cgs=c('cg09212014','cg23522522','cg06504820','cg16544010','cg16103098','cg18392604','cg10113715','cg15831296','cg09865386','cg13684115','cg18864857','cg14092276','cg18373878','cg10528576')

## make sure last char of paths are '/' to help create output file names
input_dir=clean_path(input_dir)
output_dir=clean_path(output_dir)

# set.seed(24601)
set.seed(1)
model<-readRDS(model_path)

## grab sample sheet
cat("\n**GRABBING SAMPLE SHEET**\n\n")
targets = read.metharray.sheet(base=input_dir, pattern = sample_sheet, ignore.case = TRUE, recursive = TRUE, verbose = TRUE)

## load methylation data
cat("\n**LOADING IDAT FILES**\n\n")
RGset<-read.metharray.exp(targets=targets,recursive=TRUE,force=TRUE)

## normalize data
cat("\n**NORMALIZING DATA (SWAN)**\n\n")
grSet = preprocessSWAN(RGset)

## QC check
cat("\n**CHECKING QC**\n\n")
qc <- getQC(grSet)

pdf(paste(output_dir,run_name,'_intensity_plot.pdf',sep=''))
plotQC(qc)
dev.off()

qcReport(RGset,pdf=paste(output_dir,run_name,'_qc_report.pdf',sep=''))


## get beta values
BetaValues=getBeta(grSet)

array_type=''

## sperm age, dlk1, intensities, snps
if(analyze){
  
  # sperm age
  cat("\n**CALCULATING SPERM AGE**\n\n")
  sperm_age_calculations=calc_sperm_age(BetaValues,model)
  colnames(sperm_age_calculations)=c("Sperm_Age")

  # dlk1 means
  BetaValues=data.frame(t(BetaValues))
  if(ncol(BetaValues)<=500000){
    cat("\n**THIS LOOKS LIKE 450K DATA**\n\n")
    array_type='450k'
  }
  else{
    cat("\n**THIS LOOKS LIKE EPIC DATA**\n\n")
    array_type='EPIC'
  }

  dlk1_means=get_dlk1(BetaValues,dlk1_cgs) # DLK1 CGs are same for both 450k and EPIC array data

  ###########################
  # signal intensity values #
  ###########################
  cat("\n**GRABBING SIGNAL INTENSITY VALUES**\n\n")
  unmethylated.points=getUnmeth(grSet) # unmethylated values
  unmethylated.points=t(unmethylated.points)
  methylated.points=getMeth(grSet) # methylated values
  methylated.points=t(methylated.points)
  
  meth_mean=data.frame(rowMeans(methylated.points,na.rm=T))
  colnames(meth_mean)=c("Methylated_Signal_Mean")
  unmeth_mean=data.frame(rowMeans(unmethylated.points,na.rm=T))
  colnames(unmeth_mean)=c("Unmethylated_Signal_Mean")
  signal_data=cbind(meth_mean,unmeth_mean)
  signal_data$Combined_Signal_Means=signal_data$Methylated_Signal_Mean+signal_data$Unmethylated_Signal_Mean

  ## add column with array type (450k or EPIC)
  signal_data$Array=array_type
  signal_data$Batch=run_name


  ############
  # snp data #
  ############
  cat("\n**GRABBING SNP VALUES**\n\n")
  snps=getSnpBeta(RGset)
  snps=t(snps)

  ## merge all data
  combined_data=cbind(sperm_age_calculations,dlk1_means,signal_data,snps)
 
  ## write data to file 
  write.table(combined_data,file=paste(output_dir,run_name,'_minfi_results.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)
}

if(write_betas){
  cat("\n**GRABBING AND WRITING BETA VALUES TO FILE**\n\n")
  ## check orientation of data frame and correct if necessary
#  if(analyze==FALSE){
#    BetaValues=t(BetaValues)
#  }
  # I know I already extracted beta values but I couldn't get them to write to file unless I extracted them again
  BetaValues=getBeta(grSet)
  #BetaValues=t(BetaValues) # loading this in python is faster if you don't transpose data frame for some reason
  write.table(BetaValues,file=paste(output_dir,run_name,'_beta_values.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)
}

## get m-values & write to file
if(write_mvals){
  cat("\n**GRABBING AND WRITING M VALUES TO FILE**\n\n\n")
  MValues=getM(grSet)
  #MValues=t(MValues) # loading this in python is faster if you don't transpose data frame for some reason
  write.table(MValues,file=paste(output_dir,run_name,'_m_values.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)
}



