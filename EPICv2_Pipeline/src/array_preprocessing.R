
cat("\n**LOADING PACKAGES**\n\n")
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(BiocManager))

## set up parser
parser = ArgumentParser(description='script to normalize microarray data with minfi for EPICv2 arrays')
parser$add_argument('--input_dir', type='character', help='full path input directory',required=TRUE)
parser$add_argument('--output_dir', type='character', help='full path output directory',required=TRUE)
parser$add_argument('--sample_sheet', type='character', help='name of sample sheet',required=TRUE)
parser$add_argument('--name', type='character', help='name you want to give to this run',required=TRUE)
parser$add_argument("--betas", action='store_true', help='use this flag if you want to write beta values to .csv file', default=FALSE)
parser$add_argument("--mvals", action='store_true', help='use this flag if you want to write m-values to .csv file', default=FALSE)

args = parser$parse_args()

## get arguments
input_dir=args$input_dir
output_dir=args$output_dir
sample_sheet=args$sample_sheet
run_name=args$name
write_betas=args$betas
write_mvals=args$mvals

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

## make sure last char of paths are '/' to help create output file names
input_dir=clean_path(input_dir)
output_dir=clean_path(output_dir)

# set.seed(24601)
set.seed(1)

## grab sample sheet
cat("\n**GRABBING SAMPLE SHEET**\n\n")
targets = read.metharray.sheet(base=input_dir, pattern = sample_sheet, ignore.case = TRUE, recursive = TRUE, verbose = TRUE)

## load methylation data
cat("\n**LOADING IDAT FILES**\n\n")
RGset<-read.metharray.exp(targets=targets,recursive=TRUE,force=TRUE)
annotation(RGset)["array"] = "IlluminaHumanMethylationEPICv2"
annotation(RGset)["annotation"] = "20a1.hg38"

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

## set array type
array_type=''

BetaValues=data.frame(t(BetaValues))
if(ncol(BetaValues)>=900000){
  cat("\n**THIS LOOKS LIKE EPICv2 DATA**\n\n")
  array_type='EPICv2'
}

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

## add column with array type (EPICv2)
signal_data$Array=array_type
signal_data$Batch=run_name

if(write_betas){
  cat("\n**GRABBING AND WRITING BETA VALUES TO FILE**\n\n")
  # I know I already extracted beta values but I couldn't get them to write to file unless I extracted them again
  BetaValues=getBeta(grSet)
  write.table(BetaValues,file=paste(output_dir,run_name,'_beta_values.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)
}

## get m-values & write to file
if(write_mvals){
  cat("\n**GRABBING AND WRITING M VALUES TO FILE**\n\n\n")
  MValues=getM(grSet)
  write.table(MValues,file=paste(output_dir,run_name,'_m_values.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)
}

############
# snp data #
############
cat("\n**GRABBING AND WRITING P-Values**\n\n")
detp<-detectionP(RGset,type="m+u")
write.table(detp,file=paste(output_dir,run_name,'_p_values.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)

###get unmethylated values
cat("\n**GRABBING AND WRITING Unmethylated Points**\n\n")
unmethylated_points<-getUnmeth(grSet)
write.table(unmethylated_points,file=paste(output_dir,run_name,'_unmethylated_points.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)

###get methylated values
cat("\n**GRABBING AND WRITING Methylated Points**\n\n")
methylated_points<-getMeth(grSet)
write.table(methylated_points,file=paste(output_dir,run_name,'_methylated_points.csv',sep=''),row.names = T, col.names = T, sep = ",", quote = F)