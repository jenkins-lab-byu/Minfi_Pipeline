dir=/data/EPIC_Summary_Data/Analyses
suffix=_methylation_analysis_results-with_qc.csv
name=EPIC_Summary_Data

script=/data/sperm-qt/src/concat.py

python $script --dir $dir --suffix $suffix --name $name
