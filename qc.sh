#bad_samples=203225140082_R05C01,203219750033_R05C01,203219730092_R03C01
bad_samples=none
batch_name=Batch51_06042022
base_dir="/data/${batch_name}/" # use '/' at end

analysis="${base_dir}${batch_name}_methylation_analysis_results.csv"
script_base=/data/sperm-qt/src/ # use '/' at end
qc_script="${script_base}beta_distribution_qc.py"

python $qc_script --results $analysis --bad $bad_samples
