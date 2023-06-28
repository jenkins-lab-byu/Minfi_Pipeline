import argparse
import pandas as pd
import sys

if __name__ == "__main__" :
    parser=argparse.ArgumentParser(description='Add "Beta_Distribution" column to methylation results file to indicate which samples have good beta distributions or not')
    parser.add_argument('--results', help='Full path to results file',type=str, required=True)
    parser.add_argument('--bad', help='Sentrix info of bad samples separated by commas (e.g. "200325500002_R01C01,200325500084_R02C01"). Don\'t include flag if no samples with bad beta distributions',type=str, required=False, default='none')

    args=parser.parse_args()

    results_path=args.results
    bad_samples=args.bad

    results_file=pd.read_csv(results_path,header=0,index_col=None)
    results_file['Beta_Distribution']='pass'

    # add 'fail' to 'Beta_Distribution' column for samples that have weird beta distributions
    if bad_samples != 'none':
        bad_samples=bad_samples.strip('\n').split(',')

        for sample in bad_samples:
            results_file.loc[results_file.Sentrix_Info==sample,'Beta_Distribution']='fail'

        # make sure number of failed samples actually made it into the file
        marked_as_failed=sorted(list(results_file.loc[results_file.Beta_Distribution=='fail','Sentrix_Info']))
        if marked_as_failed != sorted(bad_samples):
            print('\n###########')
            print('## ERROR ##')
            print('## Samples given: {}##'.format(sorted(bad_samples)))
            print('don\'t equal')
            print('## Samples in results file: {}'.format(marked_as_failed))
            print('## EXITING ##')
            print('###########\n')
            sys.exit()
        else:
            print('\n** Beta QC: Looks like all files were marked correctly **\n')
    else:
        print('\n** Beta QC: All were marked as passing QC **\n')

    # write results to file
    new_path=results_path.strip('\n').split('.csv')[0]
    new_path+='-with_qc.csv'
    results_file.to_csv(new_path,header=True,index=None)

