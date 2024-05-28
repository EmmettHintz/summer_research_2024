from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np

# Load the data
file_path = '/Users/emmetthintz/Documents/Computational-Biology/Data/GSE97154_Cleaned.csv'
data = pd.read_csv(file_path)

# Separate groups
treatment_group = data[data['Treatment'] == 1]
control_group = data[data['Treatment'] == 0]

# Further separate treatment group into responders and nonresponders
treatment_responders = treatment_group[treatment_group['Response'] == 1]
treatment_nonresponders = treatment_group[treatment_group['Response'] == 0]

# Further separate control group into responders and nonresponders (if applicable)
control_responders = control_group[control_group['Response'] == 1]
control_nonresponders = control_group[control_group['Response'] == 0]

def perform_bootstrap_t_tests(group, n_bootstraps=5):
    miRNA_columns = group.columns[5:]  # miRNA expression columns
    bootstrap_results = pd.DataFrame(index=miRNA_columns, columns=range(n_bootstraps))
    directionality_results = pd.DataFrame(index=miRNA_columns, columns=range(n_bootstraps))
    
    for i in range(n_bootstraps):
        # Create bootstrap sample for each group
        bootstrap_sample_t0 = group[group['Timepoint'] == 'T0'].sample(frac=1, replace=True)
        bootstrap_sample_t8 = group[group['Timepoint'] == 'T8'].sample(frac=1, replace=True)
        
        for miRNA in miRNA_columns:
            t0_values = bootstrap_sample_t0[miRNA].values
            t8_values = bootstrap_sample_t8[miRNA].values
            
            # Perform t-test
            t_stat, p_val = ttest_ind(t0_values, t8_values, equal_var=False, nan_policy='omit')
            bootstrap_results.loc[miRNA, i] = p_val
            directionality_results.loc[miRNA, i] = 'Increase' if np.mean(t8_values) > np.mean(t0_values) else 'Decrease'
    
    # Calculate the median p-value across all bootstraps for each miRNA
    bootstrap_results['median_p_value'] = bootstrap_results.median(axis=1)
    
    # Adjust median p-values using Benjamini-Hochberg
    corrected_p_values = multipletests(bootstrap_results['median_p_value'], alpha=0.05, method='fdr_bh')[1]
    bootstrap_results['corrected_p_value'] = corrected_p_values
    
    # Add directionality information
    bootstrap_results['Direction'] = directionality_results.mode(axis=1)[0]  # Get the most frequent direction
    
    return bootstrap_results[bootstrap_results['corrected_p_value'] < 0.05]

# Perform the analysis for each subgroup
significant_treatment_responders = perform_bootstrap_t_tests(treatment_responders)
significant_treatment_nonresponders = perform_bootstrap_t_tests(treatment_nonresponders)
significant_control_responders = perform_bootstrap_t_tests(control_responders)
significant_control_nonresponders = perform_bootstrap_t_tests(control_nonresponders)

# Add response information to the results
significant_treatment_responders['Response'] = 1
significant_treatment_nonresponders['Response'] = 0
significant_control_responders['Response'] = 1
significant_control_nonresponders['Response'] = 0

# Save the results to CSV files
output_dir = './data/'

# Save each DataFrame to a CSV file
significant_treatment_responders.to_csv(f'{output_dir}significant_treatment_responders.csv', index=True)
significant_treatment_nonresponders.to_csv(f'{output_dir}significant_treatment_nonresponders.csv', index=True)
significant_control_responders.to_csv(f'{output_dir}significant_control_responders.csv', index=True)
significant_control_nonresponders.to_csv(f'{output_dir}significant_control_nonresponders.csv', index=True)

print("Differential expression analysis complete. Results saved to CSV files.")