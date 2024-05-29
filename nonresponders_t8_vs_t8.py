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

# Filter data to only include T8 time point
t8_treatment_responders = treatment_responders[treatment_responders['Timepoint'] == 'T8']
t8_control_responders = control_responders[control_responders['Timepoint'] == 'T8']
t8_treatment_nonresponders = treatment_nonresponders[treatment_nonresponders['Timepoint'] == 'T8']
t8_control_nonresponders = control_nonresponders[control_nonresponders['Timepoint'] == 'T8']

# Print lengths of groups
print(f"Total samples: {len(data)}")
print(f"Total treatment group: {len(treatment_group)}")
print(f"Total control group: {len(control_group)}")
print(f"Total treatment responders: {len(treatment_responders)}")
print(f"Total treatment nonresponders: {len(treatment_nonresponders)}")
print(f"Total control responders: {len(control_responders)}")
print(f"Total control nonresponders: {len(control_nonresponders)}")
print(f"T8 treatment responders: {len(t8_treatment_responders)}")
print(f"T8 control responders: {len(t8_control_responders)}")
print(f"T8 treatment nonresponders: {len(t8_treatment_nonresponders)}")
print(f"T8 control nonresponders: {len(t8_control_nonresponders)}")

def perform_t_tests(group1, group2):
    miRNA_columns = group1.columns[5:]  # miRNA expression columns
    results = pd.DataFrame(index=miRNA_columns, columns=['p_value', 'corrected_p_value', 'Direction'])

    for miRNA in miRNA_columns:
        group1_values = group1[miRNA].values
        group2_values = group2[miRNA].values
        
        # Check for NaN values
        if np.isnan(group1_values).any() or np.isnan(group2_values).any():
            print(f"NaN values found in miRNA: {miRNA}")
            continue
        
        # Perform t-test
        t_stat, p_val = ttest_ind(group1_values, group2_values, equal_var=False, nan_policy='omit')
        results.loc[miRNA, 'p_value'] = p_val
        results.loc[miRNA, 'Direction'] = 'Increase' if np.mean(group2_values) > np.mean(group1_values) else 'Decrease'
    
    # Print raw p-values before correction
    print("Raw p-values:")
    print(results['p_value'].dropna().sort_values())
    
    # Adjust p-values using Benjamini-Hochberg
    corrected_p_values = multipletests(results['p_value'].dropna(), alpha=0.05, method='fdr_bh')[1]
    results.loc[results['p_value'].dropna().index, 'corrected_p_value'] = corrected_p_values
    
    return results[results['p_value'] < 0.05]

# Perform the analysis between T8 control responders and T8 treatment responders
significant_differential_expression_responders = perform_t_tests(t8_control_responders, t8_treatment_responders)

# Perform the analysis between T8 control nonresponders and T8 treatment nonresponders
significant_differential_expression_nonresponders = perform_t_tests(t8_control_nonresponders, t8_treatment_nonresponders)

# Add response information to the results
significant_differential_expression_responders['Response'] = 1  # since we are comparing responders
significant_differential_expression_nonresponders['Response'] = 0  # since we are comparing nonresponders

# Sort by corrected p-value
significant_differential_expression_responders = significant_differential_expression_responders.sort_values('corrected_p_value')
significant_differential_expression_nonresponders = significant_differential_expression_nonresponders.sort_values('corrected_p_value')

# Save the results to a CSV file
output_dir = './'
significant_differential_expression_responders.to_csv(f'{output_dir}t8_responders.csv', index=True)
significant_differential_expression_nonresponders.to_csv(f'{output_dir}t8_nonresponders.csv', index=True)

print("Differential expression analysis complete. Results saved to CSV files.")