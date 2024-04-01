import os
import pandas as pd
import cobra
from cobra.io import load_json_model
from dgap import m9
from concurrent.futures import ProcessPoolExecutor

def run_fva(model_file):
    # Load model
    model = load_json_model(model_file)
    # Apply media formulation
    m9(model)
    # Perform FVA
    fva_result = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=1.0)
    # Add model ID as prefix to column names
    model_id = os.path.basename(model_file)[:-5]  # Extract model ID from file name
    fva_result.columns = [f"{model_id}_{col}" for col in fva_result.columns]
    return fva_result

def run_parsimonious_fba(model_file):
    # Load model
    model = load_json_model(model_file)
    # Apply media formulation
    m9(model)
    # Perform parsimonious FBA
    pars_fba_result = model.slim_optimize()
    return pars_fba_result

def save_results(result, filename):
    result.to_csv(filename)

if __name__ == '__main__':
    # Define directory containing the GEMs
    gem_dir = '/home/omidard/PanGEMs'
    gem_files = [os.path.join(gem_dir, f) for f in os.listdir(gem_dir) if f.endswith('.json')]
    
    # Run FVA and parsimonious FBA for each model using concurrent.futures
    with ProcessPoolExecutor(max_workers=64) as executor:
        fva_results = list(executor.map(run_fva, gem_files))
        pars_fba_results = list(executor.map(run_parsimonious_fba, gem_files))
    
    # Concatenate FVA results
    fva_df = pd.concat(fva_results, axis=1)
    
    # Create DataFrame for parsimonious FBA results
    pars_fba_df = pd.DataFrame(pars_fba_results, index=[os.path.basename(f)[:-5] for f in gem_files])
    
    # Save results to CSV files
    save_results(fva_df, 'lactopan_fva.csv')
    save_results(pars_fba_df, 'lactopan_parsFBA.csv')
