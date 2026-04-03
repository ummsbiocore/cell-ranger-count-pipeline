#!/usr/bin/env python3

import scvelo as scv
import scanpy as sc
import pandas as pd
import argparse
import sys
from pathlib import Path

def parse_arguments():
    """Parse command-line arguments for Nextflow integration."""
    parser = argparse.ArgumentParser(
        description='Preprocess AnnData with velocity data from loom files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python preprocess_anndata.py \\
    --h5ad Final_Analysis.h5ad \\
    --loom sample1.loom sample2.loom sample3.loom \\
    --output processed_adata.h5ad
    
Sample names are automatically extracted from filenames by removing '_output.loom' suffix.
For example: 'DF_Def_run1_S2_L004_output.loom' -> 'DF_Def_run1_S2_L004'
        """
    )
    
    parser.add_argument(
        '--h5ad',
        required=True,
        help='Path to the input h5ad file (combined AnnData object)'
    )
    
    parser.add_argument(
        '--loom',
        nargs='+',
        required=True,
        help='List of loom file paths (sample names extracted automatically from filenames)'
    )
    
    parser.add_argument(
        '--output',
        default='processed_adata.h5ad',
        help='Output path for the processed h5ad file (default: processed_adata.h5ad)'
    )
    
    parser.add_argument(
        '--n-pcs',
        type=int,
        default=30,
        help='Number of principal components for moments calculation (default: 30)'
    )
    
    parser.add_argument(
        '--n-neighbors',
        type=int,
        default=30,
        help='Number of neighbors for moments calculation (default: 30)'
    )
    
    parser.add_argument(
        '--skip-inspection',
        action='store_true',
        help='Skip the barcode inspection step (useful for production runs)'
    )
    
    return parser.parse_args()


def extract_sample_name(loom_path):
    """
    Extract sample name from loom filename.
    Removes '_output.loom' suffix from the filename.
    
    Args:
        loom_path: Path to loom file
    
    Returns:
        Sample name string
    """
    filename = Path(loom_path).stem  # Get filename without extension
    # Remove '_output' suffix if present
    if filename.endswith('_output'):
        sample_name = filename[:-7]  # Remove last 7 characters ('_output')
    else:
        sample_name = filename
    
    return sample_name


def parse_loom_inputs(loom_paths):
    """
    Parse loom file paths and extract sample names.
    
    Args:
        loom_paths: List of loom file paths
    
    Returns:
        Dictionary mapping sample names to loom file paths
    """
    samples = {}
    
    for loom_path in loom_paths:
        if not Path(loom_path).exists():
            print(f"Error: Loom file not found: {loom_path}")
            sys.exit(1)
        
        sample_name = extract_sample_name(loom_path)
        samples[sample_name] = loom_path
    
    return samples


def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    # --- 1. Validate input files ---
    path_to_h5ad = args.h5ad
    if not Path(path_to_h5ad).exists():
        print(f"Error: H5ad file not found: {path_to_h5ad}")
        sys.exit(1)
    
    # Parse loom inputs and extract sample names automatically
    samples = parse_loom_inputs(args.loom)

    # --- 2. Load the main AnnData object ---
    adata = sc.read(path_to_h5ad)

    # --- 3. Load and concatenate loom files ---
    loom_adatas = []
    for sample_name, loom_path in samples.items():
        # Load a single loom file
        ldata = sc.read(loom_path, cache=True)
        
        # Clean and prefix the loom barcodes
        # Original format: 'local_input_JEYHR:ACTGAGTAGTACACCTx'
        # New format:      'DF_01_S2_L004-ACTGAGTAGTACACCT'
        barcodes = [bc.split(':')[1][:-1] for bc in ldata.obs_names]
        ldata.obs_names = [f"{sample_name}-{bc}" for bc in barcodes]
        
        # Make gene names unique to prevent concatenation errors
        ldata.var_names_make_unique()
        
        loom_adatas.append(ldata)

    # Concatenate all loom-derived AnnData objects into one
    if loom_adatas:
        if len(loom_adatas) > 1:
            # We set index_unique=None to prevent concatenate from adding its own suffix.
            ldata_merged = loom_adatas[0].concatenate(
                loom_adatas[1:], 
                batch_key='sample', 
                batch_categories=samples.keys(),
                index_unique=None
            )
        else:
            # Only one sample, no need to concatenate
            ldata_merged = loom_adatas[0]
    else:
        raise ValueError("No loom files were loaded. Please check the 'samples' dictionary.")

    # --- 4. Merge the data ---
    # The rest of the notebook will continue from here.
    # The original notebook had a manual barcode cleaning step.
    # The merging logic is now handled in the next cell.
    # We will rename the merged loom data to 'ldata' to match the next cell.
    ldata = ldata_merged



    # --- 5. Unify barcode formats and Merge ---

    # The goal is to create a common barcode format: 'sample_name-barcode_sequence'

    # 1. Transform adata barcodes
    # Original format: 'AAAGATGTCATAACCG-1'
    # We use the 'orig.ident' column to get the sample name.
    # New format:      'DF_01_S2_L004-AAAGATGTCATAACCG'
    adata.obs_names = adata.obs['orig.ident'].astype(str) + '-' + adata.obs_names.str.split('-').str[0]

    # Now, find the common barcodes between the two datasets
    common_barcodes = adata.obs_names.intersection(ldata.obs_names)
    if len(common_barcodes) == 0:
        print("Error: No common barcodes found.")
        sys.exit(1)
    
    # Subset both anndata objects to keep only the common cells
    adata = adata[common_barcodes].copy()
    ldata = ldata[common_barcodes].copy()

    # Finally, merge the two objects. 
    # scv.utils.merge will add the 'spliced' and 'unspliced' layers from ldata to adata.
    adata = scv.utils.merge(adata, ldata)

    # --- 6. Velocity analysis ---
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)
    scv.tl.score_genes_cell_cycle(adata)
    scv.tl.velocity_confidence(adata)
    
    # --- 7. Save output ---
    adata.write(args.output)


if __name__ == "__main__":
    main()

