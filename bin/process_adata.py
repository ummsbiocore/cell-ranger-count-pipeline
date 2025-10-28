"""Process and merge AnnData / Loom files for scVelo.

This script is Nextflow-friendly: pass input file paths and an output path.

Usage example:
    python process_adata.py --ldata LDATA.loom --adata ANALYSIS.h5ad --output processed_adata.h5ad

Inputs:
  --ldata   Path to loom file (or readable file for sc.read)
  --adata   Path to anndata h5ad file (converted from Seurat)
  --output  Path where processed h5ad will be written

The script validates the inputs, selects common barcodes, runs scVelo preprocessing
and velocities, then writes the merged/processed AnnData to --output.
"""

import argparse
import logging
import sys
from pathlib import Path

import scvelo as scv
import scanpy as sc


def setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("process_adata")


def main(ldata_path: Path, adata_path: Path, out_path: Path) -> int:
    log = setup_logger()

    log.info("Reading loom / AnnData inputs")
    if not ldata_path.exists():
        log.error("Loom/input file not found: %s", ldata_path)
        return 2
    if not adata_path.exists():
        log.error("AnnData file not found: %s", adata_path)
        return 2

    # sc.read works for loom and h5ad; keep cache=False for reproducible behavior in pipelines
    try:
        ldata = sc.read(str(ldata_path), cache=False)
    except Exception as e:
        log.exception("Failed to read ldata (%s): %s", ldata_path, e)
        return 3

    try:
        adata = sc.read(str(adata_path))
    except Exception as e:
        log.exception("Failed to read adata (%s): %s", adata_path, e)
        return 3

    # Normalize / fix barcodes in loom obs_names if needed (original pipeline used name.split(':')[1])
    try:
        # Only transform if the pattern with ':' appears in names
        sample_names = list(ldata.obs_names)
        if any(':' in n for n in sample_names):
            log.info("Normalizing loom obs_names to match h5ad barcodes")
            ldata.obs_names = [name.split(':')[1].replace('x', '-1') for name in ldata.obs_names]
    except Exception:
        log.exception("Failed to normalize ldata.obs_names; continuing with original names")

    # Intersect barcodes
    common_barcodes = ldata.obs_names.intersection(adata.obs_names)
    if len(common_barcodes) == 0:
        log.error("No common barcodes found between ldata (%s) and adata (%s)", ldata_path, adata_path)
        return 4

    log.info("Found %d common barcodes; subsetting and merging", len(common_barcodes))

    adata = adata[common_barcodes].copy()
    ldata = ldata[common_barcodes].copy()

    # Merge loom into adata using scVelo helper
    try:
        adata = scv.utils.merge(adata, ldata)
    except Exception:
        log.exception("Failed to merge ldata into adata")
        return 5

    # Preprocessing and velocity pipeline
    try:
        log.info("Running scVelo preprocessing: filter_and_normalize -> moments")
        scv.pp.filter_and_normalize(adata)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

        log.info("Computing velocities and velocity graph")
        scv.tl.velocity(adata)
        scv.tl.velocity_graph(adata)

        log.info("Scoring cell cycle and velocity confidence")
        scv.tl.score_genes_cell_cycle(adata)
        scv.tl.velocity_confidence(adata)
    except Exception:
        log.exception("Processing failed during scVelo steps")
        return 6

    # Ensure output directory exists
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        log.info("Writing processed AnnData to %s", out_path)
        adata.write(str(out_path))
    except Exception:
        log.exception("Failed to write output file %s", out_path)
        return 7

    log.info("Processing completed successfully")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and merge Loom + AnnData for scVelo (Nextflow friendly)")
    parser.add_argument("--ldata", required=True, help="Path to loom file (input)")
    parser.add_argument("--adata", required=True, help="Path to AnnData h5ad file (input)")
    parser.add_argument("--output", required=True, help="Path to output processed h5ad")
    args = parser.parse_args()

    exit_code = main(Path(args.ldata), Path(args.adata), Path(args.output))
    sys.exit(exit_code)