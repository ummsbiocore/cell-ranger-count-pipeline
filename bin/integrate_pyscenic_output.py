#!/usr/bin/env python3
"""
Integrate pySCENIC outputs into a SCope-compatible loom file
This script combines pySCENIC results with the original expression data
"""

import argparse
import json
import zlib
import base64
import sys
import numpy as np
import pandas as pd
import loompy as lp
import umap
from MulticoreTSNE import MulticoreTSNE as TSNE


def dfToNamedMatrix(df):
    """Convert DataFrame to numpy structured array"""
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr


def main():
    parser = argparse.ArgumentParser(description='Integrate pySCENIC outputs into SCope-compatible loom')
    parser.add_argument('-i', '--input_loom', required=True, 
                        help='Input loom file with expression data (e.g., Data.loom)')
    parser.add_argument('-p', '--pyscenic_output', required=True,
                        help='pySCENIC output loom file with AUC matrix (from aucell step)')
    parser.add_argument('-o', '--output_loom', required=True,
                        help='Output SCope-compatible loom file')
    parser.add_argument('--export_auc_csv', type=str, default=None,
                        help='Optional: Export AUC matrix to CSV file')
    parser.add_argument('-t', '--threads', type=int, default=10,
                        help='Number of threads for tSNE (default: 10)')
    parser.add_argument('--title', default='pySCENIC_analysis',
                        help='Title for the SCope viewer (default: pySCENIC_analysis)')
    parser.add_argument('--genome', default='hg38',
                        help='Genome build (default: hg38)')
    
    args = parser.parse_args()
    
    # Load input loom file
    lf_input = lp.connect(args.input_loom, mode='r', validate=False)
    
    # Load pySCENIC output loom file
    lf_scenic = lp.connect(args.pyscenic_output, mode='r', validate=False)
    
    # Extract metadata and AUC matrix
    try:
        meta = json.loads(zlib.decompress(base64.b64decode(lf_scenic.attrs.MetaData)))
    except:
        meta = {'regulonThresholds': []}
    
    auc_mtx = pd.DataFrame(lf_scenic.ca.RegulonsAUC, index=lf_scenic.ca.CellID)
    
    # Export AUC matrix to CSV if requested
    if args.export_auc_csv:
        auc_mtx.to_csv(args.export_auc_csv, sep='\t')
    
    regulons = lf_scenic.ra.Regulons
    
    # Fix regulon names to display properly in SCope
    auc_mtx.columns = auc_mtx.columns.str.replace('\(', '_(')
    regulons.dtype.names = tuple([x.replace("(", "_(") for x in regulons.dtype.names])
    
    # Fix regulon thresholds
    rt = meta.get('regulonThresholds', [])
    for i, x in enumerate(rt):
        tmp = x.get('regulon', '').replace("(", "_(")
        x.update({'regulon': tmp})
    
    # Compute embeddings
    # UMAP
    runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
    dr_umap = runUmap(auc_mtx)
    dr_umap_df = pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index)
    
    # tSNE
    tsne = TSNE(n_jobs=args.threads)
    dr_tsne = tsne.fit_transform(auc_mtx)
    dr_tsne_df = pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index)
    
    # Create embeddings matrices
    tsneDF = pd.DataFrame(dr_tsne, columns=['_X', '_Y'], index=auc_mtx.index)
    
    Embeddings_X = pd.DataFrame(index=lf_scenic.ca.CellID)
    Embeddings_X = pd.concat([
        dr_tsne_df['X'],
        dr_umap_df['X']
    ], sort=False, axis=1, join='outer')
    Embeddings_X.columns = ['1', '2']
    
    Embeddings_Y = pd.DataFrame(index=lf_scenic.ca.CellID)
    Embeddings_Y = pd.concat([
        dr_tsne_df['Y'],
        dr_umap_df['Y']
    ], sort=False, axis=1, join='outer')
    Embeddings_Y.columns = ['1', '2']
    
    # Create metadata JSON
    metaJson = {}
    
    metaJson['embeddings'] = [
        {
            "id": -1,
            "name": "SCENIC AUC t-SNE"
        },
        {
            "id": 1,
            "name": "SCENIC AUC t-SNE"
        },
        {
            "id": 2,
            "name": "SCENIC AUC UMAP"
        }
    ]
    
    metaJson["clusterings"] = []
    
    metaJson["metrics"] = [
        {
            "name": "nUMI"
        },
        {
            "name": "nGene"
        }
    ]
    
    metaJson["annotations"] = []
    
    # SCENIC regulon thresholds
    metaJson["regulonThresholds"] = rt
    
    # Get basic metrics from input loom
    try:
        nGene = np.array(lf_input.ca.nGene if 'nGene' in lf_input.ca.keys() else np.sum(lf_input[:, :] > 0, axis=0)).flatten()
        nUMI = np.array(lf_input.ca.nUMI if 'nUMI' in lf_input.ca.keys() else np.sum(lf_input[:, :], axis=0)).flatten()
    except:
        nGene = np.array(np.sum(lf_input[:, :] > 0, axis=0)).flatten()
        nUMI = np.array(np.sum(lf_input[:, :], axis=0)).flatten()
    
    # Assemble column attributes
    col_attrs = {
        "CellID": np.array(lf_scenic.ca.CellID),
        "nUMI": nUMI,
        "nGene": nGene,
        "Embedding": dfToNamedMatrix(tsneDF),
        "Embeddings_X": dfToNamedMatrix(Embeddings_X),
        "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
        "RegulonsAUC": dfToNamedMatrix(auc_mtx),
    }
    
    row_attrs = {
        "Gene": lf_input.ra.Gene,
        "Regulons": regulons,
    }
    
    attrs = {
        "title": args.title,
        "MetaData": json.dumps(metaJson),
        "Genome": args.genome,
        "SCopeTreeL1": "",
        "SCopeTreeL2": "",
        "SCopeTreeL3": ""
    }
    
    # Compress the metadata field
    attrs['MetaData'] = base64.b64encode(
        zlib.compress(json.dumps(metaJson).encode('ascii'))
    ).decode('ascii')
    
    # Create output loom file
    lp.create(
        filename=args.output_loom,
        layers=lf_input[:, :],
        row_attrs=row_attrs,
        col_attrs=col_attrs,
        file_attrs=attrs
    )
    
    # Close connections
    lf_input.close()
    lf_scenic.close()


if __name__ == "__main__":
    main()
