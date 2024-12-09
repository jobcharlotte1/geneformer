"""
This is a boilerplate pipeline 'read_anndata'
generated using Kedro 0.19.10
"""
from pathlib import Path
from typing import Dict, List, Optional, Union
from typing import Literal
import pandas as pd
import scipy.sparse as sp_sparse
import tables
import collections
from anndata import AnnData
from anndata import read_h5ad
import numpy as np 
import scanpy as sc


def read_h5ad_file(file_path):
    """
    Reads an .h5ad file and returns an AnnData object.

    Args:
        filepath: Path to the .h5ad file.

    Returns:
        AnnData object containing the loaded data.
    """
    adata = sc.read_h5ad(file_path)

    return adata

def get_matrix_from_h5(filename):
    CountMatrix = collections.namedtuple('CountMatrix', ['feature_ref', 'barcodes', 'matrix'])
    with tables.open_file(filename, 'r') as f:
        mat_group = f.get_node(f.root, 'matrix')
        barcodes = f.get_node(mat_group, 'barcodes').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
         
        feature_ref = {}
        feature_group = f.get_node(mat_group, 'features')
        feature_ref['id'] = getattr(feature_group, 'id').read()
        feature_ref['name'] = getattr(feature_group, 'name').read()
        feature_ref['feature_type'] = getattr(feature_group, 'feature_type').read()
         
        return CountMatrix(feature_ref, barcodes, matrix)
    
    

def read_tisch2_data(file_path, metainfo_file, output_filepath):

    """
    """
    metainfo = pd.read_csv(metainfo_file, delimiter='\t')

    bc_matrix = get_matrix_from_h5(file_path)
    gene_names = [v.decode('UTF-8') for v in bc_matrix.feature_ref['id']]
    cell_names = [v.decode('UTF-8') for v in bc_matrix.barcodes]

    df_data = pd.DataFrame(bc_matrix.matrix.toarray(), index=gene_names, columns=cell_names).T

    adata = sc.AnnData(X=df_data.values)  # Use the values from the DataFrame
    adata.obs = metainfo  # Add the metadata to obs (observations)
    adata.var_names = df_data.columns  # Assign gene names to var_names (variables)
    adata.obs_names = df_data.index 
    adata.obs.rename(columns={'n_count': 'n_counts'}, inplace=True)

    adata.write(Path(output_filepath))

    return adata

def read_data(file_path, format, metainfo_file, output_filepath):

    if format == 'tisch2':
        adata = read_tisch2_data(file_path, metainfo_file, output_filepath)
    else:
        adata = sc.read_h5ad(file_path)
    return adata 