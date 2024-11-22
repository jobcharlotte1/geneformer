
from pathlib import Path
from typing import Dict, List, Optional, Union
import pandas as pd
from anndata import AnnData
from anndata import read_h5ad
import numpy as np
import seaborn as sns
import scipy.sparse as sp_sparse
from typing import Literal
from datasets import Dataset, load_from_disk
import scanpy as sc
import logging
logger = logging.getLogger(__name__)


def read_h5ad_file(filepath):
    """
    Reads an .h5ad file and returns an AnnData object.

    Args:
        filepath: Path to the .h5ad file.

    Returns:
        AnnData object containing the loaded data.
    """
    return sc.read_h5ad(filepath)



def select_highly_variable_genes(adata, n_top_genes, filepath):
    """
    select top highly variable genes using scanpy function 

    Args:
        adata file and number of top genes

    Returns:
        adata filtered by number of top genes
    """

    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=n_top_genes)
    new_adata = adata[:, adata.var['highly_variable']]
    new_adata.write(Path(filepath))
    
    return new_adata


