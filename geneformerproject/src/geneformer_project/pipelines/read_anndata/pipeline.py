"""
This is a boilerplate pipeline 'read_anndata'
generated using Kedro 0.19.10
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import read_h5ad_file, read_tisch2_data, read_data


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
        
        node(
                func=read_data,
                inputs=["cardiomyocytes_anndata_raw", "h5ad", "None", "None"],  # Input defined in the catalog
                outputs="cardiomyocytes_anndata",  # Output name in memory
                name="read_anndata_cardiomyocytes_node",
            ),

        node(
                    func=read_data,
                    inputs=["PAAD_CRA001160_expression_raw", "tisch2", "PAAD_CRA001160_metainfo_raw", "param: path_PAAD_CRA001160_anndata"],
                    outputs="PAAD_CRA001160_anndata",
                    name="read_anndata_PAAD_CRA001160_node",
                ),


    ])
