from kedro.pipeline import Pipeline, node, pipeline

from .nodes import read_h5ad_file, select_highly_variable_genes


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [

            node(
                func=select_highly_variable_genes,
                inputs=["cardiomyocytes_anndata", "params:select_hvg_genes", "params:output_filepath"],
                outputs="cardiomyocytes_anndata_preprocessed",
                name="select_hvg_cardiomyocytes_node",
            ),
            
            node(
                func=select_highly_variable_genes,
                inputs=["PAAD_CRA001160_anndata", "params:select_hvg_genes", "params:output_filepath"],
                outputs="PAAD_CRA001160_anndata_preprocessed",
                name="select_hvg_PAAD_CRA001160_node",
            ),

        ]
    )
