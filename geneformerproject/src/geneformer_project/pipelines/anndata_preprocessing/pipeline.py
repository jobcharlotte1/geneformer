from kedro.pipeline import Pipeline, node, pipeline

from .nodes import read_h5ad_file, select_highly_variable_genes


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=read_h5ad_file,
                inputs="cardiomyocytes_anndata_raw",  # Input defined in the catalog
                outputs="cardiomyocytes_anndata",  # Output name in memory
                name="read_anndata_node",
            ),
            node(
                func=select_highly_variable_genes,
                inputs=["cardiomyocytes_anndata", "params:select_hvg_genes", "params:output_filepath"],
                outputs="cardiomyocytes_anndata_preprocessed",
                name="select_hvg_node",
            ),
            

        ]
    )
