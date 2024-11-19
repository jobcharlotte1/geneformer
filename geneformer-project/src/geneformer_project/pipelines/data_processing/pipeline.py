from kedro.pipeline import Pipeline, node, pipeline

from .nodes import (
    read_h5ad_file,
    read_csv_file,
    get_genes_names,
    get_cells_names,
    build_dataframe_from_expressionmatrix,
    select_highly_variable_genes,
    get_matrix_from_h5,
    rank_genes,
    tokenize_cell,
    sum_ensembl_ids,
    TranscriptomeTokenizer,
    transcriptome_tokenizer_function,
    load_and_filter
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=read_h5ad_file,
                inputs="cardiomyocytes_anndata_preprocessed",  # Input defined in the catalog
                outputs="anndata",  # Output name in memory
                name="read_anndata_node",
            ),
            node(
                func=select_highly_variable_genes,
                inputs=["anndata", "params:select_hvg_genes"],
                outputs="adata_hvg",
                name="select_hvg_node",
            ),
            node(
                func=transcriptome_tokenizer_function,
                inputs=["params:tokenize_data_params", "adata_hvg", "params:paths.tokenized_data_dir", "params: tokenized_format_params"],
                outputs=None,
                name="transcriptome_tokenizer_node",
            ),

        ]
    )
