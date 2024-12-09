"""
This is a boilerplate pipeline 'data_tokenization'
generated using Kedro 0.19.9
"""
from kedro.pipeline import Pipeline, node, pipeline

from .nodes import (
    rank_genes,
    tokenize_cell,
    sum_ensembl_ids,
    transcriptome_tokenizer_function,
    load_and_filter
)


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [

            node(
                    func=transcriptome_tokenizer_function,
                    inputs=["cardiomyocytes_anndata_preprocessed", "params:output_path_cardiomyocytes_tokenized", "params:output_prefix", "params:file_format", "params:use_generator"],
                    outputs="cardiomyocytes_tokenized",
                    name="transcriptome_tokenizer_cardiomyocytes_node",
                ),

            node(
                    func=transcriptome_tokenizer_function,
                    inputs=["PAAD_CRA001160_anndata_preprocessed", "params:output_path_PAAD_CRA001160_tokenized", "params:output_prefix", "params:file_format", "params:use_generator"],
                    outputs="PAAD_CRA001160_tokenized",
                    name="transcriptome_tokenizer_PAAD_CRA001160_node",
                ),
        ]
        )
