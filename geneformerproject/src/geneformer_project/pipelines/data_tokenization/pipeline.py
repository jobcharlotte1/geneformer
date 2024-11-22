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
                    inputs=["cardiomyocytes_anndata_preprocessed", "cardiomyocytes_tokenized", "params:output_prefix", "params:file_format", "params:use_generator"],
                    outputs=None,
                    name="transcriptome_tokenizer_node",
                ),
        ]
        )
