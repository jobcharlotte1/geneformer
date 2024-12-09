"""
This is a boilerplate pipeline 'data_classification'
generated using Kedro 0.19.9
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import create_classifier_node, prepare_data_node, split_anndata


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([

        node(
            func=split_anndata,
            inputs=["PAAD_CRA001160_anndata_preprocessed", "params:id_split.PAAD_CRA001160", "params:length_test_size", "params:random_state"],
            outputs=["PAAD_CRA001160_adata_splitted", "train_test_id_split_dict_PAAD_CRA001160_Cell"],
            name="split_anndata_PAAD_CRA001160_node",
        ),

        node(
                func=create_classifier_node,
                inputs=["params:classifier_args_PAAD_CRA001160"],
                outputs="classifier_PAAD_CRA001160",
                name="initiate_classifier_PAAD_CRA001160_node",
            ),

        node(
                func=prepare_data_node,
                inputs=["PAAD_CRA001160_tokenized", "params:output_path_PAAD_CRA001160_tokenized", "params:output_prefix_classifier"],
                outputs=["PAAD_CRA001160_tokenized_trainset","PAAD_CRA001160_tokenized_testset"],
                name='prepare_data_PAAD_CRA001160_node',
            ),

    ])
