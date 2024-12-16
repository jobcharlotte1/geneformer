"""
This is a boilerplate pipeline 'data_classification'
generated using Kedro 0.19.9
"""

from kedro.pipeline import Pipeline, pipeline, node
from .nodes import create_classifier_node, prepare_data_node, split_anndata, validate_node, evaluate_saved_model_node


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([

        node(
            func=split_anndata,
            inputs=["PAAD_CRA001160_anndata_preprocessed", "params:id_split.PAAD_CRA001160", "params:length_test_size", "params:random_state"],
            outputs=["PAAD_CRA001160_adata_splitted", "train_test_id_split_dict_PAAD_CRA001160_Cell", "train_valid_id_split_dict_PAAD_CRA001160_Cell"],
            name="split_anndata_PAAD_CRA001160_node",
        ),

        node(
                func=create_classifier_node,
                inputs=['Cell', "['Ductal cell type 1', 'Ductal cell type 2']", "None", "params:training_args.PAAD_CRA001160", "2", "1", "100", "96"],
                outputs="classifier_PAAD_CRA001160",
                name="initiate_classifier_PAAD_CRA001160_node",
            ),

        node(
                func=prepare_data_node,
                inputs=["classifier_PAAD_CRA001160", "PAAD_CRA001160_tokenized", "params:output_path_PAAD_CRA001160_tokenized", "None", "None", "None", "0.2"],
                outputs=["PAAD_CRA001160_tokenized_trainset","PAAD_CRA001160_tokenized_testset", "PAAD_CRA001160_cm_classifier_test_id_class_dict"],
                name='prepare_data_PAAD_CRA001160_node',
            ),
        
        node(
                func=validate_node,
                inputs=["classifier_PAAD_CRA001160", "params:model_directory", "PAAD_CRA001160_tokenized_trainset", "PAAD_CRA001160_cm_classifier_test_id_class_dict", 
                       "params:validate_args_PAAD_CRA001160.output_directory", "params:validate_args_PAAD_CRA001160.output_prefix", "train_valid_id_split_dict_PAAD_CRA001160_Cell",
                       "params:validate_args_PAAD_CRA001160.save_eval_output", "params:validate_args_PAAD_CRA001160.predict_eval", "params:validate_args_PAAD_CRA001160.predict_trainer",
                       "params:validate_args_PAAD_CRA001160.n_hyperopt_trials"],
                outputs=["geneformer_model_PAAD_CRA001160", "all_metrics_PAAD_CRA001160"],
                name="validate_PAAD_CRA001160_node",
        ),
        
        node(
                func=evaluate_saved_model_node,
                inputs=["classifier_PAAD_CRA001160", "geneformer_model_PAAD_CRA001160", "PAAD_CRA001160_tokenized_testset", "PAAD_CRA001160_cm_classifier_test_id_class_dict", "params:validate_args_PAAD_CRA001160.output_directory", "params:validate_args_PAAD_CRA001160.output_prefix", "params:evaluate_model_PAAD_CRA001160.batch_size"],
                outputs="all_metrics_test_PAAD_CRA001160",
                name="evaluate_saved_model_PAAD_CRA001160_node",
        ),

    ])
