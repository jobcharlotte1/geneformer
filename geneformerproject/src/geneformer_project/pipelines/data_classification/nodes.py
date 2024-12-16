"""
This is a boilerplate pipeline 'data_classification'
generated using Kedro 0.19.9
"""
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from datasets import Dataset, load_from_disk
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    auc,
    confusion_matrix,
    f1_score,
    roc_curve,
)
import torch
import transformers
import sys
sys.path.append("/home/BS94_SUR/kedro_geneformer/geneformerproject")
from Geneformer.geneformer.classifier import Classifier
from Geneformer.geneformer.perturber_utils import filter_by_dict

def split_anndata(adata, id_split, length_test_size, random_state_nb):

    train_ids, test_val_ids = train_test_split(adata.obs_names, test_size=length_test_size, random_state=random_state_nb)
    eval_ids, test_ids = train_test_split(test_val_ids, test_size=0.5, random_state=42)
    
    train_test_id_split_dict = {"attr_key": id_split ,
                            "train": list(train_ids.values) + list(eval_ids.values),
                            "test": list(test_ids.values)}
    
    train_valid_id_split_dict = {"attr_key": "Cell",
                            "train": train_ids,
                            "eval": eval_ids}
    
    conditions = [
        (adata.obs.Cell.isin(list(train_ids.values))),
        (adata.obs.Cell.isin(list(eval_ids.values))),
        (adata.obs.Cell.isin(list(test_ids.values)))
    ]
    choices = ['Train', 'Eval', 'Test']
    adata.obs['set'] = np.select(conditions, choices, default='Other')
    adata.obs.groupby('set')[id_split].value_counts()

    return adata, train_test_id_split_dict, train_valid_id_split_dict



def create_classifier_node(classifier_type,
    cell_state_dict,
    filter_data,
    training_args,
    freeze_layers,
    num_crossval_splits,
    forward_batch_size,
    nproc):
    """
    Node to create an instance of the Classifier.
    """

    classifier = Classifier(
        classifier=classifier_type,
        cell_state_dict=cell_state_dict,
        filter_data=filter_data,
        training_args=training_args,
        freeze_layers=freeze_layers,
        num_crossval_splits=num_crossval_splits,
        forward_batch_size=forward_batch_size,
        nproc=nproc
    )

    return classifier 

def load_and_filter(filter_data, nproc, input_data_file):
    data = load_from_disk(input_data_file)
    if filter_data is not None:
        data = filter_by_dict(data, filter_data, nproc)
    return data


def prepare_data_node(
    classifier: Classifier,
    input_data: any,
    output_directory: str,
    num_classes: int = None,
    stratify_column: str = None,
    predict_col: str = None,
    test_split_size: float = 0.2,
    ):
    
    prepared_data = classifier.prepare_data(
        data=input_data,
        output_directory=output_directory,
        num_classes=num_classes,
        stratify_column=stratify_column,
        predict_col=predict_col,
        test_split_size=test_split_size,
    )

    return prepared_data

def validate_node(
        classifier: Classifier,
        model_directory: str,
        input_data: any,
        id_class_dict: dict,
        output_directory: str,
        output_prefix: str,
        split_id_dict: dict,
        save_eval_output: bool,
        predict_eval: bool,
        predict_trainer: bool,
        n_hyperopt_trials: int
):

    all_metrics = classifier.validate(
        model_directory=model_directory,
        prepared_input_data_file=input_data,
        id_class_dict_file=id_class_dict,
        output_directory=output_directory,
        output_prefix=output_prefix,
        split_id_dict=split_id_dict,
        save_eval_output=save_eval_output,
        predict_eval=predict_eval,
        predict_trainer=predict_trainer,
        n_hyperopt_trials=n_hyperopt_trials
)

    return all_metrics


def evaluate_saved_model_node(
    classifier: Classifier,
    model_directory: str,
    evaluation_data: dict,
    id_class_dict: dict,
    output_directory: str,
    output_prefix: str,
    batch_size: int = 32,
    evaluation_metrics: list = None,
):

    evaluation_results = classifier.evaluate_saved_model(
        model_directory=model_directory,
        prepared_evaluation_data_file=evaluation_data,
        id_class_dict_file=id_class_dict,
        output_directory=output_directory,
        output_prefix=output_prefix,
        batch_size=batch_size,
        evaluation_metrics=evaluation_metrics,
    )

    return evaluation_results
