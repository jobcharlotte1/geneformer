"""Project pipelines."""

from kedro.framework.project import find_pipelines
from kedro.pipeline import Pipeline
import geneformer_project.pipelines.anndata_preprocessing as annp
import geneformer_project.pipelines.data_tokenization as dt


def register_pipelines() -> dict[str, Pipeline]:
    """Register the project's pipelines.

    Returns:
        A mapping from pipeline names to ``Pipeline`` objects.
    """

    anndata_preprocessing_pipeline = annp.create_pipeline()
    data_tokenization_pipeline = dt.create_pipeline()

    merged_pipeline = anndata_preprocessing_pipeline + data_tokenization_pipeline

    return {
        "__default__": merged_pipeline,
        "data_preprocessing": anndata_preprocessing_pipeline,
        "data_tokenization": data_tokenization_pipeline ,
    }
