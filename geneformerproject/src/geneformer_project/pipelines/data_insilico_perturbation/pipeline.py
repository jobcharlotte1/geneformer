"""
This is a boilerplate pipeline 'data_insilico_perturbation'
generated using Kedro 0.19.10
"""

from kedro.pipeline import Pipeline, pipeline
from .nodes import read_pkl_json_file, emb_extractornode


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline([
    
        node(func=read_pkl_json_file,
             inputs=["params:gene_info_path", "pkl"],
             outputs=["gene_info"],
             name="read_pkl_json_file_node",
        ),
        node(func=read_pkl_json_file,
             inputs=["params:token_info_path", "pkl"],
             outputs=["token_info"],
             name="read_pkl_json_file_node",
        ),
        node(func=get_infos_gene,
            inputs="params:gene",
            outputs=["gene_ids", "genes_to_perturb", "token_id"],
            name="get_infos_gene_node",
        ),
        node(func=emb_extractor_node,
             inputs=["params:model_type",
                    "params.dataset.PAAD_CRA001160.nb_classes",
                    "params.states_parameters.PAAD_CRA001160.Source.normal_to_tumor.filter_data_dict",
                    "params.states_parameters.PAAD_CRA001160.Source.normal_to_tumor.max_n_cells",
                    "params.emb_layer",
                    "params.summary_stat",
                    "params.forward_batch_size",
                    "params.nproc"]
             outputs="embex_PAAD_CRA001160",
             name="embextractor_PAAD_CRA001160_node",
        ),
        node(func=get_state_embs_node,
             inputs=["embex_PAAD_CRA001160","params.states_parameters.PAAD_CRA001160.Source.normal_to_tumor.cell_states_to_model", "geneformer_model_PAAD_CRA001160", "params.dataset.PAAD_CRA001160.saved_file", "params.dataset.PAAD_CRA001160_saved_prefix"],
             outputs="state_embs_dict_PAAD_CRA001160_Source_normaltotumor",
             name="get_state_embs_PAAD_CRA001160_Source_normaltotumor_node",
        ),
        node(func=In_Silico_perturber_node,
             inputs=["params:perturbation",
                        "params:perturb_rank_shift",
                        "genes_to_perturb",
                        "params:combos",
                        "params:anchor_gene",
                        "params:model_type",
                        "params.dataset.PAAD_CRA001160.nb_classes",
                        "params:emb_mode",
                        "params:cell_emb_style",
                        "params.states_parameters.PAAD_CRA001160.Source.normal_to_tumor.filter_data_dict",
                        "params.states_parameters.PAAD_CRA001160.Source.normal_to_tumor.cell_states_to_model",
                        "state_embs_dict_PAAD_CRA001160_Source_normaltotumor",
                        "params.states_parameters.PAAD_CRA001160.Source.normal_to_tumor.max_n_cells",
                        "params:emb_layer",
                        "params:forward_batch_size",
                        "params:nproc"],
             outputs="isp_PAAD_CRA001160_Source_normaltotumor",
             name="in_silico_perturber_PAAD_CRA001160_Source_normaltotumor_node",
        ),
        node(func=perturb_data,
             inputs=["isp_PAAD_CRA001160_Source_normaltotumor",
                    "params.dataset.PAAD_CRA001160.path_model",
                    "params.dataset.PAAD_CRA001160.path_data_tokenized",
                    "params.dataset.PAAD_CRA001160.output_dir",],
        ),
    ])
