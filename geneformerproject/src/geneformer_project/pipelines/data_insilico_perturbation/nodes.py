"""
This is a boilerplate pipeline 'data_insilico_perturbation'
generated using Kedro 0.19.10
"""
import sys
import os
from pathlib import Path
from typing import Literal
import json
import pickle
import logging
import math
import pandas as pd
from anndata import AnnData
from anndata import read_h5ad
import datetime
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import anndata as ad
from typing import Dict, List, Optional, Union
import scipy.sparse as sp
import numpy as np
import umap
import pacmap
from matplotlib.patches import Circle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from sklearn.manifold import TSNE, Isomap, LocallyLinearEmbedding
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import binned_statistic_2d
from scipy.interpolate import griddata
from scipy.stats import gaussian_kde
from anndata import AnnData
from anndata import read_h5ad
import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
import anndata as ad
from typing import Dict, List, Optional, Union
from collections import Counter, defaultdict
import scipy.sparse as sp
import numpy as np
import datasets
from datasets import Dataset
import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    accuracy_score,
    auc,
    confusion_matrix,
    f1_score,
    roc_curve,
)
from sklearn import preprocessing
from tqdm.auto import tqdm, trange
from transformers import Trainer
from transformers.training_args import TrainingArguments
import umap
from transformers import (
    BertForMaskedLM,
    BertForSequenceClassification,
    BertForTokenClassification,
)
from torch.nn.parallel import DistributedDataParallel as DDP
import time
import gseapy as gp
import umap


sys.path.append("/home/BS94_SUR/kedro_geneformer/geneformerproject")
from Geneformer.geneformer import perturber_utils as pu
from Geneformer.geneformer.emb_extractor import get_embs
from Geneformer.geneformer.tokenizer import TOKEN_DICTIONARY_FILE
from Geneformer.geneformer import EmbExtractor


def read_pkl_json_file(path_file, type_file):
    if type_file == 'json':
        with open(path_file, 'r') as file:
            data = json.load(file)
    if type_file == 'pkl':
        with open(path_file, "rb") as f:
            data = pickle.load(f)
    return data 

def emb_extractor_node(model_type,
    nb_classes,
    filter_data_dict,
    max_n_cells,
    emb_layer,
    summary_stat,
    forward_batch_size,
    nproc):
    
    embex = EmbExtractor(model_type=model_type,
                     num_classes=nb_classes,
                     filter_data=filter_data_dict,
                     max_ncells=max_n_cells ,
                     emb_layer=emb_layer,
                     summary_stat=summary_stat,
                     forward_batch_size=forward_batch_size,
                     nproc=nproc)

    return embex


def get_state_embs_node(embex,
                    cell_states_to_model,
                       path_model,
                       saved_file,
                       saved_prefix)

    state_embs_dict = embex.get_state_embs(cell_states_to_model,
                            path_model,
                              path_dataset,
                              saved_file,
                              saved_prefix)
    
    return state_embs_dict

def get_infos_gene(gene):
    gene_ids = [gene]
    genes_to_perturb = [gene_info[gene] for gene in gene_ids]
    token_id = token_info[genes_to_perturb[0]]
    return gene_ids, genes_to_perturb, token_id



class InSilicoPerturber_custom:
    valid_option_dict = {
        "perturb_type": {"delete", "overexpress", "inhibit", "activate"},
        "perturb_rank_shift": {None, 1, 2, 3},
        "genes_to_perturb": {"all", list},
        "combos": {0, 1},
        "anchor_gene": {None, str},
        "model_type": {"Pretrained", "GeneClassifier", "CellClassifier"},
        "num_classes": {int},
        "emb_mode": {"cell", "cell_and_gene"},
        "cell_emb_style": {"mean_pool"},
        "filter_data": {None, dict},
        "cell_states_to_model": {None, dict},
        "state_embs_dict": {None, dict},
        "max_ncells": {None, int},
        "cell_inds_to_perturb": {"all", dict},
        "emb_layer": {-1, 0},
        "forward_batch_size": {int},
        "nproc": {int},
    }

    def __init__(
        self,
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb="all",
        combos=0,
        anchor_gene=None,
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell",
        cell_emb_style="mean_pool",
        filter_data=None,
        cell_states_to_model=None,
        state_embs_dict=None,
        max_ncells=None,
        cell_inds_to_perturb="all",
        emb_layer=-1,
        forward_batch_size=100,
        nproc=4,
        token_dictionary_file=TOKEN_DICTIONARY_FILE,
    ):

        self.perturb_type = perturb_type
        self.perturb_rank_shift = perturb_rank_shift
        self.genes_to_perturb = genes_to_perturb
        self.combos = combos
        self.anchor_gene = anchor_gene
        if self.genes_to_perturb == "all":
            self.perturb_group = False
        else:
            self.perturb_group = True
            if (self.anchor_gene is not None) or (self.combos != 0):
                self.anchor_gene = None
                self.combos = 0
                logger.warning(
                    "anchor_gene set to None and combos set to 0. "
                    "If providing list of genes to perturb, "
                    "list of genes_to_perturb will be perturbed together, "
                    "without anchor gene or combinations."
                )
        self.model_type = model_type
        self.num_classes = num_classes
        self.emb_mode = emb_mode
        self.cell_emb_style = cell_emb_style
        self.filter_data = filter_data
        self.cell_states_to_model = cell_states_to_model
        self.state_embs_dict = state_embs_dict
        self.max_ncells = max_ncells
        self.cell_inds_to_perturb = cell_inds_to_perturb
        self.emb_layer = emb_layer
        self.forward_batch_size = forward_batch_size
        self.nproc = nproc

        self.validate_options()

        # load token dictionary (Ensembl IDs:token)
        with open(token_dictionary_file, "rb") as f:
            self.gene_token_dict = pickle.load(f)

        self.pad_token_id = self.gene_token_dict.get("<pad>")

        if self.anchor_gene is None:
            self.anchor_token = None
        else:
            try:
                self.anchor_token = [self.gene_token_dict[self.anchor_gene]]
            except KeyError:
                logger.error(f"Anchor gene {self.anchor_gene} not in token dictionary.")
                raise

        if self.genes_to_perturb == "all":
            self.tokens_to_perturb = "all"
        else:
            missing_genes = [
                gene
                for gene in self.genes_to_perturb
                if gene not in self.gene_token_dict.keys()
            ]
            if len(missing_genes) == len(self.genes_to_perturb):
                logger.error(
                    "None of the provided genes to perturb are in token dictionary."
                )
                raise
            elif len(missing_genes) > 0:
                logger.warning(
                    f"Genes to perturb {missing_genes} are not in token dictionary."
                )
            self.tokens_to_perturb = [
                self.gene_token_dict.get(gene) for gene in self.genes_to_perturb
            ]

    def validate_options(self):
        # first disallow options under development
        if self.perturb_type in ["inhibit", "activate"]:
            logger.error(
                "In silico inhibition and activation currently under development. "
                "Current valid options for 'perturb_type': 'delete' or 'overexpress'"
            )
            raise
        if (self.combos > 0) and (self.anchor_gene is None):
            logger.error(
                "Combination perturbation without anchor gene is currently under development. "
                "Currently, must provide anchor gene for combination perturbation."
            )
            raise

        # confirm arguments are within valid options and compatible with each other
        for attr_name, valid_options in self.valid_option_dict.items():
            attr_value = self.__dict__[attr_name]
            if type(attr_value) not in {list, dict}:
                if attr_value in valid_options:
                    continue
                if attr_name in ["anchor_gene"]:
                    if type(attr_name) in {str}:
                        continue
            valid_type = False
            for option in valid_options:
                if (option in [bool, int, list, dict]) and isinstance(
                    attr_value, option
                ):
                    valid_type = True
                    break
            if valid_type:
                continue
            logger.error(
                f"Invalid option for {attr_name}. "
                f"Valid options for {attr_name}: {valid_options}"
            )
            raise

        if self.perturb_type in ["delete", "overexpress"]:
            if self.perturb_rank_shift is not None:
                if self.perturb_type == "delete":
                    logger.warning(
                        "perturb_rank_shift set to None. "
                        "If perturb type is delete then gene is deleted entirely "
                        "rather than shifted by quartile"
                    )
                elif self.perturb_type == "overexpress":
                    logger.warning(
                        "perturb_rank_shift set to None. "
                        "If perturb type is overexpress then gene is moved to front "
                        "of rank value encoding rather than shifted by quartile"
                    )
            self.perturb_rank_shift = None

        if (self.anchor_gene is not None) and (self.emb_mode == "cell_and_gene"):
            self.emb_mode = "cell"
            logger.warning(
                "emb_mode set to 'cell'. "
                "Currently, analysis with anchor gene "
                "only outputs effect on cell embeddings."
            )

        if self.cell_states_to_model is not None:
            pu.validate_cell_states_to_model(self.cell_states_to_model)

            if self.anchor_gene is not None:
                self.anchor_gene = None
                logger.warning(
                    "anchor_gene set to None. "
                    "Currently, anchor gene not available "
                    "when modeling multiple cell states."
                )

            if self.state_embs_dict is None:
                logger.error(
                    "state_embs_dict must be provided for mode with cell_states_to_model. "
                    "Format is dictionary with keys specifying each possible cell state to model. "
                    "Values are target embedding positions as torch.tensor."
                )
                raise

            for state_emb in self.state_embs_dict.values():
                if not torch.is_tensor(state_emb):
                    logger.error(
                        "state_embs_dict must be dictionary with values being torch.tensor."
                    )
                    raise

            keys_absent = []
            for k, v in self.cell_states_to_model.items():
                if (k == "start_state") or (k == "goal_state"):
                    if v not in self.state_embs_dict.keys():
                        keys_absent.append(v)
                if k == "alt_states":
                    for state in v:
                        if state not in self.state_embs_dict.keys():
                            keys_absent.append(state)
            if len(keys_absent) > 0:
                logger.error(
                    "Each start_state, goal_state, and alt_states in cell_states_to_model "
                    "must be a key in state_embs_dict with the value being "
                    "the state's embedding position as torch.tensor. "
                    f"Missing keys: {keys_absent}"
                )
                raise

        if self.perturb_type in ["inhibit", "activate"]:
            if self.perturb_rank_shift is None:
                logger.error(
                    "If perturb_type is inhibit or activate then "
                    "quartile to shift by must be specified."
                )
                raise

        if self.filter_data is not None:
            for key, value in self.filter_data.items():
                if not isinstance(value, list):
                    self.filter_data[key] = [value]
                    logger.warning(
                        "Values in filter_data dict must be lists. "
                        f"Changing {key} value to list ([{value}])."
                    )

        if self.cell_inds_to_perturb != "all":
            if set(self.cell_inds_to_perturb.keys()) != {"start", "end"}:
                logger.error(
                    "If cell_inds_to_perturb is a dictionary, keys must be 'start' and 'end'."
                )
                raise
            if (
                self.cell_inds_to_perturb["start"] < 0
                or self.cell_inds_to_perturb["end"] < 0
            ):
                logger.error("cell_inds_to_perturb must be positive.")
                raise

    def perturb_data(
        self, model_directory, input_data_file, output_directory, output_prefix
    ):
        """
        Perturb genes in input data and save as results in output_directory.

        **Parameters:**

        model_directory : Path
            | Path to directory containing model
        input_data_file : Path
            | Path to directory containing .dataset inputs
        output_directory : Path
            | Path to directory where perturbation data will be saved as batched pickle files
        output_prefix : str
            | Prefix for output files
        """

        ### format output path ###
        output_path_prefix = os.path.join(
            output_directory, f"in_silico_{self.perturb_type}_{output_prefix}"
        )

        ### load model and define parameters ###
        model = pu.load_model(
            self.model_type, self.num_classes, model_directory, mode="eval"
        )
        self.max_len = pu.get_model_input_size(model)
        layer_to_quant = pu.quant_layers(model) + self.emb_layer

        ### filter input data ###
        # general filtering of input data based on filter_data argument
        filtered_input_data = pu.load_and_filter(
            self.filter_data, self.nproc, input_data_file
        )
        filtered_input_data = self.apply_additional_filters(filtered_input_data)

        if self.perturb_group is True:
            self.isp_perturb_set(
                model, filtered_input_data, layer_to_quant, output_path_prefix
            )
        else:
            #print('do perturbation embeddings')
            self.isp_perturb_all(
                model, filtered_input_data, layer_to_quant, output_path_prefix
            )

    def apply_additional_filters(self, filtered_input_data):
        # additional filtering of input data dependent on isp mode
        if self.cell_states_to_model is not None:
            # filter for cells with start_state and log result
            filtered_input_data = pu.filter_data_by_start_state(
                filtered_input_data, self.cell_states_to_model, self.nproc
            )

        if (self.tokens_to_perturb != "all") and (self.perturb_type != "overexpress"):
            # filter for cells with tokens_to_perturb and log result
            filtered_input_data = pu.filter_data_by_tokens_and_log(
                filtered_input_data,
                self.tokens_to_perturb,
                self.nproc,
                "genes_to_perturb",
            )

        if self.anchor_token is not None:
            # filter for cells with anchor gene and log result
            filtered_input_data = pu.filter_data_by_tokens_and_log(
                filtered_input_data, self.anchor_token, self.nproc, "anchor_gene"
            )

        # downsample and sort largest to smallest to encounter memory constraints earlier
        filtered_input_data = pu.downsample_and_sort(
            filtered_input_data, self.max_ncells
        )

        # slice dataset if cells_inds_to_perturb is not "all"
        if self.cell_inds_to_perturb != "all":
            filtered_input_data = pu.slice_by_inds_to_perturb(
                filtered_input_data, self.cell_inds_to_perturb
            )

        return filtered_input_data

    def isp_perturb_set(
        self,
        model,
        filtered_input_data: Dataset,
        layer_to_quant: int,
        output_path_prefix: str,
    ):
        def make_group_perturbation_batch(example):
            example_input_ids = example["input_ids"]
            example["tokens_to_perturb"] = self.tokens_to_perturb
            indices_to_perturb = [
                example_input_ids.index(token) if token in example_input_ids else None
                for token in self.tokens_to_perturb
            ]
            indices_to_perturb = [
                item for item in indices_to_perturb if item is not None
            ]
            if len(indices_to_perturb) > 0:
                example["perturb_index"] = indices_to_perturb
            else:
                # -100 indicates tokens to overexpress are not present in rank value encoding
                example["perturb_index"] = [-100]
            if self.perturb_type == "delete":
                example = pu.delete_indices(example)
            elif self.perturb_type == "overexpress":
                example = pu.overexpress_tokens(example, self.max_len)
                example["n_overflow"] = pu.calc_n_overflow(
                    self.max_len,
                    example["length"],
                    self.tokens_to_perturb,
                    indices_to_perturb,
                )
            return example

        total_batch_length = len(filtered_input_data)
        if self.cell_states_to_model is None:
            cos_sims_dict = defaultdict(list)
        else:
            cos_sims_dict = {
                state: defaultdict(list)
                for state in pu.get_possible_states(self.cell_states_to_model)
            }

        perturbed_data = filtered_input_data.map(
            make_group_perturbation_batch, num_proc=self.nproc
        )
        if self.perturb_type == "overexpress":
            filtered_input_data = filtered_input_data.add_column(
                "n_overflow", perturbed_data["n_overflow"]
            )
            # remove overflow genes from original data so that embeddings are comparable
            # i.e. if original cell has genes 0:2047 and you want to overexpress new gene 2048,
            # then the perturbed cell will be 2048+0:2046 so we compare it to an original cell 0:2046.
            # (otherwise we will be modeling the effect of both deleting 2047 and adding 2048,
            # rather than only adding 2048)
            filtered_input_data = filtered_input_data.map(
                pu.truncate_by_n_overflow, num_proc=self.nproc
            )

        if self.emb_mode == "cell_and_gene":
            stored_gene_embs_dict = defaultdict(list)
        
        original_emb_dict = {}
        perturbed_emb_dict = {}
        cell_cossim_dict = {}

        # iterate through batches
        for i in trange(0, total_batch_length, self.forward_batch_size):
            max_range = min(i + self.forward_batch_size, total_batch_length)
            inds_select = [i for i in range(i, max_range)]

            minibatch = filtered_input_data.select(inds_select)
            perturbation_batch = perturbed_data.select(inds_select)

            if self.cell_emb_style == "mean_pool":
                full_original_emb = get_embs(
                    model,
                    minibatch,
                    "gene",
                    layer_to_quant,
                    self.pad_token_id,
                    self.forward_batch_size,
                    summary_stat=None,
                    silent=True,
                )
                indices_to_perturb = perturbation_batch["perturb_index"]
                # remove indices that were perturbed
                original_emb = pu.remove_perturbed_indices_set(
                    full_original_emb,
                    self.perturb_type,
                    indices_to_perturb,
                    self.tokens_to_perturb,
                    minibatch["length"],
                )
                full_perturbation_emb = get_embs(
                    model,
                    perturbation_batch,
                    "gene",
                    layer_to_quant,
                    self.pad_token_id,
                    self.forward_batch_size,
                    summary_stat=None,
                    silent=True,
                )

                # remove overexpressed genes
                if self.perturb_type == "overexpress":
                    perturbation_emb = full_perturbation_emb[
                        :, len(self.tokens_to_perturb) :, :
                    ]

                elif self.perturb_type == "delete":
                    perturbation_emb = full_perturbation_emb[
                        :, : max(perturbation_batch["length"]), :
                    ]
                    

                n_perturbation_genes = perturbation_emb.size()[1]

                # if no goal states, the cosine similarties are the mean of gene cosine similarities
                if (
                    self.cell_states_to_model is None
                    or self.emb_mode == "cell_and_gene"
                ):
                    gene_cos_sims = pu.quant_cos_sims(
                        perturbation_emb,
                        original_emb,
                        self.cell_states_to_model,
                        self.state_embs_dict,
                        emb_mode="gene",
                    )

                # if there are goal states, the cosine similarities are the cell cosine similarities
                if self.cell_states_to_model is not None:
                    original_cell_emb = pu.mean_nonpadding_embs(
                        full_original_emb,
                        torch.tensor(minibatch["length"], device="cuda"),
                        dim=1,
                    )
                    original_emb_dict[i] = original_cell_emb
                    
                    perturbation_cell_emb = pu.mean_nonpadding_embs(
                        full_perturbation_emb,
                        torch.tensor(perturbation_batch["length"], device="cuda"),
                        dim=1,
                    )
                    
                    perturbed_emb_dict[i] = perturbation_cell_emb
                    
                    cell_cos_sims = pu.quant_cos_sims(
                        perturbation_cell_emb,
                        original_cell_emb,
                        self.cell_states_to_model,
                        self.state_embs_dict,
                        emb_mode="cell",
                    )
                    
                    cell_cossim_dict[i] = cell_cos_sims

                # get cosine similarities in gene embeddings
                # if getting gene embeddings, need gene names
                if self.emb_mode == "cell_and_gene":
                    gene_list = minibatch["input_ids"]
                    # need to truncate gene_list
                    gene_list = [
                        [g for g in genes if g not in self.tokens_to_perturb][
                            :n_perturbation_genes
                        ]
                        for genes in gene_list
                    ]

                    for cell_i, genes in enumerate(gene_list):
                        for gene_j, affected_gene in enumerate(genes):
                            if len(self.genes_to_perturb) > 1:
                                tokens_to_perturb = tuple(self.tokens_to_perturb)
                            else:
                                tokens_to_perturb = self.tokens_to_perturb[0]

                            # fill in the gene cosine similarities
                            try:
                                stored_gene_embs_dict[
                                    (tokens_to_perturb, affected_gene)
                                ].append(gene_cos_sims[cell_i, gene_j].item())
                            except KeyError:
                                stored_gene_embs_dict[
                                    (tokens_to_perturb, affected_gene)
                                ] = gene_cos_sims[cell_i, gene_j].item()
                else:
                    gene_list = None

            if self.cell_states_to_model is None:
                # calculate the mean of the gene cosine similarities for cell shift
                # tensor of nonpadding lengths for each cell
                if self.perturb_type == "overexpress":
                    # subtract number of genes that were overexpressed
                    # since they are removed before getting cos sims
                    n_overexpressed = len(self.tokens_to_perturb)
                    nonpadding_lens = [
                        x - n_overexpressed for x in perturbation_batch["length"]
                    ]
                else:
                    nonpadding_lens = perturbation_batch["length"]
                cos_sims_data = pu.mean_nonpadding_embs(
                    gene_cos_sims, torch.tensor(nonpadding_lens, device="cuda")
                )
                cos_sims_dict = self.update_perturbation_dictionary(
                    cos_sims_dict,
                    cos_sims_data,
                    filtered_input_data,
                    indices_to_perturb,
                    gene_list,
                )
            else:
                cos_sims_data = cell_cos_sims
                for state in cos_sims_dict.keys():
                    cos_sims_dict[state] = self.update_perturbation_dictionary(
                        cos_sims_dict[state],
                        cos_sims_data[state],
                        filtered_input_data,
                        indices_to_perturb,
                        gene_list,
                    )
            del minibatch
            del perturbation_batch
            del original_emb
            del perturbation_emb
            del cos_sims_data

            torch.cuda.empty_cache()
            time.sleep(5)

        pu.write_perturbation_dictionary(
            cos_sims_dict,
            f"{output_path_prefix}_cell_embs_dict_{self.tokens_to_perturb}",
        )
        
        pu.write_perturbation_dictionary(
            original_emb_dict,
            f"{output_path_prefix}_cell_embs_original_{self.tokens_to_perturb}",
        )
        
        pu.write_perturbation_dictionary(
            perturbed_emb_dict,
            f"{output_path_prefix}_cell_embs_perturbed_{self.tokens_to_perturb}",
        )
        
        pu.write_perturbation_dictionary(
            cell_cossim_dict,
            f"{output_path_prefix}_cell_cossim_dict_{self.tokens_to_perturb}",
        )

        if self.emb_mode == "cell_and_gene":
            pu.write_perturbation_dictionary(
                stored_gene_embs_dict,
                f"{output_path_prefix}_gene_embs_dict_{self.tokens_to_perturb}",
            )
        torch.cuda.empty_cache()
        time.sleep(5)
        
        return perturbed_emb_dict, cell_cossim_dict

    def isp_perturb_all(
        self,
        model,
        filtered_input_data: Dataset,
        layer_to_quant: int,
        output_path_prefix: str,
    ):
        pickle_batch = -1
        if self.cell_states_to_model is None:
            cos_sims_dict = defaultdict(list)
        else:
            cos_sims_dict = {
                state: defaultdict(list)
                for state in pu.get_possible_states(self.cell_states_to_model)
            }

        if self.emb_mode == "cell_and_gene":
            stored_gene_embs_dict = defaultdict(list)
        for i in trange(len(filtered_input_data)):
            example_cell = filtered_input_data.select([i])
            full_original_emb = get_embs(
                model,
                example_cell,
                "gene",
                layer_to_quant,
                self.pad_token_id,
                self.forward_batch_size,
                summary_stat=None,
                silent=True,
            )

            # gene_list is used to assign cos sims back to genes
            # need to remove the anchor gene
            gene_list = example_cell["input_ids"][0][:]
            if self.anchor_token is not None:
                for token in self.anchor_token:
                    gene_list.remove(token)

            perturbation_batch, indices_to_perturb = pu.make_perturbation_batch(
                example_cell,
                self.perturb_type,
                self.tokens_to_perturb,
                self.anchor_token,
                self.combos,
                self.nproc,
            )

            full_perturbation_emb = get_embs(
                model,
                perturbation_batch,
                "gene",
                layer_to_quant,
                self.pad_token_id,
                self.forward_batch_size,
                summary_stat=None,
                silent=True,
            )
            

            num_inds_perturbed = 1 + self.combos
            # need to remove overexpressed gene to quantify cosine shifts
            if self.perturb_type == "overexpress":
                perturbation_emb = full_perturbation_emb[:, num_inds_perturbed:, :]
                gene_list = gene_list[
                    num_inds_perturbed:
                ]  # index 0 is not overexpressed

            elif self.perturb_type == "delete":
                perturbation_emb = full_perturbation_emb

            

            if self.cell_states_to_model is None or self.emb_mode == "cell_and_gene":
                original_batch = pu.make_comparison_batch(
                full_original_emb, indices_to_perturb, perturb_group=False
                )
                
                gene_cos_sims = pu.quant_cos_sims(
                    perturbation_emb,
                    original_batch,
                    self.cell_states_to_model,
                    self.state_embs_dict,
                    emb_mode="gene",
                )
                
            if self.cell_states_to_model is not None:
                original_cell_emb = pu.compute_nonpadded_cell_embedding(
                    full_original_emb, "mean_pool"
                )
                perturbation_cell_emb = pu.compute_nonpadded_cell_embedding(
                    full_perturbation_emb, "mean_pool"
                )

                cell_cos_sims = pu.quant_cos_sims(
                    perturbation_cell_emb,
                    original_cell_emb,
                    self.cell_states_to_model,
                    self.state_embs_dict,
                    emb_mode="cell",
                )
                
                del original_cell_emb
                del perturbation_cell_emb
                del full_original_emb
                del full_perturbation_emb
            
            torch.cuda.empty_cache()
            
            if self.emb_mode == "cell_and_gene":
                # remove perturbed index for gene list
                perturbed_gene_dict = {
                    gene: gene_list[:i] + gene_list[i + 1 :]
                    for i, gene in enumerate(gene_list)
                }

                for perturbation_i, perturbed_gene in enumerate(gene_list):
                    for gene_j, affected_gene in enumerate(
                        perturbed_gene_dict[perturbed_gene]
                    ):
                        try:
                            stored_gene_embs_dict[
                                (perturbed_gene, affected_gene)
                            ].append(gene_cos_sims[perturbation_i, gene_j].item())
                        except KeyError:
                            stored_gene_embs_dict[
                                (perturbed_gene, affected_gene)
                            ] = gene_cos_sims[perturbation_i, gene_j].item()

            if self.cell_states_to_model is None:
                cos_sims_data = torch.mean(gene_cos_sims, dim=1).detach().cpu()
                cos_sims_dict = self.update_perturbation_dictionary(
                    cos_sims_dict,
                    cos_sims_data,
                    filtered_input_data,
                    indices_to_perturb,
                    gene_list,
                )
            else:
                cos_sims_data = cell_cos_sims
                del cell_cos_sims
                for state in cos_sims_dict.keys():
                    cos_sims_dict[state] = self.update_perturbation_dictionary(
                        cos_sims_dict[state],
                        cos_sims_data[state],
                        filtered_input_data,
                        indices_to_perturb,
                        gene_list,
                    )
            

            # save dict to disk every 100 cells
            if i % 10 == 0:
                print('save dict every 1 cell')
                pu.write_perturbation_dictionary(
                    cos_sims_dict,
                    f"{output_path_prefix}_dict_cell_embs_1Kbatch{pickle_batch}",
                )
                if self.emb_mode == "cell_and_gene":
                    pu.write_perturbation_dictionary(
                        stored_gene_embs_dict,
                        f"{output_path_prefix}_dict_gene_embs_1Kbatch{pickle_batch}",
                    )
            torch.cuda.empty_cache()
            
            print(i)
            print(i % 10)

            # reset and clear memory every 1000 cells
            if i % 10 == 0:
                print('reset and clear memory every  cell')
                pickle_batch += 1
                if self.cell_states_to_model is None:
                    cos_sims_dict = defaultdict(list)
                else:
                    cos_sims_dict = {
                        state: defaultdict(list)
                        for state in pu.get_possible_states(self.cell_states_to_model)
                    }
                

                if self.emb_mode == "cell_and_gene":
                    stored_gene_embs_dict = defaultdict(list)

            torch.cuda.empty_cache()
            time.sleep(5)
                
        pu.write_perturbation_dictionary(
            cos_sims_dict, f"{output_path_prefix}_dict_cell_embs_1Kbatch{pickle_batch}"
        )
        
        torch.cuda.empty_cache()
        time.sleep(5)

        if self.emb_mode == "cell_and_gene":
            pu.write_perturbation_dictionary(
                stored_gene_embs_dict,
                f"{output_path_prefix}_dict_gene_embs_1Kbatch{pickle_batch}",
            )
        
        del cos_sims_data
        del cos_sims_dict
        
        #del stored_gene_embs_dict

    def update_perturbation_dictionary(
        self,
        cos_sims_dict: defaultdict,
        cos_sims_data: torch.Tensor,
        filtered_input_data: Dataset,
        indices_to_perturb: List[List[int]],
        gene_list=None,
    ):
        if gene_list is not None and cos_sims_data.shape[0] != len(gene_list):
            logger.error(
                f"len(cos_sims_data.shape[0]) != len(gene_list). \n \
                            cos_sims_data.shape[0] = {cos_sims_data.shape[0]}.\n \
                            len(gene_list) = {len(gene_list)}."
            )
            raise

        if self.perturb_group is True:
            if len(self.tokens_to_perturb) > 1:
                perturbed_genes = tuple(self.tokens_to_perturb)
            else:
                perturbed_genes = self.tokens_to_perturb[0]

            # if cell embeddings, can just append
            # shape will be (batch size, 1)
            cos_sims_data = torch.squeeze(cos_sims_data).tolist()

            # handle case of single cell left
            if not isinstance(cos_sims_data, list):
                cos_sims_data = [cos_sims_data]

            cos_sims_dict[(perturbed_genes, "cell_emb")] += cos_sims_data

        else:
            for i, cos in enumerate(cos_sims_data.tolist()):
                cos_sims_dict[(gene_list[i], "cell_emb")].append(cos)
        
        del cos_sims_data
        torch.cuda.empty_cache()
        time.sleep(5)
        
        print('write in update perturbation dictionary')

        return cos_sims_dict


def In_Silico_perturber_node(perturbation,
                        perturb_rank_shift,
                        genes_to_perturb=genes_to_perturb,
                        combos,
                        anchor_gene,
                        model_type,
                        num_classes,
                        emb_mode,
                        cell_emb_style,
                        filter_data_dict,
                        cell_states_to_model,
                        state_embs_dict,
                        n_cells,
                        emb_layer,
                        forward_batch_size,
                        nproc):


    isp = InSilicoPerturber_custom(perturb_type=perturbation,
                        perturb_rank_shift=perturb_rank_shift,
                        genes_to_perturb=genes_to_perturb,
                        combos=combos,
                        anchor_gene=anchor_gene,
                        model_type=model_type,
                        num_classes=nb_classes,
                        emb_mode=emb_mode,
                        cell_emb_style=cell_emb_style,
                        filter_data=filter_data_dict,
                        cell_states_to_model=cell_states_to_model,
                        state_embs_dict=state_embs_dict,
                        max_ncells=n_cells,
                        emb_layer=emb_layer,
                        forward_batch_size=forward_batch_size,
                        nproc=nproc)
    
    return isp 


def perturb_data(isp,
                path_model,
                path_data_tokenized,
                path_output,
                prefix_output):
    
    perturbations = isp.perturb_data(path_model,
                 path_data_tokenized,
                 path_output,
                 prefix_output)
    
    return perturbations


def in_silico_perturber_stats(mode,
                             genes,
                             combos,
                             anchor_gene,
                             cell_states_to_model):
    
    ispstats = InSilicoPerturberStats(mode=mode,
                                  genes_perturbed=genes,
                                  combos=combos,
                                  anchor_gene=anchor_gene,
                                  cell_states_to_model=cell_states_to_model)
    
    return ispstats


def get_statistics(ipsstats,
                  path,
                    path_output,
                   prefix_output
                  ):
    
    statistics = ispstats.get_stats(path,
                   None,
                   path_output,
                   prefix_output)
    
    return statistics