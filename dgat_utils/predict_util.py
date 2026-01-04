import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)
import requests
import numpy as np

import scanpy as sc
#import torch
import warnings
#import pandas as pd
import torch
#from torch_geometric.loader import DataLoader
#from Model.dgat import GATEncoder, Decoder_Protein, Decoder_mRNA
from utils.Preprocessing import qc_control_cytassist, normalize, clean_protein_names, fill_genes, preprocess_ST
from Model.Train_and_Predict import train_and_evaluate_fold, protein_predict,get_activity

import random
import os
#from muon import prot as pt

#from utils.idk_utils import leiden_plot, find_edges,leiden_plot_scatter, leiden_plot_eva, plot_spatial_expression,merge_celltypes, infer_celltype_activity, plot_heatmap

warnings.filterwarnings('ignore')
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

dataset_save_dir = './DGAT_training_datasets'
pyg_save_dir = './pyg_data' # Building a graph might take 5 mins or more (depends on the scale of the dataset).This directory is used to save the graph data in PyG format, which can be reused for prediction without rebuilding the graph.
model_save_dir = './DGAT_pretrained_models'# Directory to save the trained DGAT model. If you have a pre-trained model , please set this to the directory where the model is saved.
pred_result_path = './DGAT_results' # Not necessary to set this directory, but if you want to save the prediction results, please set this to the directory where you want to save the results and uncomment the last line of each section.

seed = 2025
torch.manual_seed(seed)
np.random.seed(seed)
random.seed(seed)


def web_predict(url,adata):
    common_gene = requests.get(f"{url}/common_gene_11535.txt").text.strip().splitlines()
    common_protein = requests.get(f"{url}/common_protein_31.txt").text.strip().splitlines()
    adata.var_names_make_unique()
    # Performs quality control and normalization on a Spatial Transcriptomics dataset.
    adata,missing_genes = fill_genes(adata, common_gene)
    preprocess_ST(adata)
    # Ensures that the test dataset contains the same genes as the training dataset. Fills missing genes in the test data with zeros.
    pred_adata = protein_predict(adata, common_gene, common_protein, url, pyg_save_dir)

    return pred_adata,missing_genes
