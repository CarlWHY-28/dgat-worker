import pandas as pd
from scipy.stats import wilcoxon
from scipy import sparse
from typing import Optional, Tuple, List, Dict, Any


import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import tempfile
import os
import sys
import psutil
import squidpy as sq
import random
from scipy import sparse as sp



FIGSIZE = (4.8, 4.8)

IMAGE_NA_PATH = "./logo/no_available_icon.png"

def compute_moran_for_pair(pre_adata,
                           rna_adata,
                           tissue_name=None,
                           coord_type="grid",
                           n_perms=10):
    """
    Moran's I for a single pair of AnnData objects (pre_adata for protein, rna_adata for RNA).
    """

    def return_nan_df(marker_list):
        return pd.DataFrame({
            "tissue": [tissue_name] * len(marker_list),
            "marker": marker_list,
            "moran_mrna": [np.nan] * len(marker_list),
            "moran_protein": [np.nan] * len(marker_list),
        })

    try:
        rna = rna_adata.copy()
        rna = rna[rna.obs_names.intersection(pre_adata.obs_names), :].copy()
        prot = pre_adata.copy()

        if "spatial_connectivities" not in rna.obsp:
            if "spatial" not in rna.obsm:
                print(f"[{tissue_name}] ERROR: rna_adata.obsm['spatial'] 不存在 → 返回 NaN")
                return return_nan_df(prot.var_names)

            sq.gr.spatial_neighbors(rna, coord_type=coord_type)

        prot.obsp["spatial_connectivities"] = rna.obsp["spatial_connectivities"].copy()
        if "spatial_distances" in rna.obsp:
            prot.obsp["spatial_distances"] = rna.obsp["spatial_distances"].copy()

        protein_all_names = np.array(prot.var_names)
        gene_names = np.array(rna.var_names)
        common_markers = np.intersect1d(protein_all_names, gene_names)

        print(f"[{tissue_name}] Protein total: {len(protein_all_names)}, Intersection with RNA: {len(common_markers)}")


        try:
            sq.gr.spatial_autocorr(
                prot,
                mode="moran",
                genes=protein_all_names.tolist(),
                n_perms=max(n_perms, 1)
            )
            moran_prot_series = prot.uns["moranI"]["I"]
        except Exception as e:
            print(f"[{tissue_name}] Protein spatial_autocorr error: {e}")
            moran_prot_series = pd.Series(np.nan, index=protein_all_names)

        try:
            if len(common_markers) > 0:
                rna_sub = rna[:, common_markers].copy()
                sq.gr.spatial_autocorr(
                    rna_sub,
                    mode="moran",
                    genes=common_markers.tolist(),
                    n_perms=max(n_perms, 1)
                )
                moran_rna_series = rna_sub.uns["moranI"]["I"]
            else:
                moran_rna_series = pd.Series(dtype=float)
        except Exception as e:
            print(f"[{tissue_name}] RNA spatial_autocorr error: {e}")
            moran_rna_series = pd.Series(dtype=float)


        df = pd.DataFrame({
            "tissue": tissue_name,
            "marker": protein_all_names
        })

        df["moran_protein"] = moran_prot_series.reindex(protein_all_names).values.astype(float)

        df["moran_mrna"] = moran_rna_series.reindex(protein_all_names).values.astype(float)

        return df

    except Exception as e:
        print(f"[{tissue_name}] Spatial Moran's I calculation error:")
        print(e)
        try:
            cols = pre_adata.var_names
        except:
            cols = []
        return return_nan_df(cols)




def compute_moran_single(pre_adata,
                         adata,
                         tissue_name="Sample",
                         coord_type="grid",
                         n_perms=10):
    """

    """

    moran_df = compute_moran_for_pair(
        pre_adata,
        adata,
        tissue_name=tissue_name,
        coord_type=coord_type,
        n_perms=n_perms,
    )

    stat, pval = np.nan, np.nan
    valid_len = 0

    try:
        valid_rows = moran_df.dropna(subset=["moran_protein", "moran_mrna"])
        valid_len = len(valid_rows)

        if valid_len > 0:
            kwargs = {}
            if "method" in wilcoxon.__code__.co_varnames:
                kwargs["method"] = "approx"

            stat, pval = wilcoxon(
                valid_rows["moran_protein"].values,
                valid_rows["moran_mrna"].values,
                alternative="greater",
                **kwargs
            )
    except Exception as e:
        print(f"Wilcoxon test failed: {e}")

    print(f"Paired Wilcoxon (on {valid_len} common markers for {tissue_name}): Stat={stat}, P-val={pval}")

    return moran_df, {"stat": stat, "pvalue": pval}




def compute_bivariate_moran_single(
        adata,
        tissue_name="Sample",
        output_dir=None,
        coord_type="grid"
):


    print(f"Processing Bivariate Moran's I: {tissue_name} ...")

    try:
        adata_tmp = adata.copy()

        if "spatial_connectivities" not in adata_tmp.obsp:
            sq.gr.spatial_neighbors(adata_tmp, coord_type=coord_type)

        W = adata_tmp.obsp["spatial_connectivities"].copy()
        if not sparse.isspmatrix_csr(W):
            W = sparse.csr_matrix(W)

        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1.0
        norm_diag = sparse.diags(1.0 / row_sums)
        W_norm = norm_diag @ W

        if sparse.issparse(adata_tmp.X):
            X = adata_tmp.X.toarray()
        else:
            X = adata_tmp.X

        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0)
        std[std == 0] = 1e-9
        Z = (X - mean) / std

        N = Z.shape[0]
        spatial_lag = W_norm @ Z
        moran_matrix = (Z.T @ spatial_lag) / N

        markers = adata_tmp.var_names
        df_matrix = pd.DataFrame(moran_matrix, index=markers, columns=markers)

        mask = np.triu(np.ones(df_matrix.shape), k=1).astype(bool)
        df_masked = df_matrix.where(mask)

        result_melted = df_masked.stack().reset_index()
        result_melted.columns = ['Marker_A', 'Marker_B', 'Bivariate_Moran_I']

        result_melted.insert(0, 'Tissue', tissue_name)

        result_melted['Abs_Value'] = result_melted['Bivariate_Moran_I'].abs()
        result_melted = result_melted.sort_values(by='Abs_Value', ascending=False).drop(columns=['Abs_Value'])

        return result_melted

    except Exception as e:
        print(f"  - [{tissue_name}] Error: {e}")
        return pd.DataFrame(columns=['Tissue', 'Marker_A', 'Marker_B', 'Bivariate_Moran_I'])



def _varnames(adata) -> List[str]:
    return [str(x) for x in getattr(adata, "var_names", [])]

def _to_1d_vals(adata, gene) -> np.ndarray:
    X = adata[:, str(gene)].X
    if sp is not None and sp.issparse(X):
        return np.ravel(X.A)
    return np.asarray(X).ravel()

def _has_spatial_coords(adata) -> bool:
    try:
        return "spatial" in adata.obsm and getattr(adata.obsm["spatial"], "shape", (0, 0))[1] >= 2
    except Exception:
        return False

def _probe_spatial_meta(adata) -> Tuple[bool, Optional[str], Optional[str], Dict[str, Any]]:

    meta = {}
    try:
        if "spatial" not in adata.uns or not isinstance(adata.uns["spatial"], dict):
            return False, None, None, meta

        spatial_uns = adata.uns["spatial"]
        libs = list(spatial_uns.keys())
        if len(libs) == 0:
            return False, None, None, meta
        lib = libs[0]
        lib_dict = spatial_uns.get(lib, {})
        images_dict = lib_dict.get("images", {})

        candidates = ["hires", "image", "lowres"]
        img_key = None
        if isinstance(images_dict, dict) and len(images_dict) > 0:
            for k in candidates:
                if k in images_dict:
                    img_key = k
                    break
            if img_key is None:
                img_key = list(images_dict.keys())[0]
            return True, lib, img_key, {"libs": libs, "img_keys": list(images_dict.keys())}
        for k in candidates:
            if k in lib_dict:
                return True, lib, k, {"libs": libs, "img_keys": [k]}

        return False, lib, None, {"libs": libs, "img_keys": []}
    except Exception:
        return False, None, None, meta


def _plot_spatial_tissue_scanpy(adata, library_id: Optional[str], img_key: Optional[str]) -> Optional[plt.Figure]:

    if library_id is not None and img_key is not None:
        try:
            fig_obj = sc.pl.spatial(
                adata,
                color=None,
                library_id=library_id,
                img_key=img_key,
                show=False,
                return_fig=True,
                figsize=FIGSIZE,

            )
            return fig_obj if fig_obj is not None else plt.gcf()
        except Exception:
            pass

    try:
        fig_obj = sc.pl.spatial(
            adata,
            color=None,
            show=False,
            return_fig=True,
            figsize=FIGSIZE
        )
        return fig_obj if fig_obj is not None else plt.gcf()
    except Exception:
        return None


def _plot_spatial_expr_scanpy(adata, gene: str, library_id: Optional[str], img_key: Optional[str]) -> Optional[plt.Figure]:
    if library_id is None or img_key is None:
        return None
    try:
        fig_obj = sc.pl.spatial(
            adata,
            color=str(gene),
            library_id=library_id,
            img_key=img_key,
            show=False,
            return_fig=True,
            figsize=FIGSIZE
        )
        fig = fig_obj if fig_obj is not None else plt.gcf()
        return fig
    except Exception:
        return None

def _plot_scatter_expr(adata, gene: Optional[str]) -> plt.Figure:
    fig = plt.figure(figsize=FIGSIZE)
    if not _has_spatial_coords(adata) or gene is None or str(gene) not in _varnames(adata):
        plt.axis("off")
        plt.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16)
        return fig

    coords = adata.obsm["spatial"]
    vals = _to_1d_vals(adata, gene)
    sca = plt.scatter(coords[:, 0], coords[:, 1], s=20, c=vals)
    #plt.gca().invert_yaxis()
    plt.xticks([]); plt.yticks([])
    plt.colorbar(sca, shrink=0.75).set_label(str(gene))
    return fig

def _plot_image_placeholder(img_path: str) -> plt.Figure:
    fig = plt.figure(figsize=FIGSIZE)
    try:
        img = plt.imread(img_path)
        plt.imshow(img)
        plt.axis("off")
    except Exception:
        plt.axis("off")
        plt.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16)
    return fig


def _plot_tissue_only(adata, library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    try:
        if library_id and img_key:
            sc.pl.spatial(
                adata,
                color=None,
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax
            )
        else:
            sc.pl.spatial(adata, color=None, show=False, ax=ax)
        ax.set_xlim(ax.get_xlim())
        ax.set_ylim(ax.get_ylim())
        ax.set_aspect('equal', adjustable='box')
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        return fig
    except Exception:
        plt.close(fig)

    # fallback: NA
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    ax.axis("off")
    ax.text(0.5, 0.5, "No tissue image", ha="center", va="center", fontsize=14, transform=ax.transAxes)
    return fig

def _plot_spatial_expr_mrna(adata, gene: Optional[str], library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    if gene is None or str(gene) not in _varnames(adata):
        ax.axis("off")
        ax.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16, transform=ax.transAxes)
        return fig

    if library_id and img_key:
        try:
            sc.pl.spatial(
                adata,
                color=str(gene),
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax,
                size=1.7,
                cmap="viridis"   # ⭐ 新增：mRNA 用 viridis
            )
            ax.set_xlim(ax.get_xlim())
            ax.set_ylim(ax.get_ylim())
            ax.set_aspect('equal', adjustable='box')
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            return fig
        except Exception:
            plt.close(fig)

    return _plot_scatter_expr(adata, gene)


def _plot_spatial_expr(adata, gene: Optional[str], library_id: Optional[str], img_key: Optional[str]) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=100)
    if gene is None or str(gene) not in _varnames(adata):
        ax.axis("off")
        ax.text(0.5, 0.5, "NA", ha="center", va="center", fontsize=16, transform=ax.transAxes)
        return fig

    if library_id and img_key:
        try:
            sc.pl.spatial(
                adata,
                color=str(gene),
                library_id=library_id,
                img_key=img_key,
                show=False,
                ax=ax,
                size=1.7,
                cmap="plasma"   # ⭐ 新增：protein 用 plasma
            )
            ax.set_xlim(ax.get_xlim())
            ax.set_ylim(ax.get_ylim())
            ax.set_aspect('equal', adjustable='box')
            plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
            return fig
        except Exception:
            plt.close(fig)

    # fallback scatter
    return _plot_scatter_expr(adata, gene)

def _plot_leiden_clustering(
    adata,
    ax,
    n_neighbors = 10,
    resolution = 0.5,
    title = None,
    seed = 0
    ):
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, use_rep = 'X', random_state = seed)
    sc.tl.leiden(adata, resolution = resolution, random_state = seed)
    if 'leiden_colors' in adata.uns.keys():
        adata.uns.pop('leiden_colors')
    sq.pl.spatial_scatter(adata, color = "leiden", title = title, ax = ax)

    ax.get_legend().set_title("Leiden cluster")