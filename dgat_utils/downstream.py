import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import os
from scipy import sparse
import squidpy as sq


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