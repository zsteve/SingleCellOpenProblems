# Copyright (C) 2020 Krishnaswamy Lab, Yale University

import numpy as np
import pandas as pd
import scipy
import sklearn

import scanpy as sc


def _preprocess(adata):
    """Library size normalize and sqrt transform adata"""
    result_dict = sc.pp.normalize_total(adata, target_sum=1e5, inplace=False)
    adata.layers["X_norm"] = result_dict["X"]
    adata.obs["norm_factor"] = result_dict["norm_factor"]
    adata.layers["X_norm_sqrt"] = adata.layers["X_norm"].sqrt()

    # Do PCA
    adata.obsm["X_pca"] = sc.pp.pca(adata.layers["X_norm_sqrt"])
    return adata


def _create_pdf(data_embedding):
    # Given a data_embedding, sample a simplex to weight each dimension to get a
    # new PDF for each condition

    # Create an array of values that sums to 1
    n_components = data_embedding.shape[1]
    data_simplex = np.sort(np.random.uniform(size=(n_components - 1)))
    data_simplex = np.hstack([0, data_simplex, 1])
    data_simplex = np.diff(data_simplex)
    np.random.shuffle(data_simplex)  # operates inplace

    # Weight each embedding component by the simplex weights
    sort_axis = np.sum(data_embedding * data_simplex, axis=1)

    # Pass the weighted components through a logit
    pdf = scipy.special.expit(sort_axis)
    if np.random.choice([True, False]):
        pdf = 1 - pdf
    return pdf


def simulate_treatment(
    adata: AnnData,
    n_conditions: Optional[int] = 2,
    n_replicates: Optional[int] = 2,
    embedding_name: Optional[str] = "X_pca",
    n_components: Optional[int] = 10,
    seed: Optional[int, np.random.RandomState] = None,
) -> np.ndarray:
    """Creates random differentially expressed regions over a dataset for benchmarking.

    Parameters
    ----------
    seed : integer or numpy.RandomState, optional, default: None
        Random state. Defaults to the global `numpy` random number generator
    adata : AnnData
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_conditions : int, optional, default: 2
        Number of conditions to simulate
    n_replicates : int, optional, default: 2
        Number of replicates to simulate
    embedding_name : str, optional, default: "X_pca"
        Name of embedding in adata.obsm to use to simulate differential abundance
    n_components : int, optional, default: 10
        Number of dimensions of adata.obsm['embedding_name'] to use for simulation.
        For embeddings sorted by explained variance, like PCA, more components results
        in more uniform enrichment / depletion
    seed : [int, np.RandomState], default: None


    Returns
    ----------
    condition : array-like, shape=[n_obs,]
        Condition assiment for each cell
    replicate : array-like, shape=[n_obs,]
        Replicate assiment for each cell
    condition_probability : pandas.DataFrame, shape=[n_obs, n_conditions]
        DataFrame with the corresponding probabiltiy for each condition
    """

    np.random.seed(seed)

    data_embedding = adata.obsm[embedding_name]
    if not np.isclose(data_embedding.mean(), 0):
        # embedding data must be mean-centered
        data_embedding = scipy.stats.zscore(data_embedding, axis=0)

    # Randomly flip sign of each embedding dimension
    data_embedding *= np.random.choice([-1, 1], size=data_embedding.shape[1])

    # Create information about each condition and replicate
    conditions = ["condition{}".format(i) for i in range(1, n_conditions + 1)]
    replicates = ["replicate{}".format(i) for i in range(1, n_replicates + 1)]

    # Create one PDF for each condition
    condition_probability = []
    for condition in conditions:
        pdf = _create_pdf(data_embedding[:, :n_components]).reshape(-1, 1)
        condition_probability.append(pdf)
    condition_probability = np.concatenate(condition_probability, axis=1)

    # Normalize PDF for each condition to sum to 1
    condition_probability = sklearn.preprocessing.normalize(
        condition_probability,
        norm="l1",
    )
    condition_probability = pd.DataFrame(
        condition_probability,
        columns=conditions,
        index=adata.obs_names,
    )

    condition = []
    for ix, prob in condition_probability.iterrows():
        condition.append(np.random.choice(condition_probability.columns, p=prob))

    replicate = np.random.choice(replicates, size=adata.n_obs)

    # Assign attributes to adata
    adata.obs["condition"] = condition
    adata.obs["replicate"] = replicate
    adata.obs["sample"] = [
        ".".join(row) for _, row in adata.obs[["condition", "replicate"]].iterrows()
    ]
    adata.obsm["ground_truth_probability"] = condition_probability
    adata.uns["conditions"] = conditions
    adata.uns["replicates"] = replicates
    adata.uns["samples"] = np.unique(adata.obs["sample"])
    adata.uns["n_conditions"] = n_conditions
    adata.uns["n_replicates"] = n_replicates
    adata.uns["n_samples"] = n_conditions * n_replicates

    return adata
