# from ....tools.normalize import log_cpm
from .....tools.decorators import method
from .....tools.utils import check_version


@method(
    method_name="Saucie embedding output",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("saucie"),
    image="openproblems-r-scib",  # only if required
)
def saucie_embed_full_unscaled(adata):
    from scIB.integration import runSaucie
    from scIB.preprocessing import reduce_data

    adata = runSaucie(adata, "batch")
    reduce_data(adata, use_rep="X_emb")
    # Complete the result in-place
    return adata


@method(
    method_name="Saucie embedding output (hvg/unscaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("saucie"),
    image="openproblems-r-scib",  # only if required
)
def saucie_embed_hvg_unscaled(adata):
<<<<<<< HEAD
    from scIB.integration import runSaucie
    from ._hvg import hvg_batch
=======
    from _hvg import hvg_batch
    from scIB.integration import runSaucie
>>>>>>> e8516528a5b0d05e063df8d4a219be965fdc52a0
    from scIB.preprocessing import reduce_data

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = runSaucie(adata, "batch")
    reduce_data(adata, use_rep="X_emb")
    return adata


@method(
    method_name="Saucie embedding output (hvg/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("saucie"),
    image="openproblems-r-scib",  # only if required
)
def saucie_embed_hvg_scaled(adata):
<<<<<<< HEAD
    from scIB.integration import runSaucie
    from ._hvg import hvg_batch
=======
    from _hvg import hvg_batch
    from scIB.integration import runSaucie
>>>>>>> e8516528a5b0d05e063df8d4a219be965fdc52a0
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = hvg_batch(adata, "batch", target_genes=2000, adataOut=True)
    adata = scale_batch(adata, "batch")
    adata = runSaucie(adata, "batch")
    reduce_data(adata, use_rep="X_emb")
    return adata


@method(
    method_name="Saucie embedding output (full/scaled)",
    paper_name="Sc",
    paper_url="temp",
    paper_year=2020,
    code_url="",
    code_version=check_version("saucie"),
    image="openproblems-r-scib",  # only if required
)
def saucie_embed_full_scaled(adata):
    from scIB.integration import runSaucie
    from scIB.preprocessing import reduce_data
    from scIB.preprocessing import scale_batch

    adata = scale_batch(adata, "batch")
    adata = runSaucie(adata, "batch")
    reduce_data(adata, use_rep="X_emb")
    return adata
