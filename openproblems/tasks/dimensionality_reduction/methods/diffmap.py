from ....tools.decorators import method
from ....tools.utils import check_version

import scanpy as sc


@method(
    method_name="Diffusion map embeddings",
    paper_name="Diffusion maps for high-dimensional single-cell analysis of "
    "differentiation data",
    paper_url="https://academic.oup.com/bioinformatics/article/31/18/2989/" "241305",
    paper_year=2015,
    code_url="https://scanpy.readthedocs.io/en/stable/api/" "scanpy.tl.diffmap.html",
    code_version=check_version("scanpy"),
)
def diffmap(adata):
    sc.tl.diffmap(adata, n_comps=2)
    adata.obsm["X_emb"] = adata.obsm["X_diffmap"]
    return adata
