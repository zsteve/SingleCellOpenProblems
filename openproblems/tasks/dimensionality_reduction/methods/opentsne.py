from ....tools.decorators import method
from ....tools.utils import check_version
from openTSNE import TSNE


@method(
    method_name="“t-Distributed Stochastic Neighbor Embedding (t-SNE)",
    paper_name="Visualizing Data using t-SNE",
    paper_url="https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf",
    paper_year=2008,
    code_url="https://github.com/pavlin-policar/openTSNE",
    code_version=check_version("opentsne"),
    image="openproblems-python-extras",
)
def opentsne(adata):
    sc.pp.pca(adata)
    embedding = TSNE().fit(adata.obsm["X_pca"])
    adata.obsm["X_emb"] = embedding.embedding
    return adata
