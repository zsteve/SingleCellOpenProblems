from .....tools.decorators import metric

import numpy as np


@metric(
    metric_name="Graph connectivity",
    maximize=True,
    image="openproblems-python-batch-integration" # only if required
)
def graph_connectivity(adata):
    import scIB.metrics
    return scIB.metrics.graph_connectivity(adata, "labels")