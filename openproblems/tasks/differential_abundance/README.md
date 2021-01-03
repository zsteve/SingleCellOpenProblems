<!--- TODO: update --->

# Template Dataset

Here's a brief task description, maybe link to some seminal papers.

## API

Datasets should contain the following attributes:

* `adata.obs["condition"]`
* `adata.obs["replicate"]`

Methods should assign output to `adata.obs['probability_estimate']`.

Metrics should compare `adata['ground_truth_probability']` to `adata.obs['probability_estimate']`.
