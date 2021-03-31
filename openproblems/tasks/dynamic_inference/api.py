from ...data.sample import load_sample_data

import scanpy as sc


def check_dataset(adata):
    """Check that dataset output fits expected API."""
    return True


def check_method(adata):
    """Check that method output fits expected API."""
    # bunch of asserts here
    return True


def sample_dataset():
    """Create a simple dataset to use for testing methods in this task."""
    return load_sample_data()


def sample_method(adata):
    """Create sample method output for testing metrics in this task."""
    # TODO
    return adata
