from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Dynamic inference"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
