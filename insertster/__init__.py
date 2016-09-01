from ._io import available_parsers
from ._measures import available_score_functions
from ._best import available_best_functions
from ._core import decorate, propagate, insert, best, score


def query_parsers(method=None):
    if method is None:
        return list(available_parsers)
    else:
        return available_parsers[method]

def score_functions(method=None):
    if method is None:
        return list(available_score_functions)
    else:
        return available_score_functions[method]

def best_functions(method=None):
    if method is None:
        return list(available_best_functions)
    else:
        return available_best_functions[method]


__all__ = ['query_parsers', 'score_functions', 'best_functions'
           'decorate', 'propagate', 'insert', 'best', 'score']
