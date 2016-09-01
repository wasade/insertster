from insertster._io import available_parsers
from insertster._measures import available_score_functions
from insertster._best import available_best_functions


def query_parsers():
    return list(available_parsers)


def score_functions():
    return list(available_score_functions)


def best_functions():
    return list(available_best_functions)


__all__ = ['query_parsers', 'score_functions', 'best_functions']
