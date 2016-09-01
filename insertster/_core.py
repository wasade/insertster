from operator import add
from functools import reduce
from collections import defaultdict

import skbio
import numpy as np


def decorate(tree, hits, threshold=None):
    """Decorate hit details on to the tips of the tree

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on
    hits : Iterable of tuple
        The tuples are expected to contain (str, [dict]) where str corresponds
        to the queries name, and [dict] contains the subject hits and
        sequence score for each hit. It is expected that the subject
        identifiers correspond to tips in tree
    threshold : float, optional
        A minimim threshold at which a subject sequence score must be in order
        to map to the tree

    Notes
    -----
    The tree is modified in place.

    Every node will have a .hits attribute set. This attribute will be a dict
    keyed by query ID and valued by its sequence scores.

    This method will only set the .hits attribute of the tips.

    Returns
    -------
    tree
        For convenience as the structure is modified in place

    """
    if threshold is None:
        threshold = 0.0

    for node in tree.traverse():
        node.hits = defaultdict(list)

    for name, details in hits:
        for detail in details:
            node = tree.find(detail['subject'])
            if detail['seq_score'] >= threshold:
                node.hits[name].append(detail['seq_score'])

    return tree


def propagate(tree):
    """Propagate the .hits attribute up the tree

    Notes
    -----
    The tree is modified in place.

    Returns
    -------
    tree
        For convenience as the structure is modified in place
    """
    for node in tree.postorder():
        if node.is_tip():
            continue

        for child in node.children:
            for name, identities in child.hits.items():
                node.hits[name].extend(identities)

    for node in tree.traverse():
        for name, values in node.hits.items():
            node.hits[name] = np.asarray(values)

    return tree


def _set_number_of_tips(tree):
    """Sets .ntips

    Notes
    -----
    Operates in place
    """
    for node in tree.postorder():
        if node.is_tip():
            node.ntips = 1
        else:
            node.ntips = reduce(add, [child.ntips for child in node.children])
    return tree


def _set_minimum_distance_to_tips(tree):
    """Sets .min_tip_dist which is the minimum number of branches to a tip

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on.
    """
    for node in tree.postorder():
        if node.is_tip():
            node.min_tip_dist = 0
        else:
            node.min_tip_dist = min(c.min_tip_dist for c in node.children) + 1
    return tree


def _set_distance_to_root(tree):
    """Sets .root_dist which is the number of branches to the root

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on.
    """
    tree.root_dist = 0
    for node in tree.preorder(include_self=False):
        node.root_dist = node.parent.root_dist + 1
    return tree

def _all_seq_scores(tree):
    """Aggregate all of the sequence scores

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on. It is expected that this tree has .hits set

    Returns
    -------
    dict
        A dict keyed by query name, valued by a list of the observed sequence
        scores
    """
    # root has all that we need
    return {name: np.asarray(values) for name, values in tree.hits.items()}


def score(tree, score_f):
    """Set the .score attribute for all nodes and queries

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on. It is expected that this tree has a .hits set
    score_f : function
        This method to compute a score for a query sequence. This method must
        conform to the following signature:

            f(np.array, np.array, skbio.TreeNode) -> float

        The first argument are the total seq_scores (see _all_seq_scores()) for
        a given query, the second argument are the local scores for the current
        node for a given query, and the last argument is the node itself. The
        return value is expected to be a floating point value.

    Notes
    -----
    The tree is modified in place

    Returns
    -------
    tree
        For convenience as the structure is modified in place
    """
    seq_scores = _all_seq_scores(tree)
    _set_number_of_tips(tree)

    for node in tree.traverse():
        node.scores = {query: score_f(seq_scores[query], local_scores, node)
                       for query, local_scores in node.hits.items()}
    return tree


def best(tree, battle_f=None):
    """Get the best nodes for a query, ties prefer nodes at or near tips

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to operate on. It is expected that this tree has a .scores set
    battle_f : function
        A method to determine what node is best. This method is expected to
        resolve ties. The method signature is:

            f(TreeNode, TreeNode, float, float) -> bool

        Where the first TreeNode is the current node, the second argument and
        second TreeNode is the existing node, the first float is the current
        score, and the second float is the contenting score. A True result
        indicates that the contender (i.e., the current node) "wins."

        This function has, at its disposal, the minimum distance to a tip as
        well as the distance to the root of the tree via attributes set on
        the node.

        The default method will take into set the victor as the one with the
        highest score, and in the event of a tie, select the node that is
        closest to a tip.

    Returns
    -------
    dict
        A dict keyed by the query, and valued by a dict of the node and score
    """
    if battle_f is None:
        def battle_f(cur_node, existing_node, cur_score, existing_score):
            if cur_score > existing_score:
                return True
            if cur_score == existing_score:
                return cur_node.min_tip_dist < existing_node.min_tip_dist
            return False

    _set_minimum_distance_to_tips(tree)
    _set_distance_to_root(tree)

    results = {}

    def store(query, node, score):
        results[query] = {'node': node, 'score': score}

    for node in tree.postorder():
        for query, score in node.scores.items():
            if query not in results:
                store(query, node, score)
            else:
                existing_node = results[query]['node']
                existing_score = results[query]['score']

                if battle_f(node, existing_node, score, existing_score):
                    store(query, node, score)

    return results


def insert(scores, threshold, length_f=None):
    """Insert nodes that pass the threshold

    Parameters
    ----------
    scores : dict
        The result of best()
    threshold : float
        A minimum scoring threshold which must be exceeded for a node to be
        inserted
    length_f : function, optional
        A function which can provide a branch length to set. This function
        must conform to the following signature:

            f(TreeNode, str, float) -> float

        Where TreeNode is the node to insert at, str is the name of the node
        being inserted, and float is the score of the query at the node.

        The default is to set a branch length of 0.0
    """
    if length_f is None:
        def length_f(a, b, c):
            return 0.0

    for query, result in scores.items():
        if result['score'] >= threshold:
            length = length_f(result['node'], query, result['score'])
            result['node'].append(skbio.TreeNode(name=query, length=length))
