from unittest import TestCase, main

import skbio
import numpy.testing as npt

from insertster._core import (decorate, propagate, _set_number_of_tips,
                              _all_seq_scores, score, best, insert,
                              _set_minimum_distance_to_tips,
                              _set_distance_to_root)


class CoreTests(TestCase):
    def setUp(self):
        self.tree = skbio.TreeNode.read(['(((a,b)c,d)e,(f,g)h)root;'])
        self.mock_hits = [('query1', [{'subject': 'a', 'seq_score': 95.6},
                                      {'subject': 'd', 'seq_score': 10.7},
                                      {'subject': 'f', 'seq_score': 15.7}]),
                          ('query2', [{'subject': 'f', 'seq_score': 90.6},
                                      {'subject': 'g', 'seq_score': 12.7}])]

    def test_decorate(self):
        exp = self.tree.copy()

        exp.find('a').hits = {'query1': [95.6]}
        exp.find('b').hits = {}
        exp.find('d').hits = {}
        exp.find('f').hits = {'query1': [15.7], 'query2': [90.6]}
        exp.find('g').hits = {'query2': [12.7]}

        obs = decorate(self.tree, self.mock_hits, threshold=11.0)
        for o, e in zip(obs.tips(), exp.tips()):
            self.assertEqual(o.hits, e.hits)
        for o in obs.non_tips():
            self.assertEqual(o.hits, {})

    def test_propagate(self):
        exp = self.tree.copy()
        exp.find('a').hits = {'query1': [95.6]}
        exp.find('b').hits = {}
        exp.find('c').hits = {'query1': [95.6]}
        exp.find('d').hits = {}
        exp.find('e').hits = {'query1': [95.6]}
        exp.find('f').hits = {'query1': [15.7], 'query2': [90.6]}
        exp.find('g').hits = {'query2': [12.7]}
        exp.find('h').hits = {'query1': [15.7], 'query2': [90.6, 12.7]}
        exp.hits = {'query1': [95.6, 15.7], 'query2': [90.6, 12.7]}

        obs = propagate(decorate(self.tree, self.mock_hits, threshold=11.0))

        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.hits.keys(), e.hits.keys())
            for o_key in o.hits:
                npt.assert_equal(o.hits[o_key], e.hits[o_key])

    def test_set_number_of_tips(self):
        exp = self.tree.copy()
        exp.find('a').ntips = 1
        exp.find('b').ntips = 1
        exp.find('c').ntips = 2
        exp.find('d').ntips = 1
        exp.find('e').ntips = 3
        exp.find('f').ntips = 1
        exp.find('g').ntips = 1
        exp.find('h').ntips = 2
        exp.ntips = 5

        obs = _set_number_of_tips(self.tree)
        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.ntips, e.ntips)

    def test_all_seq_scores(self):
        exp = {'query1': [95.6, 10.7, 15.7], 'query2': [90.6, 12.7]}
        obs = _all_seq_scores(propagate(decorate(self.tree, self.mock_hits)))
        self.assertEqual(obs.keys(), exp.keys())
        for k in obs:
            npt.assert_equal(obs[k], exp[k])

    def test_score(self):
        def score_f(global_seq_scores, local_seq_scores, node):
            return len(local_seq_scores)

        exp = self.tree.copy()
        exp.find('a').scores = {'query1': 1}
        exp.find('b').scores = {}
        exp.find('c').scores = {'query1': 1}
        exp.find('d').scores = {}
        exp.find('e').scores = {'query1': 1}
        exp.find('f').scores = {'query1': 1, 'query2': 1}
        exp.find('g').scores = {'query2': 1}
        exp.find('h').scores = {'query1': 1, 'query2': 2}
        exp.scores = {'query1': 2, 'query2': 2}

        obs = propagate(decorate(self.tree, self.mock_hits, threshold=11.0))
        obs = score(obs, score_f)

        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.scores, e.scores)

    def test_best_allow_root(self):
        def score_f(global_seq_scores, local_seq_scores, node):
            return len(local_seq_scores)

        exp = {'query1': {'node': self.tree.find('root'), 'score': 2},
               'query2': {'node': self.tree.find('h'), 'score': 2}}

        def score_f(global_seq_scores, local_seq_scores, node):
            return len(local_seq_scores)

        obs = best(score(propagate(decorate(self.tree, self.mock_hits,
                                            threshold=11.0)),
                         score_f))
        self.assertEqual(obs, exp)

    def test_best_no_root(self):
        def score_f(global_seq_scores, local_seq_scores, node):
            return len(local_seq_scores)

        def battle_f(cur_node, existing_node, cur_score, existing_score):
            if cur_node.is_root():
                return False
            if cur_score > existing_score:
                return True
            if cur_score == existing_score:
                return cur_node.min_tip_dist < existing_node.min_tip_dist
            return False

        exp = {'query1': {'node': self.tree.find('a'), 'score': 1},
               'query2': {'node': self.tree.find('h'), 'score': 2}}

        obs = best(score(propagate(decorate(self.tree, self.mock_hits,
                                            threshold=11.0)),
                         score_f),
                   battle_f=battle_f)
        self.assertEqual(obs, exp)


    def test_insert(self):
        def battle_f(cur_node, existing_node, cur_score, existing_score):
            if cur_node.is_root():
                return False
            if cur_score > existing_score:
                return True
            if cur_score == existing_score:
                return cur_node.min_tip_dist < existing_node.min_tip_dist
            return False

        def score_f(global_seq_scores, local_seq_scores, node):
            return len(local_seq_scores)

        insert(best(score(propagate(decorate(self.tree, self.mock_hits,
                                             threshold=11.0)),
                          score_f),
                    battle_f=battle_f),
               threshold=2)

        node = self.tree.find('query2')
        self.assertTrue(node.parent is self.tree.find('h'))
        with self.assertRaises(skbio.tree.MissingNodeError):
            self.tree.find('query1')

    def test_set_minimum_distance_to_tips(self):
        exp = self.tree.copy()

        exp.find('a').min_tip_dist = 0
        exp.find('b').min_tip_dist = 0
        exp.find('c').min_tip_dist = 1
        exp.find('d').min_tip_dist = 0
        exp.find('e').min_tip_dist = 1
        exp.find('f').min_tip_dist = 0
        exp.find('g').min_tip_dist = 0
        exp.find('h').min_tip_dist = 1
        exp.find('root').min_tip_dist = 2

        obs = _set_minimum_distance_to_tips(self.tree)

        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.min_tip_dist, e.min_tip_dist)

    def test_set_distance_to_root(self):
        self.tree = skbio.TreeNode.read(['(((a,b)c,d)e,(f,g)h)root;'])
        exp = self.tree.copy()

        exp.find('a').root_dist = 3
        exp.find('b').root_dist = 3
        exp.find('c').root_dist = 2
        exp.find('d').root_dist = 2
        exp.find('e').root_dist = 1
        exp.find('f').root_dist = 2
        exp.find('g').root_dist = 2
        exp.find('h').root_dist = 1
        exp.find('root').root_dist = 0

        obs = _set_distance_to_root(self.tree)

        for o, e in zip(obs.traverse(), exp.traverse()):
            self.assertEqual(o.root_dist, e.root_dist)


if __name__ == '__main__':
    main()
