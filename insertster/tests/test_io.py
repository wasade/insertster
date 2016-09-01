from unittest import TestCase, main
import os

import insertster

class IOTests(TestCase):
    def setUp(self):
        self.sam_fp = os.path.join(os.path.abspath(insertster.__file__),
                                   'data/test.sam')

    def test_from_sam(self):
        exp = [('query1', [{'subject': '10', 'identity': 95.6},
                           {'subject': '123', 'identity': 10.7}]),
               ('query2', [{'subject': '15', 'identity': 90.6},
                           {'subject': '13', 'identity': 12.7}])]

        obs = list(from_sam(self.sam_fp))
        self.assertEqual(obs, exp)
