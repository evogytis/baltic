import unittest
import imp

class test_parsers(unittest.TestCase):
    
    def test_newick(self):
        bt = imp.load_source('baltic', '../baltic/baltic.py')
        tree = bt.loadNewick('./data/zika.nwk')
        expected_num_nodes = 564
        assert len(tree.Objects) == expected_num_nodes, 'Newick tree does not contain correct number of nodes. Expected: {}. Observed: {}'.format(expected_num_nodes, len(tree.Objects))
        max_height = round(max([i.height for i in tree.Objects]), 4)
        expected_height = 0.0058
        assert max_height == expected_height, 'Newick tree height is not correct. Expected: {}. Observed: {}'.format(expected_height, max_height)

    def test_nexus(self):
        bt = imp.load_source('baltic', '../baltic/baltic.py')
        tree = bt.loadNewick('./data/2020-04-13_treetime/divergence_tree.nexus')
        tree.treeStats()
        pass

if __name__ == '__main__':
    unittest.main()
