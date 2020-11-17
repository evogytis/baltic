import unittest
import baltic as bt

class test_parsers(unittest.TestCase):

    def test_beast1_traits(self):
        print('Testing BEAST v1 trait parsing')

        ll=bt.loadNexus('./data/miniFluB.mcc.tree',tip_regex='_([0-9-]+)')

        ll.traverse_tree()
        print('Test if branches have correct number of traits')
        assert len(ll.root.traits)==77
        print('Root has correct number of traits')

        for k in ll.Objects:
            if k.is_node() and k!=ll.root:
                assert len(k.traits) in [77,74,71,41]
            elif k.is_leaf():
                assert len(k.traits) in [76,73,70]

        print('Branches have correct number of traits')

    def test_beast2_traits(self):
        print('Testing BEAST v2 trait parsing')
        ll=bt.loadNexus('./data/MERS.mcc.tree')

        ll.traverse_tree()
        print('Test if branches have correct number of traits')
        assert len(ll.root.traits)==10
        print('Root has correct number of traits')

        for k in ll.Objects:
            if k.is_node() and k!=ll.root:
                assert len(k.traits)==13
            elif k.is_leaf():
                assert len(k.traits) in [9,12]

        print('Branches have correct number of traits')

    def test_newick(self):
        tree = bt.loadNewick('./data/zika.nwk')
        expected_num_nodes = 564
        assert len(tree.Objects) == expected_num_nodes, 'Newick tree does not contain correct number of nodes. Expected: {}. Observed: {}'.format(expected_num_nodes, len(tree.Objects))
        max_height = round(max([i.height for i in tree.Objects]), 4)
        expected_height = 0.0058
        assert max_height == expected_height, 'Newick tree height is not correct. Expected: {}. Observed: {}'.format(expected_height, max_height)

    def test_nexus(self):
        tree = bt.loadNexus('./data/2020-04-13_treetime/divergence_tree.nexus')
        tree.treeStats()
        pass

if __name__ == '__main__':
    unittest.main()
