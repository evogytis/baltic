import unittest

import baltic as bt


class TestReduceTree(unittest.TestCase):

    @staticmethod
    def make_target_tree():
        tree = bt.make_tree(
            "(((A:1.0,B:1.0):1.0,C:1.0):1.0,D:1.0);",
            treeType="divergence",
        )
        tree.traverse_tree()
        return tree

    @staticmethod
    def make_foreign_tree():
        tree = bt.make_tree(
            "((W:1.0,X:1.0):1.0,(Y:1.0,Z:1.0):1.0);",
            treeType="divergence",
        )
        tree.traverse_tree()
        return tree

    def test_reduces_tree_using_owned_leaf_objects(self):
        tree = self.make_target_tree()

        reduced = tree.reduce_tree([
            tree.get_leaf("A"),
            tree.get_leaf("D"),
        ])

        self.assertEqual(
            sorted(tip.name for tip in reduced.get_external()),
            ["A", "D"],
        )
        self.assertEqual(
            sorted(tip.name for tip in tree.get_external()),
            ["A", "B", "C", "D"],
        )
        self.assertFalse(
            any(node.is_node() and not node.children for node in reduced.Objects)
        )

    def test_rejects_foreign_leaf_objects_regardless_of_index_collision(self):
        tree = self.make_target_tree()
        foreign = self.make_foreign_tree()

        cases = {
            "internal-node collision": foreign.get_leaf("W"),
            "leaf collision": foreign.get_leaf("Y"),
            "no collision": foreign.get_leaf("X"),
        }

        for label, tip in cases.items():
            with self.subTest(label=label):
                with self.assertRaisesRegex(ValueError, "do not belong"):
                    tree.reduce_tree([tip])

    def test_rejects_empty_and_non_leaf_inputs(self):
        tree = self.make_target_tree()

        with self.assertRaisesRegex(ValueError, "No tips"):
            tree.reduce_tree([])

        with self.assertRaisesRegex(ValueError, "not leaf-like"):
            tree.reduce_tree([tree.root])


if __name__ == "__main__":
    unittest.main()
