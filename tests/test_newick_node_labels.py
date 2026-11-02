import io
import unittest

import baltic as bt


class TestNewickNodeLabels(unittest.TestCase):
    def test_labels_without_branch_lengths(self):
        for label in ("95", "95.5", "95/100", "inner_node", "inner-node"):
            for template in ("((A,B){label},C);", "(((A,B){label}),C);", "(A,B){label};"):
                with self.subTest(label=label, template=template):
                    tree = bt.make_tree(template.format(label=label), treeType="divergence")
                    baseline = bt.make_tree(template.format(label=""), treeType="divergence")
                    node = tree.get_leaf("A").parent

                    self.assertIs(node, tree.get_leaf("B").parent)
                    self.assertEqual(node.traits, {"label": label})
                    self.assertEqual(
                        [tip.name for tip in tree.get_external()],
                        [tip.name for tip in baseline.get_external()],
                    )
                    self.assertEqual(
                        tree.get_parameter_list("length"),
                        baseline.get_parameter_list("length"),
                    )

    def test_labels_with_branch_lengths(self):
        tree = bt.make_tree("((A:1,B:2)95:3,C:4)root:5;", treeType="divergence")

        self.assertEqual(tree.get_leaf("A").parent.traits, {"label": "95"})
        self.assertEqual(tree.get_leaf("A").parent.length, 3.0)
        self.assertEqual(tree.root.traits, {"label": "root"})
        self.assertEqual(tree.root.length, 5.0)

    def test_label_before_annotation_without_branch_lengths(self):
        tree = bt.make_tree('((A,B)inner[&state="X"],C);', treeType="divergence")

        self.assertEqual(tree.get_leaf("A").parent.traits, {"label": "inner", "state": "X"})

    def test_file_loaders_preserve_labels_without_branch_lengths(self):
        newick = "((A,B)95,C)root;"
        trees = (
            bt.load_newick(io.StringIO(newick), treeType="divergence"),
            bt.load_nexus(
                io.StringIO(f"#NEXUS\nBegin trees;\ntree TREE1 = {newick}\nEnd;\n"),
                treeType="divergence",
                absoluteTime=False,
            ),
        )
        for tree in trees:
            with self.subTest(loader=type(tree).__name__):
                self.assertEqual(tree.get_leaf("A").parent.traits, {"label": "95"})
                self.assertEqual(tree.root.traits, {"label": "root"})


if __name__ == "__main__":
    unittest.main()
