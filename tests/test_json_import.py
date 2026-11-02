import unittest

import baltic as bt


class TestJsonImport(unittest.TestCase):
    def test_load_json_accepts_unnamed_internal_nodes(self):
        auspice = {
            "meta": {"colorings": []},
            "tree": {
                "node_attrs": {"div": 0.0},
                "children": [
                    {"name": "A", "node_attrs": {"div": 1.0}},
                    {
                        "node_attrs": {"div": 0.5},
                        "children": [
                            {"name": "B", "node_attrs": {"div": 1.2}},
                            {"name": "C", "node_attrs": {"div": 1.3}},
                        ],
                    },
                ],
            },
        }

        tree, meta = bt.load_JSON(auspice, "divergence", stats=False)

        self.assertEqual(meta["colorings"], [])
        self.assertEqual(tree.root.index, "NODE_0000001")
        self.assertEqual([node.index for node in tree.get_internal()], ["NODE_0000001", "NODE_0000003"])
        self.assertEqual(sorted(leaf.name for leaf in tree.get_external()), ["A", "B", "C"])
        self.assertEqual(tree.get_leaf("A").length, 1.0)

    def test_load_json_reconstructs_divergence_from_branch_lengths(self):
        auspice = {
            "meta": {"colorings": []},
            "tree": {
                "node_attrs": {"length": {"value": 0.0}},
                "children": [
                    {"name": "A", "node_attrs": {"length": {"value": 1.0}}},
                    {
                        "node_attrs": {"length": {"value": 0.5}},
                        "children": [
                            {"name": "B", "node_attrs": {"length": {"value": 0.7}}},
                            {"name": "C", "node_attrs": {"length": {"value": 0.8}}},
                        ],
                    },
                ],
            },
        }

        tree, _ = bt.load_JSON(auspice, "divergence", stats=False)

        self.assertEqual(tree.root.height, 0.0)
        self.assertEqual(tree.get_leaf("A").height, 1.0)
        self.assertEqual(tree.get_leaf("B").height, 1.2)
        self.assertEqual(tree.get_leaf("C").length, 0.8)


if __name__ == "__main__":
    unittest.main()
