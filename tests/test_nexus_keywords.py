import io
import unittest

import baltic as bt


def nexus(*lines):
    return io.StringIO("\n".join(lines) + "\n")


class TestNexusKeywordCase(unittest.TestCase):
    def test_treetime_nexus_loads(self):
        ## shortened from a TreeTime divergence_tree.nexus: taxa block, no translate block, "Tree tree1="
        lines = (
            "#NEXUS",
            "Begin Taxa;",
            " Dimensions NTax=3;",
            " TaxLabels MF073359|Human|2015-03-01|Brazil MF073358|Human|2015-06-01|Brazil KU729217|Human|2015|Brazil",
            "End;",
            "Begin Trees;",
            ' Tree tree1=((MF073359|Human|2015-03-01|Brazil:0.001073[&mutations="A1367G,T3051C",date=2015.16],'
            'MF073358|Human|2015-06-01|Brazil:0.000878[&mutations="C3527T",date=2015.42]):0.0002[&mutations="",date=2014.9],'
            'KU729217|Human|2015|Brazil:0.00147[&mutations="G8424A",date=2015.5])NODE_0000000:0.0[&mutations="",date=2014.5];',
            "End;",
        )
        from_nexus = bt.load_nexus(nexus(*lines), treeType="divergence", absoluteTime=False)
        from_newick = bt.load_newick(nexus(*lines), treeType="divergence")

        self.assertEqual(len(from_nexus.get_external()), 3)
        self.assertEqual(
            sorted((tip.name, tip.length) for tip in from_nexus.get_external()),
            sorted((tip.name, tip.length) for tip in from_newick.get_external()),
        )
        self.assertEqual(
            from_nexus.get_leaf("MF073359|Human|2015-03-01|Brazil").traits,
            {"mutations": "A1367G,T3051C", "date": 2015.16},
        )

    def test_lowercase_translate_block(self):
        handle = nexus(
            "#NEXUS",
            "begin trees;",
            "\ttranslate",
            "\t\t1\tBicolor,",
            "\t\t2\tBoliviana,",
            "\t\t3\tMexicana",
            "\t\t;",
            "   tree con_50_majrule = [&U] (1:0.1,2:0.2,3:0.3);",
            "end;",
        )
        tree = bt.load_nexus(handle, treeType="divergence", absoluteTime=False)

        self.assertEqual(sorted(tip.name for tip in tree.get_external()), ["Bicolor", "Boliviana", "Mexicana"])
        self.assertEqual(tree.tipMap, {"1": "Bicolor", "2": "Boliviana", "3": "Mexicana"})

    def test_uppercase_keywords_and_any_tree_name(self):
        handle = nexus(
            "#NEXUS",
            "BEGIN TREES;",
            "TRANSLATE",
            "    1 A,",
            "    2 B;",
            "TREE best = (1:0.1,2:0.2);",
            "END;",
        )
        tree = bt.load_nexus(handle, treeType="divergence", absoluteTime=False)

        self.assertEqual(sorted(tip.name for tip in tree.get_external()), ["A", "B"])

    def test_no_tree_line_still_fails(self):
        with self.assertRaises(AssertionError):
            bt.load_nexus(nexus("#NEXUS", "begin trees;", "end;"), treeType="divergence", absoluteTime=False)


if __name__ == "__main__":
    unittest.main()
