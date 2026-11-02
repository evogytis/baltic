import io
import unittest

import baltic as bt

## shortened from a MrBayes 3.2 .con.tre file: node block before the branch length, branch block after it
MRBAYES_NEXUS = """#NEXUS
begin trees;
\ttranslate
\t\t1\tBicolor,
\t\t2\tBoliviana,
\t\t3\tMexicana,
\t\t4\tTenella
\t\t;
   tree con_50_majrule = [&U] (1[&prob=1.00000000e+00,prob_stddev=0.00000000e+00,prob_range={1.00000000e+00,1.00000000e+00},prob(percent)="100",prob+-sd="100+-0"]:1.653468e-03[&length_mean=1.65558861e-03,length_median=1.65346800e-03,length_95%HPD={1.44636200e-03,1.89358800e-03}],2[&prob=1.00000000e+00]:2.1e-02,(3[&prob=1.00000000e+00]:3.0e-02,4[&prob=1.00000000e+00]:1.2e+00)[&prob=9.75000000e-01,prob(percent)="98",prob+-sd="98+-1"]:4.0e-02[&length_mean=4.1e-02]);
end;
"""


class TestMrBayesConsensusTrees(unittest.TestCase):
    def setUp(self):
        self.tree = bt.load_nexus(io.StringIO(MRBAYES_NEXUS), treeType="divergence", absoluteTime=False)

    def test_tips_are_translated(self):
        self.assertEqual(sorted(tip.name for tip in self.tree.get_external()), ["Bicolor", "Boliviana", "Mexicana", "Tenella"])

    def test_node_and_branch_annotations_are_both_read(self):
        tip = self.tree.get_leaf("Bicolor")

        self.assertEqual(tip.traits, {
            "prob": 1.0,
            "prob_stddev": 0.0,
            "prob_range": [1.0, 1.0],
            "prob(percent)": "100",
            "prob+-sd": "100+-0",
            "length_mean": 1.65558861e-03,
            "length_median": 1.65346800e-03,
            "length_95%HPD": [1.44636200e-03, 1.89358800e-03],
        })
        self.assertAlmostEqual(tip.length, 1.653468e-03)

    def test_clade_support_and_positive_exponent_lengths(self):
        clade = self.tree.get_leaf("Mexicana").parent

        self.assertEqual(clade.traits["prob"], 0.975)
        self.assertEqual(clade.traits["prob+-sd"], "98+-1")
        self.assertAlmostEqual(self.tree.get_leaf("Tenella").length, 1.2)


class TestNumbers(unittest.TestCase):
    def test_scientific_notation(self):
        for text, expected in (("1.2E+3", 1200.0), ("1e-05", 1e-05), ("-1.5e-3", -1.5e-3), (".5", 0.5), ("+2", 2.0)):
            with self.subTest(text=text):
                tree = bt.make_tree(f"(A[&rate={text}]:{text},B:1);", treeType="divergence")
                leaf = tree.get_leaf("A")

                self.assertEqual(leaf.traits["rate"], expected)
                self.assertEqual(leaf.length, expected)

    def test_malformed_branch_length_raises(self):
        with self.assertRaises(bt.TreeParseError):
            bt.make_tree("(A:1.2e,B:1);", treeType="divergence")

    def test_non_numeric_unquoted_value_is_a_string(self):
        tree = bt.make_tree("(A[&rate=NaN,host=human]:1,B:1);", treeType="divergence")

        self.assertEqual(tree.get_leaf("A").traits, {"rate": "NaN", "host": "human"})

    def test_unsupported_unquoted_value_raises(self):
        with self.assertRaises(bt.TreeParseError):
            bt.make_tree("(A[&host=abc#1]:1,B:1);", treeType="divergence")


class TestTiedStates(unittest.TestCase):
    def test_tied_states_keep_first_and_warn_once(self):
        data = '((A[&location="Asia+Europe"]:1,B[&location=Asia+Africa]:1):1,C[&location="Asia"]:1);'
        with self.assertLogs("baltic", level="WARNING") as logs:
            tree = bt.make_tree(data, treeType="divergence")

        self.assertEqual(tree.get_leaf("A").traits["location"], "Asia")
        self.assertEqual(tree.get_leaf("B").traits["location"], "Asia")
        self.assertEqual(len(logs.records), 1)
        self.assertIn("2 annotation values held tied states", logs.output[0])

    def test_numeric_plus_values_are_not_split(self):
        tree = bt.make_tree('(A[&prob+-sd="100+-0",code="1+2"]:1,B:1);', treeType="divergence")

        self.assertEqual(tree.get_leaf("A").traits, {"prob+-sd": "100+-0", "code": "1+2"})


if __name__ == "__main__":
    unittest.main()
