import io
import unittest

import baltic as bt


def clade_label(tree):
    """Label of the (A,B) clade."""
    return tree.get_leaf("A").parent.traits.get("label")


class TestQuotedNodeLabels(unittest.TestCase):
    def test_quoted_and_unquoted_support_give_the_same_label(self):
        for data in ("((A:0.1,B:0.2)95.5/100:0.1,C:0.3);", "((A:0.1,B:0.2)'95.5/100':0.1,C:0.3);", '((A:0.1,B:0.2)"95.5/100":0.1,C:0.3);'):
            with self.subTest(data=data):
                tree = bt.make_tree(data, treeType="divergence")

                self.assertEqual(clade_label(tree), "95.5/100")
                self.assertAlmostEqual(tree.get_leaf("A").parent.length, 0.1)

    def test_concordance_factor_label(self):
        tree = bt.make_tree("((A:0.1,B:0.2)'100/45.2/50.1':0.1,C:0.3);", treeType="divergence")

        self.assertEqual(clade_label(tree), "100/45.2/50.1")

    def test_astral_label_with_brackets_and_semicolons(self):
        tree = bt.make_tree("((A,B)'[pp1=0.95;q1=0.6]':0.4,C);", treeType="divergence")

        self.assertEqual(clade_label(tree), "[pp1=0.95;q1=0.6]")
        self.assertEqual(sorted(tip.name for tip in tree.get_external()), ["A", "B", "C"])

    def test_doubled_quote_is_an_escape(self):
        tree = bt.make_tree("((A:0.1,B:0.2)'O''Brien':0.1,C:0.3);", treeType="divergence")

        self.assertEqual(clade_label(tree), "O'Brien")

    def test_parentheses_inside_quoted_label(self):
        tree = bt.make_tree("((A:0.1,B:0.2)'support (SH-aLRT':0.1,C:0.3);", treeType="divergence")

        self.assertEqual(clade_label(tree), "support (SH-aLRT")

    def test_quoted_label_before_annotation(self):
        tree = bt.make_tree("((A:0.1,B:0.2)'95/100'[&rate=0.5]:0.1,C:0.3);", treeType="divergence")

        self.assertEqual(tree.get_leaf("A").parent.traits, {"label": "95/100", "rate": 0.5})

    def test_unterminated_quoted_label_raises(self):
        with self.assertRaises(bt.TreeParseError) as cm:
            bt.make_tree("((A:0.1,B:0.2)'95/100:0.1,C:0.3);", treeType="divergence")

        self.assertIn("Unterminated quoted node label at position 14", str(cm.exception))


class TestConcordanceFactorNexus(unittest.TestCase):
    def test_slash_keys_are_read(self):
        ## shortened from an IQ-TREE 2 .cf.tree.nex file
        handle = io.StringIO(
            "#NEXUS\n"
            "[ See turtle.fa.treefile.cf.stat for branch annotation meanings. ]\n"
            "begin trees;\n"
            '  tree tree_1 = (C:0.33,(A:0.09,B:0.06):0.0045[&gCF="75",gCF/gDF1/gDF2/gDFP="75/0/25/0",gN="4",'
            'label="39",sCF="38.06",sCF/sDF1/sDF2="38.06/28.37/33.57"],D:0.10);\n'
            "end;\n"
        )
        tree = bt.load_nexus(handle, treeType="divergence", absoluteTime=False)

        self.assertEqual(tree.get_leaf("A").parent.traits, {
            "gCF": "75",
            "gCF/gDF1/gDF2/gDFP": "75/0/25/0",
            "gN": "4",
            "label": "39",
            "sCF": "38.06",
            "sCF/sDF1/sDF2": "38.06/28.37/33.57",
        })


if __name__ == "__main__":
    unittest.main()
