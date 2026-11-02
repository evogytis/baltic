import io
import unittest

import baltic as bt

## from an ASTRAL-scored species tree (Mikkelsen & Weir 2023, Syst. Biol.), shortened; the file ends in "; "
SKUA = ("(antarcticus:0.45,(maccormicki:0.23,((skua:0.39,(pomarinus:0.017,(parasiticus:8.0,longicaudus)"
        "N5'[pp1=1;pp2=0;pp3=0]':1.72)N7'[pp1=1;pp2=0;pp3=0]':0.69)N8'[pp1=0.74;pp2=0;pp3=0.26]':0.05,"
        "chilensis)N9:0.03)N10)N11; \n")


class TestAstralLabels(unittest.TestCase):
    def test_name_and_quoted_label_are_kept_apart(self):
        tree = bt.load_newick(io.StringIO(SKUA), treeType="divergence")
        clade = tree.get_leaf("parasiticus").parent

        self.assertEqual(clade.traits, {"label": "N5", "quoted_label": "[pp1=1;pp2=0;pp3=0]"})
        self.assertAlmostEqual(clade.length, 1.72)
        self.assertEqual(tree.get_leaf("chilensis").parent.traits, {"label": "N9"})
        self.assertEqual(len(tree.get_external()), 7)

    def test_plain_astral_label(self):
        tree = bt.make_tree("((A,B)'[pp1=0.95;q1=0.6;f1=95.0]':0.4,C);", treeType="divergence")

        self.assertEqual(tree.get_leaf("A").parent.traits, {"label": "[pp1=0.95;q1=0.6;f1=95.0]"})

    def test_unterminated_quote_after_name_raises(self):
        with self.assertRaises(bt.TreeParseError):
            bt.make_tree("((A,B)N5'[pp1=1:0.4,C);", treeType="divergence")


class TestTrailingWhitespace(unittest.TestCase):
    def test_whitespace_after_semicolon_is_ignored(self):
        for ending in (";", "; ", ";\r", ";\t \r\n"):
            with self.subTest(ending=repr(ending)):
                tree = bt.make_tree(f"((A:1,B:2):1,C:3){ending}", treeType="divergence")

                self.assertEqual(sorted(tip.name for tip in tree.get_external()), ["A", "B", "C"])

    def test_text_after_semicolon_still_raises(self):
        with self.assertRaises(bt.TreeParseError):
            bt.make_tree("((A:1,B:2):1,C:3); x", treeType="divergence")


class TestRoundTrip(unittest.TestCase):
    def test_astral_labels_survive_to_string(self):
        tree = bt.load_newick(io.StringIO(SKUA), treeType="divergence")
        again = bt.make_tree(tree.to_string(), treeType="divergence")

        self.assertEqual(again.get_leaf("parasiticus").parent.traits, {"label": "N5", "quoted_label": "[pp1=1;pp2=0;pp3=0]"})
        self.assertEqual(again.get_leaf("chilensis").parent.traits, {"label": "N9"})


if __name__ == "__main__":
    unittest.main()
