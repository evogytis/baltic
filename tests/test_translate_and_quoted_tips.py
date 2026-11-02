import io
import os
import tempfile
import unittest

import baltic as bt
from baltic import samogitia


def nexus(translate, tree="((1,2),3)"):
    return io.StringIO(f"#NEXUS\nbegin trees;\n{translate}\ntree t = {tree};\nend;\n")


def tip_names(tree):
    return sorted(tip.name for tip in tree.get_external())


class TestTranslateBlock(unittest.TestCase):
    def load(self, handle):
        return bt.load_nexus(handle, treeType="divergence", absoluteTime=False)

    def test_several_entries_per_line(self):
        tree = self.load(nexus("translate 1 A, 2 B,\n 3 C;"))

        self.assertEqual(tip_names(tree), ["A", "B", "C"])

    def test_quoted_names_hold_any_character(self):
        tree = self.load(nexus("translate\n1 'Mus (house), adult',\n2 'O''Brien',\n3 \"Homo sapiens\";"))

        self.assertEqual(tip_names(tree), ["Homo sapiens", "Mus (house), adult", "O'Brien"])

    def test_unquoted_names_keep_every_character(self):
        tree = self.load(nexus("Translate\n\t1 Côte_d'Ivoire|2020-01-01,\n\t2 A/B+C|x,\n\t3 Homo sapiens\n;"))

        self.assertEqual(tip_names(tree), ["A/B+C|x", "Côte_d'Ivoire|2020-01-01", "Homo sapiens"])

    def test_utf8_file_on_disk(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "tree.nex")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("#NEXUS\nbegin trees;\ntranslate 1 São_Paulo, 2 B;\ntree t = (1,2);\nend;\n")
            tree = self.load(path)

        self.assertEqual(tip_names(tree), ["B", "São_Paulo"])

    def test_unreadable_entries_raise(self):
        for translate in ("translate 1 'A' junk, 2 B, 3 C;", "translate 1 A, 1 B, 3 C;", "translate 1 , 2 B, 3 C;", "translate 1 A, 2 B, 3 C"):
            with self.subTest(translate=translate):
                with self.assertRaises(bt.TreeParseError):
                    self.load(nexus(translate))

    def test_posterior_iterator_uses_the_same_reader(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "posterior.trees")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("#NEXUS\nBegin trees;\n\tTranslate\n\t\t1 'São_Paulo|2020-01-01', 2 B|2019-01-01;\n"
                             "tree STATE_0 = [&R] (1:0.1,2:0.2);\ntree STATE_10 = [&R] (1:0.2,2:0.1);\nEnd;\n")
            records = list(samogitia.posterior_tree_iterator(path, 0, None, r"\|([0-9]+\-[0-9]+\-[0-9]+)", "%Y-%m-%d", r"tree STATE_([0-9]+)"))

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0][3], {"1": "São_Paulo|2020-01-01", "2": "B|2019-01-01"})


class TestSeveralTrees(unittest.TestCase):
    def test_newick_keeps_last_and_warns(self):
        with self.assertLogs("baltic.io", level="WARNING") as logs:
            tree = bt.load_newick(io.StringIO("((A:1,B:1):1,C:1);\n((A:1,C:1):1,B:1);\n"), treeType="divergence")

        self.assertEqual(tree.get_leaf("A").parent, tree.get_leaf("C").parent)
        self.assertIn("2 trees found; only the last one was kept", logs.output[0])

    def test_nexus_keeps_last_and_warns(self):
        handle = io.StringIO("#NEXUS\nbegin trees;\ntree a = ((A:1,B:1):1,C:1);\ntree b = ((A:1,C:1):1,B:1);\nend;\n")
        with self.assertLogs("baltic.io", level="WARNING") as logs:
            tree = bt.load_nexus(handle, treeType="divergence", absoluteTime=False)

        self.assertEqual(tree.get_leaf("A").parent, tree.get_leaf("C").parent)
        self.assertIn("2 trees found; only the last one was kept", logs.output[0])


class TestQuotedTips(unittest.TestCase):
    def test_quoted_tip_names(self):
        tree = bt.make_tree("(('A, b':1,'C:(x)':2):1,'O''Brien':3,\"D e\":4);", treeType="divergence")

        self.assertEqual(tip_names(tree), ["A, b", "C:(x)", "D e", "O'Brien"])
        self.assertEqual(tree.get_leaf("C:(x)").length, 2.0)

    def test_quoted_tip_with_annotation(self):
        tree = bt.make_tree("('A b'[&host=\"bat\"]:1,B:1);", treeType="divergence")

        self.assertEqual(tree.get_leaf("A b").traits, {"host": "bat"})

    def test_unterminated_quoted_tip_raises(self):
        with self.assertRaises(bt.TreeParseError) as cm:
            bt.make_tree("('A b:1,B:1);", treeType="divergence")

        self.assertIn("Unterminated quoted tip name at position 1", str(cm.exception))

    def test_quotes_in_names_survive_to_string(self):
        tree = bt.make_tree("('O''Brien':1,'A, b':1);", treeType="divergence")
        tree.traverse_tree()
        again = bt.make_tree(tree.to_string(), treeType="divergence")

        self.assertEqual(tip_names(again), ["A, b", "O'Brien"])


if __name__ == "__main__":
    unittest.main()
