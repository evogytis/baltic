import unittest

import baltic as bt


class TestLongNames(unittest.TestCase):
    def test_long_tip_names_are_not_truncated(self):
        name = "sample_" + "x" * 300
        for template in ("({name}:1,B:2);", "({name},B);", "(B:1,{name}:2);", "(B:1,'{name}':2);"):
            with self.subTest(template=template):
                tree = bt.make_tree(template.format(name=name), treeType="divergence")
                tips = {tip.name: tip for tip in tree.get_external()}

                self.assertEqual(sorted(tips), sorted([name, "B"]))
                self.assertNotIn("label", tips[name].traits)

    def test_long_annotations_are_read_in_full(self):
        value = "y" * 500
        tree = bt.make_tree(f'(A[&host="{value}",rate=0.5]:1,B:2);', treeType="divergence")

        self.assertEqual(tree.get_leaf("A").traits, {"host": value, "rate": 0.5})


class TestEmptySets(unittest.TestCase):
    def test_empty_sets_are_read_as_empty_lists(self):
        for field in ("segments={}", "segments={ }", "location.set={}", "height_95%_HPD={}"):
            with self.subTest(field=field):
                tree = bt.make_tree(f"(A[&{field},rate=0.5]:1,B:2);", treeType="divergence")
                key = field.split("=")[0]

                self.assertEqual(tree.get_leaf("A").traits, {key: [], "rate": 0.5})

    def test_empty_segments_on_reassortment_network(self):
        tree = bt.make_tree("((A:1,B:1)#H1[&segments={}]:1,(#H1[&segments={0,3}]:1,C:1):1);", treeType="time")
        hybrid = [k for k in tree.Objects if k.traits.get("label") == "#H1"][0]
        retic = [k for k in tree.Objects if isinstance(k, bt.Reticulation)][0]

        self.assertEqual(hybrid.traits["segments"], [])
        self.assertEqual(retic.traits["segments"], [0.0, 3.0])


class TestParseErrors(unittest.TestCase):
    def assertParseError(self, data, message, position):
        with self.assertRaises(bt.TreeParseError) as cm:
            bt.make_tree(data, treeType="divergence")
        text = str(cm.exception)
        self.assertIn(message, text)
        self.assertIn(f"at position {position}", text)
        return text

    def test_error_type_is_backwards_compatible(self):
        with self.assertRaises(ValueError):
            bt.make_tree("(A,B)", treeType="divergence")
        with self.assertRaises(AssertionError):
            bt.make_tree("(A,B)", treeType="divergence")

    def test_missing_semicolon(self):
        self.assertParseError("(A:1,B:2)", "must end in semicolon", 8)

    def test_unbalanced_parentheses(self):
        self.assertParseError("(A:1,B:2));", "unmatched ')'", 9)
        self.assertParseError("((A:1,B:2);", "2 '(' but 1 ')'", 0)

    def test_unsupported_annotation_character(self):
        data = '((A[&host=bat@cave]:1,B:2):1,C:3);'
        text = self.assertParseError(data, "Unsupported character '@' in annotation", data.index("@"))

        # the caret sits under the offending character in the snippet
        snippet, caret = text.splitlines()[1:3]
        self.assertEqual(snippet[caret.index("^")], "@")
        self.assertIn("Check the tree file", text)

    def test_quoted_values_may_hold_any_character(self):
        tree = bt.make_tree('((A[&host="bat@cave; [roost]"]:1,B:2):1,C:3);', treeType="divergence")

        self.assertEqual(tree.get_leaf("A").traits, {"host": "bat@cave; [roost]"})

    def test_unterminated_quote_in_annotation(self):
        self.assertParseError('(A[&host="bat]:1,B:2);', "Unterminated quote in annotation", 9)

    def test_unterminated_annotation(self):
        self.assertParseError("(A[&x=1,B:2);", 'Unterminated annotation (no closing "]")', 2)

    def test_semicolon_before_end(self):
        self.assertParseError("(A:1,B:2);(C:1,D:2);", 'Unexpected ";"', 9)

    def test_unreadable_annotation_field_raises(self):
        data = '(A[&x={"a+b"}]:1,B:2);'
        self.assertParseError(data, "Could not read annotation field 'x={\"a+b\"}'", 2)

    def test_non_numeric_set_raises(self):
        self.assertParseError("(A[&range={a,b}]:1,B:2);", "expected a set of numbers", 2)

    def test_error_snippet_is_short(self):
        data = "(" + ",".join(f"T{n}:1" for n in range(2000)) + ",X[&a=@]:1);"
        with self.assertRaises(bt.TreeParseError) as cm:
            bt.make_tree(data, treeType="divergence")

        self.assertLess(len(str(cm.exception)), 400)


if __name__ == "__main__":
    unittest.main()
