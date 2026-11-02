import unittest

import baltic as bt


class TestAnnotationStrings(unittest.TestCase):
    def test_issue_49_preserves_mutations_and_empty_internal_trait(self):
        tree = bt.make_tree(
            '((A:1[&mutations="C1T,-5A,G7R",date=2000.5],'
            'B:1[&mutations="G9A"]):1[&mutations="",date=1999.0]);',
            treeType="time",
        )

        self.assertEqual(
            tree.get_leaf("A").traits,
            {"mutations": "C1T,-5A,G7R", "date": 2000.5},
        )
        self.assertEqual(tree.get_leaf("B").traits, {"mutations": "G9A"})
        self.assertEqual(
            tree.get_leaf("A").parent.traits,
            {"mutations": "", "date": 1999.0},
        )

    def test_both_quote_styles_and_annotation_positions(self):
        for quote in ('"', "'"):
            for before_length in (True, False):
                with self.subTest(quote=quote, before_length=before_length):
                    comment = f"[&label={quote}A,B{quote},empty={quote}{quote}]"
                    tip = f"A{comment}:1" if before_length else f"A:1{comment}"
                    tree = bt.make_tree(f"({tip},B:1);", treeType="divergence")

                    self.assertEqual(
                        tree.get_leaf("A").traits, {"label": "A,B", "empty": ""}
                    )

    def test_quoted_contents_do_not_create_traits(self):
        label = "A,date=2000,flag=true,state.set={X,Y}"
        tree = bt.make_tree(
            f'(A[&label="{label}",date=1999]:1,B:1);',
            treeType="divergence",
        )

        self.assertEqual(
            tree.get_leaf("A").traits, {"label": label, "date": 1999.0}
        )

    def test_quoted_values_keep_spaces_and_numeric_text(self):
        tree = bt.make_tree(
            '(A[&label=" A,B ",code="123",other="O\'Brien"]:1,B:1);',
            treeType="divergence",
        )

        self.assertEqual(
            tree.get_leaf("A").traits,
            {"label": " A,B ", "code": "123", "other": "O'Brien"},
        )

    def test_existing_conversions_sets_and_histories_are_preserved(self):
        tree = bt.make_tree(
            '(A[&label="A,B",state="X+Y",enabled=true,disabled="false",'
            'score=-0.001,range={0.1,0.9},location.set={"X","Y"},'
            'location.set.prob={0.4,0.6},'
            'all_N={{259,9.053121200761384E-4,ACT,ACC}}]:1,B:1);',
            treeType="divergence",
        )

        self.assertEqual(
            tree.get_leaf("A").traits,
            {
                "label": "A,B",
                "state": "X",
                "enabled": True,
                "disabled": False,
                "score": -0.001,
                "range": [0.1, 0.9],
                "location.set": ["X", "Y"],
                "location.set.prob": [0.4, 0.6],
                "all_N": [[259, 9.053121200761384e-4, "ACT", "ACC"]],
            },
        )

    def test_quoted_scalar_traits_survive_serialization(self):
        tree = bt.make_tree(
            '(A[&label="A,B",empty=""]:1,B:1);', treeType="divergence"
        )

        reparsed = bt.make_tree(
            tree.to_string(traits=["label", "empty"]), treeType="divergence"
        )

        self.assertEqual(
            reparsed.get_leaf("A").traits, {"label": "A,B", "empty": ""}
        )


if __name__ == "__main__":
    unittest.main()
