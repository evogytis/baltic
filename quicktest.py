import baltic as bt

print("Baltic package information")
print(bt)
print(type(bt))
print(dir(bt))

print("\nbaltic.io package information")
print(bt.io)
print(type(bt.io))
print(dir(bt.io))

print("\n\n")

print('Testing BEAST v1 trait parsing')

ll=bt.io.loadNexus('./tests/data/miniFluB.mcc.tree',tip_regex='_([0-9-]+)')

ll.traverse_tree()
print('Test if branches have correct number of traits')
assert len(ll.root.traits)==77
print('Root has correct number of traits')

for k in ll.Objects:
    if k.is_node() and k!=ll.root:
        assert len(k.traits) in [77,74,71,41]
    elif k.is_leaf():
        assert len(k.traits) in [76,73,70]

print('Branches have correct number of traits')