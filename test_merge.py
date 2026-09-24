"""Checks for the evidence merge and the poster. Run: python test_merge.py (or pytest)."""
from pathlib import Path
import tempfile

from merge import merge_dino_trees
from tree_utils import MergeRules


def parents(tree, parent=None, out=None):
    out = {} if out is None else out
    for name, children in tree.items():
        out[name] = parent
        parents(children, name, out)
    return out


def flat(*names):
    return {name: [] for name in names}


def test_merge_rules():
    data = {
        # Flat lists: first name contains the rest; species collapse to their genus.
        "Adasaurus": flat("Velociraptorinae", "Linheraptor", "Adasaurus mongoliensis"),
        "Linheraptor": flat("Dromaeosauridae", "Velociraptorinae", "Linheraptor", "Microraptor"),
        # Outgroup-rooted list (starts with a genus) contributes no containment.
        "Microraptor": flat("Linheraptor", "Microraptor", "Tyrannosauridae"),
        # Rank inversion (a tribe cannot hold a subfamily) discredits the whole list.
        "Tarbosaurus": flat("Tyrannosaurini", "Albertosaurinae", "Tarbosaurus", "Adasaurus"),
        # Stems: Tyrannosaurus -> Tyrannosaurini -> Tyrannosaurinae; junk labels are dropped.
        "Tyrannosaurus": flat("Tyrannosauridae", "Tyrannosaurinae", "Tyrannosaurini", "Tyrannosaurus", "IGM 100/42", "Aguja Fm."),
        # Nested lineage, and a list that tries to move a backbone clade.
        "Herrerasaurus": {"Dinosauria": {"Saurischia": {"Herrerasauridae": {"Herrerasaurus": []}}}},
        "Eoraptor": flat("Theropoda", "Eoraptor", "Saurischia"),
        "Coelophysis": flat("Theropoda", "Coelophysis", "Dromaeosauridae", "Tyrannosauridae"),
    }
    result = merge_dino_trees(data, MergeRules(rename={}, collapse=set(), remove=set(), drop_leaves=set()))
    up = parents(result.tree)

    assert up["Adasaurus"] == "Velociraptorinae"
    assert up["Velociraptorinae"] == "Dromaeosauridae"
    assert up["Microraptor"] == "Dromaeosauridae"
    assert "Adasaurus mongoliensis" not in up
    assert up["Tyrannosaurus"] == "Tyrannosaurini"
    assert up["Tyrannosaurini"] == "Tyrannosaurinae"
    assert up["Tyrannosaurinae"] == "Tyrannosauridae"
    assert up["Herrerasaurus"] == "Herrerasauridae"
    assert up["Theropoda"] == "Saurischia", "backbone must beat the Eoraptor list"
    assert "IGM 100/42" not in up and "Aguja Fm." not in up
    assert result.report["dropped_names"]["IGM 100/42"] == "specimen number"
    rooted = result.report["outgroup_rooted_lists"]
    assert rooted["Microraptor"].startswith("starts with the genus")
    assert rooted["Tarbosaurus"].startswith("rank:")
    assert rooted["Eoraptor"].startswith("backbone:")
    # With no usable claim left, Tarbosaurus is placed by the company it keeps in its list.
    assert up["Tarbosaurus"] == "Tyrannosauridae"
    assert result.report["placed_by_company"]["Tarbosaurus"] == "Tyrannosauridae"

    removed = merge_dino_trees(data, MergeRules(rename={}, collapse={"Velociraptorinae"}, remove={"Tyrannosauridae"}, drop_leaves=set()))
    up = parents(removed.tree)
    assert up["Adasaurus"] == "Dromaeosauridae" and "Velociraptorinae" not in up
    assert "Tyrannosaurus" not in up


def test_tied_rivals_without_common_ancestor():
    # Rival claims with no shared container leave the taxon at the root instead of crashing.
    data = {
        "Chilesaurus": flat("Ornithischia", "Chilesaurus", "Lesothosaurus"),
        "Allosaurus": flat("Tetanurae", "Chilesaurus", "Allosaurus"),
        "Megalosaurus": flat("Tetanurae", "Chilesaurus", "Megalosaurus"),
    }
    result = merge_dino_trees(data)
    assert "Chilesaurus" in result.tree and result.report["unplaced_taxa"] == ["Chilesaurus"]
    assert set(result.tree["Tetanurae"]) == {"Allosaurus", "Megalosaurus"}


def test_poster_renders():
    try:
        from poster import render_poster
    except ImportError:
        print("skipped poster check: matplotlib/numpy not installed")
        return

    tree = {"Dinosauria": {"Theropoda": {"Tyrannosauridae": {"Tyrannosaurus": {}, "Tarbosaurus": {}}, "Coelophysis": {}},
                           "Ornithischia": {"Stegosauria": {"Stegosaurus": {}}, "Ceratopsia": {"Triceratops": {}}}}}
    with tempfile.TemporaryDirectory() as folder:
        for suffix in (".png", ".svg"):
            path = Path(folder) / f"poster{suffix}"
            ages = {"Tyrannosaurus": [72.2, 68.0, 68.0, 66.0], "Stegosaurus": [155.0, 152.0, 150.0, 145.0]}
            render_poster(tree, path, source_count=3, ages=ages)
            assert path.stat().st_size > 1000


if __name__ == "__main__":
    test_merge_rules()
    test_tied_rivals_without_common_ancestor()
    test_poster_renders()
    print("ok")
