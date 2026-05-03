import argparse
import json
import logging
from pathlib import Path
import sys

# Fix relative import when running from root
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from tree_utils import (
    DATA_FILE,
    TREES_DIR,
    FINAL_TREE_FILE,
    load_dino_dict,
    sync_tree_files,
    clean_generated_tree_files,
    list_clade_names,
    load_merge_rules,
    merge_dino_trees,
    prune_clades,
    write_merged_tree_files,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed

def cmd_generate(args):
    data_file = args.data or DATA_FILE
    dino_dict = load_dino_dict(data_file)

    stats = sync_tree_files(dino_dict, args.trees_dir, prune_orphans=not args.keep_orphans)
    logger.info(
        "Generated %s files, skipped %s unchanged files, removed %s orphaned files",
        stats.written,
        stats.skipped,
        stats.removed,
    )
    return 0

def cmd_clean(args):
    removed = clean_generated_tree_files(args.trees_dir)
    print(f"Removed {removed} file(s)")
    return 0

def cmd_list_clades(args):
    names = list_clade_names(args.input)
    if args.unique:
        names = list(dict.fromkeys(names))
    for name in names:
        print(name)
    return 0

def cmd_prune(args):
    output = args.output or args.input.with_name(f"{args.input.stem}_pruned.xml")
    removed = prune_clades(args.input, args.names, output)
    print(f"Processed {removed} clade(s)")
    return 0

def cmd_merge(args):
    data_file = args.data or DATA_FILE
    dino_dict = load_dino_dict(data_file)
    rules = load_merge_rules(args.rules)
    result = merge_dino_trees(dino_dict, rules=rules, include_entry_names=args.include_entry_names)

    newick_output = args.newick_output
    if newick_output is None and not args.no_newick:
        newick_output = args.output.with_suffix(".nwk")

    conflicts_output = args.conflicts_output
    if conflicts_output is None and not args.no_conflicts:
        conflicts_output = args.output.with_name(f"{args.output.stem}_conflicts.json")

    dot_output = args.dot_output
    if dot_output is None and not args.no_dot:
        dot_output = args.output.with_suffix(".dot")

    try:
        stats = write_merged_tree_files(
            result,
            root_name=args.root_name,
            xml_output=args.output,
            newick_output=newick_output,
            conflicts_output=conflicts_output,
            dot_output=dot_output,
            png_output=args.png_output,
        )
    except RuntimeError as exc:
        logger.error("%s", exc)
        return 1

    logger.info(
        "Merged %s source tree(s); wrote %s file(s), skipped %s unchanged file(s), found %s conflict(s)",
        result.source_count,
        stats.written,
        stats.skipped,
        len(result.conflicts),
    )

    if result.conflicts and args.strict_conflicts:
        logger.error("Merge conflicts found. See %s", conflicts_output or "conflict output")
        return 1
    return 0

def cmd_fetch_wiki(args):
    try:
        from wiki_parser import fetch_wiki_trees
    except ImportError as exc:
        missing = exc.name or "optional dependency"
        logger.error(
            "fetch-wiki requires optional dependency '%s'. Install requests and beautifulsoup4.",
            missing,
        )
        return 1

    logger.info("Fetching taxonomy data for %d genera from Wikipedia...", len(args.genera))
    trees = fetch_wiki_trees(args.genera, num_workers=args.workers)
    
    if not trees:
        logger.warning("No genera trees were successfully fetched.")
        return 1
        
    logger.info("Successfully fetched %d trees. Unknown ancestors were handled recursively.", len(trees))
    
    if args.update_dict:
        data_file = args.data or DATA_FILE
        existing_data = {}
        if data_file.exists():
            try:
                with data_file.open("r", encoding="utf-8") as f:
                    existing_data = json.load(f)
            except Exception as e:
                logger.warning("Could not read existing data file. Starting fresh. Error: %s", e)
        
        for name, tree in trees.items():
            existing_data[name] = tree
            
        with data_file.open("w", encoding="utf-8") as f:
            json.dump(existing_data, f, indent=2, ensure_ascii=False)
            f.write("\n")
        logger.info("Updated %s with new genera.", data_file)
    else:
        print(json.dumps(trees, indent=2, ensure_ascii=False))
        
    return 0

def main():
    parser = argparse.ArgumentParser(description="Dinosaur Trees processing toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Generate
    parser_gen = subparsers.add_parser("generate", help="Generate per-dinosaur tree files from dino_dict.json")
    parser_gen.add_argument("--data", type=Path, default=DATA_FILE, help="Path to dino_dict.json")
    parser_gen.add_argument("--trees-dir", type=Path, default=TREES_DIR, help="Output directory for tree files")
    parser_gen.add_argument("--keep-orphans", action="store_true", help="Keep old files that are not present in the JSON")
    parser_gen.set_defaults(func=cmd_generate)

    # Clean
    parser_clean = subparsers.add_parser("clean", help="Remove generated tree files")
    parser_clean.add_argument("--trees-dir", type=Path, default=TREES_DIR, help="Directory to clean")
    parser_clean.set_defaults(func=cmd_clean)

    # List Clades
    parser_list = subparsers.add_parser("list-clades", help="List clade names from a phyloxml file")
    parser_list.add_argument("--input", type=Path, default=FINAL_TREE_FILE, help="Input phyloxml file")
    parser_list.add_argument("--unique", action="store_true", help="Print each clade name only once")
    parser_list.set_defaults(func=cmd_list_clades)

    # Prune
    parser_prune = subparsers.add_parser("prune", help="Remove one or more clades from a phyloxml tree")
    parser_prune.add_argument("--input", type=Path, default=FINAL_TREE_FILE, help="Input phyloxml file")
    parser_prune.add_argument("--output", type=Path, default=None, help="Output phyloxml file")
    parser_prune.add_argument("names", nargs="+", help="Clade names to remove or collapse")
    parser_prune.set_defaults(func=cmd_prune)

    # Merge
    parser_merge = subparsers.add_parser("merge", help="Merge all dino_dict.json entries into one tree")
    parser_merge.add_argument("--data", type=Path, default=DATA_FILE, help="Path to dino_dict.json")
    parser_merge.add_argument("--output", type=Path, default=Path("merged_tree.xml"), help="Output merged phyloxml file")
    parser_merge.add_argument("--newick-output", type=Path, default=None, help="Output merged Newick file")
    parser_merge.add_argument("--no-newick", action="store_true", help="Do not write a merged Newick file")
    parser_merge.add_argument("--dot-output", type=Path, default=None, help="Output merged Graphviz DOT file")
    parser_merge.add_argument("--no-dot", action="store_true", help="Do not write a merged Graphviz DOT file")
    parser_merge.add_argument("--png-output", type=Path, default=None, help="Render the merged tree to a PNG file; uses Graphviz dot or matplotlib")
    parser_merge.add_argument("--conflicts-output", type=Path, default=None, help="Output merge conflict report JSON")
    parser_merge.add_argument("--no-conflicts", action="store_true", help="Do not write a merge conflict report")
    parser_merge.add_argument("--rules", type=Path, default=None, help="Optional JSON rules for rename, collapse, remove, and drop_leaves")
    parser_merge.add_argument("--root-name", default="DinoTree", help="Root name for the merged output tree")
    parser_merge.add_argument("--include-entry-names", action="store_true", help="Add top-level dino_dict keys as leaves when absent from their own tree data")
    parser_merge.add_argument("--strict-conflicts", action="store_true", help="Return a non-zero exit code if duplicate taxa appear under multiple parents")
    parser_merge.set_defaults(func=cmd_merge)

    # Fetch Wiki
    parser_wiki = subparsers.add_parser("fetch-wiki", help="Fetch specific dinosaur genera taxonomy from Wikipedia")
    parser_wiki.add_argument("genera", nargs="+", help="One or more dinosaur genera to fetch (e.g. Alpkarakush Tyrannosaurus)")
    parser_wiki.add_argument("--update-dict", action="store_true", help="Update dino_dict.json with the generated trees")
    parser_wiki.add_argument("--data", type=Path, default=DATA_FILE, help="Path to dino_dict.json")
    parser_wiki.add_argument("--workers", type=positive_int, default=5, help="Number of concurrent Wikipedia requests")
    parser_wiki.set_defaults(func=cmd_fetch_wiki)

    args = parser.parse_args()
    return args.func(args)

if __name__ == "__main__":
    sys.exit(main())
