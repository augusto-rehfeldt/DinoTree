import argparse
import sys
import logging
from pathlib import Path

# Fix relative import when running from root
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from tree_utils import (
    DATA_FILE,
    TREES_DIR,
    FINAL_TREE_FILE,
    load_dino_dict,
    sync_tree_files,
    clean_generated_tree_files,
    list_clade_names,
    prune_clades,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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

    args = parser.parse_args()
    return args.func(args)

if __name__ == "__main__":
    sys.exit(main())
