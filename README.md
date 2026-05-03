# DinoTree

> A lightning-fast Python CLI toolkit for generating, pruning, and managing dinosaur phylogeny trees in phyloXML and Newick standard formats.

## Overview

DinoTree is designed to simplify the processing of large taxonomic data sets for dinosaurs. It can ingest a JSON source (`dino_dict.json`) containing nested phylogenetic dictionaries and generate standard tree representation formats suitable for biological computing and visualization.

The tool operates via a centralized CLI controller (`main.py`) and is pure Python without reliance on external biological packages. It utilizes parallel processing to ensure rapid tree generation and syncs files intelligently by only overwriting files that have changed.

## Features

- **Tree Generation**: Converts JSON taxonomy data into individual standard `phyloXML` and `Newick` tree files for each dinosaur entry.
- **Smart Synchronization**: Tracks changes and only writes updated trees. It automatically cleans up "orphaned" files that are no longer present in your active taxonomy dataset.
- **Parallel Processing**: Uses CPU-based concurrent multi-processing to scale efficiently over thousands of records and I/O tasks.
- **Clade Extraction**: Quickly inspects a compiled `phyloXML` tree and lists all parsed clades, with options for checking uniqueness.
- **Pruning**: Programmatically removes specific clades from a tree by name. Rather than breaking the tree, it intelligently patches the hierarchy by shifting sub-clades upwards into the removed clade's position.

## Files

- `dino_dict.json`: Source taxonomy data containing dictionaries of dinosaur names mapped to their taxonomic structures.
- `trees/`: Default output directory for generated `phyloXML` and `Newick` files, one pair per dinosaur.
- `final_tree.xml`: An example `phyloXML` tree used for downstream pruning and inspection.
- `src/`: Source code folder containing the fast XML and Newick building, parsing, and pruning utilities.
- `main.py`: The main entrypoint file and CLI.

## Usage

You can use the toolkit through `main.py` which provides several subcommands.

### Generate Trees

To generate or refresh the tree files from the main source JSON:

```powershell
python main.py generate --data dino_dict.json
```
*Options:*
- `--trees-dir`: Set a custom output folder for generated trees.
- `--keep-orphans`: Prevent the deletion of existing files in the output directory that are not found in the current JSON.

### Clean Generated Files

To quickly delete all generated `.xml`, `.nwk`, and `.tre` files from the trees directory:

```powershell
python main.py clean
```

### Merge Trees

Merge all entries from `dino_dict.json` into one combined tree. By default this writes `merged_tree.xml`, `merged_tree.nwk`, `merged_tree.dot`, and `merged_tree_conflicts.json`.

```powershell
python main.py merge --data dino_dict.json --output merged_tree.xml
```

Render a PNG cladistic tree image:

```powershell
python main.py merge --data dino_dict.json --output merged_tree.xml --png-output merged_tree.png
```

PNG rendering uses Graphviz `dot` when it is installed and available on `PATH`. If Graphviz is not installed, DinoTree falls back to matplotlib when available. The command also writes the `.dot` file, which can be rendered later:

```powershell
dot -Tpng merged_tree.dot -o merged_tree.png
```

Use `--rules` to normalize and prune while merging:

```powershell
python main.py merge --rules merge_rules.example.json --output merged_tree.xml
```

Rule fields:
- `rename`: maps one clade name to another before merge.
- `collapse`: removes a clade but lifts its children into the parent.
- `remove`: removes a clade and all descendants.
- `drop_leaves`: removes matching leaf-only clades.

Conflict handling:
- Flat source cladograms are interpreted as an earliest ancestor followed by contained taxa.
- Explicit nested lineages are used as canonical anchors when the same clade also appears at root.
- Root copies of clades are merged into the earliest non-root placement when no explicit anchor exists.
- Any remaining duplicate clade names are resolved to one path by preferring non-root placements, stronger source support, and deeper paths.
- Taxonomic hindsight moves binomial species under their matching genus and type genera under matching `-idae` families when those nodes exist.
- The conflicts JSON is empty when every clade name has one final placement.
- Add `--strict-conflicts` to make those conflicts fail the command.
- Add `--include-entry-names` to insert top-level `dino_dict.json` keys as leaves when they are absent from their own tree data.

Rendering options:
- `--dot-output`: choose the Graphviz `.dot` output path.
- `--no-dot`: skip writing the Graphviz `.dot` file.
- `--png-output`: render a PNG via Graphviz `dot`, or matplotlib as a fallback.

### List Clades

Extract and print all clade names present in a specific phyloXML file:

```powershell
python main.py list-clades --input final_tree.xml
```
*Options:*
- `--unique`: Print each clade name only once, ignoring duplicates.

### Fetch from Wikipedia

Fetch the full hierarchy of ancestors for multiple specific valid dinosaur genera from Wikipedia:

This command requires the optional scraping dependencies:

```powershell
python -m pip install -r requirements-wiki.txt
```

```powershell
python main.py fetch-wiki Alpkarakush Tyrannosaurus
```
*Options:*
- `--update-dict`: Automatically updates and merges the result into `dino_dict.json`.
- `--workers`: Specify the number of concurrent web requests to use (defaults to 5).

### Prune Clades

Remove specific clade elements from a phyloXML tree while keeping their sub-clades. 

```powershell
python main.py prune --input final_tree.xml Coelophysoidea
```
*Options:*
- `--input`: The input `.xml` tree to prune.
- `--output`: File path to save the pruned tree. Defaults to `<original_name>_pruned.xml`.

## Notes

- Core commands use only the Python standard library.
- `fetch-wiki` requires the optional `requests` and `beautifulsoup4` packages.
- `merge --png-output` requires either Graphviz `dot` or matplotlib.
