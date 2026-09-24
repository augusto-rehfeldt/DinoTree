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
- `tree_utils.py`: XML and Newick building, parsing and pruning utilities.
- `main.py`: The main entrypoint file and CLI.
- `merge.py`: Evidence-weighted merge of all source cladograms, with the curated backbone and name cleaning.
- `poster.py` / `silhouettes.py`: Poster renderer and the code-drawn dinosaur silhouettes.
- `ages.py` / `genus_ages.json`: Genus stratigraphic ranges from the Paleobiology Database (cached).

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

Merge all entries from `dino_dict.json` into one tree. By default this writes `merged_tree.xml`, `merged_tree.nwk`, `merged_tree.dot` and `merged_tree_report.json`.

```powershell
python main.py merge --data dino_dict.json --output merged_tree.xml
```

Render the poster: a circular cladogram with clade-coloured branches, italic genus names, a colour band and a silhouette per major group, a legend and a short "how to read" note. The extension picks the format (`.png` at 150 dpi, or vector `.svg` / `.pdf` for print):

```powershell
python -m pip install -r requirements-poster.txt
python main.py merge --poster merged_tree.png
```

The poster also carries a geologic-time ring between the genus names and the colour bands: Triassic, Jurassic and Cretaceous are shaded rings (time runs outward), and each genus has a radial bar from its first to its last appearance (solid: central 80% of well-dated fossils; faint: their full dating windows). Ranges are built from securely identified Paleobiology Database fossil occurrences (well-dated ones, central 80%), with the taxon summary as fallback and a few curated overrides (`AGE_OVERRIDES` in `ages.py`); they are cached in `genus_ages.json`; refresh them with:

```powershell
python main.py fetch-ages
```

`--png-output` is kept as an alias of `--poster`. The silhouettes are drawn from code in `silhouettes.py` (`python silhouettes.py` writes a contact sheet to `.tmp_trees/`).

#### How the merge decides

Most sources are flattened Wikipedia cladograms: a list whose first name is usually the clade that contains the rest. Every source is turned into weighted "clade C contains taxon X" claims (logic in `merge.py`):

| Evidence | Weight |
|---|---|
| Curated backbone of well-supported clades (outgroups down to subfamilies), plus a few pinned genera (`PINS`) | overrides everything |
| Name stem: a taxon sits in the next rank up of its name family (Tyrannosaurus in Tyrannosaurini in Tyrannosaurinae) | 10 |
| The genus's own article cladogram | 3 |
| Nested lineage (e.g. from `fetch-wiki`) | 2 |
| Another genus's cladogram | 1 |
| Every `dino_dict` entry is a dinosaur | 0.1, and a hard limit: entries never leave Dinosauria |

Names are cleaned first: species collapse to their genus, `P. tubicen` expands when one genus in the list starts with P, and specimen numbers, formations and footnotes are dropped. A few typos and competing names are unified (`SYNONYMS` in `merge.py`).

A whole list is ignored when it is rooted on an outgroup: it starts with a genus, the backbone says its first clade cannot contain another clade in it, or ranks are inverted (a tribe cannot contain a subfamily). Claims are then inserted strongest first. A claim is rejected when a genus would contain something, when ranks are inverted, when it would move a backbone clade, or when it contradicts stronger claims (a cycle).

Placement walks down from the root and follows a branch only while its support is at least 3× its best rival (`DOMINANCE`). When rivals are close, the taxon stays at their common ancestor, so disagreement shows as a polytomy and not as a wrong split. Clades are placed first, then genera on the resulting clade tree, so no support is counted twice. `KNOWN_CLADES` lists unranked clade names that would otherwise look like genera. A family-group clade known only from its name stem (e.g. Shamosaurinae) borrows its members' containers (report source: "inherited from its members"), so it cannot pull its own type genus upwards. Genera with no usable claim go to the deepest clade holding at least 60% of their list-mates (the report lists them under `placed_by_company`). Taxa whose rival claims share no ancestor stay at the root and are reported as unplaced.

List order is not used. Measured on this data, "taxa after a clade name belong to it" held only ~39% of the time.

The report JSON records what the tree disagrees with:
- `contradictions`: every rejected or outvoted claim, with its support, reason, sources and where the taxon went instead.
- `outgroup_rooted_lists`, `placed_by_company`, `unplaced_taxa`, `empty_clades` and `dropped_names` (with reasons).

Use `--rules` to normalise and prune:

```powershell
python main.py merge --rules merge_rules.example.json --output merged_tree.xml
```

Rule fields:
- `rename`: maps one name to another before merging.
- `collapse`: removes a clade and lifts its children into its parent.
- `remove`: removes a clade and all descendants.
- `drop_leaves`: removes matching leaves.

Other options: `--report-output` / `--no-report`, `--dot-output` / `--no-dot`, `--newick-output` / `--no-newick`, `--root-name`.

Run the checks with `python test_merge.py` (or `pytest`).

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
python main.py prune --input merged_tree.xml Riojasauridae
```
*Options:*
- `--input`: The input `.xml` tree to prune.
- `--output`: File path to save the pruned tree. Defaults to `<original_name>_pruned.xml`.

## Notes

- Core commands use only the Python standard library.
- `fetch-wiki` requires the optional `requests` and `beautifulsoup4` packages.
- `merge --poster` requires the optional `matplotlib` and `numpy` packages (`requirements-poster.txt`).
