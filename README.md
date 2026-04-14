# DinoTree

> A lightning-fast Python CLI toolkit for generating, pruning, and managing dinosaur phylogeny trees in phyloXML and Newick standard formats.

## What changed

- Removed the hard dependency on Biopython, matplotlib, requests, and BeautifulSoup for the local tree-generation flow.
- Replaced repeated full rewrites with a sync step that only writes changed files.
- Consolidated the helper logic so the API uses one `main.py` controller.
- Parallelized IO processing and updated translations from Spanish legacy code.
- Added a centralized `main.py` CLI interface.

## Files

- `dino_dict.json`: source taxonomy data.
- `trees/`: generated phyloXML and Newick files, one pair per dinosaur.
- `final_tree.xml`: example phyloXML tree for downstream pruning and inspection.
- `src/`: source code folder containing tree utilities.

## Usage

Generate or refresh the tree files:

```powershell
python main.py generate --data dino_dict.json
```

Remove generated tree files:

```powershell
python main.py clean
```

List clades from `final_tree.xml`:

```powershell
python main.py list-clades --input final_tree.xml
```

Prune clades from `final_tree.xml`:

```powershell
python main.py prune --input final_tree.xml Coelophysoidea
```

## Notes

- The scripts only require the Python standard library for the refactored flow.
- If you want to keep stale files in `trees/`, pass `--keep-orphans` to `main.py generate`.
