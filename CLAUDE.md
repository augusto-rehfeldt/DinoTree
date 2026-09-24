# DinoTree

## What This Is
CLI toolkit for dinosaur phylogeny trees: generate/prune/merge phyloXML + Newick from `dino_dict.json`.

## Commands
- Generate: `python main.py generate --data dino_dict.json`
- Clean: `python main.py clean`
- Merge: `python main.py merge --data dino_dict.json --output merged_tree.xml` (writes `merged_tree_report.json`)
- Poster: `python main.py merge --poster merged_tree.png` (needs `requirements-poster.txt`; `.svg`/`.pdf` also work)
- Checks: `python test_merge.py`
- Help: `python main.py --help`

## Merge notes
- Merge logic lives in `merge.py`: weighted containment claims, curated `BACKBONE`, `SYNONYMS`, rank/stem rules. Read its docstring and the README section before changing weights.
- Validate merge changes by spot-checking placements of well-known genera and reading `merged_tree_report.json` contradictions.
