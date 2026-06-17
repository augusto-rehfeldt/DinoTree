from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
from xml.sax.saxutils import escape
import json
import re
import shutil
import subprocess


ROOT_DIR = Path(__file__).resolve().parent
DATA_FILE = ROOT_DIR / "dino_dict.json"
TREES_DIR = ROOT_DIR / "trees"
FINAL_TREE_FILE = ROOT_DIR / "final_tree.xml"

XML_NAMESPACE = "http://www.phyloxml.org"

_XML_HEADER = (
    '<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
    f'xmlns="{XML_NAMESPACE}" '
    f'xsi:schemaLocation="{XML_NAMESPACE} {XML_NAMESPACE}/1.10/phyloxml.xsd">'
)


@dataclass
class SyncStats:
    written: int = 0
    skipped: int = 0
    removed: int = 0


@dataclass
class MergeConflict:
    name: str
    parents: list[str]


@dataclass
class MergeResult:
    tree: dict[str, Any]
    conflicts: list[MergeConflict]
    source_count: int


@dataclass
class MergeRules:
    rename: dict[str, str]
    collapse: set[str]
    remove: set[str]
    drop_leaves: set[str]


def load_dino_dict(path: Path = DATA_FILE) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError("dino_dict.json must contain a JSON object")

    return data


def load_merge_rules(path: Path | None) -> MergeRules:
    if path is None:
        return MergeRules(rename={}, collapse=set(), remove=set(), drop_leaves=set())

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError("Merge rules must contain a JSON object")

    rename = data.get("rename", {})
    if not isinstance(rename, dict):
        raise ValueError("Merge rule 'rename' must be an object")

    return MergeRules(
        rename={str(key): str(value) for key, value in rename.items()},
        collapse=_load_name_set(data, "collapse"),
        remove=_load_name_set(data, "remove"),
        drop_leaves=_load_name_set(data, "drop_leaves"),
    )


def _load_name_set(data: dict[str, Any], key: str) -> set[str]:
    value = data.get(key, [])
    if not isinstance(value, list):
        raise ValueError(f"Merge rule '{key}' must be a list")
    return {str(item).strip() for item in value if str(item).strip()}


def write_text_if_changed(path: Path, content: str) -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") == content:
        return False
    path.write_text(content, encoding="utf-8", newline="\n")
    return True


def clean_generated_tree_files(trees_dir: Path = TREES_DIR) -> int:
    if not trees_dir.exists():
        return 0

    removed = 0
    for path in trees_dir.iterdir():
        if path.is_file() and path.suffix.lower() in {".xml", ".nwk", ".tre"}:
            path.unlink()
            removed += 1
    return removed


def _process_single_tree(args_tuple):
    dino_name, tree_data, trees_dir = args_tuple
    xml_path = trees_dir / f"{dino_name}.xml"
    nwk_path = trees_dir / f"{dino_name}.nwk"

    xml_content = build_phyloxml(dino_name, tree_data)
    nwk_content = build_newick(dino_name, tree_data)

    xml_written = write_text_if_changed(xml_path, xml_content)
    nwk_written = write_text_if_changed(nwk_path, nwk_content)
    return xml_written, nwk_written


def sync_tree_files(
    dino_dict: dict[str, Any],
    trees_dir: Path = TREES_DIR,
    prune_orphans: bool = True,
) -> SyncStats:
    trees_dir.mkdir(parents=True, exist_ok=True)
    stats = SyncStats()

    if prune_orphans:
        expected = set(dino_dict)
        for path in trees_dir.iterdir():
            if path.is_file() and path.suffix.lower() in {".xml", ".nwk", ".tre"} and path.stem not in expected:
                path.unlink()
                stats.removed += 1

    tasks = [(dino_name, tree_data, trees_dir) for dino_name, tree_data in dino_dict.items()]

    # ponytail: process pool was I/O-bound, a loop is faster here
    for task in tasks:
        xml_written, nwk_written = _process_single_tree(task)
        if xml_written:
            stats.written += 1
        else:
            stats.skipped += 1
                
        if nwk_written:
            stats.written += 1
        else:
            stats.skipped += 1

    return stats


def merge_dino_trees(
    dino_dict: dict[str, Any],
    rules: MergeRules | None = None,
    include_entry_names: bool = False,
) -> MergeResult:
    rules = rules or MergeRules(rename={}, collapse=set(), remove=set(), drop_leaves=set())
    merged: dict[str, Any] = {}
    anchor_paths = _collect_explicit_anchor_paths(dino_dict, rules)
    parent_support = _collect_parent_support(dino_dict, rules)

    for entry_name, tree_data in dino_dict.items():
        normalized = normalize_tree_data(tree_data, rules)
        if include_entry_names:
            normalized = _with_entry_name(entry_name, normalized, rules)
        _merge_tree_children(merged, normalized)

    _canonicalize_duplicate_clades(merged, anchor_paths)
    _deduplicate_clade_names(merged, anchor_paths, parent_support)
    _apply_taxonomic_hindsight(merged)
    conflicts = find_parent_conflicts(merged)
    return MergeResult(tree=merged, conflicts=conflicts, source_count=len(dino_dict))


def normalize_tree_data(data: Any, rules: MergeRules | None = None) -> dict[str, Any]:
    rules = rules or MergeRules(rename={}, collapse=set(), remove=set(), drop_leaves=set())
    return _normalize_tree_children(data, rules, infer_flat=True)


def write_merged_tree_files(
    result: MergeResult,
    root_name: str,
    xml_output: Path,
    newick_output: Path | None = None,
    conflicts_output: Path | None = None,
    dot_output: Path | None = None,
    png_output: Path | None = None,
) -> SyncStats:
    stats = SyncStats()
    if write_text_if_changed(xml_output, build_phyloxml(root_name, result.tree)):
        stats.written += 1
    else:
        stats.skipped += 1

    if newick_output is not None:
        if write_text_if_changed(newick_output, build_newick(root_name, result.tree)):
            stats.written += 1
        else:
            stats.skipped += 1

    if conflicts_output is not None:
        content = json.dumps(serialize_merge_conflicts(result.conflicts), indent=2, ensure_ascii=False)
        if write_text_if_changed(conflicts_output, f"{content}\n"):
            stats.written += 1
        else:
            stats.skipped += 1

    if dot_output is not None:
        if write_text_if_changed(dot_output, build_graphviz_dot(root_name, result.tree)):
            stats.written += 1
        else:
            stats.skipped += 1

    if png_output is not None:
        if dot_output is None:
            dot_output = png_output.with_suffix(".dot")
            write_text_if_changed(dot_output, build_graphviz_dot(root_name, result.tree))
        render_graphviz_png(dot_output, png_output)
        stats.written += 1

    return stats


def build_graphviz_dot(root_name: str, data: Any) -> str:
    lines = [
        "digraph DinoTree {",
        "  graph [rankdir=LR, overlap=false, splines=true];",
        "  node [shape=box, style=\"rounded,filled\", fillcolor=\"#f7f4e8\", color=\"#8a7a52\", fontname=\"Arial\", fontsize=10];",
        "  edge [color=\"#8a7a52\", arrowsize=0.5];",
    ]
    next_id = 0

    def add_node(label: str) -> str:
        nonlocal next_id
        node_id = f"n{next_id}"
        next_id += 1
        lines.append(f"  {node_id} [label={_quote_dot_label(label)}];")
        return node_id

    def add_children(parent_id: str, children: Any) -> None:
        for child_name, child_data in _iter_raw_children(children):
            child_id = add_node(str(child_name))
            lines.append(f"  {parent_id} -> {child_id};")
            add_children(child_id, child_data)

    root_id = add_node(root_name)
    add_children(root_id, data)
    lines.append("}")
    return "\n".join(lines) + "\n"


def render_graphviz_png(dot_path: Path, png_path: Path, dot_command: str = "dot") -> None:
    if shutil.which(dot_command) is None:
        raise RuntimeError(
            "Graphviz 'dot' was not found on PATH. Install Graphviz to render PNG output."
        )

    png_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [dot_command, "-Tpng", str(dot_path), "-o", str(png_path)],
        check=True,
        capture_output=True,
        text=True,
    )


def _quote_dot_label(label: str) -> str:
    return json.dumps(label, ensure_ascii=False)


def serialize_merge_conflicts(conflicts: Iterable[MergeConflict]) -> list[dict[str, Any]]:
    return [
        {
            "name": conflict.name,
            "parents": [parent or "<root>" for parent in conflict.parents],
        }
        for conflict in conflicts
    ]


def find_parent_conflicts(tree: dict[str, Any]) -> list[MergeConflict]:
    parents_by_name: dict[str, set[str]] = {}

    def visit(children: dict[str, Any], parent_path: str) -> None:
        for name, grandchildren in children.items():
            parents_by_name.setdefault(name, set()).add(parent_path)
            if isinstance(grandchildren, dict):
                child_path = f"{parent_path}/{name}" if parent_path else name
                visit(grandchildren, child_path)

    visit(tree, "")
    conflicts = [
        MergeConflict(name=name, parents=sorted(parents))
        for name, parents in parents_by_name.items()
        if len(parents) > 1
    ]
    return sorted(conflicts, key=lambda conflict: conflict.name.casefold())


def _with_entry_name(entry_name: str, children: dict[str, Any], rules: MergeRules) -> dict[str, Any]:
    normalized_name = rules.rename.get(entry_name, entry_name).strip()
    if not normalized_name or normalized_name in rules.remove:
        return children
    if normalized_name in rules.collapse:
        return children
    if normalized_name in rules.drop_leaves and not children:
        return children
    if _tree_contains_name(children, normalized_name):
        return children
    return {**children, normalized_name: {}}


def _tree_contains_name(children: dict[str, Any], target_name: str) -> bool:
    for name, grandchildren in children.items():
        if name == target_name:
            return True
        if isinstance(grandchildren, dict) and _tree_contains_name(grandchildren, target_name):
            return True
    return False


def _normalize_tree_children(data: Any, rules: MergeRules, infer_flat: bool = False) -> dict[str, Any]:
    normalized: dict[str, Any] = {}
    raw_items = _iter_raw_children(data)
    if infer_flat and _looks_like_flat_cladogram(raw_items):
        root_name, root_children = raw_items[0]
        root_name = rules.rename.get(str(root_name), str(root_name)).strip()
        if not root_name or root_name in rules.remove:
            raw_items = raw_items[1:]
        elif root_name in rules.collapse:
            raw_items = raw_items[1:]
        else:
            children = _normalize_tree_children(dict(raw_items[1:]), rules, infer_flat=False)
            if root_name in rules.drop_leaves and not children:
                return normalized
            return {root_name: children}

    for raw_name, raw_children in raw_items:
        name = rules.rename.get(str(raw_name), str(raw_name)).strip()
        if not name or name in rules.remove:
            continue

        children = _normalize_tree_children(raw_children, rules, infer_flat=False)
        if name in rules.drop_leaves and not children:
            continue
        if name in rules.collapse:
            _merge_tree_children(normalized, children)
            continue

        _merge_tree_children(normalized, {name: children})
    return normalized


def _looks_like_flat_cladogram(raw_items: list[tuple[str, Any]]) -> bool:
    if len(raw_items) < 2:
        return False
    return all(not _iter_raw_children(children) for _name, children in raw_items)


def _collect_explicit_anchor_paths(
    dino_dict: dict[str, Any],
    rules: MergeRules,
) -> dict[str, tuple[str, ...]]:
    anchors: dict[str, tuple[str, ...]] = {}

    for tree_data in dino_dict.values():
        normalized = normalize_tree_data(tree_data, rules)
        if _looks_like_flat_cladogram(_iter_raw_children(tree_data)):
            continue
        for path in _iter_tree_paths(normalized):
            name = path[-1]
            current = anchors.get(name)
            if current is None or len(path) > len(current):
                anchors[name] = path

    return anchors


def _collect_parent_support(
    dino_dict: dict[str, Any],
    rules: MergeRules,
) -> dict[tuple[str, tuple[str, ...]], int]:
    support: dict[tuple[str, tuple[str, ...]], int] = {}

    for tree_data in dino_dict.values():
        normalized = normalize_tree_data(tree_data, rules)
        for path in _iter_tree_paths(normalized):
            name = path[-1]
            parent_path = path[:-1]
            key = (name, parent_path)
            support[key] = support.get(key, 0) + 1

    return support


def _iter_tree_paths(tree: dict[str, Any]) -> Iterable[tuple[str, ...]]:
    def visit(children: dict[str, Any], path: tuple[str, ...]) -> Iterable[tuple[str, ...]]:
        for name, grandchildren in children.items():
            child_path = (*path, name)
            yield child_path
            if isinstance(grandchildren, dict):
                yield from visit(grandchildren, child_path)

    yield from visit(tree, ())


def _iter_raw_children(data: Any) -> list[tuple[str, Any]]:
    if isinstance(data, dict):
        return [(str(key), value) for key, value in data.items()]
    if isinstance(data, list):
        return _normalize_sequence(data)
    return []


def _merge_tree_children(target: dict[str, Any], source: dict[str, Any]) -> None:
    for name, children in source.items():
        if name not in target:
            target[name] = {}
        if isinstance(children, dict):
            _merge_tree_children(target[name], children)


def _canonicalize_duplicate_clades(
    tree: dict[str, Any],
    anchor_paths: dict[str, tuple[str, ...]],
) -> None:
    while True:
        occurrences = _collect_clade_occurrences(tree)
        changed = False

        for name, paths in occurrences.items():
            if len(paths) < 2:
                continue

            target_path = anchor_paths.get(name)
            if target_path not in paths:
                root_path = (name,)
                non_root_paths = [path for path in paths if path != root_path]
                if root_path not in paths or not non_root_paths:
                    continue
                target_path = min(non_root_paths, key=lambda path: (len(path), path))

            target = _get_subtree(tree, target_path)
            if target is None:
                continue

            for path in sorted(paths, key=len, reverse=True):
                if path == target_path or _is_prefix(path, target_path):
                    continue
                if name not in anchor_paths and len(path) > 1:
                    continue

                source = _get_subtree(tree, path)
                if source is not None:
                    _merge_tree_children(target, source)
                if _remove_subtree(tree, path):
                    changed = True

        if not changed:
            return


def _deduplicate_clade_names(
    tree: dict[str, Any],
    anchor_paths: dict[str, tuple[str, ...]],
    parent_support: dict[tuple[str, tuple[str, ...]], int],
) -> None:
    while True:
        occurrences = _collect_clade_occurrences(tree)
        duplicates = {name: paths for name, paths in occurrences.items() if len(paths) > 1}
        if not duplicates:
            return

        changed = False
        for name, paths in duplicates.items():
            target_path = _choose_duplicate_target(name, paths, anchor_paths, parent_support)
            target = _get_subtree(tree, target_path)
            if target is None:
                continue

            for path in sorted(paths, key=len, reverse=True):
                if path == target_path:
                    continue

                source = _get_subtree(tree, path)
                if source is not None:
                    _merge_tree_children(target, source)
                if _remove_subtree(tree, path):
                    changed = True

        if not changed:
            return


def _choose_duplicate_target(
    name: str,
    paths: list[tuple[str, ...]],
    anchor_paths: dict[str, tuple[str, ...]],
    parent_support: dict[tuple[str, tuple[str, ...]], int],
) -> tuple[str, ...]:
    anchor_path = anchor_paths.get(name)
    if anchor_path in paths and name not in anchor_path[:-1]:
        return anchor_path

    candidates = [path for path in paths if name not in path[:-1]]
    if not candidates:
        candidates = paths

    def score(path: tuple[str, ...]) -> tuple[int, int, int, tuple[str, ...]]:
        support = parent_support.get((name, path[:-1]), 0)
        is_non_root = 1 if len(path) > 1 else 0
        return (is_non_root, support, len(path), path)

    return max(candidates, key=score)


def _apply_taxonomic_hindsight(tree: dict[str, Any]) -> None:
    while True:
        changed = False
        occurrences = _collect_clade_occurrences(tree)

        for species_name in sorted(occurrences):
            genus_name = _binomial_genus(species_name)
            if not genus_name or genus_name not in occurrences:
                continue

            species_path = occurrences[species_name][0]
            genus_path = occurrences[genus_name][0]
            if _is_prefix(genus_path, species_path):
                continue
            if _move_subtree(tree, species_path, genus_path):
                changed = True
                break

        if changed:
            continue

        occurrences = _collect_clade_occurrences(tree)
        family_by_stem = {
            _family_stem(name): paths[0]
            for name, paths in occurrences.items()
            if _family_stem(name)
        }

        for genus_name in sorted(occurrences):
            if not _looks_like_genus_name(genus_name):
                continue

            family_path = _matching_family_path(genus_name, family_by_stem)
            if family_path is None:
                continue

            genus_path = occurrences[genus_name][0]
            if _is_prefix(family_path, genus_path):
                continue
            if _is_prefix(genus_path, family_path):
                continue
            if _move_subtree(tree, genus_path, family_path):
                changed = True
                break

        if not changed:
            return


def _binomial_genus(name: str) -> str | None:
    match = re.fullmatch(r"([A-Z][A-Za-z-]+) [a-z][a-z-]+", name)
    if match is None:
        return None
    return match.group(1)


def _looks_like_genus_name(name: str) -> bool:
    if name.endswith("idae"):
        return False
    return re.fullmatch(r"[A-Z][A-Za-z-]+", name) is not None


def _family_stem(name: str) -> str | None:
    if not name.endswith("idae") or re.fullmatch(r"[A-Z][A-Za-z-]+", name) is None:
        return None
    stem = name[:-4]
    return stem if len(stem) >= 4 else None


def _matching_family_path(
    genus_name: str,
    family_by_stem: dict[str, tuple[str, ...]],
) -> tuple[str, ...] | None:
    matches = [
        (stem, path)
        for stem, path in family_by_stem.items()
        if genus_name.startswith(stem)
    ]
    if not matches:
        return None
    return max(matches, key=lambda item: len(item[0]))[1]


def _move_subtree(
    tree: dict[str, Any],
    source_path: tuple[str, ...],
    target_parent_path: tuple[str, ...],
) -> bool:
    if not source_path or not target_parent_path:
        return False
    if source_path == target_parent_path:
        return False
    if _is_prefix(source_path, target_parent_path):
        return False

    source = _get_subtree(tree, source_path)
    target_parent = _get_subtree(tree, target_parent_path)
    source_parent = _get_subtree(tree, source_path[:-1])
    if source is None or target_parent is None or source_parent is None:
        return False

    name = source_path[-1]
    if name in target_parent:
        _merge_tree_children(target_parent[name], source)
    else:
        target_parent[name] = source
    source_parent.pop(name, None)
    return True


def _collect_clade_occurrences(tree: dict[str, Any]) -> dict[str, list[tuple[str, ...]]]:
    occurrences: dict[str, list[tuple[str, ...]]] = {}

    def visit(children: dict[str, Any], path: tuple[str, ...]) -> None:
        for name, grandchildren in children.items():
            child_path = (*path, name)
            occurrences.setdefault(name, []).append(child_path)
            if isinstance(grandchildren, dict):
                visit(grandchildren, child_path)

    visit(tree, ())
    return occurrences


def _get_subtree(tree: dict[str, Any], path: tuple[str, ...]) -> dict[str, Any] | None:
    current = tree
    for name in path:
        child = current.get(name)
        if not isinstance(child, dict):
            return None
        current = child
    return current


def _remove_subtree(tree: dict[str, Any], path: tuple[str, ...]) -> bool:
    if not path:
        return False

    parent = _get_subtree(tree, path[:-1])
    if parent is None:
        return False
    return parent.pop(path[-1], None) is not None


def _is_prefix(prefix: tuple[str, ...], path: tuple[str, ...]) -> bool:
    return len(prefix) < len(path) and path[: len(prefix)] == prefix


def build_phyloxml(name: str, data: Any) -> str:
    return (
        f"{_XML_HEADER}<phylogeny rooted=\"false\">"
        f"<name>{escape(name)}</name>"
        f"<clade>{_render_phyloxml_children(data)}</clade>"
        f"</phylogeny></phyloxml>"
    )


def build_newick(name: str, data: Any) -> str:
    return f"{_render_newick_children(data)}{_quote_newick_label(name)};\n"


def _render_phyloxml_children(data: Any) -> str:
    if isinstance(data, dict):
        return "".join(
            _render_phyloxml_clade(child_name, child_data)
            for child_name, child_data in data.items()
        )
    if isinstance(data, list):
        return "".join(
            _render_phyloxml_clade(child_name, child_data)
            for child_name, child_data in _normalize_sequence(data)
        )
    return ""


def _render_phyloxml_clade(name: str, data: Any) -> str:
    return (
        "<clade>"
        f"<name>{escape(str(name))}</name>"
        f"{_render_phyloxml_children(data)}"
        "</clade>"
    )


def _render_newick_children(data: Any) -> str:
    if isinstance(data, dict):
        children = [
            _render_newick_clade(child_name, child_data)
            for child_name, child_data in data.items()
        ]
    elif isinstance(data, list):
        children = [
            _render_newick_clade(child_name, child_data)
            for child_name, child_data in _normalize_sequence(data)
        ]
    else:
        children = []

    if not children:
        return ""
    return f"({','.join(children)})"


def _render_newick_clade(name: str, data: Any) -> str:
    branch = _render_newick_children(data)
    return f"{branch}{_quote_newick_label(str(name))}"


def _normalize_sequence(items: Iterable[Any]) -> list[tuple[str, Any]]:
    normalized: list[tuple[str, Any]] = []
    for item in items:
        if isinstance(item, dict):
            normalized.extend((str(key), value) for key, value in item.items())
        elif isinstance(item, str):
            normalized.append((item, []))
    return normalized


def _quote_newick_label(label: str) -> str:
    if not label:
        return ""
    escaped = label.replace("'", "''")
    needs_quotes = any(ch.isspace() or ch in "():;,[]" for ch in escaped)
    return f"'{escaped}'" if needs_quotes else escaped


def list_clade_names(xml_path: Path) -> list[str]:
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_path)
    root = tree.getroot()
    ns = {"p": XML_NAMESPACE}
    names: list[str] = []
    for clade in root.findall(".//p:clade", ns):
        name = clade.findtext("p:name", default="", namespaces=ns).strip()
        if name:
            names.append(name)
    return names


def prune_clades(xml_path: Path, target_names: Iterable[str], output_path: Path | None = None) -> int:
    import xml.etree.ElementTree as ET

    target_set = {name.strip() for name in target_names if name.strip()}
    if not target_set:
        raise ValueError("At least one clade name is required")

    tree = ET.parse(xml_path)
    root = tree.getroot()
    ns = {"p": XML_NAMESPACE}
    removed = _prune_clade_container(root, target_set, ns)

    if output_path is None:
        output_path = xml_path
    tree.write(output_path, encoding="utf-8", xml_declaration=False)
    return removed


def _prune_clade_container(element, target_names: set[str], ns: dict[str, str]) -> int:
    removed = 0
    i = 0
    while i < len(element):
        child = element[i]
        if _local_name(child.tag) == "clade":
            child_name = child.findtext("p:name", default="", namespaces=ns).strip()
            if child_name in target_names:
                grandchildren = [gc for gc in child if _local_name(gc.tag) == "clade"]
                element.remove(child)
                for offset, grandchild in enumerate(grandchildren):
                    element.insert(i + offset, grandchild)
                removed += 1
                # Re-evaluate the elements shifted into the current index i
                continue
            removed += _prune_clade_container(child, target_names, ns)
        i += 1
    return removed


def _local_name(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", 1)[1]
    return tag
