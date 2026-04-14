from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
from xml.sax.saxutils import escape
import json
import concurrent.futures


ROOT_DIR = Path(__file__).resolve().parent.parent
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


def load_dino_dict(path: Path = DATA_FILE) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError("dino_dict.json must contain a JSON object")

    return data


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

    # Use ProcessPoolExecutor to speed up CPU-bound tree generation strings and I/O writes
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for xml_written, nwk_written in executor.map(_process_single_tree, tasks):
            if xml_written:
                stats.written += 1
            else:
                stats.skipped += 1
                
            if nwk_written:
                stats.written += 1
            else:
                stats.skipped += 1

    return stats


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
