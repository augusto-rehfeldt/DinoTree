"""Render a merged tree as a circular, poster-style cladogram with clade colours and silhouettes."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any
import math
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
from matplotlib.patches import Polygon, Rectangle, Wedge
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath
import numpy as np

from silhouettes import silhouette

PAPER = "#f4eedf"
INK = "#2b2823"
MUTED = "#8b8374"
NEUTRAL = "#a39a8a"
MIN_BAND_CAPTION = 6.5  # points; smaller band names move beneath the silhouette
LEGEND_ROW = 1.2  # inches per legend entry
MIN_ART_GENERA = 4  # smaller groups get a radial name instead of a silhouette
TIME_RING = 2.4  # inches of radius for the geologic-time ring
TIME_START, TIME_END = 252.0, 66.0
# name, start, end (Ma), ICS colour
PERIODS = [("Triassic", 252.0, 201.4, "#812b92"), ("Jurassic", 201.4, 145.0, "#34b2c9"), ("Cretaceous", 145.0, 66.0, "#7fc64e")]
EPOCH_BOUNDARIES = [252.0, 247.2, 237.0, 201.4, 174.7, 161.5, 145.0, 100.5, 66.0]
SERIF = ["Georgia", "Palatino Linotype", "Book Antiqua", "DejaVu Serif"]

# clade: (poster label, colour, silhouette, lineage, description).
# The deepest listed clade colours a genus; "Early ..." groups hold genera outside the subgroups.
GROUPS: dict[str, tuple[str, str, str, str, str]] = {
    "Dinosauria": ("Early dinosaurs", "#6f685d", "eoraptor", "Dinosauria", "Dinosaurs of uncertain position near the root"),
    "Saurischia": ("Early saurischians", "#aaa08b", "herrerasaur", "Dinosauria", "Herrerasaurids and other basal saurischians"),
    "Theropoda": ("Early theropods", "#bf9678", "coelophysoid", "Theropoda", "Basal theropods and early neotheropods"),
    "Coelophysoidea": ("Coelophysoids", "#c2703f", "coelophysoid", "Theropoda", "Slender Late Triassic–Early Jurassic predators"),
    "Ceratosauria": ("Ceratosaurs", "#9b3d2c", "abelisaur", "Theropoda", "Abelisaurids, noasaurids and kin; reduced forelimbs"),
    "Tetanurae": ("Early tetanurans", "#b79a7e", "megalosaur", "Theropoda", "Tetanurans outside the major named clades"),
    "Megalosauroidea": ("Megalosauroids", "#cf8a3a", "spinosaur", "Theropoda", "Megalosaurids and piscivorous spinosaurids"),
    "Allosauroidea": ("Allosauroids", "#a8662d", "allosaur", "Theropoda", "Large Jurassic–Cretaceous apex predators"),
    "Coelurosauria": ("Other coelurosaurs", "#dcbd86", "compsognathid", "Theropoda", "Coelurosaurs outside the named subclades"),
    "Tyrannosauroidea": ("Tyrannosauroids", "#7e2a26", "tyrannosaur", "Theropoda", "From small proceratosaurids to Tyrannosaurus"),
    "Ornithomimosauria": ("Ornithomimosaurs", "#d8a54f", "ornithomimid", "Theropoda", "Mostly toothless, cursorial “ostrich dinosaurs”"),
    "Alvarezsauria": ("Alvarezsaurs", "#b5733f", "alvarezsaur", "Theropoda", "Small insectivores with one enlarged hand claw"),
    "Therizinosauria": ("Therizinosaurs", "#8f4f33", "therizinosaur", "Theropoda", "Herbivorous maniraptorans with scythe-like claws"),
    "Oviraptorosauria": ("Oviraptorosaurs", "#c86a4c", "oviraptor", "Theropoda", "Beaked, often crested; some brooded their nests"),
    "Paraves": ("Other paravians", "#8f6a5e", "paravian", "Theropoda", "Paravians outside the named subclades"),
    "Troodontidae": ("Troodontids", "#a4503c", "troodontid", "Theropoda", "Lightly built, large-brained paravians"),
    "Dromaeosauridae": ("Dromaeosaurs", "#6e2f28", "raptor", "Theropoda", "Feathered predators with a hyperextensible claw"),
    "Sauropodomorpha": ("Early sauropodomorphs", "#9aa162", "prosauropod", "Sauropodomorpha", "Early, partly bipedal sauropodomorphs"),
    "Sauropoda": ("Early sauropods", "#6d8a52", "sauropod", "Sauropodomorpha", "Early obligate quadrupeds with columnar limbs"),
    "Diplodocoidea": ("Diplodocoids", "#3f7d5c", "diplodocid", "Sauropodomorpha", "Long-tailed sauropods with peg-like teeth"),
    "Macronaria": ("Macronarians", "#79a383", "brachiosaur", "Sauropodomorpha", "Large-nostrilled sauropods: camarasaurs, brachiosaurs"),
    "Titanosauria": ("Titanosaurs", "#2d5e55", "titanosaur", "Sauropodomorpha", "Cretaceous sauropods, including the largest land animals"),
    "Ornithischia": ("Early ornithischians", "#8d95b0", "small_ornithischian", "Ornithischia", "Basal beaked herbivores"),
    "Thyreophora": ("Early thyreophorans", "#7f93ad", "thyreophoran", "Ornithischia", "Early armoured forms such as Scelidosaurus"),
    "Stegosauria": ("Stegosaurs", "#4d6f96", "stegosaur", "Ornithischia", "Plated thyreophorans with tail spikes"),
    "Ankylosauria": ("Ankylosaurs", "#2f4f74", "ankylosaur", "Ornithischia", "Heavily armoured; some with a tail club"),
    "Neornithischia": ("Early neornithischians", "#9b9cc3", "small_ornithischian", "Ornithischia", "Small bipedal neornithischians"),
    "Ornithopoda": ("Ornithopods", "#6f7ebb", "iguanodont", "Ornithischia", "Bipedal and facultatively quadrupedal browsers"),
    "Hadrosauridae": ("Hadrosaurs", "#454b8c", "hadrosaur", "Ornithischia", "Duck-billed herbivores with dental batteries"),
    "Pachycephalosauria": ("Pachycephalosaurs", "#8c6aa6", "pachycephalosaur", "Ornithischia", "Marginocephalians with thickened skull domes"),
    "Ceratopsia": ("Ceratopsians", "#5f437f", "ceratopsian", "Ornithischia", "Rostral-beaked; later forms frilled and horned"),
}
# Legend columns: lineage -> (heading, one-line characterisation).
LINEAGES = {
    "Dinosauria": ("Early dinosaurs", "Taxa near the base of Dinosauria and Saurischia"),
    "Theropoda": ("Theropoda", "Bipedal, predominantly carnivorous saurischians; birds are living theropods"),
    "Sauropodomorpha": ("Sauropodomorpha", "Long-necked herbivorous saurischians, ultimately quadrupedal"),
    "Ornithischia": ("Ornithischia", "Herbivores united by the predentary bone at the tip of the lower jaw"),
}
MAJOR_CLADES = {"Saurischia", "Ornithischia", "Theropoda", "Sauropodomorpha", "Dinosauria"}


@dataclass
class _Node:
    name: str
    depth: int
    children: list["_Node"] = field(default_factory=list)
    theta: float = 0.0
    leaves: int = 0
    group: str | None = None  # nearest enclosing clade listed in GROUPS


def render_poster(
    tree: dict[str, Any], path: Path, source_count: int, root: str = "Dinosauria", ages: dict[str, list[float]] | None = None
) -> None:
    """Draw the poster. `ages` maps genus -> [first_max, first_min, last_max, last_min] in Ma."""
    ages = ages or {}
    subtree = _find(tree, root)
    if subtree is None:
        root, subtree = "Life", tree
    top = _build(root, subtree, 0, None)
    _group_order(top)
    leaves = _leaves(top)
    if len(leaves) < 2:
        raise ValueError("Nothing to draw: the tree has fewer than two genera")

    # Angular layout: one slot per genus, a small gap wherever the colour group changes, and an
    # opening at twelve o'clock that carries the time scale.
    gap = max(12, len(leaves) // 70) if ages else 0
    slots, previous = [None] * gap, None
    for leaf in leaves:
        if previous is not None and leaf.group != previous:
            slots.append(None)
            slots.append(None)
        slots.append(leaf)
        previous = leaf.group
    step = 2 * math.pi / len(slots)
    for i, leaf in enumerate(slots):
        if leaf is not None:
            leaf.theta = math.pi / 2 + gap * step / 2 - (i + 0.5) * step
    _assign_internal_angles(top)

    label_pt = 4.6
    r_tip = max(8.0, len(slots) * label_pt * 1.22 / 72 / (2 * math.pi))
    r_hub = 0.9  # keeps the first splits from piling up on the centre
    max_depth = max(_depths(top))
    label_room = max(len(leaf.name) for leaf in leaves) * label_pt * 0.5 / 72 + 0.12
    time_in = r_tip + label_room + 0.15
    time_out = time_in + (TIME_RING if ages else 0.0)
    band_in = time_out + 0.15
    band_out = band_in + 0.45
    art_size = 1.7

    # Colour bands: one per contiguous run of a group. Its art and name go on the group's longest run.
    runs: list[tuple[str, float, float, int]] = []
    for leaf in leaves:
        if leaf.group is None:
            continue
        if runs and runs[-1][0] == leaf.group and abs(runs[-1][2] - leaf.theta - step) < step * 0.5:
            name, first, _, count = runs[-1]
            runs[-1] = (name, first, leaf.theta, count + 1)
        else:
            runs.append((leaf.group, leaf.theta, leaf.theta, 1))
    best_run: dict[str, tuple[str, float, float, int]] = {}
    for run in runs:
        if run[3] > best_run.get(run[0], ("", 0.0, 0.0, 0))[3]:
            best_run[run[0]] = run

    # Silhouettes sit just outside the band; crowded neighbours are pushed further out.
    # A group whose band is too short for its name gets the name beneath its silhouette instead.
    def band_caption(first: float, last: float, label: str) -> float:
        return min(15.0, (first - last + step) * (band_in + band_out) / 2 * 72 / max(len(label), 1))

    def text_width(label: str, points: float) -> float:
        return len(label) * points * 0.66 / 72

    # Radial names of tiny groups are placed first and act as obstacles for the art.
    obstacles: list[tuple[float, float, float]] = []  # r, theta, diameter
    for name, first, last, count in best_run.values():
        if count < MIN_ART_GENERA:
            length = text_width(GROUPS[name][0], 10.5)
            for k in range(3):
                obstacles.append((band_out + 0.12 + length * (k + 0.5) / 3, (first + last) / 2, length / 3 + 0.25))

    art: list[tuple[str, float, float, float, bool, bool]] = []  # group, r, theta, size, pushed, captioned
    for name, first, last, count in best_run.values():
        if count < MIN_ART_GENERA:
            continue
        mid = (first + last) / 2
        size = art_size * min(1.0, max(0.5, (first - last + step) * band_out / 2.0))
        captioned = band_caption(first, last, GROUPS[name][0]) < MIN_BAND_CAPTION
        footprint = max(size * 1.35, text_width(GROUPS[name][0], 10.5) + 0.1) if captioned else size
        r = band_out + 0.2 + footprint * 0.5
        taken = [(pr, pt, ps) for _, pr, pt, ps, _, _ in art] + obstacles
        while any(_polar_dist(r, mid, pr, pt) < (footprint + ps) * 0.52 for pr, pt, ps in taken):
            r += art_size * 0.75
        art.append((name, r, mid, footprint, r > band_out + 0.2 + footprint * 0.5, captioned))
    radius = max((r + size * 0.55 for _, r, _, size, _, _ in art), default=band_out)

    margin, title_h = 0.9, 5.0
    legend_h = _legend_height(leaves)
    width = 2 * (radius + margin)
    height = width + title_h + legend_h
    cx, cy = width / 2, legend_h + radius + margin

    plt.rcParams.update({"font.family": "serif", "font.serif": SERIF, "svg.fonttype": "none", "pdf.fonttype": 42})
    fig = plt.figure(figsize=(width, height), facecolor=PAPER)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(Rectangle((0.45, 0.45), width - 0.9, height - 0.9, fill=False, edgecolor=INK, linewidth=2.4))
    ax.add_patch(Rectangle((0.62, 0.62), width - 1.24, height - 1.24, fill=False, edgecolor=INK, linewidth=0.7))

    def polar(r: float, theta: float) -> tuple[float, float]:
        return cx + r * math.cos(theta), cy + r * math.sin(theta)

    def node_radius(node: _Node) -> float:
        if not node.children:
            return r_tip
        return 0.0 if node.depth == 0 else r_hub + (r_tip * 0.93 - r_hub) * node.depth / max_depth

    # Faint rings give the eye a sense of depth, like growth rings in a trunk.
    t = np.linspace(0, 2 * math.pi, 400)
    for r in np.linspace(r_hub, r_tip * 0.93, 7)[1:]:
        ax.plot(cx + r * np.cos(t), cy + r * np.sin(t), color=MUTED, linewidth=0.35, alpha=0.22, zorder=0)

    total = top.leaves

    def colour(node: _Node) -> str:
        return GROUPS[node.group][1] if node.group else NEUTRAL

    def line_width(node: _Node) -> float:
        return 0.35 + 3.2 * (node.leaves / total) ** 0.5

    segments, colors, widths = [], [], []
    for node in _walk(top):
        if not node.children:
            continue
        r = node_radius(node)
        thetas = [child.theta for child in node.children]
        arc = np.linspace(min(thetas), max(thetas), max(2, int(math.degrees(max(thetas) - min(thetas))) + 2))
        segments.append(np.column_stack([cx + r * np.cos(arc), cy + r * np.sin(arc)]))
        colors.append(colour(node))
        widths.append(line_width(node))
        for child in node.children:
            segments.append(np.array([polar(r, child.theta), polar(node_radius(child), child.theta)]))
            colors.append(colour(child))
            widths.append(line_width(child))
    ax.add_collection(LineCollection(segments, colors=colors, linewidths=widths, capstyle="round", zorder=2))
    ax.plot(cx, cy, marker="o", markersize=11, color=INK, zorder=5)

    for leaf in leaves:
        deg = math.degrees(leaf.theta) % 360
        flip = 90 < deg < 270
        ax.text(
            *polar(r_tip + 0.06, leaf.theta), leaf.name, fontsize=label_pt,
            fontstyle="normal" if leaf.name.startswith('"') else "italic", color=_shade(colour(leaf), 0.5),
            rotation=deg + 180 if flip else deg, rotation_mode="anchor", ha="right" if flip else "left", va="center", zorder=3,
        )

    if ages:
        # Dotted leaders run from the end of each name across the ring, so a bar can be traced to its genus.
        measure = TextToPath()
        leaders, leader_colors = [], []
        for leaf in (leaf for leaf in leaves if leaf.name in ages):
            style = "normal" if leaf.name.startswith('"') else "italic"
            name_w = measure.get_text_width_height_descent(leaf.name, FontProperties(family="serif", style=style, size=label_pt), ismath=False)[0] / 72
            c, s = math.cos(leaf.theta), math.sin(leaf.theta)
            start = r_tip + 0.06 + name_w + 0.04
            leaders.append([(cx + start * c, cy + start * s), (cx + time_out * c, cy + time_out * s)])
            leader_colors.append(colour(leaf))
        ax.add_collection(LineCollection(leaders, colors=leader_colors, linewidths=0.35, linestyles=(0, (1.2, 1.6)), alpha=0.7, zorder=1.5))
        _draw_time_ring(ax, leaves, ages, cx, cy, time_in, time_out, gap * step, colour)

    # Named clades with enough genera get a label on their node.
    # Most important first; a label that would overlap one already placed is skipped.
    halo = [patheffects.withStroke(linewidth=3.2, foreground=PAPER)]
    boxes: list[tuple[float, float, float, float]] = []
    named = [node for node in _walk(top) if node is not top and node.children and (node.leaves >= 12 or node.name in MAJOR_CLADES)]
    for node in sorted(named, key=lambda node: (node.name not in MAJOR_CLADES, -node.leaves)):
        major = node.name in MAJOR_CLADES
        size = 17 if major else 7.5 + min(5.0, node.leaves / 12)
        x, y = polar(node_radius(node), node.theta)
        half_w, half_h = len(node.name) * size * 0.29 / 72, size * 0.6 / 72
        box = (x - half_w, y - half_h, x + half_w, y + half_h)
        if any(box[0] < b[2] and b[0] < box[2] and box[1] < b[3] and b[1] < box[3] for b in boxes):
            continue
        boxes.append(box)
        ax.text(x, y, node.name, fontsize=size, color=INK, fontweight="bold" if major else "normal",
                ha="center", va="center", zorder=4, path_effects=halo)

    # Every band that can hold its name shows it, including the extra pieces of a split group.
    for name, first, last, count in runs:
        label, color = GROUPS[name][:2]
        lo, hi = math.degrees(last - step / 2), math.degrees(first + step / 2)
        ax.add_patch(Wedge((cx, cy), band_out, lo, hi, width=band_out - band_in, facecolor=color, edgecolor=PAPER, linewidth=1.0, zorder=2))
        size = band_caption(first, last, label)
        if size >= MIN_BAND_CAPTION:
            deg = math.degrees((first + last) / 2) % 360
            ax.text(*polar((band_in + band_out) / 2, (first + last) / 2), label.upper(), fontsize=size, color=PAPER, fontweight="bold",
                    rotation=deg - 90 if deg <= 180 else deg + 90, rotation_mode="anchor", ha="center", va="center", zorder=3)

    # Groups too small for art are named radially, just outside their band.
    for name, first, last, count in best_run.values():
        if count >= MIN_ART_GENERA:
            continue
        label, color = GROUPS[name][:2]
        mid = (first + last) / 2
        deg = math.degrees(mid) % 360
        flip = 90 < deg < 270
        ax.text(*polar(band_out + 0.12, mid), label.upper(), fontsize=10.5, fontweight="bold", color=_shade(color, 0.75),
                rotation=deg + 180 if flip else deg, rotation_mode="anchor", ha="right" if flip else "left", va="center", zorder=3)

    for name, r, theta, footprint, pushed, captioned in art:
        label, color, kind = GROUPS[name][:3]
        size = min(footprint / 1.35, art_size) if captioned else footprint
        x, y = polar(r, theta)
        if captioned:
            y += size * 0.12
            ax.text(x, y - size * 0.42, label.upper(), fontsize=10.5, fontweight="bold", color=_shade(color, 0.75), ha="center", va="top", zorder=3)
        if pushed:
            ax.plot(*zip(polar(band_out + 0.06, theta), polar(r - footprint * 0.45, theta)), color=color, linewidth=0.8, alpha=0.8, zorder=1)
        _draw_art(ax, kind, x, y, size, color)

    # Title block.
    ax.text(1.4, height - 1.2, root.upper(), fontsize=132, fontweight="bold", color=INK, va="top", ha="left")
    ax.text(1.48, height - 3.25, f"A consensus cladogram of {len(leaves):,} {'dinosaur ' if root == 'Dinosauria' else ''}genera", fontsize=42, fontstyle="italic", color=INK, va="top")
    ax.text(1.5, height - 4.15, f"Compiled from {source_count:,} Wikipedia article cladograms   ·   {date.today():%B %Y}", fontsize=22, color=MUTED, va="top")
    _how_to_read(ax, width - 1.4, height - 1.35, source_count)
    _legend(ax, leaves, width, legend_h, bool(ages))

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, facecolor=PAPER)
    plt.close(fig)


def _draw_time_ring(ax, leaves, ages, cx, cy, r_in, r_out, opening, colour) -> None:
    """Geologic periods as shaded rings (older inside) and each genus's range as a radial bar."""

    def radius(ma: float) -> float:
        return r_in + (TIME_START - ma) / (TIME_START - TIME_END) * (r_out - r_in)

    for name, start, end, color in PERIODS:
        ax.add_patch(Wedge((cx, cy), radius(end), 0, 360, width=radius(end) - radius(start), facecolor=color, alpha=0.14, edgecolor="none", zorder=1))
    t = np.linspace(0, 2 * math.pi, 720)
    for ma in EPOCH_BOUNDARIES:
        major = ma in {p[1] for p in PERIODS} | {PERIODS[-1][2]}
        ax.plot(cx + radius(ma) * np.cos(t), cy + radius(ma) * np.sin(t), color=MUTED,
                linewidth=0.7 if major else 0.3, alpha=0.6 if major else 0.35, zorder=1)

    envelope, core, envelope_colors, core_colors = [], [], [], []
    for leaf in leaves:
        span = ages.get(leaf.name)  # quoted names ("Velociraptor" osmolskae) are not that genus: no bar
        if not span:
            continue
        envelope_old, first, last, envelope_young = span
        c, s = math.cos(leaf.theta), math.sin(leaf.theta)
        if first - last < 0.8:  # keep very short ranges visible
            first, last = first + 0.4, last - 0.4
        for bucket, colors, (old, young) in ((envelope, envelope_colors, (envelope_old, envelope_young)), (core, core_colors, (first, last))):
            bucket.append([(cx + radius(old) * c, cy + radius(old) * s), (cx + radius(young) * c, cy + radius(young) * s)])
            colors.append(_shade(colour(leaf), 0.82))
    ax.add_collection(LineCollection(envelope, colors=envelope_colors, linewidths=1.0, alpha=0.45, zorder=2))
    ax.add_collection(LineCollection(core, colors=core_colors, linewidths=2.6, capstyle="butt", zorder=2))

    # Scale in the opening at twelve o'clock.
    for name, start, end, color in PERIODS:
        ax.text(cx, cy + (radius(start) + radius(end)) / 2, name.upper(), fontsize=8.5, fontweight="bold",
                color=_shade(color, 0.7), ha="center", va="center", zorder=3)
    for ma in (TIME_START, 201.4, 145.0, TIME_END):
        ax.text(cx + opening * radius(ma) * 0.28, cy + radius(ma) + 0.07, f"{ma:g}", fontsize=6.5, color=MUTED, ha="left", va="center", zorder=3)
    ax.text(cx, cy + r_out + 0.12, "Ma", fontsize=7, color=MUTED, ha="center", va="bottom", zorder=3)


def _time_key(ax, x: float, y: float) -> None:
    ax.text(x, y, "GEOLOGIC TIME", fontsize=20, fontweight="bold", color=INK, va="top")
    ax.text(x, y - 0.5, "Ring between names and bands; time runs outward", fontsize=13, fontstyle="italic", color=MUTED, va="top")
    y -= 1.35
    for name, start, end, color in PERIODS:
        ax.add_patch(Rectangle((x, y - 0.14), 0.7, 0.28, facecolor=color, alpha=0.22, edgecolor=_shade(color, 0.8), linewidth=0.6))
        ax.text(x + 0.9, y, f"{name}   {start:g}–{end:g} Ma", fontsize=13, color=INK, va="center")
        y -= 0.45
    y -= 0.2
    ax.plot([x, x + 0.7], [y, y], color=INK, linewidth=0.8, alpha=0.35)
    ax.plot([x + 0.18, x + 0.52], [y, y], color=INK, linewidth=3.2, solid_capstyle="butt")
    ax.text(x + 0.9, y, "Solid: central 80% of well-dated, securely identified fossils", fontsize=11, color=INK, va="center")
    ax.text(x + 0.9, y - 0.27, "Faint: their full dating windows. Source: Paleobiology Database", fontsize=11, color=MUTED, va="center")


def _polar_dist(r1: float, t1: float, r2: float, t2: float) -> float:
    return math.dist((r1 * math.cos(t1), r1 * math.sin(t1)), (r2 * math.cos(t2), r2 * math.sin(t2)))


def _how_to_read(ax, x: float, y: float, source_count: int) -> None:
    blocks = [
        ("Reading the cladogram", [
            "Genera occupy the outer ring; each internal node marks a common ancestor.",
            "The nearer the rim two branches join, the more closely related the genera.",
            "Nodes with many branches (polytomies) mark unresolved relationships.",
            "Branch lengths are not scaled to time or to morphological change.",
        ]),
        ("Method", [
            f"Topology is a weighted consensus of {source_count:,} source cladograms,",
            "constrained by a curated backbone of well-supported clades and by",
            "nomenclatural rank. Where sources conflict without a clear majority,",
            "the genus is placed at the deepest clade on which they agree.",
            "Stratigraphic ranges: Paleobiology Database (paleobiodb.org).",
        ]),
    ]
    for heading, lines in blocks:
        ax.text(x, y, heading, fontsize=24, fontweight="bold", color=INK, ha="right", va="top")
        for i, line in enumerate(lines):
            ax.text(x, y - 0.62 - i * 0.33, line, fontsize=15, color=MUTED, ha="right", va="top")
        y -= 0.62 + len(lines) * 0.33 + 0.35


def _legend_layout(leaves: list[_Node]) -> tuple[dict[str, int], dict[str, list[str]], dict[str, int], int]:
    counts: dict[str, int] = {}
    for leaf in leaves:
        if leaf.group:
            counts[leaf.group] = counts.get(leaf.group, 0) + 1
    members = {lineage: [name for name, spec in GROUPS.items() if spec[3] == lineage and counts.get(name)] for lineage in LINEAGES}
    subcolumns = {lineage: 2 if len(names) > 6 else 1 for lineage, names in members.items() if names}
    rows = max((math.ceil(len(members[lineage]) / cols) for lineage, cols in subcolumns.items()), default=0)
    return counts, members, subcolumns, rows


def _legend_height(leaves: list[_Node]) -> float:
    return 3.1 + _legend_layout(leaves)[3] * LEGEND_ROW


def _legend(ax, leaves: list[_Node], width: float, legend_h: float, time_key: bool = False) -> None:
    counts, members, subcolumns, rows = _legend_layout(leaves)
    sub_w = (width - 2.8) / max(1, sum(subcolumns.values()))
    x0 = 1.4
    for lineage, cols in subcolumns.items():
        heading, gloss = LINEAGES[lineage]
        top = legend_h - 0.55
        ax.text(x0, top, heading.upper(), fontsize=20, fontweight="bold", color=INK, va="top")
        for i, line in enumerate(textwrap.wrap(gloss, int(cols * sub_w * 72 / (13 * 0.5)) - 4)[:2]):
            ax.text(x0, top - 0.5 - i * 0.27, line, fontsize=13, fontstyle="italic", color=MUTED, va="top")
        per_column = math.ceil(len(members[lineage]) / cols)
        chars = int((sub_w - 1.35) * 72 / (10.5 * 0.5))
        for i, name in enumerate(members[lineage]):
            label, color, kind, _, description = GROUPS[name]
            x = x0 + (i // per_column) * sub_w
            yy = top - 1.75 - (i % per_column) * LEGEND_ROW
            _draw_art(ax, kind, x + 0.5, yy, 0.9, color)
            ax.text(x + 1.15, yy + 0.24, label, fontsize=15, color=INK, va="center")
            # "Early"/"Other" groups only count genera outside the named subgroups.
            detail = f"{counts[name]} other {name}" if label.startswith(("Early", "Other")) else f"{name}  ·  {counts[name]} genera"
            ax.text(x + 1.15, yy - 0.03, detail, fontsize=10.5, color=MUTED, fontstyle="italic", va="center")
            for j, line in enumerate(textwrap.wrap(description, chars)[:2]):
                ax.text(x + 1.15, yy - 0.28 - j * 0.2, line, fontsize=10.5, color=INK, alpha=0.8, va="center")
        if time_key and lineage == "Dinosauria":
            _time_key(ax, x0, top - 1.75 - len(members[lineage]) * LEGEND_ROW - 0.2)
        x0 += cols * sub_w


def _draw_art(ax, kind: str, x: float, y: float, size: float, color: str) -> None:
    """Draw a silhouette centred on (x, y), `size` inches wide."""
    polygons = silhouette(kind)
    points = np.vstack(polygons)
    lo, hi = points.min(axis=0), points.max(axis=0)
    scale = size / max(hi[0] - lo[0], hi[1] - lo[1])
    centre = (lo + hi) / 2
    for poly in polygons:
        ax.add_patch(Polygon((poly - centre) * scale + (x, y), closed=True, facecolor=color, edgecolor="none", zorder=3))


def _find(tree: dict[str, Any], name: str) -> dict[str, Any] | None:
    for child, children in tree.items():
        if child == name:
            return children
        found = _find(children, name)
        if found is not None:
            return found
    return None


def _build(name: str, children: dict[str, Any], depth: int, group: str | None) -> _Node:
    group = name if name in GROUPS else group
    node = _Node(name, depth, group=group)
    node.children = [_build(child, grandchildren, depth + 1, group) for child, grandchildren in children.items()]
    node.leaves = sum(child.leaves for child in node.children) or 1
    return node


def _group_order(node: _Node) -> str | None:
    """Draw first the children that continue this node's colour, so each group stays one band.

    Returns the colour group of the node's first genus.
    """
    firsts = {id(child): _group_order(child) for child in node.children}
    node.children.sort(key=lambda child: firsts[id(child)] != node.group)
    return firsts[id(node.children[0])] if node.children else node.group


def _assign_internal_angles(node: _Node) -> float:
    if node.children:
        thetas = [_assign_internal_angles(child) for child in node.children]
        node.theta = (min(thetas) + max(thetas)) / 2
    return node.theta


def _walk(node: _Node):
    yield node
    for child in node.children:
        yield from _walk(child)


def _leaves(node: _Node) -> list[_Node]:
    return [n for n in _walk(node) if not n.children]


def _depths(node: _Node):
    return (n.depth for n in _walk(node) if n.children)


def _shade(color: str, factor: float) -> tuple[float, float, float]:
    return tuple(channel * factor for channel in to_rgb(color))
