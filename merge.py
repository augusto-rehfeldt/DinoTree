"""Merge many small dinosaur cladograms into one tree by weighing containment evidence.

Every source says "clade C contains taxon X" in some way:

* flat lists (flattened Wikipedia cladograms): the first name contains the rest,
  unless the list is rooted on an outgroup (it starts with a genus, or the backbone
  says the first clade cannot contain another clade in the list);
* nested lineages: every ancestor contains every descendant;
* name stems: Tyrannosaurini < Tyrannosaurinae < Tyrannosauridae < Tyrannosauroidea,
  and the genus Tyrannosaurus sits inside all of them;
* a curated backbone of well-supported clades, which outranks everything;
* every dino_dict entry is a dinosaur genus (weakest fallback).

Claims are inserted strongest first into a containment graph. A claim is rejected
when a genus would contain something, when ranks are inverted (a subfamily cannot
contain a family), when it moves a backbone clade, or when it would contradict
stronger claims (a cycle). Each taxon is then attached to its most specific
well-supported container. Everything the final tree disagrees with lands in the
report, together with the sources that made the claim.

List order is not used: measured on the data, "the taxa after a clade name belong
to it" held only ~39% of the time, so it would add more errors than it fixes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from typing import Any
import re

from tree_utils import MergeRules, _iter_raw_children, load_merge_rules


# Outgroups, so evidence from wider cladograms cannot pull them inside Dinosauria.
OUTGROUPS = """
Animalia
  Chordata
    Tetrapoda
      Amniota
        Reptilia
          Euparkeria
          Archosauria
            Pseudosuchia
            Dinosauromorpha
              Silesauridae
              Dinosauria
"""

BACKBONE = """
Dinosauria
  Saurischia
    Herrerasauridae
    Eusaurischia
      Sauropodomorpha
        Plateosauria
          Plateosauridae
          Massopoda
            Massospondylidae
            Sauropodiformes
              Sauropoda
                Lessemsauridae
                Gravisauria
                  Vulcanodontidae
                  Eusauropoda
                    Mamenchisauridae
                    Turiasauria
                    Neosauropoda
                      Diplodocoidea
                        Rebbachisauridae
                          Limaysaurinae
                          Nigersaurinae
                        Flagellicaudata
                          Dicraeosauridae
                          Diplodocidae
                            Apatosaurinae
                            Diplodocinae
                      Macronaria
                        Camarasauromorpha
                          Camarasauridae
                          Titanosauriformes
                            Brachiosauridae
                            Somphospondyli
                              Euhelopodidae
                              Titanosauria
                                Lithostrotia
                                  Saltasauridae
                                    Saltasaurinae
                                    Opisthocoelicaudiinae
                                  Colossosauria
                                    Lognkosauria
                                    Rinconsauria
      Theropoda
        Neotheropoda
          Coelophysoidea
            Coelophysidae
          Averostra
            Ceratosauria
              Ceratosauridae
              Abelisauroidea
                Noasauridae
                  Elaphrosaurinae
                  Noasaurinae
                Abelisauridae
                  Majungasaurinae
                  Brachyrostra
            Tetanurae
              Orionides
                Megalosauroidea
                  Piatnitzkysauridae
                  Megalosauridae
                    Megalosaurinae
                  Spinosauridae
                    Baryonychinae
                    Spinosaurinae
                Avetheropoda
                  Carnosauria
                    Allosauroidea
                      Metriacanthosauridae
                      Allosauria
                        Allosauridae
                        Carcharodontosauriformes
                          Carcharodontosauria
                            Neovenatoridae
                            Carcharodontosauridae
                  Coelurosauria
                    Compsognathidae
                    Tyrannoraptora
                      Tyrannosauroidea
                        Proceratosauridae
                        Eutyrannosauria
                          Tyrannosauridae
                            Albertosaurinae
                            Tyrannosaurinae
                      Maniraptoromorpha
                        Maniraptoriformes
                          Ornithomimosauria
                            Deinocheiridae
                            Ornithomimidae
                          Maniraptora
                            Alvarezsauria
                              Alvarezsauroidea
                                Alvarezsauridae
                            Therizinosauria
                              Therizinosauroidea
                                Therizinosauridae
                            Pennaraptora
                              Oviraptorosauria
                                Caudipteridae
                                Caenagnathoidea
                                  Caenagnathidae
                                  Oviraptoridae
                              Paraves
                                Scansoriopterygidae
                                Dromaeosauridae
                                  Microraptorinae
                                  Unenlagiinae
                                  Halszkaraptorinae
                                  Eudromaeosauria
                                    Saurornitholestinae
                                    Dromaeosaurinae
                                    Velociraptorinae
                                Troodontidae
                                Avialae
  Ornithischia
    Heterodontosauridae
    Genasauria
      Thyreophora
        Eurypoda
          Stegosauria
            Stegosauridae
          Ankylosauria
            Ankylosauridae
              Ankylosaurinae
            Nodosauridae
              Nodosaurinae
      Neornithischia
        Thescelosauridae
          Orodrominae
        Cerapoda
          Ornithopoda
            Clypeodonta
              Iguanodontia
                Rhabdodontomorpha
                  Rhabdodontidae
                Dryomorpha
                  Dryosauridae
                  Ankylopollexia
                    Camptosauridae
                    Styracosterna
                      Hadrosauriformes
                        Hadrosauroidea
                          Hadrosauromorpha
                            Hadrosauridae
                              Saurolophinae
                              Lambeosaurinae
          Marginocephalia
            Pachycephalosauria
              Pachycephalosauridae
            Ceratopsia
              Leptoceratopsidae
              Coronosauria
                Protoceratopsidae
                Ceratopsidae
                  Centrosaurinae
                  Chasmosaurinae
"""

# Spelling fixes and synonyms seen in the source data.
SYNONYMS = {
    "Ceratopidae": "Ceratopsidae",
    "Tetrapods": "Tetrapoda",
    "Acristravus": "Acristavus",
    "Leallynasaura": "Leaellynasaura",
    "Microraptoria": "Microraptorinae",
    "Unenlagiinia": "Unenlagiinae",
    "Unenlagiidae": "Unenlagiinae",
    "Sinraptoridae": "Metriacanthosauridae",
    # Older or competing names for the same clades.
    "Hadrosaurinae": "Saurolophinae",
    "Panoplosauridae": "Nodosauridae",
    "Polacanthidae": "Polacanthinae",
    "Ingenia": "Ajancingenia",  # preoccupied name, replaced in 2013
    '"Ingenia" yanshini': "Ajancingenia",
}

# Genera whose placement the sources get wrong or never state; these outrank all claims.
PINS = {
    "Euparkeria": "Reptilia",  # outgroup genus that heads several lists
    "Hadrosaurus": "Hadrosauridae",  # outside Saurolophinae, whatever its name stem suggests
    "Saurophaganax": "Allosauridae",  # only listed in a mixed, outgroup-rooted list
    "Staurikosaurus": "Herrerasauridae",  # its lists also claim wider outgroup clades
}

# Unranked clade names without a telling suffix, so they are never mistaken for genera.
KNOWN_CLADES = {"Parapredentata", "Hypsilophodontia", "Elasmaria"}

# Rank by ending; a genus is rank 0. Unranked clades (-ia, -poda, ...) return None.
RANK_SUFFIXES = (("oidea", 4), ("idae", 3), ("inae", 2), ("ini", 1))
CLADE_SUFFIXES = ("sauria", "formes", "morpha", "morphae", "poda", "cephalia", "ptera")
GENUS_ENDINGS = (
    "saurus", "ceratops", "long", "raptor", "mimus", "don", "venator", "titan", "pelta", "onyx",
    "nykus", "saura", "suchus", "gnathus", "pteryx", "lophus", "cephale", "tholus", "dromeus", "cursor",
)
# Labels that pass the name pattern but are not taxa (formations, fields, informal groups).
NOT_TAXA = {
    "Basal", "Incertae", "Unnamed", "Indet", "Cloverly", "Kaiparowits", "Paleontology",
    "Amphibians", "Archosaurs", "Crocodilians", "Mammals", "Reptiles", "Synapsids",
}

BACKBONE_WEIGHT = 1000.0
STEM_WEIGHT = 10.0  # a family-group name is built on its type genus, so this is nomenclature
NESTED_WEIGHT = 2.0
LIST_WEIGHT = 1.0
OWN_PAGE_WEIGHT = 3.0  # a genus's own article cladogram is the most careful source about it
ENTRY_WEIGHT = 0.1
# ponytail: fixed ratio; a branch must carry 3x the support of its best rival to be followed,
# otherwise the taxon stays at their common ancestor (a polytomy). Tune per dataset.
DOMINANCE = 3.0
# Share of list-mates a clade must hold to adopt a taxon that has no direct claim.
COMPANY_SHARE = 0.6

_TAXON = re.compile(r"[A-Z][a-z]+(?:-[A-Z][a-z]+)?")
_QUOTED = re.compile(r'"[A-Z][a-z]+" [a-z]+')
_ABBREVIATED = re.compile(r"([A-Z])\. [a-z]+")


@dataclass
class MergeResult:
    tree: dict[str, Any]
    report: dict[str, Any]
    source_count: int


@dataclass
class _Claim:
    weight: float = 0.0
    sources: list[str] = field(default_factory=list)


def merge_dino_trees(dino_dict: dict[str, Any], rules: MergeRules | None = None) -> MergeResult:
    rules = rules or load_merge_rules(None)
    dropped: dict[str, str] = {}
    claims: dict[tuple[str, str], _Claim] = {}

    def claim(container: str, taxon: str, weight: float, source: str) -> None:
        if container != taxon:
            entry = claims.setdefault((container, taxon), _Claim())
            entry.weight += weight
            entry.sources.append(source)

    # 1. Clean every source into flat lists and nested lineages.
    entry_genera: set[str] = set()
    own_name: dict[str, str | None] = {}
    lists: list[tuple[str, list[str]]] = []
    lineages: list[tuple[str, list[str]]] = []
    for entry, data in dino_dict.items():
        entry_name = own_name[str(entry)] = _clean(str(entry), rules, dropped)
        if entry_name:
            entry_genera.add(entry_name)
        items = _iter_raw_children(data)
        if _is_flat(items):
            taxa = _clean_list([name for name, _ in items], rules, dropped)
            if entry_name and entry_name not in taxa:
                taxa.append(entry_name)
            lists.append((str(entry), taxa))
        else:
            for path in _iter_paths(items):
                cleaned = [name for name in (_clean(raw, rules, dropped) for raw in path) if name]
                lineages.append((str(entry), list(dict.fromkeys(cleaned))))

    names = entry_genera | {name for _, taxa in lists + lineages for name in taxa}
    backbone_parent = _present_backbone(_parse_outline(OUTGROUPS) | _parse_outline(BACKBONE), names)
    backbone_ancestors = {name: set(_walk_up(name, backbone_parent)) for name in backbone_parent}
    in_backbone = backbone_parent.keys() | set(backbone_parent.values())
    pins = {taxon: home for taxon, home in PINS.items() if taxon in names and home in names}
    genera = entry_genera | set(_species_genera(dino_dict, rules)) | pins.keys() | {
        name for name in names
        if name not in in_backbone and _suffix_rank(name) is None and name.endswith(GENUS_ENDINGS)
    }

    # 2. Collect claims.
    outgroup_rooted: dict[str, str] = {}
    for entry, taxa in lists:
        if len(taxa) < 2:
            continue
        head = taxa[0]
        if head in genera:
            outgroup_rooted[entry] = f"starts with the genus {head}"
            continue
        misplaced = [
            taxon for taxon in taxa[1:]
            if head in in_backbone and taxon in backbone_parent and taxon not in genera and head not in backbone_ancestors[taxon]
        ]
        if misplaced:
            outgroup_rooted[entry] = f"backbone: {head} does not contain {misplaced[0]}"
            continue
        head_rank = _suffix_rank(head)
        inverted = [
            taxon for taxon in taxa[1:]
            if head_rank is not None and (_suffix_rank(taxon) or 0) >= head_rank
        ]
        if inverted:
            outgroup_rooted[entry] = f"rank: {head} cannot contain {inverted[0]}"
            continue
        for taxon in taxa[1:]:
            claim(head, taxon, OWN_PAGE_WEIGHT if taxon == own_name[entry] else LIST_WEIGHT, entry)

    for entry, path in lineages:
        for i, ancestor in enumerate(path):
            for descendant in path[i + 1:]:
                claim(ancestor, descendant, NESTED_WEIGHT, entry)

    containers = {container for container, _ in claims}
    genera |= {
        name for name in names - containers - in_backbone - KNOWN_CLADES
        if _suffix_rank(name) is None and not name.endswith(CLADE_SUFFIXES)
    }

    def rank(name: str) -> int | None:
        return 0 if name in genera else _suffix_rank(name)

    # Stems link each name only to the next rank up (Saltasaurus -> Saltasaurini -> Saltasaurinae),
    # so a stray superfamily name cannot compete with the family chain below it.
    stems = {container: stem for container in names if (stem := _stem(container))}
    for taxon in names:
        taxon_rank = rank(taxon)
        if taxon_rank is None:
            continue
        matches = [
            container for container, stem in stems.items()
            if container != taxon and rank(container) > taxon_rank and _stem_matches(taxon, stem)
        ]
        lowest = min((rank(container) for container in matches), default=None)
        for container in matches:
            if rank(container) == lowest:
                claim(container, taxon, STEM_WEIGHT, "name stem")

    for child, parent in backbone_parent.items():
        claim(parent, child, BACKBONE_WEIGHT, "backbone")
    for taxon, home in pins.items():
        claim(home, taxon, BACKBONE_WEIGHT, "curated pin")

    if "Dinosauria" in names:
        for genus in entry_genera:
            claim("Dinosauria", genus, ENTRY_WEIGHT, "dino_dict entry")

    # 3. Insert claims strongest first; reject impossible or cycle-forming ones.
    children: dict[str, set[str]] = {}
    accepted: dict[str, dict[str, float]] = {}
    rejected: list[dict[str, Any]] = []

    def reaches(start: str, target: str) -> bool:
        stack, seen = [start], {start}
        while stack:
            for child in children.get(stack.pop(), ()):
                if child == target:
                    return True
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
        return False

    for (container, taxon), entry in sorted(claims.items(), key=lambda item: (-item[1].weight, item[0])):
        container_rank, taxon_rank = rank(container), rank(taxon)
        if container_rank == 0:
            reason = "a genus cannot contain other taxa"
        elif container_rank is not None and taxon_rank is not None and container_rank <= taxon_rank:
            reason = "rank inversion"
        elif taxon in backbone_parent and container not in backbone_ancestors[taxon]:
            reason = "conflicts with the curated backbone"
        elif taxon in pins and container != pins[taxon] and container not in backbone_ancestors.get(pins[taxon], ()):
            reason = "conflicts with a curated pin"
        elif reaches(taxon, container):
            reason = "contradicts stronger evidence"
        else:
            children.setdefault(container, set()).add(taxon)
            accepted.setdefault(taxon, {})[container] = entry.weight
            continue
        rejected.append(_claim_record(container, taxon, entry, reason))

    def below(start: str) -> set[str]:
        found, stack = set(), [start]
        while stack:
            for child in children.get(stack.pop(), ()):
                if child not in found:
                    found.add(child)
                    stack.append(child)
        return found

    # A clade known only from its name stem (Shamosaurinae) borrows its members' containers;
    # otherwise it would become a rival branch that drags its own type genus upwards.
    unplaced: list[str] = []
    outgroups = {"Dinosauria", *backbone_ancestors.get("Dinosauria", ())}
    for name in sorted(names - accepted.keys() - outgroups):
        members = below(name)
        borrowed: dict[str, float] = {}
        for member in children.get(name, ()) if _suffix_rank(name) else ():
            for container, weight in accepted.get(member, {}).items():
                container_rank = rank(container)
                if container != name and container not in members and not (
                    container_rank is not None and rank(name) is not None and container_rank <= rank(name)
                ):
                    borrowed[container] = borrowed.get(container, 0.0) + weight
        for container, weight in borrowed.items():
            if not reaches(name, container):
                claim(container, name, weight, "inherited from its members")
                children.setdefault(container, set()).add(name)
                accepted.setdefault(name, {})[container] = weight
        if name not in accepted and "Dinosauria" in names and members & entry_genera and "Dinosauria" not in members:
            children.setdefault("Dinosauria", set()).add(name)
            accepted[name] = {"Dinosauria": ENTRY_WEIGHT}
            unplaced.append(name)
    dinosaurs = {"Dinosauria"} | below("Dinosauria") if "Dinosauria" in names else set()

    @cache
    def descendants(name: str) -> frozenset[str]:
        found = set(children.get(name, ()))
        for child in children.get(name, ()):
            found |= descendants(child)
        return frozenset(found)

    def place(options: dict[str, float], under) -> str | None:
        """Walk down nested containers, following a branch only while it clearly dominates."""
        scope, current = set(options), None
        while scope:
            tops = [option for option in scope if not any(option in under(other) for other in scope)]
            mass = {top: sum(w for option, w in options.items() if option == top or option in under(top)) for top in tops}
            tops.sort(key=lambda top: (-mass[top], top))
            if len(tops) > 1 and mass[tops[0]] < DOMINANCE * mass[tops[1]]:
                shared = [name for name in children if all(top in under(name) for top in tops)]
                return min(shared, key=lambda name: (len(under(name)), name), default=current)
            current = tops[0]
            scope = {option for option in scope if option in under(current)}
        return current

    # 4. Place clades on the claim graph, then genera on the resulting clade tree, so that
    #    support is never counted twice through two paths to the same clade.
    parent: dict[str, str] = {}
    for taxon, options in accepted.items():
        if taxon not in genera and (home := place(options, descendants)):
            parent[taxon] = home

    clade_children: dict[str, set[str]] = {}
    for clade, home in parent.items():
        clade_children.setdefault(home, set()).add(clade)

    @cache
    def clades_under(name: str) -> frozenset[str]:
        found = set(clade_children.get(name, ()))
        for child in clade_children.get(name, ()):
            found |= clades_under(child)
        return frozenset(found)

    for taxon, options in accepted.items():
        if taxon not in genera:
            continue
        if taxon in entry_genera and dinosaurs & options.keys():
            # Every dino_dict entry is a dinosaur: containers outside Dinosauria are outgroup noise.
            options = {option: weight for option, weight in options.items() if option in dinosaurs}
        if home := place(options, clades_under):
            parent[taxon] = home

    # Genera with no real claim (only in outgroup-rooted lists) are placed by the company they keep.
    by_company: dict[str, str] = {}
    for taxon in sorted(names):
        if taxon in genera and taxon not in pins and set(accepted.get(taxon, {})) <= {"Dinosauria"}:
            home = _company_parent(taxon, lists, parent, rank)
            if home:
                parent[taxon] = by_company[taxon] = home

    outvoted = []
    for taxon, options in accepted.items():
        ancestors = set(_walk_up(taxon, parent))
        for container in options:
            if container not in ancestors:
                record = _claim_record(container, taxon, claims[(container, taxon)], "outvoted by a more specific or stronger placement")
                record["placed_under"] = parent.get(taxon)
                outvoted.append(record)

    # 5. Build, clean up and order the tree.
    nodes: dict[str, dict[str, Any]] = {name: {} for name in names}
    for taxon, container in parent.items():
        nodes[container][taxon] = nodes[taxon]
    tree = {name: nodes[name] for name in sorted(names - parent.keys())}
    _apply_rules(tree, rules)
    empty_clades = _prune_empty_clades(tree, genera)
    unplaced += [name for name, sub in tree.items() if not sub]
    _ladderize(tree)

    contradictions = sorted(rejected + outvoted, key=lambda record: (-record["support"], record["taxon"]))
    report = {
        "summary": {
            "sources": len(dino_dict),
            "outgroup_rooted_lists": len(outgroup_rooted),
            "taxa": _count(tree),
            "claims": len(claims),
            "rejected_claims": len(rejected),
            "outvoted_claims": len(outvoted),
            "placed_by_company": len(by_company),
            "unplaced_taxa": len(unplaced),
            "dropped_names": len(dropped),
        },
        "contradictions": contradictions,
        "placed_by_company": by_company,
        "unplaced_taxa": sorted(unplaced),
        "empty_clades": empty_clades,
        "outgroup_rooted_lists": dict(sorted(outgroup_rooted.items())),
        "dropped_names": dict(sorted(dropped.items())),
    }
    return MergeResult(tree=tree, report=report, source_count=len(dino_dict))


def _clean(raw: str, rules: MergeRules, dropped: dict[str, str]) -> str | None:
    """Return a taxon name (the genus for species), or None for non-taxa and record why."""
    raw = rules.rename.get(raw, raw)
    name = re.sub(r"\[[^\]]*\]", "", raw.replace("\xa0", " "))
    name = " ".join(name.replace("†", "").replace("?", "").split())
    name = rules.rename.get(name, SYNONYMS.get(name, name))
    if not name:
        dropped.setdefault(raw, "empty or footnote")
        return None
    if any(ch.isdigit() for ch in name):
        dropped.setdefault(raw, "specimen number")
        return None
    if _QUOTED.fullmatch(name):
        return name
    words = name.split()
    if len(words) >= 2 and _TAXON.fullmatch(words[0]) and all(word.islower() for word in words[1:]):
        return rules.rename.get(words[0], SYNONYMS.get(words[0], words[0]))
    if _TAXON.fullmatch(name) and name not in NOT_TAXA:
        return name
    dropped.setdefault(raw, "not a taxon name (locality, formation or informal label)")
    return None


def _clean_list(raw_names: list[str], rules: MergeRules, dropped: dict[str, str]) -> list[str]:
    """Clean one flat list; expand "P. tubicen" when one genus in the list starts with P."""
    cleaned = {raw: _clean(raw, rules, {}) for raw in raw_names}
    genera = {name for name in cleaned.values() if name and _suffix_rank(name) is None and not name.endswith(CLADE_SUFFIXES)}
    taxa: list[str] = []
    for raw, name in cleaned.items():
        if name is None:
            match = _ABBREVIATED.fullmatch(" ".join(raw.replace("\xa0", " ").split()))
            options = [genus for genus in genera if match and genus.startswith(match.group(1))]
            if len(options) != 1:
                _clean(raw, rules, dropped)
                if match:
                    dropped[raw] = "abbreviated species without a unique genus in its list"
                continue
            name = options[0]
        if name not in taxa:
            taxa.append(name)
    return taxa


def _species_genera(dino_dict: dict[str, Any], rules: MergeRules):
    """Yield the genus of every species name, so it is known to be a genus."""
    for data in dino_dict.values():
        for path in _iter_paths(_iter_raw_children(data)):
            if " " in path[-1].replace("\xa0", " ").strip():
                name = _clean(path[-1], rules, {})
                if name and '"' not in name:
                    yield name


def _is_flat(items: list[tuple[str, Any]]) -> bool:
    return len(items) >= 2 and all(not _iter_raw_children(children) for _, children in items)


def _iter_paths(items: list[tuple[str, Any]], prefix: tuple[str, ...] = ()):
    for name, children in items:
        path = (*prefix, name)
        yield path
        yield from _iter_paths(_iter_raw_children(children), path)


def _parse_outline(text: str) -> dict[str, str]:
    """Indented outline (two spaces per level) to a child -> parent map."""
    parents: dict[str, str] = {}
    stack: list[str] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        depth = (len(line) - len(line.lstrip())) // 2
        del stack[depth:]
        if stack:
            parents[line.strip()] = stack[-1]
        stack.append(line.strip())
    return parents


def _present_backbone(backbone: dict[str, str], names: set[str]) -> dict[str, str]:
    """Backbone edges between names present in the data, skipping absent levels."""
    edges = {}
    for child, parent in backbone.items():
        if child not in names:
            continue
        while parent is not None and parent not in names:
            parent = backbone.get(parent)
        if parent is not None:
            edges[child] = parent
    return edges


def _company_parent(taxon: str, lists: list[tuple[str, list[str]]], parent: dict[str, str], rank) -> str | None:
    """Deepest clade that holds at least COMPANY_SHARE of the taxon's list-mates."""
    votes: dict[str, int] = {}
    mates = 0
    for _, taxa in lists:
        if taxon not in taxa:
            continue
        for mate in taxa:
            if mate != taxon and mate in parent:
                mates += 1
                for clade in _walk_up(mate, parent):
                    votes[clade] = votes.get(clade, 0) + 1
    taxon_rank = rank(taxon)
    eligible = [
        clade for clade, count in votes.items()
        if count >= COMPANY_SHARE * mates and clade != taxon and taxon not in _walk_up(clade, parent)
        and not (taxon_rank is not None and rank(clade) is not None and rank(clade) <= taxon_rank)
    ]
    return max(eligible, key=lambda clade: (len(_walk_up(clade, parent)), clade), default=None)


def _walk_up(name: str, parent: dict[str, str]) -> list[str]:
    path = []
    while name in parent:
        name = parent[name]
        path.append(name)
    return path


def _suffix_rank(name: str) -> int | None:
    if _TAXON.fullmatch(name):
        for suffix, rank in RANK_SUFFIXES:
            if name.endswith(suffix):
                return rank
    return None


def _stem(name: str) -> str | None:
    if _suffix_rank(name) is None:
        return None
    for suffix, _ in RANK_SUFFIXES:
        if name.endswith(suffix):
            stem = name[: -len(suffix)]
            return stem if len(stem) >= 5 else None
    return None


def _stem_matches(taxon: str, stem: str) -> bool:
    # Troodontidae -> Troodon, Baryonychinae -> Baryonyx, Archaeopterygidae -> Archaeopteryx
    return (
        taxon.startswith(stem)
        or (stem.endswith("odont") and taxon == stem[:-1])
        or (stem.endswith("ch") and taxon == stem[:-2] + "x")
        or (stem.endswith("yg") and taxon == stem[:-1] + "x")
    )


def _claim_record(container: str, taxon: str, claim: _Claim, reason: str) -> dict[str, Any]:
    return {
        "taxon": taxon,
        "claimed_container": container,
        "support": round(claim.weight, 2),
        "sources": sorted(set(claim.sources))[:8],
        "reason": reason,
    }


def _apply_rules(tree: dict[str, Any], rules: MergeRules) -> None:
    for name in list(tree):
        children = tree[name]
        _apply_rules(children, rules)
        if name in rules.remove or (name in rules.drop_leaves and not children):
            del tree[name]
        elif name in rules.collapse:
            del tree[name]
            for child, grandchildren in children.items():
                tree.setdefault(child, {}).update(grandchildren)


def _prune_empty_clades(tree: dict[str, Any], genera: set[str]) -> list[str]:
    removed = []
    for name in list(tree):
        removed += _prune_empty_clades(tree[name], genera)
        if not tree[name] and name not in genera:
            del tree[name]
            removed.append(name)
    return sorted(removed)


def _ladderize(tree: dict[str, Any]) -> None:
    """Order siblings small-to-large so the drawn tree reads as a clean ladder."""
    for children in tree.values():
        _ladderize(children)
    ordered = sorted(tree.items(), key=lambda item: (_count(item[1]), item[0]))
    tree.clear()
    tree.update(ordered)


def _count(tree: dict[str, Any]) -> int:
    return sum(1 + _count(children) for children in tree.values())
