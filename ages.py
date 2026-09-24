"""Stratigraphic ranges for genera from the Paleobiology Database (PBDB), cached as JSON.

Ranges are built from individual fossil occurrences, not from PBDB's taxon summary: the
summary includes every doubtful referral ever made, which stretched e.g. Stegosaurus into the
Late Cretaceous. Only occurrences identified with certainty to genus (no "cf." or "?") count,
and of those only the well-dated ones (dating window at most WELL_DATED My) when there are any.

Each genus maps to [envelope_old, core_old, core_young, envelope_young] in Ma:
the core spans the 10th-90th percentile of occurrence midpoints, the envelope the full dating
windows of the occurrences inside that core. Genera without usable occurrences fall back to
the taxon summary.
"""
from __future__ import annotations

from pathlib import Path
from urllib.request import Request, urlopen
import json

AGES_FILE = Path(__file__).resolve().parent / "genus_ages.json"
PBDB = "https://paleobiodb.org/data1.2"
OCCURRENCES_URL = f"{PBDB}/occs/list.json?base_name=Dinosauria&idreso=genus&idqual=certain&vocab=pbdb"
TAXA_URL = f"{PBDB}/taxa/list.json?base_name=Dinosauria&rank=genus&show=app&vocab=pbdb"
WELL_DATED = 12.0  # My

# Wastebasket genera whose PBDB occurrences still include many historical referrals.
AGE_OVERRIDES = {
    "Megalosaurus": [170.3, 168.3, 166.1, 164.7],  # Bathonian of England
    "Cetiosaurus": [171.6, 170.3, 166.1, 164.7],  # Bajocian-Bathonian of England
}


def fetch_pbdb_ages(timeout: float = 300) -> dict[str, list[float]]:
    occurrences: dict[str, list[tuple[float, float]]] = {}
    for record in _get(OCCURRENCES_URL, timeout):
        name = (record.get("accepted_name") or record.get("identified_name") or "").split()
        if name and "max_ma" in record:
            occurrences.setdefault(name[0], []).append((float(record["max_ma"]), float(record["min_ma"])))
    ages = {genus: occurrence_range(windows) for genus, windows in occurrences.items()}

    for record in _get(TAXA_URL, timeout):
        if "firstapp_max_ma" not in record:
            continue
        first_max = float(record["firstapp_max_ma"])
        first_min = float(record.get("firstapp_min_ma", first_max))
        last_max = float(record.get("lastapp_max_ma", first_min))
        last_min = float(record.get("lastapp_min_ma", last_max))
        summary = [first_max, (first_max + first_min) / 2, (last_max + last_min) / 2, last_min]
        accepted = record.get("accepted_name")
        if record.get("taxon_name") == accepted:
            ages.setdefault(accepted, summary)
        for name in (accepted, record.get("taxon_name")):  # junior synonyms share the range
            if name and accepted in ages:
                ages.setdefault(name, ages[accepted])
    ages.update(AGE_OVERRIDES)
    return dict(sorted(ages.items()))


def occurrence_range(windows: list[tuple[float, float]]) -> list[float]:
    dated = [w for w in windows if w[0] - w[1] <= WELL_DATED] or windows
    mids = sorted((old + young) / 2 for old, young in dated)
    core_old, core_young = _percentile(mids, 0.9), _percentile(mids, 0.1)
    core = [w for w in dated if core_young <= (w[0] + w[1]) / 2 <= core_old] or dated
    return [max(w[0] for w in core), round(core_old, 2), round(core_young, 2), min(w[1] for w in core)]


def _percentile(sorted_values: list[float], q: float) -> float:
    k = (len(sorted_values) - 1) * q
    low = int(k)
    high = min(low + 1, len(sorted_values) - 1)
    return sorted_values[low] + (sorted_values[high] - sorted_values[low]) * (k - low)


def _get(url: str, timeout: float) -> list[dict]:
    request = Request(url, headers={"User-Agent": "DinoTree/1.0 (dinosaur cladogram poster)"})
    with urlopen(request, timeout=timeout) as response:
        return json.load(response)["records"]


def load_ages(path: Path = AGES_FILE) -> dict[str, list[float]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def save_ages(ages: dict[str, list[float]], path: Path = AGES_FILE) -> None:
    path.write_text(json.dumps(ages, indent=0, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    # The percentile core ignores stray outliers.
    assert occurrence_range([(72.2, 66.0)] * 19 + [(83.6, 72.2)])[1:3] == [69.1, 69.1]
    assert occurrence_range([(100.0, 60.0)]) == [100.0, 80.0, 80.0, 60.0]
    print("ok")
