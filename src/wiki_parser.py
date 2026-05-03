import requests
from bs4 import BeautifulSoup
import re
import concurrent.futures
from typing import Dict, Any, List

REQUEST_TIMEOUT_SECONDS = 10
HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) DinoTreeBot/1.0"}


def clean_taxon_name(name: str) -> str:
    name = re.sub(r'\[.*?\]', '', name)
    name = name.replace('†', '').replace('?', '').strip()
    return name

def fetch_dinosaur_genera_links() -> Dict[str, str]:
    url = "https://en.wikipedia.org/wiki/List_of_dinosaur_genera"
    resp = requests.get(url, headers=HEADERS, timeout=REQUEST_TIMEOUT_SECONDS)
    resp.raise_for_status()
    
    # Matching pattern: <li><i><a href="/wiki/... " title="...">Name</a></i>
    matches = re.findall(r'<li><i><a href="(/wiki/[^"]+)"[^>]*>([^<]+)</a></i>', resp.text)
    
    genera_links = {}
    for href, raw_name in matches:
        name = clean_taxon_name(raw_name)
        if name:
            genera_links[name.lower()] = "https://en.wikipedia.org" + href
            
    return genera_links

def fetch_genus_lineage(url: str) -> List[str]:
    try:
        resp = requests.get(url, headers=HEADERS, timeout=REQUEST_TIMEOUT_SECONDS)
        resp.raise_for_status()
    except Exception:
        return []
        
    soup = BeautifulSoup(resp.text, "html.parser")
    lineage = []
    
    infoboxes = soup.find_all("table", class_="infobox")
    for ib in infoboxes:
        sc_found = False
        for th in ib.find_all("th"):
            if "Scientific classification" in th.text:
                sc_found = True
                break
        
        if sc_found or "biota" in ib.get("class", []):
            in_sc = False
            for tr in ib.find_all("tr"):
                cells = tr.find_all(["th", "td"])
                
                # Check for section headers (colspan)
                if len(cells) == 1 and cells[0].get("colspan"):
                    if "Scientific classification" in cells[0].text:
                        in_sc = True
                        continue
                    elif in_sc:
                        break  # reached next section (e.g. Type species)

                if in_sc and len(cells) >= 2:
                    # Ranks are first column, taxon in second
                    rank_text = cells[0].text
                    td = cells[1]
                    taxon = ""
                    a_tag = td.find("a")
                    if a_tag:
                        taxon = clean_taxon_name(a_tag.text)
                    else:
                        taxon = clean_taxon_name(td.text)
                        
                    # Handle multiple words
                    if taxon:
                        taxon = taxon.split('\n')[0].strip()
                        first_word = taxon.split()[0]
                        lineage.append(first_word)
                        
                    if "Genus:" in rank_text:
                        break
            break
            
    return lineage

def build_nested_dict(lineage: List[str]) -> Dict[str, Any]:
    if not lineage:
        return {}
    
    # We output them as nested dictionaries. The lowest level is an empty dict or empty list depending on tree format.
    # The existing dino_dict uses empty list `[]` for clades without further taxa.
    root = {}
    current = root
    for taxon in lineage[:-1]:
        current[taxon] = {}
        current = current[taxon]
    
    if lineage:
        current[lineage[-1]] = []
        
    return root

def fetch_wiki_trees(genera_names: List[str], num_workers: int=10) -> Dict[str, Dict[str, Any]]:
    links = fetch_dinosaur_genera_links()
    results = {}
    
    def process_genus(name):
        target = name.lower()
        if target in links:
            url = links[target]
        else:
            url = f"https://en.wikipedia.org/wiki/{name.capitalize()}"
            
        lineage = fetch_genus_lineage(url)
        if lineage:
            return name, build_nested_dict(lineage)
        return name, None

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_genus, name): name for name in genera_names}
        for future in concurrent.futures.as_completed(futures):
            name, tree = future.result()
            if tree:
                results[name] = tree
                
    return results
