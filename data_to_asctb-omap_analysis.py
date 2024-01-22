from pprint import pprint
from pathlib import Path
from collections import Counter
import itertools
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import json

fname = Path.home() / "Downloads/asct+b_graph_data_selected_organs_2024.01.20_03.14.json"
with open(fname, "r") as fh:
    data = json.loads(fh.read())
    
G = nx.DiGraph()
G.add_nodes_from((node["id"], node) for node in data["nodes"])
G.add_edges_from((e["source"], e["target"]) for e in data["edges"])

for n, ntype in G.nodes(data="type"):
    if ntype == "AS":
        G.nodes[n]["layer"] = 0
        G.nodes[n]["color"] = "tab:red"
    elif ntype == "CT":
        G.nodes[n]["layer"] = 1
        G.nodes[n]["color"] = "tab:blue"
    elif ntype == "BM":
        G.nodes[n]["layer"] = 2
        G.nodes[n]["color"] = "tab:green"
    else:
        raise ValueError(f"Unrecognized node type {ntype} for node {n}")
        
pos = nx.multipartite_layout(G, subset_key="layer")
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values())
all_markers = {(G.nodes[n]["name"], G.nodes[n]["metadata"]["label"]) for n in G if G.nodes[n]["type"] == "BM"}
sorted(all_markers)

# Manual mapping
data_to_asct_marker_mapping = {"AE1": "AE1", "AQP1": "AQP1", "CALBINDIN": "CALB1", "CD31":"CD31", "COLIV": 'COL4A1', "NAKATPASE": 'ATP1A1', "PCAD": 'CDH1', "S6": "RPS6", "SMA": "ACTA2", "VIMENTIN": "VIM"}
name_to_id = {G.nodes[n]["name"]: G.nodes[n]["id"] for n in G}
markers_in_asctb = set(data_to_asct_marker_mapping.values()) & {name for _, name in G.nodes(data="name")}
all_ancestors = set(
    itertools.chain.from_iterable(  # Flatten
        nx.ancestors(G, source=name_to_id[marker]) for marker in markers_in_asctb
    )
)

detectable = nx.subgraph(G, nbunch=all_ancestors | {name_to_id[marker] for marker in markers_in_asctb})
nx.draw(detectable, pos=pos, with_labels=True)

print(Counter([ntype for _, ntype in G.nodes(data="type")]))
print(Counter([ntype for _, ntype in detectable.nodes(data="type")]))

pprint([ndata["name"] for _, ndata in detectable.nodes(data=True) if ndata["type"] == "CT"])

### 

fibroblast_markers = set(G.nodes[n]["name"] for n in G[name_to_id["Fibroblast"]])
med_fibro_markers = set(G.nodes[n]["name"] for n in G[name_to_id["Medullary Fibroblast"]])
