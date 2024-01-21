from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import itertools
import numpy as np
import pandas as pd
from pprint import pprint
plt.ion()

fname = Path.home() / "Downloads/Kidney_v1.4 - Kidney_v1.4.csv"
data = pd.read_csv(fname, skiprows=10)
data[["CT/1", "CT/2"]]
pprint(list(data.columns))

all_markers = set()
for i in range(1, 5):
    all_markers |= set(data[f"BProtein/{i}"].dropna())
    
for m in all_markers:
    print(f"==================={m}================")
    print(data[["CT/1", "CT/2"] + [f"BProtein/{i}" for i in range(1, 5)]][data["BProtein/1"] == m])
    print("\n\n\n")
    
G = nx.DiGraph()
for m in all_markers:
    celltypes = data["CT/1"][data["BProtein/1"] == m].dropna()
    G.add_node(m, layer=1, color="tab:green")
    for ct in celltypes:
        G.add_node(ct, layer=0, color="tab:blue")
        G.add_edge(ct, m)
 
pos = nx.multipartite_layout(G, subset_key="layer")
nx.draw_networkx_labels(G, pos=pos);
nx.draw_networkx_edges(G, pos=pos, min_source_margin=100, min_target_margin=50);

# Manual
data_to_asct_marker_mapping = {"AE1": "AE1", "AQP1": "AQP1", "CALBINDIN": "CALB1", "CD31":"CD31", "COLIV": 'COL4A1', "NAKATPASE": 'ATP1A1', "PCAD": 'CDH1', "S6": "RPS6", "SMA": "SMA", "VIMENTIN": "VIM"}

markers_in_asctb = all_markers & set(data_to_asct_marker_mapping.values())
all_ancestors = set(itertools.chain.from_iterable(nx.ancestors(G, source=m) for m in markers_in_asctb))
all_nodes = all_ancestors | markers_in_asctb
plt.figure()
nx.draw_networkx_labels(G.subgraph(all_nodes), pos=pos);
nx.draw_networkx_edges(G.subgraph(all_nodes), pos=pos, min_source_margin=100, min_target_margin=50);
