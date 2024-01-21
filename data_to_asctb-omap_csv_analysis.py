# IPython log file

import pandas as pd
fname = "Kidney_v1.4 - Kidney_v1.4.csv"
pd.read_csv(fname)
data = pd.read_csv(fname, skiprows=10)
data
data["CT/1"]
data[("CT/1", "CT/2")]
data["CT/1", "CT/2"]
data[["CT/1", "CT/2"]]
data["columns"]
data.columns
print(data.columns)
list(data.columns)
data[["CT/1/ID", "CT/2"]]
data.head[["CT/1", "CT/2"]]
data[["CT/1", "CT/2"]]
data[["CT/1", "CT/2", "BProtein/1"]]
set(data["BProtein/1"])
data.get("BProtein/1")
all_makers = set()
for i in range(1, 5):
    all_markers += data[f"BProtein/{i}"]
    
all_markers = set()
for i in range(1, 5):
    all_markers += data[f"BProtein/{i}"]
    
for i in range(1, 5):
    all_markers += set(data[f"BProtein/{i}"])
    
for i in range(1, 5):
    all_markers |= set(data[f"BProtein/{i}"])
all_markers
all_markers = set()
for i in range(1, 5):
    all_markers |= set(data[f"BProtein/{i}"])
    
data[["CT/1", "CT/2"]][data["BProtein/1"] == "AQP1"]
data[["CT/1", "CT/2", "BProtein/2", "BProtein3", "BProtein4"]][data["BProtein/1"] == "AQP1"]
data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == "AQP1"]
data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == "AQP2"]
all_markers
data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == "ATP1A1"]
for m in all markers:
    print(data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == m])
for m in all_markers:
    print(data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == m])
    
for m in all_markers:
    print(f"==================={m}================")
    print(data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == m])
    print("\n\n\n")
    
np.unique(data["CT/1"])
import numpy as np
np.unique(data["CT/1"])
len(np.unique(data["CT/1"]))
G = nx.DiGraph()
import networkx as nx
G = nx.DiGraph()
for m in all_markers:
    celltypes = data["CT/1"][data["BProtein/1"] == m]
    G.add_node(m, layer=1, color="tab:green")
    for ct in celltypes:
        G.add_node(ct, layer=0, color="tab:blue")
        G.add_edge(ct, m)
        
pos = nx.multipartite_layout(G, subset_key="layer")
nx.draw(G, pos=pos, with_labels=True)
plt.ion()
import matplotlib.pyplot as plt
plt.ion()
plt.show()
all_markers
all_markers - {nan}
all_markers - {float.nan}
all_markers - {float("nan")}
pd.na
pd.nan
list(all_markers)
list(all_markers)[10]
type(list(all_markers)[10])
float("nan")
all_markers - {float("nan")}
list(all_markers).pop(float("nan"))
list(all_markers).remove(float("nan"))
list(all_markers)
list(all_markers).remove(float("nan"))
for m in all_markers:
    celltypes = pd.notna(data["CT/1"][data["BProtein/1"] == m])
    G.add_node(m, layer=1, color="tab:green")
    for ct in celltypes:
        G.add_node(ct, layer=0, color="tab:blue")
        G.add_edge(ct, m)
        
G = nx.DiGraph()
G = nx.DiGraph()
for m in all_markers:
    celltypes = pd.notna(data["CT/1"][data["BProtein/1"] == m])
    G.add_node(m, layer=1, color="tab:green")
    for ct in celltypes:
        G.add_node(ct, layer=0, color="tab:blue")
        G.add_edge(ct, m)
        
nx.draw(G, pos=pos, with_labels=True, node_colors=dict(G.nodes(data="color")).values())
nx.draw(G, pos=pos, with_labels=True, node_color=dict(G.nodes(data="color")).values())
G = nx.DiGraph()
for m in all_markers:
    celltypes = data["CT/1"][data["BProtein/1"] == m].dropna()
    G.add_node(m, layer=1, color="tab:green")
    for ct in celltypes:
        G.add_node(ct, layer=0, color="tab:blue")
        G.add_edge(ct, m)
        
nx.draw(G, pos=pos, with_labels=True, node_color=dict(G.nodes(data="color")).values())
nx.draw(G, pos=pos, with_labels=True, node_color=dict(G.nodes(data="color")).values())
nx.draw_networkx_labels(G, pos=pos)
nx.draw_networkx_edges(G, pos=pos, min_source_margin=0.15, min_target_margin=0.15)
get_ipython().run_line_magic('pinfo', 'nx.draw_networkx_edges')
nx.draw_networkx_labels(G, pos=pos)
nx.draw_networkx_edges(G, pos=pos, min_source_margin=0.15, min_target_margin=100);
nx.draw_networkx_labels(G, pos=pos);
nx.draw_networkx_edges(G, pos=pos, min_source_margin=10, min_target_margin=10);
nx.draw_networkx_labels(G, pos=pos);
nx.draw_networkx_edges(G, pos=pos, min_source_margin=100, min_target_margin=50);
get_ipython().run_line_magic('pinfo', 'nx.subgraph')
all_markers
data_to_asct_marker_mapping = {"AE1": "AE1", "AQP1": "AQP1", "CALBINDIN": "CALB1", "CD31":"CD31", "COLIV": 'COL4A1', "NAKATPASE": 'ATP1A1', "PCAD": 'CDH1', "S6": "RPS6", "SMA": "SMA", "VIMENTIN": "VIM"}
all_markers & set(data_to_asct_marker_mapping.values())
marker_union = all_markers & set(data_to_asct_marker_mapping.values())
all_ancestors = set(itertools.chain.from_iterable(nx.ancestors(G, source=m) for m in markers_in_asctb))
import itertools
markers_in_asctb = all_markers & set(data_to_asct_marker_mapping.values())
all_ancestors = set(itertools.chain.from_iterable(nx.ancestors(G, source=m) for m in markers_in_asctb))
list(all_ancestors)
len(all_ancestors)
nx.draw_networkx_edges(G.subgraph(all_ancestors), pos=pos, min_source_margin=100, min_target_margin=50);
nx.draw_networkx_labels(G.subgraph(all_ancestors), pos=pos);
nx.draw_networkx_labels(G.subgraph(all_ancestors), pos=pos);
nx.draw_networkx_edges(G.subgraph(all_ancestors), pos=pos, min_source_margin=100, min_target_margin=50);
nx.draw_networkx_labels(G.subgraph(all_ancestors), pos=pos);
nx.draw_networkx_edges(G.subgraph(all_ancestors), pos=pos, min_source_margin=100, min_target_margin=50);
nx.draw_networkx_edges(G.subgraph(all_ancestors), pos=pos, min_source_margin=100, min_target_margin=50);
nx.draw_networkx_labels(G.subgraph(all_ancestors), pos=pos);
G.sugraph(all_ancestors).nodes
G.subgraph(all_ancestors).nodes
all_ancestors = set(itertools.chain.from_iterable(nx.ancestors(G, source=m) for m in markers_in_asctb))
all_nodes = all_ancestors | markers_in_asctb
all_nodes
nx.draw_networkx_edges(G.subgraph(all_ancestors), pos=pos, min_source_margin=100, min_target_margin=50);
nx.draw_networkx_labels(G.subgraph(all_ancestors), pos=pos);
nx.draw_networkx_labels(G.subgraph(all_nodes), pos=pos);
nx.draw_networkx_edges(G.subgraph(all_nodes), pos=pos, min_source_margin=100, min_target_margin=50);
for m in all_markers:
    print(f"==================={m}================")
    print(data[["CT/1", "CT/2", "BProtein/2", "BProtein/3", "BProtein/4"]][data["BProtein/1"] == m])
    print("\n\n\n")
    
get_ipython().run_line_magic('logstart', 'data_to_asctb-omap_csv_analysis.py')
exit()
