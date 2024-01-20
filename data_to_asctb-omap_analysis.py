# IPython log file

fname = "asct+b_graph_data_selected_organs_2024.01.20_11.27.json"
import json
with open(fname, "r") as fh:
    data = json.loads(fh)
    
with open(fname, "r") as fh:
    data = json.loads(fh.read())
    
data
data.keys()
data["nodes"]
import networkx as nx
data["nodes"][0]
celltypes = [node for data["nodes"] if node["type"] == "CT"]
celltypes = [node if node["type"] == "CT" for data["nodes"]]
celltypes = [node for node in data["nodes"] if node["type"] == "CT"]
celltypes
len(celltypes)
data["nodes"]
biomarkers = [node for node in data["nodes"] if node["id"] == "BM"]
biomarkers
biomarkers = [node for node in data["nodes"] if node["type"] == "BM"]
biomarkers
len(biomarkers)
assert {ct["id"] for ct in celltypes} | {bm["id"] for bm in biomarkers} == {}
assert {ct["id"] for ct in celltypes} & {bm["id"] for bm in biomarkers} == {}
{ct["id"] for ct in celltypes} | {bm["id"] for bm in biomarkers}
{ct["id"] for ct in celltypes} & {bm["id"] for bm in biomarkers}
assert {ct["id"] for ct in celltypes} & {bm["id"] for bm in biomarkers} == set()
data["edges"]
G = nx.DiGraph()
get_ipython().run_line_magic('pinfo', 'G.add_nodes_from')
celltypes
G = nx.DiGraph()
G.add_nodes_from(node["id"] for node in data["nodes"])
G.nodes()
get_ipython().run_line_magic('pinfo', 'nx.set_node_attributes')
nx.DiGraph(pe)
get_ipython().run_line_magic('pinfo', 'nx.set_node_attributes')
G = nx.DiGraph()
G.add_nodes_from(node["id"], **node for node in data["nodes"])
G.add_nodes_from(node["id"], {**node} for node in data["nodes"])
G.add_nodes_from(node["id"], **dict(node) for node in data["nodes"])
G.add_nodes_from(node["id"], dict(node) for node in data["nodes"])
get_ipython().run_line_magic('pinfo', 'G.add_nodes_from')
G.add_nodes_from((node["id"], node) for node in data["nodes"])
G.nodes
G.nodes(data=True)
G.nodes(data="id")
G.nodes(data="type")
data["edges"]
G.add_edges_from((e["source"], e["target"]) for e in data["edges"])
nx.draw(G, with_labels=True)
import matplotlib.pyplot as plt
plt.show()
for layer, nodes in enumerate(nx.topological_generations(G)):
    for node in nodes:
        G.nodes[node]["layer"] = layer
        
pos = nx.multipartite_layout(G, subset_key="layer")
nx.draw(G, pos=pos, with_labels=True)
plt.ion()
plt.show()
pos = {}
set(ntype for _, n in G.nodes(data="type"))
set(ntype for _, ntype in G.nodes(data="type"))
get_ipython().run_line_magic('pinfo', 'nx.subgraph')
nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == "AS"))
AS = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == "AS"))
BM = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == "BM"))
CT = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == "CT"))
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == ntype))
        for yval, node in enumerate(sg):
            pos[node] = np.array([xval, yval]) 
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == ntype))
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval])
        
pos
node
yval
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == ntype))
    print(sg)
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval])
        
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == ntype))
    print(xval, ntype, sg)
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval])
        
G = nx.DiGraph()
G.add_nodes_from((node["id"], node) for node in data["nodes"])
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data=type) if nt == ntype))
    print(xval, ntype, sg)
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval])
        
G.nodes(data="type")
[n for n, nt in G.nodes(data=type) if nt == "BM"]
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data="type") if nt == ntype))
    print(xval, ntype, sg)
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval])
        
import numpy as np
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data="type") if nt == ntype))
    print(xval, ntype, sg)
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval])
        
pos
nx.draw(G, pos=pos, with_labels=True)
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data="type") if nt == ntype))
    print(xval, ntype, sg)
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval*10])
        
pos = {}
for xval, ntype in enumerate(set(ntype for _, ntype in G.nodes(data="type"))):
    sg = nx.subgraph(G, (n for n, nt in G.nodes(data="type") if nt == ntype))
    for yval, node in enumerate(sg):
        pos[node] = np.array([xval, yval*10])
        
nx.draw(G, pos=pos, with_labels=True)
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
nx.draw(G, pos=pos, with_labels=True)
get_ipython().run_line_magic('pinfo', 'nx.draw_networkx_nodes')
nx.draw_networkx_nodes(G, pos=pos, node_color=G.nodes(data="color"))
nx.draw_networkx_nodes(G, pos=pos, node_color=G.nodes(data="color").values())
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values())
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values())
get_ipython().run_line_magic('pinfo', 'nx.multipartite_layout')
pos = nx.multipartite_layout(G, subset_key="layer", scale=100)
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values())
pos = nx.multipartite_layout(G, subset_key="layer", scale=1000)
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values())
pos
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values(), node_size=100)
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values(), node_size=10)
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values(), node_size=1)
pos = nx.multipartite_layout(G, subset_key="layer", scale=100)
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values(), node_size=1)
pos = nx.multipartite_layout(G, subset_key="layer")
nx.draw_networkx_nodes(G, pos=pos, node_color=dict(G.nodes(data="color")).values())
array([  6.6184508 , 868.02030457]) * [1, 2]
np.array([  6.6184508 , 868.02030457]) * [1, 2]
pos = nx.multipartite_layout(G, subset_key="layer")
biomarkers = nx.subgraph(G, (n for n, ntype in G.nodes(data="type") if ntype == "BM"))
nx.draw(G, pos=pos, with_labels=True)
nx.draw(biomarkers, pos=pos, with_labels=True)
[G.nodes[n]["name"] for G.nodes[n]]
[G.nodes[n]["name"] for n in G]
all_marker_names = sorted(G.nodes[n]["name"] for n in G if G.nodes[n]["type"] == "BM")
all_marker_names
all_marker_names = {G.nodes[n]["name"] for n in G if G.nodes[n]["type"] == "BM"}
all_marker_names
kidney_dataset_markers = {"AE1", "AQP1", "CALBINDIN", "CD31", "COLIV", "NAKATPASE", "PCAD", "S6", "SMA", "VIMENTIN"}
all_marker_names | kidney_dataset_markers
all_marker_names & kidney_dataset_markers
all_marker_names
G.nodes(data=True)
all_marker_names = {(G.nodes[n]["name"], G.nodes[n]["metadata"]["label"]) for n in G if G.nodes[n]["type"] == "BM"}
all_markers = {(G.nodes[n]["name"], G.nodes[n]["metadata"]["label"]) for n in G if G.nodes[n]["type"] == "BM"}
sorted(all_markers)
sorted((name, descr) for name, descr in all_names if "ase" in descr)
sorted((name, descr) for name, descr in all_markers if "ase" in descr)
sorted(all_markers)
sorted((name, descr) for name, descr in all_markers if "adherin" in descr)
sorted((name, descr) for name, descr in all_markers if "cad" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "SOME" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "some" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "somal" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "ribo" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "s6" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "rp" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "6" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "sma" in descr.lower())
sorted((name, descr) for name, descr in all_markers if "smooth" in descr.lower())
sorted(all_markers)
data_to_asct_marker_mapping = {"AE1": "AE1", "AQP1": "AQP1", "CALBINDIN": "CALB1", "CD31":"CD31", "COLIV": 'COL4A1', "NAKATPASE": 'ATP1A1', "PCAD": 'CDH1', "S6": "RPS6", "SMA": "ACTA2", "VIMENTIN": "VIM"}
set(data_to_asct_marker_mapping.values())
get_ipython().run_line_magic('pinfo', 'nx.ancestors')
all_ancestors = nx.ancestors(G, source=marker) for marker in set(data_to_asct_marker_mapping.values())
all_ancestors = (nx.ancestors(G, source=marker) for marker in set(data_to_asct_marker_mapping.values()))
all_ancestors
detectable = nx.subgraph(G, nbunch=all_ancestors)
set(data_to_asct_marker_mapping.values())
set(data_to_asct_marker_mapping.values()) & G.nodes
G.nodes
detectable = set(data_to_asct_marker_mapping.values()) & {name for _, name in G.nodes(data="name")}
detectable
markers_in_asctb = set(data_to_asct_marker_mapping.values()) & {name for _, name in G.nodes(data="name")}
all_ancestors = (nx.ancestors(G, source=marker) for marker in markers_in_asctb)
all_ancestors
list(all_ancestors)
sorted(all_markers)
markers_in_asctb = set(data_to_asct_marker_mapping.values()) & {name for _, name in G.nodes(data="name")}
markers_in_asctb
name_to_id = {G.nodes[n]["name"]: G.nodes[n]["id"] for n in G}
all_ancestors = (nx.ancestors(G, source=name_to_id[marker]) for marker in markers_in_asctb)
list(all_ancestors)
all_ancestors = (nx.ancestors(G, source=name_to_id[marker]) for marker in markers_in_asctb)
import itertools
all_ancestors = itertools.chain.from_iterable(nx.ancestors(G, source=name_to_id[marker]) for marker in markers_in_asctb)
all_ancestors
list(all_ancestors)
all_ancestors = set(itertools.chain.from_iterable(nx.ancestors(G, source=name_to_id[marker]) for marker in markers_in_asctb))
all_ancestors
detectable = nx.subgraph(G, nbunch=all_ancestors)
nx.draw(detectable, pos=pos, with_labels=True)
detectable = nx.subgraph(G, nbunch=all_ancestors | {name_to_id[marker] for marker in markers_in_asctb})
nx.draw(detectable, pos=pos, with_labels=True)
from collections import Counter
get_ipython().run_line_magic('pinfo', 'Counter')
Counter([ntype for _, ntype in G.nodes(data="type")])
Counter([ntype for _, ntype in detectable.nodes(data="type")])
[ndata["name"] for _, ndata in G.nodes(data=True) if data["type"] == "CT"]
[ndata["name"] for _, ndata in G.nodes(data=True) if ndata["type"] == "CT"]
len([ndata["name"] for _, ndata in G.nodes(data=True) if ndata["type"] == "CT"])
len([ndata["name"] for _, ndata in detectable.nodes(data=True) if ndata["type"] == "CT"])
[ndata["name"] for _, ndata in detectable.nodes(data=True) if ndata["type"] == "CT"]
get_ipython().run_line_magic('pinfo', 'nx.relabel_nodes')
name_to_id
id_to_name = {v: k for k, v in name_to_id.items()}
id_to_name
get_ipython().run_line_magic('pinfo', 'nx.relabel_nodes')
G = nx.relabel_nodes(G, id_to_name)
G.nodes
sorted(G)
G
G = nx.relabel_nodes(G, name_to_id)
len(G)
name_to_id
set(name_to_id) - set(G)
set(name_to_id.values()) - set(G)
set(type(k) for k in name_to_id)
name_to_id
id_to_name = {v: k for k, v in name_to_id.items()}
set(type(k) for k in id_to_name)
G = nx.relabel_nodes(G, id_to_name)
set(type(n) for n in G)
G.nodes
id_to_name[23]
name_to_id["23"]
G = nx.relabel_nodes(G, name_to_id)
G.nodes
name_to_id["Fibroblast"]
id_to_name = {str(v): k for k, v in name_to_id.items()}
G = nx.relabel_nodes(G, id_to_name)
sorted(G)
G = nx.relabel_nodes(G, name_to_ide)
G = nx.relabel_nodes(G, name_to_id)
sorted(G)
G = nx.relabel_nodes(G, id_to_name)
id_to_name
id_to_name = {v: k for k, v in name_to_id.items()}
id_to_name
G.nodes
G = nx.relabel_nodes(G, name_to_ide)
G = nx.relabel_nodes(G, id_to_name)
G.nodes
for n in G:
    G.nodes[n] = str(n)
    
for n in G:
    G[n] = str(n)
    
G.nodes
n_to_str = {n: str(n) for n in G}
H = nx.relabel_nodes(G, n_to_str)
H
H.nodes
sorted(H)
H["Fibroblast"]
[n for _, ntype in detectable.nodes(data="type") if ntype == "CT"]
G
G.nodes
G = nx.DiGraph()
G.add_nodes_from((node["id"], node) for node in data["nodes"])
G.add_edges_from((e["source"], e["target"]) for e in data["edges"])
G.nodes
set(G) - set(itertools.chain.from_iterable(G.edges))
set(itertools.chain.from_iterable(G.edges)) - set(G)
G
G[id_to_name["Fibroblast"]]
G[name_to_id["Fibroblast"]]
set(G.nodes[n]["name"] for n in G[name_to_id["Fibroblast"]])
[ndata["name"] for _, ndata in detectable.nodes(data=True) if ndata["type"] == "CT"]
fibroblast_markers = set(G.nodes[n]["name"] for n in G[name_to_id["Fibroblast"]])
med_fibro_markers = set(G.nodes[n]["name"] for n in G[name_to_id["Medullary Fibroblast"]])
fibroblast_markers
med_fibro_markers
data["nodes"]
set(n["node"]["name"] for n in data["nodes"])
set(n["name"] for n in data["nodes"])
get_ipython().run_line_magic('logstart', 'data_to_asctb-omap_analysis.py')
exit()
