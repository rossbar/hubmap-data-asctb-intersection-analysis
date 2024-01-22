---
jupytext:
  notebook_metadata_filter: all
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
language_info:
  codemirror_mode:
    name: ipython
    version: 3
  file_extension: .py
  mimetype: text/x-python
  name: python
  nbconvert_exporter: python
  pygments_lexer: ipython3
  version: 3.11.6
---

# ASCT+B tables, OMAPs, and portal data

This notebook presents a preliminary investigation of the information in the
ASCT+B tables and OMAPs.

## Data annotation

We use the ASCT+B and OMAPs to inform cell-type annotations for the HubMAP data;
i.e. determining the cell-type categories and a preliminary mapping of protein
markers to cell types.

The overlap between protein markers in the tables and the markers present in the
actual multiplex imaging data impacts what we would expect to be able to
annotate.

### Example: Kidney data

The MIBI images of the kidney available on the hubmap data portal (not yet
published, still in QA) currently contain the following marker panel:

```{code-cell}
# Recorded manually since data is not publicly available
kidney_data_markers = {
    "AE1", "AQP1", "CALBINDIN", "CD31", "COLIV", "NAKATPASE", "PCAD", "S6",
    "SMA", "VIMENTIN", "DAPI_S001", "DAPI_S009"
}
print(f"Number of markers present in data: {len(kidney_data_markers)}")
```


```{code-cell}
:tags: [remove-cell]

from pprint import pprint
from pathlib import Path
from collections import Counter
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import json
```

Load data from the ASCT+B exporter into a graph data structure...

```{code-cell}
fname = Path.home() / "Downloads/asct+b_graph_data_selected_organs_2024.01.20_03.14.json"
with open(fname, "r") as fh:
    data = json.loads(fh.read())

G = nx.DiGraph()
G.add_nodes_from((node["id"], node) for node in data["nodes"])
G.add_edges_from((e["source"], e["target"]) for e in data["edges"])
```

... and extract the markers found in the table:

```{code-cell}
asctb_table_markers = {
    G.nodes[n]["name"] for n in G if G.nodes[n]["type"] == "BM"
}
print(f"Number of markers present in table: {len(asctb_table_markers)}")
```

```{code-cell}
:tags: [remove-cell]

from myst_nb import glue
glue("data_num_markers", len(kidney_data_markers), display=False)
glue("asctb_num_markers", len(asctb_table_markers), display=False)
```

```{code-cell}
:tags: [hide-output]

np.array(sorted(asctb_table_markers))  # Use array for shorter repr
```

There are {glue:text}`asctb_num_markers` total markers in the ASCT+B table,
and only {glue:text}`data_num_markers` in the actual multiplexed dataset.
We'd like to map the dataset markers to those in the ASCT+B table, so that we
can trace them back to the expected cell types (and anatomical structures).

Unfortunately, the marker names in the datasets are not standardized.
In practice, manually mapping marker names from the data to those in the table
is the only option.

```{code-cell}
:tags: [hide-cell]

# Markers with metadata - listed here to assist with manual mapping of marker
# names

sorted(
    (G.nodes[n]["name"], G.nodes[n]["metadata"]["label"]) for n in G
    if G.nodes[n]["type"] == "BM"
)
```

Here's the mapping I was able to come up with (I'm no biologist) with the help
of wikipedia + chatgpt:

```{code-cell}
data_to_asct_marker_mapping = {
    "AE1": "None",  # No candidate found
    "AQP1": "AQP1",
    "CALBINDIN": "CALB1",
    "CD31":"None",  # CD31 is not in the exported ASCTB table (see below)
    "COLIV": 'COL4A1',
    "NAKATPASE": 'ATP1A1',
    "PCAD": 'CDH1',
    "S6": "None",  # No candidate found in the table
    "SMA": "ACTA2",  # The "SMA" in the table is "aortic smooth muscle
    "VIMENTIN": "VIM"
}
```

So, in summary: there are `{glue:text}`{data_num_markers} present int he dataset,
of which 2 are DAPI stains which are not relevant for cell-type determination.
Of the remaining markers, 2 (`AE1` and `S6`) don't have obvious candidates from
the table (feedback welcome!)

```{note}
CD31 is not present in the data from the ASCT+B exporter, which is odd. It *is*
present in the published ASCT+B table. The discrepancy between these two data
sources is explored futher below.
```

This leaves us with `7` channels of protein markers that may be relevant for
cell-type determination.
Even in the case where we assume each marker provides sufficient information to
correctly identify all associated cell types, this leaves annotators with only
a subset of all possible cell-types:

```{code-cell}
name_to_id = {G.nodes[n]["name"]: G.nodes[n]["id"] for n in G}
markers_in_asctb = set(data_to_asct_marker_mapping.values()) & \
                   {name for _, name in G.nodes(data="name")}
all_ancestors = set(
    itertools.chain.from_iterable(
        nx.ancestors(G, source=name_to_id[marker]) for marker in markers_in_asctb
    )
)

detectable = nx.subgraph(G, nbunch=all_ancestors)
```

```{code-cell}
table_summary = Counter([ntype for _, ntype in G.nodes(data="type")])
subset_summary = Counter([ntype for _, ntype in detectable.nodes(data="type")])
print(f"Total number of celltypes listed in ASCT+B table for kidney: {table_summary['CT']}")
print(f"Maximum detectable number of celltypes with marker panel: {subset_summary['CT']}")
sorted([ndata["name"] for _, ndata in detectable.nodes(data=True) if ndata["type"] == "CT"])
```

## Difference between ASCT+B reporter and master data tables

The ASCT+B information differs significantly between the reporter and the "master"
data tables.
I assume that the reporter draws from all of the tables combined, so it's not
necessarily surprising that it has more information (i.e. more AS, CT, and BM
nodes) than the master tables.
However, though the master table contains fewer data, it also includes some
which are absent from the summary provided by the reporter.
For example, for the Kidney ASCT+B table, `CD31` appears as a biomarker in the
master table but is not found in the data from the reporter.

Since the master tables differ significantly from the reporter, it's worthwhile
to repeat the above analysis (i.e. determining cell-type subsets that
correspond to available markers) using the master table as the ASCT+B reference.

Start by loading the data...

```{code-cell}
:tags: [hide-output]

fname = Path.home() / "Downloads/Kidney_v1.4 - Kidney_v1.4.csv"
data = pd.read_csv(fname, skiprows=10)
list(data.columns)  # All columns contained in the master table
```

... and extracting all of the protein markers from the table

```{code-cell}
all_markers = set()
for i in range(1, 5):
    all_markers |= set(data[f"BProtein/{i}"].dropna())
```

We can then visually inspect the subsets of cells that correspond to each
marker included in the master table:

```{code-cell}
:tags: [hide-output]

for m in all_markers:
    print(f"==================={m}================")
    print(
        data[["CT/1", "CT/2"] + \
        [f"BProtein/{i}" for i in range(1, 5)]][data["BProtein/1"] == m]
    )
    print("\n\n")
```

Similar to the reporter data, there are roughly 60-70 cell types in the master
table:

```{code-cell}
print(f"Number of unique primary cell types: {len(np.unique(data['CT/1']))}")
```

We can perform the same analysis on these data as was done on the reporter
data to get a rough estimate of the greatest number of cell types we'd expect
to be able to distinguish given the marker panel of the portal data.

```{code-cell}
# Create graph from tabular data
# NOTE: only focusing on primary protein markers here
G = nx.DiGraph()
for m in all_markers:
    celltypes = data["CT/1"][data["BProtein/1"] == m].dropna()
    G.add_node(m, layer=1)
    for ct in celltypes:
        G.add_node(ct, layer=0)
        G.add_edge(ct, m)

# Another manual mapping of non-standard marker names from the data to those
# in the ASCT+B tables
data_to_asct_marker_mapping = {
    "AE1": "None",
    "AQP1": "AQP1",
    "CALBINDIN": "CALB1",
    "CD31":"CD31",
    "COLIV": 'COL4A1',
    "NAKATPASE": 'ATP1A1',
    "PCAD": 'CDH1',
    "S6": "None",
    "SMA": "SMA",  # Is this SMA different than the reporter SMA?
    "VIMENTIN": "VIM"
}

markers_in_asctb = all_markers & set(data_to_asct_marker_mapping.values())
all_ancestors = set(
    itertools.chain.from_iterable(nx.ancestors(G, source=m) for m in markers_in_asctb)
)
all_nodes = all_ancestors | markers_in_asctb

# Visualize subgraph comprising the data marker panel intersection with the
# ASCT+B master table markers
H = G.subgraph(all_nodes)
pos = nx.multipartite_layout(G, subset_key="layer")
fig, ax = plt.subplots(figsize=(12, 9));
nx.draw_networkx_edges(H, pos=pos, ax=ax, min_source_margin=100, min_target_margin=50);
nx.draw_networkx_labels(H, pos=pos, ax=ax);

# Prettify
fig.tight_layout()
ax.set_xlim(ax.get_xlim() * np.array([1.2, 1]));
```
