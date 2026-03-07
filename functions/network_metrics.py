import numpy as np
import networkx as nx
from collections import Counter

def component_number(G):
    """
    Return:
    - number of connected components
    - size of largest component
    - distribution of component sizes (size -> count)
    """
    components = list(nx.connected_components(G))
    size_dist = Counter(len(c) for c in components)

    num_components = len(components)
    largest_size = max(size_dist.keys(), default=0)

    return num_components, largest_size, dict(size_dist)

def largest_components(G, k=2):
    """Return up to k largest connected component subgraphs"""
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    subgraphs = []
    for nodes in components[:k]:
        subgraphs.append(G.subgraph(nodes).copy())
    return subgraphs

def shannon_entropy(counter):
    total = sum(counter.values())
    if total == 0:
        return 0.0
    probs = np.array(list(counter.values())) / total
    probs = probs[probs > 0]
    return float(-(probs * np.log2(probs)).sum())

def entropy_profile(G):
    """Return Shannon entropy per sequence position for graph G"""
    if G.number_of_nodes() == 0:
        return np.array([])

    seq_len = len(next(iter(G.nodes())))
    ent = np.zeros(seq_len)

    for p in range(seq_len):
        ent[p] = shannon_entropy(
            Counter(seq[p] for seq in G.nodes())
        )

    return ent