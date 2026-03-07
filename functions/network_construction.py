import os
import pandas as pd
import networkx as nx

L = 566
ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")

def build_sequence_graph(seqs, alphabet, L):
    """
    Build graph where nodes are sequences and edges connect
    sequences differing by one mutation.
    """
    seq_set = set(seqs)
    G = nx.Graph()
    G.add_nodes_from(seq_set)

    for seq in seq_set:
        for pos in range(L):
            original = seq[pos]
            for aa in alphabet:
                if aa == original:
                    continue
                neighbor = seq[:pos] + aa + seq[pos + 1 :]
                if neighbor in seq_set:
                    G.add_edge(seq, neighbor)

    return G


def build_graphs_from_neighbourhood(
    df_long,
    strain_dirs,
    seed_dir,
    alphabet=ALPHABET,
    L=L,
):
    """
    Construct neutral mutation graphs from neighbourhood enumeration output.

    Inputs
    ------
    df_long : DataFrame
        Must contain columns:
            - strain
            - sequence (run_id)
            - position (0-based)
            - neutral_neighbours
    strain_dirs : list[str]
        e.g. ["AAD17229_20", "ACP41934_20", ...]
    seed_dir : str
        Directory containing {strain}_20.txt seed files

    Returns
    -------
    seq_dfs : dict[strain -> DataFrame]
        Columns:
            - strain_id
            - run_id
            - seed_seq
            - mutant_seq
    graphs : dict[strain -> nx.Graph]
    """

    seq_dfs = {}
    graphs = {}

    for strain_dir in strain_dirs:
        strain = strain_dir.split("_")[0]
        print(f"\nProcessing strain {strain}...")

        # --- Load seed sequences ---
        seed_path = os.path.join(seed_dir, f"{strain_dir}.txt")
        if not os.path.exists(seed_path):
            print(f"Seed file not found: {seed_path}")
            continue

        with open(seed_path) as f:
            seed_sequences = [line.strip() for line in f if line.strip()]

        print(f"Loaded {len(seed_sequences)} seed sequences")

        # --- Subset neighbourhood data ---
        strain_df = df_long[df_long["strain"] == strain]

        records = []

        for run_id, run_df in strain_df.groupby("sequence"):
            run_id = int(run_id)

            if run_id >= len(seed_sequences):
                print(f"Run {run_id} exceeds seed count ({len(seed_sequences)})")
                continue

            seed = seed_sequences[run_id]

            for _, row in run_df.iterrows():
                pos = int(row["position"])
                neighbours = row["neutral_neighbours"]

                # Always include seed itself
                if pd.isna(neighbours) or not str(neighbours).strip():
                    records.append({
                        "strain_id": strain,
                        "run_id": run_id,
                        "seed_seq": seed,
                        "mutant_seq": seed,
                    })
                    continue

                for aa in str(neighbours).split(","):
                    aa = aa.strip()
                    if not aa:
                        continue

                    new_seq = list(seed)
                    new_seq[pos] = aa

                    records.append({
                        "strain_id": strain,
                        "run_id": run_id,
                        "seed_seq": seed,
                        "mutant_seq": "".join(new_seq),
                    })

        seq_df = pd.DataFrame(records)
        print(f"Generated {len(seq_df)} sequence records")

        # --- Build graph ---
        seqs = set(seq_df["seed_seq"]) | set(seq_df["mutant_seq"])
        G = build_sequence_graph(seqs, alphabet, L)

        seq_dfs[strain] = seq_df
        graphs[strain] = G

        print(
            f"Graph for {strain}: "
            f"{G.number_of_nodes()} nodes, "
            f"{G.number_of_edges()} edges"
        )

    return seq_dfs, graphs
