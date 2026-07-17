"""Covering algorithm for frame clustering from an RMSD all-vs-all matrix.

This module is intentionally side-effect free:
- it does not write PNG/CSV outputs on import or function calls,
- it builds clean in-memory objects you can inspect or plot later.

Core logic implemented here:
1. Load RMSD matrix and relabel rows/columns to integer frame IDs (1..N).
2. Build classes using your rule:
    - Start with representative frame 1.
    - Class_k contains rows with RMSD <= threshold to representative frame.
    - A future representative cannot be equal to any row number already in prior classes.
3. Build a true RMSD threshold graph (edges from matrix values <= threshold).
4. Connect representative nodes to each other to keep class progression visible.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


DEFAULT_INPUT = Path(
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\csv\rmsd_matrix_v2_all_vs_all_biopython.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k"
)


def load_frames_matrix(csv_path: Path, max_frames: int | None = None) -> pd.DataFrame:
    """Load RMSD matrix and relabel both axes to integer frame IDs."""
    matrix = pd.read_csv(csv_path, index_col=0)
    n_cols = matrix.shape[1]
    matrix.columns = list(range(1, n_cols + 1))
    matrix.index = list(range(1, matrix.shape[0] + 1))

    if max_frames is not None:
        keep = list(range(1, min(max_frames, n_cols) + 1))
        matrix = matrix.loc[keep, keep]

    return matrix


def build_classes_recursive(
    frames: pd.DataFrame,
    threshold: float,
    remaining_rows: Set[int],
    blocked_representatives: Set[int],
    class_id: int,
    representative: int,
) -> List[Dict[str, Any]]:
    """Build covering classes while honoring representative-blocking rule.

    `blocked_representatives` contains row IDs that appeared in earlier classes.
    A next representative cannot be in that set.
    """
    if not remaining_rows:
        return []

    if representative not in remaining_rows or representative in blocked_representatives:
        available = sorted(remaining_rows - blocked_representatives)
        if not available:
            # Fallback: if all remaining rows are blocked, continue with smallest remaining.
            available = sorted(remaining_rows)
        representative = available[0]

    current_rows = sorted(remaining_rows)
    eligible_mask = frames.loc[current_rows, representative] <= threshold
    class_members = [row for row, keep in zip(current_rows, eligible_mask) if bool(keep)]

    # Diagonal should already include representative, but keep a safe guard.
    if representative not in class_members:
        class_members.append(representative)
        class_members = sorted(set(class_members))

    class_info: Dict[str, Any] = {
        "class_name": f"Class_{class_id}",
        "representative": representative,
        "members": class_members,
        "member_count": len(class_members),
    }

    new_remaining = set(remaining_rows) - set(class_members)
    new_blocked = set(blocked_representatives) | set(class_members)
    if not new_remaining:
        return [class_info]

    candidate_next = sorted(new_remaining - new_blocked)
    if candidate_next:
        next_representative = candidate_next[0]
    else:
        next_representative = min(new_remaining)

    return [class_info] + build_classes_recursive(
        frames=frames,
        threshold=threshold,
        remaining_rows=new_remaining,
        blocked_representatives=new_blocked,
        class_id=class_id + 1,
        representative=next_representative,
    )


def classes_to_table(classes: List[Dict[str, Any]]) -> pd.DataFrame:
    """Return class summary as a dataframe (for display/export elsewhere)."""
    rows = []
    for info in classes:
        members = list(info["members"])
        members_text = ", ".join(str(x) for x in members)
        rows.append(
            {
                "Class": info["class_name"],
                "Representative Frame": info["representative"],
                "N Members": info["member_count"],
                "Members": members_text,
            }
        )

    return pd.DataFrame(rows)


def build_class_edge_map(
    frames: pd.DataFrame,
    classes: List[Dict[str, Any]],
    threshold: float,
) -> Dict[str, List[Tuple[int, int, float]]]:
    """Build explicit representative-to-member edges for each class.

    These edges are the class-local cluster connections you asked for, and each
    edge stores the angstrom distance as its label/weight.
    """
    class_edges: Dict[str, List[Tuple[int, int, float]]] = {}
    for info in classes:
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        members = [int(x) for x in info["members"]]
        edges: List[Tuple[int, int, float]] = []

        for node in members:
            if node == rep:
                continue
            rmsd = float(frames.loc[rep, node])
            if rmsd <= threshold:
                edges.append((rep, node, rmsd))

        class_edges[class_name] = edges

    return class_edges


def build_display_metadata(
    frames: pd.DataFrame,
    classes: List[Dict[str, Any]],
    class_edges: Dict[str, List[Tuple[int, int, float]]],
) -> Dict[str, Any]:
    """Build node-size, label, and edge-label metadata for later plotting.

    - Representative nodes get larger sizes.
    - All cluster edges carry angstrom labels.
    - Representative-to-representative links are also included.
    """
    node_sizes: Dict[int, int] = {int(node): 90 for node in frames.index.tolist()}
    node_labels: Dict[int, str] = {int(node): f"F{int(node)}" for node in frames.index.tolist()}
    edge_labels: Dict[Tuple[int, int], str] = {}

    representative_nodes: List[int] = []
    for info in classes:
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        representative_nodes.append(rep)

        node_sizes[rep] = 420
        node_labels[rep] = f"F{rep}\n{class_name}"

        for u, v, w in class_edges[class_name]:
            edge_labels[(u, v)] = f"{w:.2f} Å"

    representative_links: List[Tuple[int, int, str]] = []
    if len(representative_nodes) > 1:
        for a, b in zip(representative_nodes[:-1], representative_nodes[1:]):
            representative_links.append((a, b, "rep_link"))

    return {
        "node_sizes": node_sizes,
        "node_labels": node_labels,
        "edge_labels": edge_labels,
        "representative_nodes": representative_nodes,
        "representative_links": representative_links,
    }


def build_threshold_graph(
    frames: pd.DataFrame,
    classes: List[Dict[str, Any]],
    threshold: float,
    connect_representatives: bool = True,
) -> nx.Graph:
    """Build a true RMSD graph from matrix values.

    - Nodes are frames.
    - Edges are added when RMSD(i,j) <= threshold (i != j), weighted by RMSD.
    - Class and representative metadata are added to nodes.
    - Optional representative-to-representative chain edges are added.
    """
    graph = nx.Graph()

    # Add nodes with default metadata.
    for node in frames.index.tolist():
        graph.add_node(int(node), class_name=None, representative=False)

    # RMSD threshold edges from the matrix itself.
    idx = frames.index.tolist()
    for i, a in enumerate(idx):
        for b in idx[i + 1 :]:
            w = float(frames.loc[a, b])
            if w <= threshold:
                graph.add_edge(int(a), int(b), weight=w, edge_type="rmsd")

    reps: List[int] = []
    for info in classes:
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        reps.append(rep)
        for node in info["members"]:
            n = int(node)
            if n in graph.nodes:
                graph.nodes[n]["class_name"] = class_name
            if n == rep:
                graph.nodes[n]["representative"] = True

    if connect_representatives and len(reps) > 1:
        for a, b in zip(reps[:-1], reps[1:]):
            graph.add_edge(a, b, weight=0.0, edge_type="rep_link", label="rep_link")

    return graph


def run_covering_algorithm(
    input_csv: Path = DEFAULT_INPUT,
    threshold: float = 2.5,
    max_frames: int | None = None,
) -> Dict[str, Any]:
    """Convenience wrapper that returns all major in-memory outputs.

    Returns:
      - frames: cleaned RMSD dataframe with numeric labels
      - classes: list of class dictionaries
      - class_table: pandas dataframe summary
      - graph: networkx graph with threshold and representative edges
    """
    frames = load_frames_matrix(input_csv, max_frames=max_frames)
    all_rows = set(int(x) for x in frames.index.tolist())
    start_rep = 1 if 1 in all_rows else min(all_rows)

    classes = build_classes_recursive(
        frames=frames,
        threshold=float(threshold),
        remaining_rows=all_rows,
        blocked_representatives=set(),
        class_id=1,
        representative=start_rep,
    )

    class_table = classes_to_table(classes)
    graph = build_threshold_graph(frames, classes, threshold=float(threshold), connect_representatives=True)
    class_edges = build_class_edge_map(frames, classes, threshold=float(threshold))
    display = build_display_metadata(frames, classes, class_edges)

    return {
        "frames": frames,
        "classes": classes,
        "class_table": class_table,
        "graph": graph,
        "class_edges": class_edges,
        "display": display,
    }


def plot_class_cluster_graph(
    result: Dict[str, Any],
    output_path: Path,
    title: str = "Covering Graph by Class Cluster",
) -> None:
    """Plot a class-cluster graph with the representative centered in each cluster.

    Layout strategy:
    - Each class is placed on a large outer ring so classes are visually separated.
    - The class representative is fixed at the center of that class.
    - Members are placed around the representative in a small ring.
    - Edges from representative to members are labeled with angstrom distances.
    - Representative-to-representative links are drawn as dashed connectors.
    """
    frames: pd.DataFrame = result["frames"]
    classes: List[Dict[str, Any]] = result["classes"]
    graph: nx.Graph = result["graph"]
    class_edges: Dict[str, List[Tuple[int, int, float]]] = result["class_edges"]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    positions: Dict[int, Tuple[float, float]] = {}
    num_classes = max(1, len(classes))
    outer_radius = max(18.0, 3.5 * num_classes)

    for idx, info in enumerate(classes):
        rep = int(info["representative"])
        members = [int(x) for x in info["members"]]

        angle = 2.0 * np.pi * (idx / num_classes)
        center_x = outer_radius * np.cos(angle)
        center_y = outer_radius * np.sin(angle)
        positions[rep] = (center_x, center_y)

        other_members = [m for m in members if m != rep]
        ring_radius = 2.1 + 0.02 * len(other_members)
        for j, node in enumerate(other_members):
            theta = 2.0 * np.pi * (j / max(1, len(other_members)))
            positions[node] = (
                center_x + ring_radius * np.cos(theta),
                center_y + ring_radius * np.sin(theta),
            )

    # Safety fallback for any node not assigned a cluster position.
    for node in frames.index.tolist():
        positions.setdefault(int(node), (0.0, 0.0))

    fig, ax = plt.subplots(figsize=(20, 14))
    ax.set_axis_off()

    # Draw circles around each class cluster.
    class_colors = plt.get_cmap("tab20", num_classes)
    for idx, info in enumerate(classes):
        members = [int(x) for x in info["members"]]
        xy = np.array([positions[m] for m in members], dtype=float)
        center = xy.mean(axis=0)
        radius = float(np.max(np.linalg.norm(xy - center, axis=1)) + 1.0)
        circle = plt.Circle(
            (center[0], center[1]),
            radius,
            fill=False,
            linewidth=2.0,
            alpha=0.65,
            edgecolor=class_colors(idx),
        )
        ax.add_patch(circle)

    # Draw per-class representative-to-member edges with RMSD labels.
    edge_label_map: Dict[Tuple[int, int], str] = {}
    for idx, info in enumerate(classes):
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        edge_list = [(u, v) for u, v, _ in class_edges.get(class_name, [])]
        if edge_list:
            nx.draw_networkx_edges(
                graph,
                positions,
                edgelist=edge_list,
                edge_color="#7f8c8d",
                width=0.9,
                alpha=0.35,
                ax=ax,
            )
        for u, v, w in class_edges.get(class_name, []):
            edge_label_map[(u, v)] = f"{w:.2f} Å"

    # Draw representative-to-representative connectors.
    rep_nodes = [int(info["representative"]) for info in classes]
    rep_rep_edges = list(zip(rep_nodes[:-1], rep_nodes[1:]))
    if rep_rep_edges:
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=rep_rep_edges,
            edge_color="black",
            width=2.2,
            style="dashed",
            alpha=0.9,
            ax=ax,
        )
        for u, v in rep_rep_edges:
            edge_label_map[(u, v)] = "rep link"

    # Draw nodes, making representatives larger.
    node_sizes = []
    node_colors = []
    for node in graph.nodes():
        if graph.nodes[node].get("representative"):
            node_sizes.append(560)
            node_colors.append("#e63946")
        else:
            node_sizes.append(95)
            node_colors.append("#a8dadc")

    nx.draw_networkx_nodes(
        graph,
        positions,
        node_size=node_sizes,
        node_color=node_colors,
        edgecolors="black",
        linewidths=0.5,
        ax=ax,
    )

    labels = {
        node: (f"F{node}" if graph.nodes[node].get("representative") else str(node))
        for node in graph.nodes()
    }
    nx.draw_networkx_labels(graph, positions, labels=labels, font_size=5, ax=ax)

    if edge_label_map:
        nx.draw_networkx_edge_labels(
            graph,
            positions,
            edge_labels=edge_label_map,
            font_size=4,
            rotate=False,
            ax=ax,
        )

    ax.set_title(title, fontsize=13)
    plt.tight_layout()
    plt.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


# Intentionally no main() execution block here.
# This file is meant for side-by-side review and in-memory debugging.
