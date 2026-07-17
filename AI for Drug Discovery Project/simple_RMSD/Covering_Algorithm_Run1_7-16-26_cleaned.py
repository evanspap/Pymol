"""Covering algorithm for frame clustering from an RMSD all-vs-all matrix.

Workflow:
1. Load RMSD matrix CSV and rename columns/index to integer frame IDs (1..N).
2. Build classes recursively with threshold filtering (default 4.5 A):
   - Start from representative frame 1.
   - Class_k = all remaining rows where RMSD to representative <= threshold.
   - Next representative must come from rows not already covered by previous class.
3. Save a PNG table of classes.
4. Save a NetworkX cluster graph PNG with:
   - nodes as frames,
   - representative nodes highlighted,
   - circles around each cluster,
   - representative nodes connected to each other.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.patches import Circle


DEFAULT_INPUT = Path(
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\csv\rmsd_matrix_v2_all_vs_all_biopython.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Recursive frame covering from RMSD matrix.")
    parser.add_argument("--input", "-i", type=Path, default=DEFAULT_INPUT, help="Input RMSD CSV")
    parser.add_argument("--output-dir", "-o", type=Path, default=DEFAULT_OUTPUT_DIR, help="Output directory")
    parser.add_argument("--threshold", type=float, default=4.5, help="RMSD threshold in Angstrom")
    parser.add_argument("--table-name", type=str, default="test_table_7-16-26.png", help="Output table PNG")
    parser.add_argument("--graph-name", type=str, default="test_graph_7-16-26.png", help="Output graph PNG")
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional debug limit: use only first N frames.",
    )
    return parser.parse_args()


def load_frames_matrix(csv_path: Path, max_frames: int | None = None) -> pd.DataFrame:
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
    class_id: int,
    representative: int,
) -> List[Dict[str, object]]:
    if not remaining_rows:
        return []

    if representative not in remaining_rows:
        representative = min(remaining_rows)

    current_rows = sorted(remaining_rows)
    eligible_mask = frames.loc[current_rows, representative] <= threshold
    class_members = [row for row, keep in zip(current_rows, eligible_mask) if bool(keep)]

    # Diagonal should already include representative, but keep a safe guard.
    if representative not in class_members:
        class_members.append(representative)
        class_members = sorted(set(class_members))

    class_info = {
        "class_name": f"Class_{class_id}",
        "representative": representative,
        "members": class_members,
        "member_count": len(class_members),
    }

    new_remaining = set(remaining_rows) - set(class_members)
    if not new_remaining:
        return [class_info]

    next_representative = min(new_remaining)
    return [class_info] + build_classes_recursive(
        frames=frames,
        threshold=threshold,
        remaining_rows=new_remaining,
        class_id=class_id + 1,
        representative=next_representative,
    )


def build_class_table_png(classes: List[Dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for info in classes:
        members = info["members"]
        members_text = ", ".join(str(x) for x in members)
        rows.append([
            info["class_name"],
            info["representative"],
            info["member_count"],
            members_text,
        ])

    headers = ["Class", "Representative Frame", "N Members", "Members"]
    fig_h = max(4.0, 1.0 + 0.48 * len(rows))
    fig, ax = plt.subplots(figsize=(18, fig_h))
    ax.axis("off")

    tbl = ax.table(
        cellText=rows,
        colLabels=headers,
        cellLoc="left",
        loc="center",
        colColours=["#4f81bd"] * len(headers),
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1.0, 1.35)

    for (row, col), cell in tbl.get_celld().items():
        if row == 0:
            cell.set_text_props(color="white", weight="bold")
        else:
            cell.set_facecolor("#f9f9f9" if row % 2 == 0 else "#ffffff")

    ax.set_title("Recursive Covering Classes (Threshold <= 4.5 A)", fontsize=12, pad=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_cluster_graph(classes: List[Dict[str, object]]) -> Tuple[nx.Graph, Dict[int, Tuple[float, float]], Dict[str, Tuple[float, float, float, float]]]:
    graph = nx.Graph()
    positions: Dict[int, Tuple[float, float]] = {}
    class_colors: Dict[str, Tuple[float, float, float, float]] = {}

    cmap = plt.get_cmap("tab20", max(1, len(classes)))

    representatives: List[int] = []
    for idx, info in enumerate(classes):
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        members = [int(x) for x in info["members"]]
        representatives.append(rep)

        color = cmap(idx)
        class_colors[class_name] = color

        center_x = idx * 8.0
        center_y = 0.0
        positions[rep] = (center_x, center_y)

        graph.add_node(rep, class_name=class_name, representative=True)

        non_rep_members = [m for m in members if m != rep]
        n_non_rep = max(1, len(non_rep_members))
        radius = 1.8 + 0.03 * len(members)

        for j, node in enumerate(non_rep_members):
            angle = 2.0 * np.pi * (j / n_non_rep)
            x = center_x + radius * np.cos(angle)
            y = center_y + radius * np.sin(angle)
            positions[node] = (x, y)
            graph.add_node(node, class_name=class_name, representative=False)
            graph.add_edge(rep, node, edge_type="member")

    for a, b in zip(representatives[:-1], representatives[1:]):
        graph.add_edge(a, b, edge_type="rep_link")

    return graph, positions, class_colors


def draw_cluster_graph_png(
    classes: List[Dict[str, object]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    graph, positions, class_colors = build_cluster_graph(classes)

    fig_w = max(12.0, 3.2 * len(classes))
    fig, ax = plt.subplots(figsize=(fig_w, 8.5))

    rep_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get("edge_type") == "rep_link"]
    member_edges = [(u, v) for u, v, d in graph.edges(data=True) if d.get("edge_type") == "member"]

    nx.draw_networkx_edges(graph, positions, edgelist=member_edges, width=1.0, alpha=0.35, edge_color="#777777", ax=ax)
    nx.draw_networkx_edges(graph, positions, edgelist=rep_edges, width=2.2, alpha=0.95, edge_color="black", style="dashed", ax=ax)

    # Draw cluster circles around each class.
    for info in classes:
        class_name = str(info["class_name"])
        members = [int(x) for x in info["members"]]
        xy = np.array([positions[m] for m in members], dtype=float)
        center = xy.mean(axis=0)
        radius = float(np.max(np.linalg.norm(xy - center, axis=1)) + 0.9)
        circle = Circle(
            (center[0], center[1]),
            radius,
            fill=False,
            edgecolor=class_colors[class_name],
            linewidth=2.0,
            alpha=0.8,
        )
        ax.add_patch(circle)

    non_rep_nodes = [n for n, d in graph.nodes(data=True) if not d.get("representative")]
    rep_nodes = [n for n, d in graph.nodes(data=True) if d.get("representative")]

    non_rep_colors = [class_colors[graph.nodes[n]["class_name"]] for n in non_rep_nodes]
    rep_colors = [class_colors[graph.nodes[n]["class_name"]] for n in rep_nodes]

    nx.draw_networkx_nodes(
        graph,
        positions,
        nodelist=non_rep_nodes,
        node_size=260,
        node_color=non_rep_colors,
        edgecolors="#2a2a2a",
        linewidths=0.8,
        ax=ax,
    )
    nx.draw_networkx_nodes(
        graph,
        positions,
        nodelist=rep_nodes,
        node_size=480,
        node_color=rep_colors,
        edgecolors="black",
        linewidths=2.0,
        ax=ax,
    )

    labels = {}
    for node, data in graph.nodes(data=True):
        cls = data["class_name"].replace("Class_", "C")
        labels[node] = f"F{node}\n{cls}" if data.get("representative") else f"{node}"

    nx.draw_networkx_labels(graph, positions, labels=labels, font_size=6, ax=ax)

    ax.set_title("Frame Covering Graph: Clusters and Representative Links", fontsize=13)
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_csv = args.input
    output_dir = args.output_dir
    threshold = float(args.threshold)

    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    output_dir.mkdir(parents=True, exist_ok=True)

    frames = load_frames_matrix(input_csv, max_frames=args.max_frames)
    all_rows = set(int(x) for x in frames.index.tolist())
    start_rep = 1 if 1 in all_rows else min(all_rows)

    classes = build_classes_recursive(
        frames=frames,
        threshold=threshold,
        remaining_rows=all_rows,
        class_id=1,
        representative=start_rep,
    )

    table_png = output_dir / args.table_name
    graph_png = output_dir / args.graph_name

    build_class_table_png(classes, table_png)
    draw_cluster_graph_png(classes, graph_png)

    print(f"Loaded matrix shape: {frames.shape}")
    print(f"Threshold: {threshold:.3f} A")
    print(f"Generated {len(classes)} classes")
    print(f"Saved table PNG: {table_png}")
    print(f"Saved graph PNG: {graph_png}")


if __name__ == "__main__":
    main()
