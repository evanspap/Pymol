"""Covering algorithm for frame clustering from an RMSD all-vs-all matrix.

PowerShell usage (hypothetical examples for any institution):

1) Basic run (outputs are written beside the input CSV):
    python .\\Covering_Algorithm_Run1_7-16-26_cleaned.py --input-csv "D:\\Project\\rmsd_matrix.csv"

2) Custom output folder and threshold:
    python .\\Covering_Algorithm_Run1_7-16-26_cleaned.py --input-csv "D:\\Project\\rmsd_matrix.csv" --output-dir "D:\\Project\\outputs" --threshold 2.5

3) Limit analysis to first N frames:
    python .\\Covering_Algorithm_Run1_7-16-26_cleaned.py --input-csv "D:\\Project\\rmsd_matrix.csv" --max-frames 500

4) Custom output filenames:
    python .\\Covering_Algorithm_Run1_7-16-26_cleaned.py --input-csv "D:\\Project\\rmsd_matrix.csv" --class1-members-name "team_class1_members.csv" --class1-matrix-name "team_class1_matrix.csv" --unclassified-matrix-name "team_unclassified_matrix.csv"

Important path note:
    If another team clones/downloads this repository, their local paths will differ.
    Replace all example paths with machine-specific paths before running commands.

This module is intentionally side-effect free: it builds clean in-memory objects you can inspect or plot later.

Core logic implemented here:
1. Load RMSD matrix and relabel rows/columns to integer frame IDs (1..N).
2. Build classes using your rule:
    - Start with representative frame 1.
    - Class_k contains rows with RMSD <= threshold to representative frame.
    - A future representative cannot be equal to any row number already in prior classes.
3. Build a true RMSD threshold graph (edges from matrix values <= threshold).
4. Connect representative nodes to each other to keep class progression visible.
"""
#### Disclaimer: Generative AI was used to assist in 
# writing this code. 
# The author has reviewed and 
# edited the content as needed and 
# takes full responsibility for the final code.
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


# These defaults are kept for convenience and backward compatibility.
# Cross-institution users should pass --input-csv and --output-dir explicitly.
DEFAULT_INPUT = Path(
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\csv\rmsd_matrix_v2_all_vs_all_biopython.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\csv"
)


def load_frames_matrix(csv_path: Path, max_frames: int | None = None) -> pd.DataFrame:
    """Load RMSD matrix and relabel both axes to integer frame IDs.

    Expected input shape:
    - Square RMSD matrix CSV with frame labels in first column and header row.
    - Values represent pairwise RMSD distances between frames.

    Normalization performed here:
    - Row/column labels are reset to 1..N integer frame IDs so downstream logic
      is deterministic across differently labeled source CSV files.
    - Optional truncation keeps top-left NxN block for max_frames experiments.
    """
    # Read the CSV into a DataFrame.
    # `index_col=0` means the first column in the CSV becomes the row index.
    matrix = pd.read_csv(csv_path, index_col=0)

    # Count how many columns exist. For a valid all-vs-all matrix, this should
    # match the number of rows (square matrix), but we do not hard-fail here.
    n_cols = matrix.shape[1]

    # Standardize labels to 1..N for BOTH columns and rows.
    # This avoids dependence on any original frame naming style.
    matrix.columns = list(range(1, n_cols + 1))
    matrix.index = list(range(1, matrix.shape[0] + 1))

    if max_frames is not None:
        # Keep only the top-left NxN section when running a smaller test.
        # Example: if max_frames=500, keep rows/columns 1..500.
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

    Recursive process overview:
    1) Choose valid representative (or recover to smallest available row).
    2) Select all currently remaining rows within threshold of representative.
    3) Record class metadata.
    4) Remove class members from remaining pool.
    5) Recurse until no rows remain.
    """
    if not remaining_rows:
        # Base case: no frames left to classify.
        return []

    if representative not in remaining_rows or representative in blocked_representatives:
        # If representative is invalid (already removed or blocked), choose a new one.
        # Preferred choice: smallest remaining frame not blocked.
        available = sorted(remaining_rows - blocked_representatives)
        if not available:
            # Fallback: if all remaining rows are blocked, continue with smallest remaining.
            available = sorted(remaining_rows)
        representative = available[0]

    # Work on a sorted list so behavior is deterministic and reproducible.
    current_rows = sorted(remaining_rows)

    # Core class rule: frame belongs to current class if its RMSD to the
    # representative frame is <= threshold.
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

    # Remove everything assigned to this class from the remaining pool.
    new_remaining = set(remaining_rows) - set(class_members)

    # Block all class members from becoming future representatives.
    new_blocked = set(blocked_representatives) | set(class_members)
    if not new_remaining:
        # No unassigned frames left: recursion ends here.
        return [class_info]

    # Pick the next representative, preferring one not blocked.
    candidate_next = sorted(new_remaining - new_blocked)
    if candidate_next:
        next_representative = candidate_next[0]
    else:
        # Emergency fallback if every remaining frame is blocked.
        next_representative = min(new_remaining)

    # Recursive step: build later classes and append to current class.
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
        # Convert member list to a human-readable comma-separated string.
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


def save_class_table_png(
    class_table: pd.DataFrame,
    output_path: Path,
    title: str = "Covering Class Table",
) -> None:
    """Render a class summary dataframe to a PNG image."""
    # Create output directory if needed.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Dynamically scale figure height by number of classes to avoid clipping.
    fig_height = max(6.0, 0.35 * (len(class_table) + 2))
    fig, ax = plt.subplots(figsize=(20, fig_height))
    ax.axis("off")
    ax.set_title(title, fontsize=14, pad=14)

    table = ax.table(
        cellText=class_table.values,
        colLabels=class_table.columns,
        loc="center",
        cellLoc="left",
        colLoc="left",
    )
    # Tune typography so larger tables remain legible.
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.35)

    plt.tight_layout()
    plt.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)


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
        # One edge list per class.
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        members = [int(x) for x in info["members"]]
        edges: List[Tuple[int, int, float]] = []

        for node in members:
            if node == rep:
                # Skip self-edge representative -> representative.
                continue
            rmsd = float(frames.loc[rep, node])
            if rmsd <= threshold:
                # Keep only edges that satisfy the threshold criterion.
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
    # Start with default visual settings for all frames.
    node_sizes: Dict[int, int] = {int(node): 90 for node in frames.index.tolist()}
    node_labels: Dict[int, str] = {int(node): f"F{int(node)}" for node in frames.index.tolist()}
    edge_labels: Dict[Tuple[int, int], str] = {}

    representative_nodes: List[int] = []
    for info in classes:
        class_name = str(info["class_name"])
        rep = int(info["representative"])
        representative_nodes.append(rep)

        # Representatives are emphasized with larger size + richer label.
        node_sizes[rep] = 420
        node_labels[rep] = f"F{rep}\n{class_name}"

        for u, v, w in class_edges[class_name]:
            edge_labels[(u, v)] = f"{w:.2f} Å"

    # Also prepare a chain linking class representatives in order.
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
    # Undirected graph: RMSD(i,j) equals RMSD(j,i), so one edge is enough.
    graph = nx.Graph()

    # Add nodes with default metadata.
    for node in frames.index.tolist():
        graph.add_node(int(node), class_name=None, representative=False)

    # RMSD threshold edges from the matrix itself.
    idx = frames.index.tolist()
    for i, a in enumerate(idx):
        # Only visit upper triangle (j > i) to avoid duplicate edges.
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
        # Draw progression link between representative of Class_k and Class_{k+1}.
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
    # 1) Load and standardize RMSD matrix labels.
    frames = load_frames_matrix(input_csv, max_frames=max_frames)

    # 2) Determine all frame IDs currently under consideration.
    all_rows = set(int(x) for x in frames.index.tolist())

    # 3) Start representative defaults to frame 1 when available.
    start_rep = 1 if 1 in all_rows else min(all_rows)

    # 4) Build classes recursively according to threshold rule.
    classes = build_classes_recursive(
        frames=frames,
        threshold=float(threshold),
        remaining_rows=all_rows,
        blocked_representatives=set(),
        class_id=1,
        representative=start_rep,
    )

    # TABLE OUTPUT BLOCK DISABLED FOR NOW (kept for later use):
    # class_table = classes_to_table(classes)

    # GRAPH/PLOT OUTPUT BLOCK DISABLED FOR NOW (kept for later use):
    # graph = build_threshold_graph(frames, classes, threshold=float(threshold), connect_representatives=True)
    # class_edges = build_class_edge_map(frames, classes, threshold=float(threshold))
    # display = build_display_metadata(frames, classes, class_edges)

    # Keep return contract stable while optional outputs are disabled.
    class_table = None
    graph = None
    class_edges = None
    display = None

    return {
        "frames": frames,
        "classes": classes,
        "class_table": class_table,
        "graph": graph,
        "class_edges": class_edges,
        "display": display,
    }


def export_class1_and_unclassified_csvs(
    frames: pd.DataFrame,
    classes: List[Dict[str, Any]],
    output_dir: Path,
    class1_members_name: str = "class1_members.csv",
    class1_matrix_name: str = "class1_rmsd_matrix.csv",
    unclassified_matrix_name: str = "unclassified_rmsd_matrix.csv",
) -> Dict[str, Path]:
    """Export Class 1 membership + Class 1 and unclassified RMSD submatrices.

    Outputs:
    - class1_members.csv: one row per Class 1 member frame ID.
    - class1_rmsd_matrix.csv: RMSD submatrix restricted to Class 1 members.
    - unclassified_rmsd_matrix.csv: RMSD submatrix for frames not in Class 1.
    """
    # Ensure the destination exists before writing any CSV artifacts.
    output_dir.mkdir(parents=True, exist_ok=True)

    # Locate Class_1 exactly; if missing, produce empty CSV outputs safely.
    class1_info = next((c for c in classes if str(c.get("class_name")) == "Class_1"), None)
    class1_members: List[int] = []
    if class1_info is not None:
        class1_members = sorted(int(x) for x in class1_info.get("members", []))

    # Unclassified is defined here as: all frames not belonging to Class_1.
    all_frames = [int(x) for x in frames.index.tolist()]
    class1_set = set(class1_members)
    unclassified_members = [x for x in all_frames if x not in class1_set]

    # Resolve all file destinations once, then write each artifact.
    class1_members_path = output_dir / class1_members_name
    class1_matrix_path = output_dir / class1_matrix_name
    unclassified_matrix_path = output_dir / unclassified_matrix_name

    # 1) Flat listing of Class 1 members.
    class1_members_df = pd.DataFrame({"Frame": class1_members})
    class1_members_df.to_csv(class1_members_path, index=False)

    # 2) Class 1 RMSD matrix.
    # Matrix slicing keeps only rows+columns corresponding to Class 1 frames.
    class1_matrix = frames.loc[class1_members, class1_members] if class1_members else pd.DataFrame()
    class1_matrix.to_csv(class1_matrix_path)

    # 3) Unclassified RMSD matrix (not in Class 1).
    # Same operation for everything that was not in Class 1.
    unclassified_matrix = (
        frames.loc[unclassified_members, unclassified_members]
        if unclassified_members
        else pd.DataFrame()
    )
    unclassified_matrix.to_csv(unclassified_matrix_path)

    return {
        "class1_members_csv": class1_members_path,
        "class1_matrix_csv": class1_matrix_path,
        "unclassified_matrix_csv": unclassified_matrix_path,
    }


def run_and_export_class1_outputs(
    input_csv: Path = DEFAULT_INPUT,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    threshold: float = 2.5,
    max_frames: int | None = None,
) -> Dict[str, Any]:
    """Run covering algorithm, then export requested Class 1/unclassified CSV outputs.

    This function is the main programmatic entry point for automation pipelines.
    It returns both in-memory analysis objects and concrete output file paths.
    """
    # Run analysis first.
    result = run_covering_algorithm(
        input_csv=input_csv,
        threshold=threshold,
        max_frames=max_frames,
    )
    # Then materialize on-disk CSV outputs from that in-memory result.
    exports = export_class1_and_unclassified_csvs(
        frames=result["frames"],
        classes=result["classes"],
        output_dir=output_dir,
    )
    result["exports"] = exports
    return result


def parse_cli_args() -> argparse.Namespace:
    """Parse command-line arguments for portable cross-team execution.

    Design choices:
    - --input-csv is required to force explicit user intent.
    - --output-dir is optional; defaults to input CSV parent folder.
    - Output filenames are overridable to support lab-specific naming standards.
    """
    # CLI parser defines user-facing contract for command-line execution.
    parser = argparse.ArgumentParser(
        description="Run covering algorithm and export Class 1 / unclassified RMSD CSV outputs.",
        epilog=(
            "Examples:\n"
            "  python Covering_Algorithm_Run1_7-16-26_cleaned.py --input-csv D:/Project/rmsd_matrix.csv\n"
            "  python Covering_Algorithm_Run1_7-16-26_cleaned.py --input-csv D:/Project/rmsd_matrix.csv --output-dir D:/Project/outputs --threshold 2.5"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        required=True,
        help="Path to all-vs-all RMSD matrix CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for generated CSV files (defaults to input CSV directory).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=2.5,
        help="RMSD threshold used by the covering algorithm.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional cap: only analyze first N frames from the matrix.",
    )
    parser.add_argument(
        "--class1-members-name",
        type=str,
        default="class1_members.csv",
        help="Output filename for Class 1 member listing.",
    )
    parser.add_argument(
        "--class1-matrix-name",
        type=str,
        default="class1_rmsd_matrix.csv",
        help="Output filename for Class 1 RMSD matrix.",
    )
    parser.add_argument(
        "--unclassified-matrix-name",
        type=str,
        default="unclassified_rmsd_matrix.csv",
        help="Output filename for unclassified RMSD matrix.",
    )
    return parser.parse_args()


def main() -> None:
    """CLI orchestration entry point.

    Workflow:
    1) Parse CLI values.
    2) Resolve and validate paths.
    3) Run covering analysis.
    4) Export requested CSV artifacts.
    5) Print artifact paths for reproducible audit/logging.
    """
    # Read all command-line values.
    args = parse_cli_args()

    # Expand potential ~ in paths and validate input file existence.
    input_csv = args.input_csv.expanduser()
    if not input_csv.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_csv}")

    # If output folder is omitted, write outputs beside input CSV.
    output_dir = args.output_dir.expanduser() if args.output_dir else input_csv.parent

    # Execute the class-building algorithm.
    result = run_covering_algorithm(
        input_csv=input_csv,
        threshold=float(args.threshold),
        max_frames=args.max_frames,
    )

    # Export requested CSV files using names provided by CLI flags.
    exports = export_class1_and_unclassified_csvs(
        frames=result["frames"],
        classes=result["classes"],
        output_dir=output_dir,
        class1_members_name=args.class1_members_name,
        class1_matrix_name=args.class1_matrix_name,
        unclassified_matrix_name=args.unclassified_matrix_name,
    )

    # Print all output paths so user can confirm where files were written.
    print("Generated CSV outputs:")
    for name, path in exports.items():
        print(f" - {name}: {path}")


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
    # Pull precomputed objects from result dictionary.
    frames: pd.DataFrame = result["frames"]
    classes: List[Dict[str, Any]] = result["classes"]
    graph: nx.Graph = result["graph"]
    class_edges: Dict[str, List[Tuple[int, int, float]]] = result["class_edges"]

    # Ensure destination folder exists.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # `positions` will store explicit 2D coordinates for every node.
    positions: Dict[int, Tuple[float, float]] = {}
    num_classes = max(1, len(classes))
    outer_radius = max(18.0, 3.5 * num_classes)

    for idx, info in enumerate(classes):
        rep = int(info["representative"])
        members = [int(x) for x in info["members"]]

        # Place each class center around a large circle to reduce overlap.
        angle = 2.0 * np.pi * (idx / num_classes)
        center_x = outer_radius * np.cos(angle)
        center_y = outer_radius * np.sin(angle)
        positions[rep] = (center_x, center_y)

        # Place non-representative members around that center in a small ring.
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

    # Build a large figure because dense graphs become unreadable quickly.
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
    # Keep edge labels in a dictionary so we can draw them all at once.
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
    # Choose size and color per node based on representative status.
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


if __name__ == "__main__":
    main()
