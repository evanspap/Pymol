
r"""
Curve_Fit.py
Version: 2.2.0

Purpose:
    PDB-driven RMSD histogram generator with Gaussian curve fitting.

What this script can do:
    1) Full-atom all-frames-vs-all-frames RMSD histogram
    2) Backbone-only all-frames-vs-all-frames RMSD histogram
    3) Custom atom-selection all-frames-vs-all-frames RMSD histograms
    4) Adaptive statistical fit overlay on top of probability histograms
       - Gaussian for symmetric unimodal distributions
       - Skew-Gaussian for asymmetric unimodal distributions
       - Two-Gaussian mixture for bimodal-capable distributions

Key requirements implemented:
    - Input is a PDB file (multi-frame trajectory style PDB)
    - Bin width defaults to 0.2 Angstrom
    - Y-axis is probability per bin
    - CLI controls comparisons, no hardcoded paths

AI usage disclaimer:
    Generative AI may have been used to assist with this code and its comments.
    The script is provided as a research helper, not a scientific authority.

    Before using results in analysis, publication, or decision-making:
    1) Validate parsing against known reference datasets.
    2) Verify atom selections reflect your scientific intent.
    3) Confirm Gaussian fit behavior on your distribution shape.
    4) Review outputs manually and with independent scripts/tools.

    The user and research team are responsible for quality control,
    interpretation, and scientific conclusions.

Implementation notes:
    - This script computes RMSD directly from coordinates as provided.
    - It does not apply alignment/superposition transformations.
    - All-vs-all is computed from atoms present in every frame for a selection.

Sample PowerShell commands (replace paths with your own):

     1) Full only:
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons full --bin-width 0.2

     2) Backbone only:
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons backbone --bin-width 0.2

     3) Full + backbone:
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons full backbone --bin-width 0.2

     4) Custom only (one atom-group):
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons custom --custom-selection SideChain=CB,CG,CD --bin-width 0.2

     5) Custom only (multiple atom-groups):
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons custom --custom-selection SideChain=CB,CG,CD --custom-selection Polar=N,O,S --bin-width 0.2

     6) Full + backbone + custom:
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons full backbone custom --custom-selection SideChain=CB,CG,CD --bin-width 0.2

     7) Exclude diagonal and write pairwise CSV tables:
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons full backbone --bin-width 0.2 --include-diagonal N --write-pairwise-csv Y

     8) Force a standard x-axis maximum across runs:
         py -3 ".\Curve_Fit .py" --pdb "C:\path\to\frames.pdb" --output-dir "C:\path\to\out" --comparisons full --bin-width 0.2 --x-max-angstrom 12.0
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np

# SciPy is optional. If present, we use robust non-linear fitting and
# probability density functions. If absent, we keep graceful fallbacks so
# the script remains runnable on minimal Python environments.
try:
    from scipy.optimize import OptimizeWarning, curve_fit as scipy_curve_fit
    from scipy.stats import norm as scipy_norm
    from scipy.stats import skewnorm as scipy_skewnorm
except Exception:
    scipy_curve_fit = None
    OptimizeWarning = Warning
    scipy_norm = None
    scipy_skewnorm = None


BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O", "OXT"}


# ============================================================================
# Data containers (simple structured records used throughout the pipeline)
# ============================================================================


@dataclass(frozen=True)
class AtomRecord:
    """Stores immutable atom metadata and coordinates for one frame."""

    atom_id: Tuple[str, str, str, str, str]
    atom_name: str
    coord: np.ndarray


@dataclass
class Frame:
    """Stores one frame worth of atoms parsed from the PDB."""

    atoms: List[AtomRecord]


@dataclass
class FitSummary:
    """Stores selected-fit values for reporting and CSV export."""

    comparison: str
    n_frames: int
    n_pairs: int
    n_atoms_common: int
    mean_rmsd: float
    std_rmsd: float
    fit_model: str
    fit_bic: float
    fit_amplitude: float
    fit_mu: float
    fit_sigma: float


@dataclass
class ComparisonResult:
    """Stores all computed RMSD artifacts for one comparison label."""

    label: str
    rmsd_values: np.ndarray
    pair_rows: List[Tuple[int, int, float]]
    n_common_atoms: int


@dataclass
class FitResult:
    """Stores the selected statistical fit model and its parameters.

    model_name values:
        - gaussian
        - skew_gaussian
        - gaussian_mixture_2
    """

    model_name: str
    params: Tuple[float, ...]
    bic: float


def parse_args() -> argparse.Namespace:
    """Build and parse command-line arguments for the workflow.

    The CLI is intentionally explicit and beginner-friendly:
    - `--pdb` is required so input origin is always clear.
    - `--comparisons` lets users choose built-in and custom groups.
    - `--bin-width` defaults to 0.2 Angstrom per user requirement.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate RMSD probability histograms and adaptive statistical "
            "fits from a PDB file for full/backbone/custom comparisons."
        )
    )
    # Input/Output paths
        # Histogram controls
        # Comparison controls
        # Pairwise distribution controls
    parser.add_argument(
        "--pdb",
        required=True,
        help="Path to multi-frame PDB file.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).parent / "curve_fit_output"),
        help="Directory for output CSV and PNG files.",
    )
    parser.add_argument(
        "--bin-width",
        type=float,
        default=0.2,
        help="Histogram bin width in Angstrom (default: 0.2).",
    )
    parser.add_argument(
        "--x-max-angstrom",
        type=float,
        default=None,
        help=(
            "Optional fixed x-axis maximum in Angstrom. If omitted, the script "
            "uses the maximum RMSD observed across all selected comparisons."
        ),
    )
    parser.add_argument(
        "--comparisons",
        nargs="+",
        choices=["full", "backbone", "custom"],
        default=["full", "backbone"],
        help=(
            "Which comparison groups to run. Use one or more of: full backbone custom. "
            "Default: full backbone"
        ),
    )
    parser.add_argument(
        "--custom-selection",
        action="append",
        default=[],
        help=(
            "Custom atom-name selection in the form Label=ATOM1,ATOM2,... "
            "Example: SideChain=CB,CG,CD"
        ),
    )
    parser.add_argument(
        "--include-diagonal",
        choices=["Y", "N"],
        default="Y",
        help=(
            "Include frame-vs-same-frame RMSD values (zeros) in all-vs-all distribution. "
            "Y or N (default: Y)."
        ),
    )
    parser.add_argument(
        "--write-pairwise-csv",
        choices=["Y", "N"],
        default="N",
        help="Write per-comparison pairwise RMSD tables as CSV. Y or N (default: N).",
    )
    parser.add_argument(
        "--fit-mode",
        choices=["adaptive", "gaussian", "none"],
        default="adaptive",
        help=(
            "Fit behavior for histogram overlays: adaptive (default), gaussian, "
            "or none."
        ),
    )
    return parser.parse_args()


def _safe_float(text: str) -> Optional[float]:
    """Convert text to float; return None for malformed numeric fields."""
    try:
        return float(text.strip())
    except Exception:
        return None


def _build_atom_id_from_pdb_line(line: str) -> Tuple[str, str, str, str, str]:
    """
    Build a stable atom identifier using fields that should remain constant
    for a given atom across frames in a trajectory-style PDB.

    Tuple layout:
        (chain_id, residue_index, insertion_code, atom_name, alt_loc)

    Why this matters:
        RMSD requires matching "the same atom" across frames.
        Matching only by atom serial number can fail in some PDB exports,
        so this identifier uses structural descriptors instead.
    """
    atom_name = line[12:16].strip()
    chain_id = line[21:22].strip()
    res_seq = line[22:26].strip()
    insertion_code = line[26:27].strip()
    alt_loc = line[16:17].strip()
    return (chain_id, res_seq, insertion_code, atom_name, alt_loc)


def parse_pdb_frames(pdb_path: Path) -> List[Frame]:
    """
    Parse multi-frame PDB file into a list of Frame objects.

    Frame boundaries are detected by ENDMDL or END records. If MODEL appears
    while a frame is open, the previous frame is finalized first.

    Parsing strategy:
        - Read line-by-line to handle large PDB files efficiently.
        - Keep only ATOM/HETATM records with valid numeric coordinates.
        - Ignore malformed coordinate lines instead of crashing.

    Returns:
        A list where each item represents one frame in the trajectory.
    """
    frames: List[Frame] = []
    current_atoms: List[AtomRecord] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            record = line[:6].strip().upper()

            if record == "MODEL":
                # Start of a new model. If we already collected atoms, push them
                # as a complete frame before continuing.
                if current_atoms:
                    frames.append(Frame(atoms=current_atoms))
                    current_atoms = []
                continue

            if record in {"ATOM", "HETATM"}:
                # PDB fixed-width coordinate columns.
                x = _safe_float(line[30:38])
                y = _safe_float(line[38:46])
                z = _safe_float(line[46:54])
                if x is None or y is None or z is None:
                    # Skip lines that do not contain parseable coordinates.
                    continue

                atom_id = _build_atom_id_from_pdb_line(line)
                atom_name = line[12:16].strip()
                coord = np.array([x, y, z], dtype=float)
                current_atoms.append(AtomRecord(atom_id=atom_id, atom_name=atom_name, coord=coord))
                continue

            if record in {"ENDMDL", "END"}:
                # End of current frame block.
                if current_atoms:
                    frames.append(Frame(atoms=current_atoms))
                    current_atoms = []
                continue

    if current_atoms:
        # Final safeguard: flush last frame if file does not end with END/ENDMDL.
        frames.append(Frame(atoms=current_atoms))

    return frames


# ============================================================================
# Atom-selection helpers (full / backbone / custom atom groups)
# ============================================================================


def _atom_names_to_set(selection_csv: str) -> Set[str]:
    """Parse comma-separated atom names and normalize to uppercase tokens."""
    return {token.strip().upper() for token in selection_csv.split(",") if token.strip()}


def parse_custom_selections(raw_items: Sequence[str]) -> Dict[str, Set[str]]:
    """Parse repeated Label=ATOM1,ATOM2 definitions from CLI.

    Example input list:
        ["SideChain=CB,CG,CD", "Polar=N,O,S"]

    Returns:
        {
            "SideChain": {"CB", "CG", "CD"},
            "Polar": {"N", "O", "S"},
        }
    """
    selections: Dict[str, Set[str]] = {}
    for item in raw_items:
        if "=" not in item:
            raise ValueError(
                f"Invalid --custom-selection value: {item!r}. Expected Label=ATOM1,ATOM2"
            )
        label, atoms_csv = item.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"Invalid custom selection label in: {item!r}")
        atom_names = _atom_names_to_set(atoms_csv)
        if not atom_names:
            raise ValueError(f"No atom names provided for custom selection: {item!r}")
        selections[label] = atom_names
    return selections


def _frame_atom_dict(frame: Frame) -> Dict[Tuple[str, str, str, str, str], AtomRecord]:
    """Build dictionary for O(1) atom lookup by atom_id within one frame."""
    return {atom.atom_id: atom for atom in frame.atoms}


def _select_atom_ids_for_frame(
    frame: Frame,
    selection_name: str,
    custom_selection_map: Dict[str, Set[str]],
) -> Set[Tuple[str, str, str, str, str]]:
    """Return atom IDs in a frame that satisfy a named selection.

    Supported selection names:
        - full
        - backbone
        - custom labels from --custom-selection
    """
    if selection_name == "full":
        return {atom.atom_id for atom in frame.atoms}
    if selection_name == "backbone":
        return {atom.atom_id for atom in frame.atoms if atom.atom_name.upper() in BACKBONE_ATOM_NAMES}

    # Custom selection by atom name membership.
    allowed_atom_names = custom_selection_map[selection_name]
    return {atom.atom_id for atom in frame.atoms if atom.atom_name.upper() in allowed_atom_names}


def _common_atom_ids_across_frames(
    frames: Sequence[Frame],
    selection_name: str,
    custom_selection_map: Dict[str, Set[str]],
) -> List[Tuple[str, str, str, str, str]]:
    """
    Return a stable, sorted list of atom IDs present in every frame for the selection.

    Why intersection is required:
        All-vs-all RMSD needs shape-consistent coordinate arrays. If an atom
        exists in frame A but not frame B, that atom is excluded globally so
        every pairwise comparison uses the same atom basis.
    """
    if not frames:
        return []

    common_ids = _select_atom_ids_for_frame(frames[0], selection_name, custom_selection_map)
    for frame in frames[1:]:
        frame_ids = _select_atom_ids_for_frame(frame, selection_name, custom_selection_map)
        common_ids &= frame_ids

    return sorted(common_ids)


def _coords_for_frame_and_ids(
    frame_atom_map: Dict[Tuple[str, str, str, str, str], AtomRecord],
    atom_ids: Sequence[Tuple[str, str, str, str, str]],
) -> np.ndarray:
    """Extract coordinate array for a frame in a fixed atom-ID order."""
    coords = [frame_atom_map[atom_id].coord for atom_id in atom_ids]
    return np.asarray(coords, dtype=float)


def calculate_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Compute RMSD between two arrays of coordinates with identical shape.

    Formula:
        RMSD = sqrt(mean(sum((A_i - B_i)^2 over xyz) over all atoms i))
    """
    if coords_a.shape != coords_b.shape:
        raise ValueError("Coordinate arrays must have the same shape for RMSD.")
    delta = coords_a - coords_b
    msd = np.mean(np.sum(delta * delta, axis=1))
    return float(np.sqrt(msd))


# ============================================================================
# All-vs-all RMSD generation
# ============================================================================


def compute_all_vs_all_rmsd(
    frames: Sequence[Frame],
    selection_name: str,
    custom_selection_map: Dict[str, Set[str]],
    include_diagonal: bool,
) -> Tuple[np.ndarray, List[Tuple[int, int, float]], int]:
    """
    Compute all-vs-all RMSD values for a given selection.

    Returns:
        - 1D numpy array of RMSD values
        - table rows (frame_i, frame_j, rmsd)
        - number of atoms used (common across all frames)

    Notes for beginners:
        - Frame numbering in outputs is 1-based for readability.
        - If diagonal is included, frame i vs i contributes RMSD = 0 values.
        - Without diagonal, this function still writes symmetric pairs
          (i,j) and (j,i) so the distribution reflects a full matrix style.
    """
    atom_ids = _common_atom_ids_across_frames(frames, selection_name, custom_selection_map)
    if not atom_ids:
        raise ValueError(
            f"No shared atoms found across frames for selection '{selection_name}'."
        )

    # Precompute dictionary lookups and aligned coordinate arrays once to
    # avoid repeated atom matching inside the nested pair loops.
    frame_maps = [_frame_atom_dict(frame) for frame in frames]
    coord_arrays = [_coords_for_frame_and_ids(fmap, atom_ids) for fmap in frame_maps]

    rmsd_values: List[float] = []
    pair_rows: List[Tuple[int, int, float]] = []

    n_frames = len(coord_arrays)
    for i in range(n_frames):
        j_start = 0 if include_diagonal else i + 1
        for j in range(j_start, n_frames):
            rmsd = calculate_rmsd(coord_arrays[i], coord_arrays[j])
            rmsd_values.append(rmsd)
            pair_rows.append((i + 1, j + 1, rmsd))

            # For a strict all-vs-all matrix interpretation without diagonal,
            # include both directions (i,j) and (j,i) to keep pair symmetry.
            if not include_diagonal and i != j:
                rmsd_values.append(rmsd)
                pair_rows.append((j + 1, i + 1, rmsd))

    return np.asarray(rmsd_values, dtype=float), pair_rows, len(atom_ids)


def _normal_probability_curve(
    x: np.ndarray,
    mu: float,
    sigma: float,
    bin_width: float,
) -> np.ndarray:
    """Normal model in probability-per-bin scale.

    This converts PDF to bin probability using:
        probability_per_bin ~= bin_width * pdf(x)
    """
    sigma = max(float(sigma), 1e-8)
    if scipy_norm is not None:
        return bin_width * scipy_norm.pdf(x, loc=float(mu), scale=sigma)
    # Fallback formula if scipy.stats is unavailable.
    coeff = bin_width / (sigma * math.sqrt(2.0 * math.pi))
    z = (x - float(mu)) / sigma
    return coeff * np.exp(-0.5 * z * z)


# ============================================================================
# Fit model definitions
# ============================================================================


def _skew_normal_probability_curve(
    x: np.ndarray,
    alpha: float,
    mu: float,
    sigma: float,
    bin_width: float,
) -> np.ndarray:
    """Skew-normal model in probability-per-bin scale."""
    sigma = max(float(sigma), 1e-8)
    alpha = float(alpha)
    if scipy_skewnorm is not None:
        return bin_width * scipy_skewnorm.pdf(x, a=alpha, loc=float(mu), scale=sigma)

    # If skew-normal is unavailable, degrade to normal while preserving loc/scale.
    return _normal_probability_curve(x=x, mu=mu, sigma=sigma, bin_width=bin_width)


def _gaussian_mixture2_probability_curve(
    x: np.ndarray,
    weight: float,
    mu1: float,
    sigma1: float,
    mu2: float,
    sigma2: float,
    bin_width: float,
) -> np.ndarray:
    """Two-component Gaussian mixture in probability-per-bin scale."""
    w = min(max(float(weight), 0.0), 1.0)
    y1 = _normal_probability_curve(x=x, mu=mu1, sigma=sigma1, bin_width=bin_width)
    y2 = _normal_probability_curve(x=x, mu=mu2, sigma=sigma2, bin_width=bin_width)
    return (w * y1) + ((1.0 - w) * y2)


def _compute_bic(observed: np.ndarray, predicted: np.ndarray, k_params: int) -> float:
    """Compute Bayesian Information Criterion from residual error.

    Lower BIC is better. BIC penalizes extra parameters, reducing unnecessary
    preference for complex models.
    """
    n = max(int(observed.size), 1)
    rss = float(np.sum((observed - predicted) ** 2))
    rss = max(rss, 1e-15)
    return float(n * math.log(rss / n) + k_params * math.log(n))


def _fit_candidate_models(
    centers: np.ndarray,
    bin_probabilities: np.ndarray,
    rmsd_values: np.ndarray,
    bin_width: float,
) -> List[FitResult]:
    """Fit all candidate models and return successful fits.

    Why multiple models?
        A single Gaussian cannot represent all realistic RMSD distributions.
        This helper attempts a small model family and lets selection logic
        choose the most appropriate curve for the observed histogram.
    """
    mu_init = float(np.mean(rmsd_values))
    sigma_init = max(float(np.std(rmsd_values, ddof=0)), 1e-8)

    successful_fits: List[FitResult] = []

    # Candidate 1: single normal (unimodal, symmetric).
    # Good default for many compact conformational ensembles.
    try:
        if scipy_curve_fit is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                params, _ = scipy_curve_fit(
                    lambda x, mu, sigma: _normal_probability_curve(x, mu, sigma, bin_width),
                    centers,
                    bin_probabilities,
                    p0=[mu_init, sigma_init],
                    bounds=([0.0, 1e-8], [np.inf, np.inf]),
                    maxfev=30000,
                )
            y_hat = _normal_probability_curve(centers, float(params[0]), float(params[1]), bin_width)
            successful_fits.append(
                FitResult(
                    model_name="gaussian",
                    params=(float(params[0]), max(float(params[1]), 1e-8)),
                    bic=_compute_bic(bin_probabilities, y_hat, k_params=2),
                )
            )
    except Exception:
        pass

    # Candidate 2: skew-normal (unimodal, asymmetric/skewed).
    # Captures right/left tails when one side is stretched.
    try:
        if scipy_curve_fit is not None and scipy_skewnorm is not None:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                params, _ = scipy_curve_fit(
                    lambda x, alpha, mu, sigma: _skew_normal_probability_curve(
                        x, alpha, mu, sigma, bin_width
                    ),
                    centers,
                    bin_probabilities,
                    p0=[0.0, mu_init, sigma_init],
                    bounds=([-40.0, 0.0, 1e-8], [40.0, np.inf, np.inf]),
                    maxfev=60000,
                )
            y_hat = _skew_normal_probability_curve(
                centers,
                float(params[0]),
                float(params[1]),
                float(params[2]),
                bin_width,
            )
            successful_fits.append(
                FitResult(
                    model_name="skew_gaussian",
                    params=(float(params[0]), float(params[1]), max(float(params[2]), 1e-8)),
                    bic=_compute_bic(bin_probabilities, y_hat, k_params=3),
                )
            )
    except Exception:
        pass

    # Candidate 3: two-Gaussian mixture (can model bimodality).
    # Useful when two conformational basins are present.
    try:
        if scipy_curve_fit is not None:
            q30 = float(np.quantile(rmsd_values, 0.30))
            q70 = float(np.quantile(rmsd_values, 0.70))
            mix_init = [0.5, min(q30, q70), sigma_init, max(q30, q70), sigma_init]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                params, _ = scipy_curve_fit(
                    lambda x, w, mu1, s1, mu2, s2: _gaussian_mixture2_probability_curve(
                        x, w, mu1, s1, mu2, s2, bin_width
                    ),
                    centers,
                    bin_probabilities,
                    p0=mix_init,
                    bounds=([0.0, 0.0, 1e-8, 0.0, 1e-8], [1.0, np.inf, np.inf, np.inf, np.inf]),
                    maxfev=90000,
                )
            y_hat = _gaussian_mixture2_probability_curve(
                centers,
                float(params[0]),
                float(params[1]),
                float(params[2]),
                float(params[3]),
                float(params[4]),
                bin_width,
            )
            successful_fits.append(
                FitResult(
                    model_name="gaussian_mixture_2",
                    params=(
                        float(params[0]),
                        float(params[1]),
                        max(float(params[2]), 1e-8),
                        float(params[3]),
                        max(float(params[4]), 1e-8),
                    ),
                    bic=_compute_bic(bin_probabilities, y_hat, k_params=5),
                )
            )
    except Exception:
        pass

    return successful_fits


def fit_gaussian_only_to_histogram(
    rmsd_values: np.ndarray,
    bin_edges: np.ndarray,
    bin_probabilities: np.ndarray,
) -> FitResult:
    """Force a single-Gaussian fit model.

    This is useful when users want strict Gaussian behavior regardless of
    skew/bimodal alternatives.
    """
    if rmsd_values.size == 0:
        return FitResult(model_name="gaussian", params=(0.0, 1.0), bic=float("inf"))

    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = float(bin_edges[1] - bin_edges[0])
    mu_init = float(np.mean(rmsd_values))
    sigma_init = max(float(np.std(rmsd_values, ddof=0)), 1e-8)

    if scipy_curve_fit is not None:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", OptimizeWarning)
                params, _ = scipy_curve_fit(
                    lambda x, mu, sigma: _normal_probability_curve(x, mu, sigma, bin_width),
                    centers,
                    bin_probabilities,
                    p0=[mu_init, sigma_init],
                    bounds=([0.0, 1e-8], [np.inf, np.inf]),
                    maxfev=30000,
                )
            mu, sigma = float(params[0]), max(float(params[1]), 1e-8)
            y_hat = _normal_probability_curve(centers, mu, sigma, bin_width)
            return FitResult(
                model_name="gaussian",
                params=(mu, sigma),
                bic=_compute_bic(bin_probabilities, y_hat, k_params=2),
            )
        except Exception:
            pass

    # Fallback Gaussian from raw moments.
    mu = mu_init
    sigma = sigma_init
    return FitResult(model_name="gaussian", params=(mu, sigma), bic=float("inf"))


def _evaluate_fit_model(
    fit_result: FitResult,
    x: np.ndarray,
    bin_width: float,
) -> np.ndarray:
    """Evaluate selected model on x grid in probability-per-bin units."""
    if fit_result.model_name == "gaussian":
        mu, sigma = fit_result.params
        return _normal_probability_curve(x, mu, sigma, bin_width)
    if fit_result.model_name == "skew_gaussian":
        alpha, mu, sigma = fit_result.params
        return _skew_normal_probability_curve(x, alpha, mu, sigma, bin_width)
    if fit_result.model_name == "gaussian_mixture_2":
        w, mu1, s1, mu2, s2 = fit_result.params
        return _gaussian_mixture2_probability_curve(x, w, mu1, s1, mu2, s2, bin_width)
    raise ValueError(f"Unknown fit model: {fit_result.model_name}")


def _fit_to_summary_stats(
    fit_result: FitResult,
    bin_edges: np.ndarray,
) -> Tuple[float, float, float]:
    """Convert a model fit to summary metrics used in CSV reporting.

    Returns:
        (peak_amplitude, center_mu_like, sigma_like)
    """
    smooth_x = np.linspace(float(bin_edges[0]), float(bin_edges[-1]), 4000)
    bin_width = float(bin_edges[1] - bin_edges[0])
    smooth_y = _evaluate_fit_model(fit_result, smooth_x, bin_width)

    peak_idx = int(np.argmax(smooth_y))
    peak_amp = float(smooth_y[peak_idx])
    mu_like = float(smooth_x[peak_idx])

    y_sum = float(np.sum(smooth_y))
    if y_sum <= 0:
        return peak_amp, mu_like, 0.0

    weights = smooth_y / y_sum
    mean_val = float(np.sum(weights * smooth_x))
    variance = float(np.sum(weights * (smooth_x - mean_val) ** 2))
    sigma_like = float(math.sqrt(max(variance, 0.0)))
    return peak_amp, mu_like, sigma_like


def fit_distribution_to_histogram(
    rmsd_values: np.ndarray,
    bin_edges: np.ndarray,
    bin_probabilities: np.ndarray,
) -> FitResult:
    """
    Fit histogram with an adaptive model selection strategy.

    Strategy:
        1) Estimate initial parameters from raw RMSD values.
        2) Fit multiple candidate models:
           - Gaussian (symmetric, unimodal)
           - Skew-Gaussian (asymmetric, unimodal)
           - Two-Gaussian mixture (bimodal-capable)
        3) Select model with lowest BIC score.
        4) If fitting libraries are unavailable, use Gaussian fallback.

    This design improves shape fidelity for bimodal and skewed distributions
    while remaining backward-compatible for simple unimodal data.

    Beginner note:
        The output fit curve is not "forcing" one distribution type. Instead,
        the script evaluates several models and picks the best supported one.
    """
    if rmsd_values.size == 0:
        return FitResult(model_name="gaussian", params=(0.0, 1.0), bic=float("inf"))

    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = float(bin_edges[1] - bin_edges[0])
    successful_fits = _fit_candidate_models(
        centers=centers,
        bin_probabilities=bin_probabilities,
        rmsd_values=rmsd_values,
        bin_width=bin_width,
    )

    if successful_fits:
        # Primary choice by BIC.
        sorted_fits = sorted(successful_fits, key=lambda item: item.bic)
        best_fit = sorted_fits[0]

        # Guardrail: a bimodal mixture should only be accepted when it is
        # plausibly bimodal, not just slightly lower RSS from extra params.
        if best_fit.model_name == "gaussian_mixture_2":
            w, mu1, s1, mu2, s2 = best_fit.params
            bin_width = float(bin_edges[1] - bin_edges[0])

            separation = abs(float(mu1) - float(mu2))
            min_required_separation = 1.5 * max(float(s1), float(s2), bin_width)
            weight_is_reasonable = 0.08 <= float(w) <= 0.92

            alternative_fits = [fit for fit in sorted_fits if fit.model_name != "gaussian_mixture_2"]
            if alternative_fits:
                bic_improvement = float(alternative_fits[0].bic - best_fit.bic)
            else:
                bic_improvement = float("inf")

            # Three checks are required before accepting bimodality:
            # 1) peaks are sufficiently separated,
            # 2) neither component has tiny/degenerate weight,
            # 3) BIC improvement is meaningfully better than simpler model.
            passes_bimodal_checks = (
                separation >= min_required_separation
                and weight_is_reasonable
                and bic_improvement >= 8.0
            )

            if not passes_bimodal_checks and alternative_fits:
                return alternative_fits[0]

        return best_fit

    # Final fallback when curve-fitting backends are unavailable.
    mu = float(np.mean(rmsd_values))
    sigma = max(float(np.std(rmsd_values, ddof=0)), 1e-8)
    return FitResult(model_name="gaussian", params=(mu, sigma), bic=float("inf"))


def _bin_edges_from_values(
    values: np.ndarray,
    bin_width: float,
    forced_max: Optional[float] = None,
) -> np.ndarray:
    """Create histogram bin edges starting at 0 with fixed bin width."""
    min_val = 0.0
    max_val = float(np.max(values)) if values.size else bin_width
    if forced_max is not None:
        max_val = max(max_val, float(forced_max))
    upper = max(bin_width, math.ceil(max_val / bin_width) * bin_width)
    edges = np.arange(min_val, upper + bin_width, bin_width, dtype=float)
    if edges.size < 2:
        edges = np.array([0.0, bin_width], dtype=float)
    return edges


def _hist_probabilities(values: np.ndarray, bin_edges: np.ndarray) -> np.ndarray:
    """Convert raw values into probability mass per histogram bin.

    Sum of returned probabilities is approximately 1.0.
    """
    if values.size == 0:
        return np.zeros(bin_edges.size - 1, dtype=float)
    weights = np.full(values.shape, 1.0 / float(values.size), dtype=float)
    probs, _ = np.histogram(values, bins=bin_edges, weights=weights)
    return probs.astype(float)


def save_pairwise_csv(output_path: Path, rows: Iterable[Tuple[int, int, float]]) -> None:
    """Write pairwise RMSD table rows to CSV."""
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["FrameI", "FrameJ", "RMSD"])
        for i, j, rmsd in rows:
            writer.writerow([i, j, f"{rmsd:.6f}"])


def plot_single_histogram_with_fit(
    label: str,
    rmsd_values: np.ndarray,
    bin_edges: np.ndarray,
    bin_probs: np.ndarray,
    fit_result: Optional[FitResult],
    output_path: Path,
) -> None:
    """Render one comparison histogram with adaptive fit overlay.

    Plot semantics:
        - bars: empirical probability mass in each RMSD bin
        - red line: selected adaptive model prediction
    """
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)

    # Bars represent probability mass in each bin interval.
    plt.figure(figsize=(10, 6))
    plt.bar(centers, bin_probs, width=widths, color="tab:blue", edgecolor="black", alpha=0.65)

    if fit_result is not None:
        # Smooth x-grid for a visually continuous fit line.
        smooth_x = np.linspace(bin_edges[0], bin_edges[-1], 1000)
        smooth_y = _evaluate_fit_model(
            fit_result=fit_result,
            x=smooth_x,
            bin_width=float(bin_edges[1] - bin_edges[0]),
        )
        fit_label = f"{fit_result.model_name} fit"
        plt.plot(smooth_x, smooth_y, color="tab:red", linewidth=2.0, label=fit_label)

    plt.title(f"RMSD Probability Histogram with Adaptive Fit ({label})")
    plt.xlabel("RMSD (Angstrom)")
    plt.ylabel("Probability per bin")
    plt.xlim(float(bin_edges[0]), float(bin_edges[-1]))
    plt.grid(alpha=0.2)
    if fit_result is not None:
        plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def plot_combined_histograms_with_fit(
    histogram_payloads: Sequence[Tuple[str, np.ndarray, np.ndarray, Optional[FitResult]]],
    output_path: Path,
) -> None:
    """Overlay multiple probability histograms + adaptive fits on one chart.

    This makes it easy to compare distribution shifts between full/backbone/
    custom selections in a single visual.
    """
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:cyan"]

    plt.figure(figsize=(12, 7))
    for idx, (label, bin_edges, bin_probs, fit_result) in enumerate(histogram_payloads):
        color = colors[idx % len(colors)]
        centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        widths = np.diff(bin_edges)

        plt.bar(
            centers,
            bin_probs,
            width=widths,
            alpha=0.25,
            color=color,
            edgecolor=color,
            label=f"{label} histogram",
        )

        if fit_result is not None:
            smooth_x = np.linspace(bin_edges[0], bin_edges[-1], 1000)
            smooth_y = _evaluate_fit_model(
                fit_result=fit_result,
                x=smooth_x,
                bin_width=float(bin_edges[1] - bin_edges[0]),
            )
            plt.plot(
                smooth_x,
                smooth_y,
                color=color,
                linewidth=2.0,
                label=f"{label} {fit_result.model_name} fit",
            )

    plt.title("RMSD Probability Histograms with Adaptive Fits")
    plt.xlabel("RMSD (Angstrom)")
    plt.ylabel("Probability per bin")
    if histogram_payloads:
        first_edges = histogram_payloads[0][1]
        plt.xlim(float(first_edges[0]), float(first_edges[-1]))
    plt.grid(alpha=0.2)
    plt.legend(fontsize=9, ncol=2)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def save_fit_summary_csv(output_path: Path, summaries: Sequence[FitSummary]) -> None:
    """Write per-comparison fit summary statistics to CSV.

    Important columns:
        - FitModel: selected model family
        - FitBIC: model quality score (lower is better)
        - FitAmplitude/FitMu/FitSigma: standardized summary descriptors of
          the selected curve for easier downstream comparison tables.
    """
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "Comparison",
                "Frames",
                "Pairs",
                "CommonAtoms",
                "MeanRMSD",
                "StdRMSD",
                "FitModel",
                "FitBIC",
                "FitAmplitude",
                "FitMu",
                "FitSigma",
            ]
        )
        for row in summaries:
            writer.writerow(
                [
                    row.comparison,
                    row.n_frames,
                    row.n_pairs,
                    row.n_atoms_common,
                    f"{row.mean_rmsd:.6f}",
                    f"{row.std_rmsd:.6f}",
                    row.fit_model,
                    f"{row.fit_bic:.6f}" if math.isfinite(row.fit_bic) else "inf",
                    f"{row.fit_amplitude:.6f}",
                    f"{row.fit_mu:.6f}",
                    f"{row.fit_sigma:.6f}",
                ]
            )


def build_comparison_list(
    base_comparisons: Sequence[str],
    custom_map: Dict[str, Set[str]],
) -> List[str]:
    """Expand comparison plan to include custom labels when requested.

    Example:
        base_comparisons = ["full", "custom"]
        custom_map keys = ["Polar", "SideChain"]
        output = ["full", "Polar", "SideChain"]
    """
    plan: List[str] = []
    for item in base_comparisons:
        if item == "custom":
            plan.extend(sorted(custom_map.keys()))
        else:
            plan.append(item)

    seen: Set[str] = set()
    ordered_unique: List[str] = []
    for name in plan:
        if name not in seen:
            seen.add(name)
            ordered_unique.append(name)
    return ordered_unique


def main() -> None:
    """Main orchestration:

    1) Parse and validate CLI options.
    2) Parse PDB into frames.
     3) For each requested comparison:
       - compute all-vs-all RMSD
       - histogram using requested bin width
         - adaptive distribution fit on probability histogram
       - write plot (and optional pairwise CSV)
    4) Write summary CSV and optional combined overlay plot.

     Practical note:
          For very large trajectories, all-vs-all RMSD can be computationally
          expensive because pair count scales approximately with N^2.
    """
    args = parse_args()

    pdb_path = Path(args.pdb)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if args.bin_width <= 0:
        raise ValueError("--bin-width must be > 0")
    if args.x_max_angstrom is not None and args.x_max_angstrom <= 0:
        raise ValueError("--x-max-angstrom must be > 0 when provided")

    custom_map = parse_custom_selections(args.custom_selection)
    if "custom" in args.comparisons and not custom_map:
        raise ValueError(
            "You selected --comparisons custom but did not provide --custom-selection."
        )

    comparison_labels = build_comparison_list(args.comparisons, custom_map)
    if not comparison_labels:
        raise ValueError("No comparisons were selected.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading PDB: {pdb_path}")
    frames = parse_pdb_frames(pdb_path)
    if len(frames) < 2:
        raise ValueError("Need at least 2 frames in the input PDB to compute all-vs-all RMSD.")
    print(f"[INFO] Parsed frames: {len(frames)}")

    include_diagonal = args.include_diagonal == "Y"
    write_pairwise = args.write_pairwise_csv == "Y"

    # `summary_rows` captures scalar metrics for downstream tables/reports.
    summary_rows: List[FitSummary] = []

    # `combined_payloads` stores plotting payloads for final overlay figure.
    combined_payloads: List[Tuple[str, np.ndarray, np.ndarray, Optional[FitResult]]] = []
    comparison_results: List[ComparisonResult] = []
    global_max_rmsd = 0.0

    for label in comparison_labels:
        print(f"[INFO] Computing comparison: {label}")
        rmsd_values, pair_rows, n_common_atoms = compute_all_vs_all_rmsd(
            frames=frames,
            selection_name=label,
            custom_selection_map=custom_map,
            include_diagonal=include_diagonal,
        )
        local_max = float(np.max(rmsd_values)) if rmsd_values.size else 0.0
        global_max_rmsd = max(global_max_rmsd, local_max)
        comparison_results.append(
            ComparisonResult(
                label=label,
                rmsd_values=rmsd_values,
                pair_rows=pair_rows,
                n_common_atoms=n_common_atoms,
            )
        )

    x_axis_max = args.x_max_angstrom if args.x_max_angstrom is not None else global_max_rmsd
    shared_bin_edges = _bin_edges_from_values(
        values=np.array([global_max_rmsd], dtype=float),
        bin_width=args.bin_width,
        forced_max=x_axis_max,
    )
    print(
        "[INFO] Histogram x-axis range: "
        f"{shared_bin_edges[0]:.3f} to {shared_bin_edges[-1]:.3f} Angstrom "
        f"(bin width {args.bin_width:.3f})."
    )

    # Second pass: fit and plot each comparison using the shared x-axis/bin grid.
    for result in comparison_results:
        label = result.label
        rmsd_values = result.rmsd_values
        pair_rows = result.pair_rows
        n_common_atoms = result.n_common_atoms

        # Histogram bins are fixed-width (default 0.2 Angstrom) and shared
        # across comparisons so no high-RMSD range is omitted in any plot.
        bin_edges = shared_bin_edges
        bin_probs = _hist_probabilities(rmsd_values, bin_edges)

        fit_result: Optional[FitResult]
        if args.fit_mode == "adaptive":
            fit_result = fit_distribution_to_histogram(
                rmsd_values=rmsd_values,
                bin_edges=bin_edges,
                bin_probabilities=bin_probs,
            )
            print(
                f"[INFO] Selected fit model for {label}: "
                f"{fit_result.model_name} (BIC={fit_result.bic:.3f})"
            )
        elif args.fit_mode == "gaussian":
            fit_result = fit_gaussian_only_to_histogram(
                rmsd_values=rmsd_values,
                bin_edges=bin_edges,
                bin_probabilities=bin_probs,
            )
            print(
                f"[INFO] Forced fit model for {label}: "
                f"{fit_result.model_name} (BIC={fit_result.bic:.3f})"
            )
        else:
            fit_result = None
            print(f"[INFO] Fit disabled for {label} (fit-mode=none).")

        if fit_result is not None:
            fit_amp, fit_mu, fit_sigma = _fit_to_summary_stats(fit_result, bin_edges)
            fit_model = fit_result.model_name
            fit_bic = fit_result.bic
        else:
            fit_amp, fit_mu, fit_sigma = 0.0, 0.0, 0.0
            fit_model = "none"
            fit_bic = float("inf")

        out_plot = output_dir / f"rmsd_hist_gaussian_{label}.png"
        plot_single_histogram_with_fit(
            label=label,
            rmsd_values=rmsd_values,
            bin_edges=bin_edges,
            bin_probs=bin_probs,
            fit_result=fit_result,
            output_path=out_plot,
        )
        print(f"[INFO] Saved plot: {out_plot}")

        if write_pairwise:
            out_pairs = output_dir / f"rmsd_pairs_{label}.csv"
            save_pairwise_csv(out_pairs, pair_rows)
            print(f"[INFO] Saved pairwise table: {out_pairs}")

        summary_rows.append(
            FitSummary(
                comparison=label,
                n_frames=len(frames),
                n_pairs=len(rmsd_values),
                n_atoms_common=n_common_atoms,
                mean_rmsd=float(np.mean(rmsd_values)),
                std_rmsd=float(np.std(rmsd_values, ddof=0)),
                fit_model=fit_model,
                fit_bic=fit_bic,
                fit_amplitude=fit_amp,
                fit_mu=fit_mu,
                fit_sigma=fit_sigma,
            )
        )
        combined_payloads.append((label, bin_edges, bin_probs, fit_result))

    summary_csv = output_dir / "rmsd_gaussian_fit_summary.csv"
    save_fit_summary_csv(summary_csv, summary_rows)
    print(f"[INFO] Saved fit summary: {summary_csv}")

    # Combined plot is meaningful only when more than one comparison exists.
    if len(combined_payloads) > 1:
        combined_plot = output_dir / "rmsd_hist_gaussian_combined.png"
        plot_combined_histograms_with_fit(combined_payloads, combined_plot)
        print(f"[INFO] Saved combined plot: {combined_plot}")

    print("[DONE] RMSD histogram + Gaussian fitting complete.")


if __name__ == "__main__":
    # Standard Python entry-point guard.
    main()

