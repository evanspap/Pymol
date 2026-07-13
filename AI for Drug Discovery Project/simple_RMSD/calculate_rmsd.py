"""
RMSD Calculator for Multi-Frame PDB Files
Calculates RMSD of each frame relative to frame 1 and frame 2
"""

import numpy as np
from pathlib import Path
import sys

def parse_pdb_frames(pdb_file):
    """
    Parse a multi-frame PDB file and extract individual frames.
    Each frame is separated by 'END' keyword.
    
    Args:
        pdb_file: Path to the PDB file
        
    Returns:
        List of frames, where each frame contains coordinates
    """
    frames = []
    current_frame = []
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        # Parse ATOM/HETATM line
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        current_frame.append([x, y, z])
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse line: {line.strip()}")
                        continue
                        
                elif line.startswith('END'):
                    if current_frame:
                        frames.append(np.array(current_frame))
                        current_frame = []
        
        # Handle case where file doesn't end with END
        if current_frame:
            frames.append(np.array(current_frame))
    
    except FileNotFoundError:
        print(f"Error: File not found: {pdb_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    return frames

def calculate_rmsd(coords1, coords2):
    """
    Calculate RMSD between two coordinate sets.
    
    Args:
        coords1: Nx3 array of coordinates
        coords2: Nx3 array of coordinates
        
    Returns:
        RMSD value
    """
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    
    diff = coords1 - coords2
    msd = np.mean(np.sum(diff**2, axis=1))
    rmsd = np.sqrt(msd)
    
    return rmsd

def get_rmsd_category(rmsd):
    """Return categorical label for a given RMSD value."""
    if rmsd <= 2.0:
        return "<=2 Å"
    if rmsd <= 3.0:
        return "2-3 Å"
    if rmsd <= 4.0:
        return "3-4 Å"
    if rmsd <= 5.0:
        return "4-5 Å"
    return ">5 Å"


def main():
    pdb_file = Path(r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mQucGbPypmGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k\apoa__run1_frames_1k.pdb")
    
    if not pdb_file.exists():
        print(f"Error: PDB file not found at {pdb_file}")
        sys.exit(1)
    
    print(f"Reading PDB file: {pdb_file}")
    print("Parsing frames...")
    
    frames = parse_pdb_frames(str(pdb_file))
    print(f"Successfully parsed {len(frames)} frames")
    
    if len(frames) < 2:
        print("Error: Need at least 2 frames to calculate RMSD")
        sys.exit(1)
    
    frame1 = frames[0]
    frame2 = frames[1]
    
    print(f"\nFrame 1 has {len(frame1)} atoms")
    print(f"Frame 2 has {len(frame2)} atoms")
    
    print("\n" + "="*60)
    print("RMSD Calculation: Frame 1 vs Frame 2")
    print("="*60)
    print(f"{'Frame':<8} {'RMSD (Å)':<15} {'Description':<30}")
    print("-"*60)
    
    results = []
    category_counts = {}
    
    for i, frame in enumerate(frames, 1):
        for target_name, target_frame in [("Frame 1", frame1), ("Frame 2", frame2)]:
            rmsd = calculate_rmsd(frame, target_frame)
            category = get_rmsd_category(rmsd)
            
            print(f"{i:<8} {target_name:<12} {rmsd:<12.4f} {category:<12}")
            
            results.append({
                'frame': i,
                'compared_to': target_name,
                'rmsd': rmsd,
                'category': category
            })
            
            category_counts[category] = category_counts.get(category, 0) + 1
    
    output_csv = Path(__file__).parent / "rmsd_results.csv"
    print(f"\nSaving categorized results to: {output_csv}")
    
    with open(output_csv, 'w') as f:
        f.write("Frame,Compared_To,RMSD,Category\n")
        for result in results:
            f.write(f"{result['frame']},{result['compared_to']},{result['rmsd']:.6f},{result['category']}\n")
    
    print("\nCategory counts:")
    category_order = ["<=2 Å", "2-3 Å", "3-4 Å", "4-5 Å", ">5 Å"]
    for category in sorted(category_counts.keys(), key=lambda x: category_order.index(x) if x in category_order else len(category_order)):
        print(f"  {category}: {category_counts[category]}")

    # Generate an all-vs-all CSV where each block is comparisons to a single target frame
    def generate_all_vs_all(frames, out_path):
        out_path = Path(out_path)
        print(f"\nGenerating all-vs-all RMSD CSV: {out_path}")
        with open(out_path, 'w') as f:
            f.write("TargetFrame,Frame,RMSD,Category\n")
            for t_idx, t_frame in enumerate(frames, 1):
                # Optional comment line to mark blocks for human readers
                f.write(f"# TargetFrame {t_idx}\n")
                # Per-target category counts
                per_cat = {}
                for i_idx, frame in enumerate(frames, 1):
                    rmsd = calculate_rmsd(frame, t_frame)
                    category = get_rmsd_category(rmsd)
                    f.write(f"{t_idx},{i_idx},{rmsd:.6f},{category}\n")
                    per_cat[category] = per_cat.get(category, 0) + 1
                # Blank line to separate blocks
                f.write("\n")
                # Also print a small summary to console
                print(f"TargetFrame {t_idx}: ", end='')
                print(", ".join([f"{c}: {per_cat.get(c,0)}" for c in category_order]))

    all_vs_all_csv = Path(__file__).parent / "rmsd_all_vs_all.csv"
    generate_all_vs_all(frames, all_vs_all_csv)

    print("\nDone!")

if __name__ == "__main__":
    main()
