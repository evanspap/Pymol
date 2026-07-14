from pathlib import Path
import os

paths = [
    r"G:\.shortcut-targets-by-id\1cfLzEn1DaVCZwnrRp5mGqBN\AI Drug Discovery for Cancer\MD simulations\Retuns_081424\Apo-A\Run1\frames_1k",
]
for p in paths:
    print('PATH:', p)
    print('  Path.exists:', Path(p).exists())
    print('  Path.is_dir:', Path(p).is_dir())
    print('  os.path.exists:', os.path.exists(p))
    print('  os.path.isdir:', os.path.isdir(p))
    try:
        print('  listdir first 5:', os.listdir(p)[:5])
    except Exception as exc:
        print('  listdir error:', exc)
    print()
