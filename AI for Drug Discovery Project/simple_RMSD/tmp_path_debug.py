from pathlib import Path
import os

def show(path):
    p = Path(path)
    print('PATH:', path)
    print('exists:', p.exists())
    print('is_dir:', p.is_dir())
    try:
        print('listdir first 5:', os.listdir(path)[:5])
    except Exception as e:
        print('listdir error:', repr(e))

show('G:/')
show('G:/.shortcut-targets-by-id')
show('G:/.shortcut-targets-by-id/1cfLzEn1DaVCZwnrRp5mGqBN')
show('G:/.shortcut-targets-by-id/1cfLzEn1DaVCZwnrRp5mGqBN/AI Drug Discovery for Cancer')
show('G:/.shortcut-targets-by-id/1cfLzEn1DaVCZwnrRp5mGqBN/AI Drug Discovery for Cancer/MD simulations')
show('G:/.shortcut-targets-by-id/1cfLzEn1DaVCZwnrRp5mGqBN/AI Drug Discovery for Cancer/MD simulations/Retuns_081424/Apo-A/Run1/frames_1k')
