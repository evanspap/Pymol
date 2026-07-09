from pymol import cmd
import numpy as np

def resi_rms(resi, atomsel="name CA"):
    coords = []

    for obj in cmd.get_object_list():
        sel = f"{obj} and resi {resi} and ({atomsel})"

        try:
            c = cmd.get_coords(sel)
            if c is not None:
                coords.extend(c)
        except:
            pass

    coords = np.array(coords)

    if len(coords) < 2:
        print("Not enough atoms found.")
        return

    centroid = coords.mean(axis=0)
    rms = np.sqrt(np.mean(np.sum((coords-centroid)**2, axis=1)))

    print(f"Residue {resi}: RMS = {rms:.3f} Å")

cmd.extend("resi_rms", resi_rms)