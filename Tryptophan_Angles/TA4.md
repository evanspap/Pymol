The best way to measure the plane angles between two tryptophan side chains is to calculate the angle between their **normal vectors**, which are vectors perpendicular to each indole ring. This process involves using the 3D atomic coordinates of the tryptophan residues to define the planes and then applying vector mathematics to find the angle between them.

Mathematical Framework

The angle (

![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

) between two planes is defined as the angle between their normal vectors (

![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

and

![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

). This is often referred to in geometry as the **dihedral angle** between the planes.

- **Formula:** The cosine of the angle is found using the dot product of the normal vectors and their magnitudes:

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)
- **Alternative (Cartesian) Form:** If the planes are defined by equations

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  and

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  , the formula is:

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

Computational Methods and Tools

To perform these measurements on protein structures (e.g., from a PDB file), several software tools and libraries can automate the coordinate extraction and calculation:

- **MDAnalysis:** This is a recommended library for analyzing Molecular Dynamics (MD) simulations. It can be used to extract the coordinates of the indole ring atoms, define vectors representing the rings, and calculate normal vectors for every frame of a simulation to track angular progression over time.
- **EDPDB:** A multi-functional tool for protein structure analysis that can determine the "best-fit plane" through a set of points (atoms) and extract various geometrical information, including bond lengths and torsion angles.
- **PyMOL:** While PyMOL's standard measurement wizard can calculate dihedral angles between atoms, measuring the angle between *planes* often requires custom scripts. There are user-provided scripts available via the PyMOL Wiki specifically for measuring angles between secondary structures or defined planes.
- **YASARA:** A visualization and analysis tool that includes simple built-in capabilities for making measurements and comparisons of biomolecular structures.

Structural Considerations for Tryptophan

- **Side Chain Geometry:** The tryptophan (Trp) side chain is an **indole ring**, a fused six-membered benzene ring and a five-membered pyrrole ring. Because it is bulky and planar, its orientation is defined by side-chain dihedral angles

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  and

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  .
- **Steric Hindrance:** Free rotation of the indole plane is restricted by the protein backbone. Sterically hindered ranges for

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  often occur around

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  ,

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  , and

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  , while favored ranges are typically between

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  –

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  and

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  –

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  .
- **Difference from Dihedral Angles:** Standard protein dihedral angles (

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  ,

  ![](data:image/gif;base64,R0lGODlhAQABAIAAAP///wAAACH5BAEAAAAALAAAAAABAAEAAAICRAEAOw==)

  ) measure rotations around backbone bonds. Measuring the plane angle of two separate tryptophans requires defining their individual indole planes independently, rather than just measuring a single torsion angle between four atoms.

Would you like me to help you find or draft a Python script using **MDAnalysis** to calculate these plane angles from a PDB file?