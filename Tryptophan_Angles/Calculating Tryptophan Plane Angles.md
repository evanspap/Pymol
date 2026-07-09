# Calculating Tryptophan Plane Angles

This guide describes how to compute the angle between two tryptophan indole ring planes in 3D space.

## I. Define Each Indole Plane

A plane is defined by three non-collinear points. For each tryptophan residue, use atom coordinates for:

- $C_{\gamma}$
- $C_{\delta 2}$
- $N_{\epsilon 1}$

Define two in-plane vectors for each residue:

$$
\vec{A} = \mathbf{r}_{C_{\delta 2}} - \mathbf{r}_{C_{\gamma}}
$$

$$
\vec{B} = \mathbf{r}_{N_{\epsilon 1}} - \mathbf{r}_{C_{\gamma}}
$$

## II. Compute Normal Vectors

The plane normal vector is:

$$
\vec{n} = \vec{A} \times \vec{B}
$$

Component form:

$$
n_x = A_yB_z - A_zB_y
$$

$$
n_y = A_zB_x - A_xB_z
$$

$$
n_z = A_xB_y - A_yB_x
$$

Compute this for both residues to obtain $\vec{n}_1$ and $\vec{n}_2$.

## III. Compute the Inter-Planar Angle

The angle between planes equals the angle between their normals:

$$
\cos(\theta) = \frac{|\vec{n}_1 \cdot \vec{n}_2|}{\|\vec{n}_1\|\,\|\vec{n}_2\|}
$$

Expanded form:

$$
\theta = \arccos\left(\frac{|n_{1x}n_{2x} + n_{1y}n_{2y} + n_{1z}n_{2z}|}{\sqrt{n_{1x}^2+n_{1y}^2+n_{1z}^2}\,\sqrt{n_{2x}^2+n_{2y}^2+n_{2z}^2}}\right)
$$

Notes:

- Use $|\vec{n}_1 \cdot \vec{n}_2|$ to report the acute angle ($0^\circ$ to $90^\circ$).
- Omit the absolute value if you need signed orientation information up to $180^\circ$.

## IV. Structural Interpretation

| Angle ($\theta$) | Geometry | Description |
|---|---|---|
| $0^\circ$ or $180^\circ$ | Parallel | Face-to-face or anti-parallel ring orientation |
| $90^\circ$ | Perpendicular | T-shaped (edge-to-face) orientation |
| $10^\circ$ to $30^\circ$ | Offset | Tilted/offset stacking, often observed in protein cores |


