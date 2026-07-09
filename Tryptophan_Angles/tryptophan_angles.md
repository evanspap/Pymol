To measure the angle between the planes of two tryptophan indole rings, the most robust method is to calculate the angle between their **normal vectors**. A normal vector is a vector sticking straight out of the plane (perpendicular to it).

Below is the complete mathematical and computational framework, including the missing formulas in both **LaTeX** and **ASCII** formats.

---

### I. Mathematical Framework (LaTeX)

To find the angle between two planes, we first define the planes using atomic coordinates, then find the normal vectors, and finally compute the angle between those normals.

#### 1. Define In-Plane Vectors
For each Tryptophan ($Trp_1$ and $Trp_2$), select three atoms that define the indole ring: $C\gamma$ (CG), $C\delta2$ (CD2), and $N\epsilon1$ (NE1).
For $Trp_n$, define two vectors $\vec{u}_n$ and $\vec{v}_n$:
$$\vec{u}_n = \vec{r}_{CD2} - \vec{r}_{CG}$$
$$\vec{v}_n = \vec{r}_{NE1} - \vec{r}_{CG}$$

#### 2. Calculate the Normal Vector ($\vec{n}$)
The normal vector is the cross product of the two in-plane vectors. This vector is perpendicular to the indole ring.
$$\vec{n}_n = \vec{u}_n \times \vec{v}_n = 
\begin{vmatrix} 
\mathbf{i} & \mathbf{j} & \mathbf{k} \\ 
u_x & u_y & u_z \\ 
v_x & v_y & v_z 
\end{vmatrix}$$

The components of $\vec{n} = (n_x, n_y, n_z)$ are:
$$n_x = (u_y v_z - u_z v_y)$$
$$n_y = (u_z v_x - u_x v_z)$$
$$n_z = (u_x v_y - u_y v_x)$$

#### 3. Calculate the Inter-Planar Angle ($\theta$)
The angle between the two planes is the angle between $\vec{n}_1$ and $\vec{n}_2$. We use the dot product formula:
$$\cos(\theta) = \frac{|\vec{n}_1 \cdot \vec{n}_2|}{\|\vec{n}_1\| \|\vec{n}_2\|}$$

The full expansion in Cartesian coordinates is:
$$\theta = \arccos \left( \frac{|n_{1x}n_{2x} + n_{1y}n_{2y} + n_{1z}n_{2z}|}{\sqrt{n_{1x}^2 + n_{1y}^2 + n_{1z}^2} \cdot \sqrt{n_{2x}^2 + n_{2y}^2 + n_{2z}^2}} \right)$$

---

### II. Computational Logic (ASCII)

This format is ideal for script comments or plain-text documentation.

```text
1. VECTOR DEFINITION (Per Tryptophan)
   Vector_U = [x(CD2)-x(CG), y(CD2)-y(CG), z(CD2)-z(CG)]
   Vector_V = [x(NE1)-x(CG), y(NE1)-y(CG), z(NE1)-z(CG)]

2. NORMAL VECTOR CALCULATION (Cross Product)
   n_x = (Uy * Vz) - (Uz * Vy)
   n_y = (Uz * Vx) - (Ux * Vz)
   n_z = (Ux * Vy) - (Uy * Vx)

3. DOT PRODUCT FORMULA
                   | (n1_x * n2_x) + (n1_y * n2_y) + (n1_z * n2_z) |
   cos(theta) = -------------------------------------------------------
                 sqrt(n1_x^2+n1_y^2+n1_z^2) * sqrt(n2_x^2+n2_y^2+n2_z^2)

4. FINAL ANGLE
   theta = arccos(cos_theta) * (180 / PI)  [to get degrees]
```

---

### III. Python Implementation (NumPy)

This script automates the process using the logic defined above.

```python
import numpy as np

def calculate_plane_angle(coords_trp1, coords_trp2):
    """
    Calculates the angle between two tryptophan indole planes.
    coords_trp: Dictionary or array containing coordinates for 'CG', 'CD2', 'NE1'
    """
    
    def get_normal(atoms):
        # 1. Define two vectors from 3 atoms
        u = atoms['CD2'] - atoms['CG']
        v = atoms['NE1'] - atoms['CG']
        
        # 2. Calculate normal vector via cross product
        n = np.cross(u, v)
        return n

    # Get normals for both rings
    n1 = get_normal(coords_trp1)
    n2 = get_normal(coords_trp2)

    # 3. Calculate the angle using the dot product formula
    dot_product = np.dot(n1, n2)
    mag1 = np.linalg.norm(n1)
    mag2 = np.linalg.norm(n2)
    
    # Use absolute value for plane-to-plane angle (0-90 degrees)
    cos_theta = np.abs(dot_product) / (mag1 * mag2)
    
    # Clip to avoid errors from floating point precision
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    
    angle_rad = np.arccos(cos_theta)
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg

# Example usage with dummy coordinates
trp1 = {'CG': np.array([1.0, 0.0, 0.0]), 'CD2': np.array([0.0, 1.0, 0.0]), 'NE1': np.array([0.0, 0.0, 0.0])}
trp2 = {'CG': np.array([1.0, 0.0, 1.0]), 'CD2': np.array([0.0, 1.0, 1.0]), 'NE1': np.array([0.0, 0.0, 1.0])}

print(f"Angle between planes: {calculate_plane_angle(trp1, trp2):.2f} degrees")
```

---

### IV. Interpretation of Results

*   **0° (or 180°):** The planes are **Parallel**. If they are close together, this indicates **$\pi$-$\pi$ stacking** (face-to-face).
*   **90°:** The planes are **Perpendicular**. This indicates a **T-shaped** edge-to-face orientation, common in stabilizing protein cores.
*   **45°:** The planes are at an oblique angle, often found in "offset" or "slipped" stacking geometries.