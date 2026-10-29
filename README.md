<p align="center">
  <img src="assets/intense-logo.png" width="150" alt="INTENSE logo" />
</p>

# RCoR — Riemannian Co-Registration

**RCoR** is a method for the registration of real symmetric positive-definite (SPD) matrices through a common eigenvector basis. This repository provides implementations in **Julia** and **Python**, using eigenvector sign synchronization followed by nuclear-norm optimization.

## Mathematical formulation

Let $C_1,\ldots,C_n\in\mathcal S_{++}^{d}$ be SPD matrices, with eigendecompositions

$$
C_i=U_i\Lambda_i U_i^{\mathsf T},
\qquad U_i\in O(d),\qquad
\Lambda_i=\operatorname{diag}(\lambda_{i1},\ldots,\lambda_{id}).
$$

The eigenvalues are ordered from largest to smallest. Each eigenvector is defined only up to its sign, so the matrix $U_i$ can be replaced by $U_i S_i$, where $S_i=\operatorname{diag}(s_{i1},\ldots,s_{id})$ and $s_{ij}\in\{-1,+1\}$, without changing $C_i$.

The aim is to find a common orthogonal basis while resolving these sign ambiguities. Define the **signed sum**

$$
A(S)=\sum_{i=1}^{n}U_i S_i,
\qquad S_1=I,
$$

where $S_1=I$ fixes the arbitrary global sign convention. The sign selection seeks to maximize the **nuclear norm** of the signed sum:

$$
\boxed{\ \max_{S_2,\ldots,S_n}\ \|A(S)\|_*
=\max_{S_2,\ldots,S_n}\sum_{k=1}^{d}\sigma_k(A(S))\ }.
$$

For fixed signs, the closest common orthogonal matrix in the Frobenius sense is the polar factor of $A(S)$. Indeed,

$$
\min_{\Gamma\in O(d)}\sum_{i=1}^{n}
\|U_i S_i-\Gamma\|_F^2
=2nd-2\|A(S)\|_*.
$$

### Sign synchronization and refinement

The implementation selects signs in two stages:

1. **Column-wise initialization.** For each eigenvector index $j$, it enumerates all $2^{n-1}$ sign patterns (holding $s_{1j}=+1$ fixed) and maximizes $\|A(S)_{:,j}\|_2^2$. As the squared Frobenius norm separates by columns, this initialization maximizes $\|A(S)\|_F^2$ exactly.
2. **Nuclear-norm refinement.** Starting from those signs, the algorithm visits one sign at a time and accepts a flip only if it increases $\|A(S)\|_*$ by more than `tol`. It stops when no improving flip remains, or when `maxiter` passes have been performed.

The second stage is a greedy coordinate-ascent procedure: it does **not** generally guarantee a global maximum of the nuclear norm. The exhaustive initialization is exponential in $n$, so it is best suited to modest collections of SPD matrices.

### Common basis and registered matrices

Let $A(S)=P\Sigma Q^{\mathsf T}$ be an SVD after sign optimization. The common orthogonal basis is

$$
\Gamma=P Q^{\mathsf T}.
$$

The registered matrices retain the original eigenvalues but share $\Gamma$:

$$
\boxed{\ G_i=\Gamma\Lambda_i\Gamma^{\mathsf T}\ }.
$$

The reference SPD matrix is obtained by taking the componentwise geometric mean of the eigenvalues:

$$
\ell_j=\exp\!\left(\frac1n\sum_{i=1}^{n}\log\lambda_{ij}\right),
\qquad
\boxed{\ M=\Gamma\operatorname{diag}(\ell_1,\ldots,\ell_d)\Gamma^{\mathsf T}\ }.
$$

The implementations return **`M`**, **`Gamma` / `Γ`**, and the collection of **registered matrices `G` / `𝐆`**.

## Installation

### Python

Requires Python 3.10 or newer. From the repository root:

```bash
python -m pip install -e ".[plot]"
```

The RCoR algorithm uses NumPy; the plots also require SciPy and Matplotlib.

### Julia

Install the dependencies in the Julia REPL:

```julia
using Pkg
Pkg.add(["PosDefManifold", "Diagonalizations", "Plots"])
```

The Julia implementation is provided as a source file that can be included directly.

## Usage

### Python

```python
import numpy as np
from rcor_nuclear import rr

# Collection of n real SPD matrices, each of size d x d:
C = np.stack([
    np.diag([5.0, 2.0, 0.5]),
    np.array([[3.0, 0.4, 0.0], [0.4, 2.0, 0.2], [0.0, 0.2, 0.9]]),
])

M, Gamma, G = rr(C)
print(M.shape)      # (3, 3)
print(Gamma.shape)  # (3, 3)
print(G.shape)      # (2, 3, 3)
```

Main source: [`src/rcor_nuclear/core.py`](src/rcor_nuclear/core.py). Optional arguments are `tol=1e-12` and `maxiter=100`.

### Julia

```julia
using LinearAlgebra, PosDefManifold, Diagonalizations
include("julia/RR_nuclear_optimized_final.jl")

C = ℍVector([
    ℍ(Matrix(Diagonal([5.0, 2.0, 0.5]))),
    ℍ([3.0 0.4 0.0; 0.4 2.0 0.2; 0.0 0.2 0.9]),
])

M, Γ, 𝐆 = rr(C)
```

Main source: [`julia/RR_nuclear_optimized_final.jl`](julia/RR_nuclear_optimized_final.jl). The Julia interface uses `ℍ` and `ℍVector` from PosDefManifold.jl.

## Example: RCoR scatterplot

The two examples read the **same set of 12 synthetic SPD matrices (8 × 8)** and plot four groups in a shared 2D representation:

- Original SPD matrices.
- Registered SPD matrices.
- Arithmetic mean of the original matrices.
- RCoR reference matrix $M$.

The visualization uses **classical multidimensional scaling (MDS)** of pairwise affine-invariant (Fisher) distances between SPD matrices. It is a 2D illustration of their relationships, not a direct plot of individual eigenvalues.

The same synthetic input is processed independently in Julia and Python. The panels show the original SPD matrices (circles), the registered matrices (squares), the original arithmetic mean (diamond), and the RCoR reference $M$ (star).

<table>
  <tr>
    <th align="center">Julia</th>
    <th align="center">Python</th>
  </tr>
  <tr>
    <td width="50%" align="center">
      <img src="assets/scatter_example_julia.png" width="100%" alt="Julia RCoR example: raw SPD matrices, registered SPD matrices, arithmetic mean, and reference M" />
    </td>
    <td width="50%" align="center">
      <img src="assets/scatter_example_python.png" width="100%" alt="Python RCoR example: raw SPD matrices, registered SPD matrices, arithmetic mean, and reference M" />
    </td>
  </tr>
</table>

*The 2D MDS axes are arbitrary up to rotation and reflection; the plots should be compared by their relative configurations, not by the absolute coordinates.*

### Generate the plots

Run from the repository root:

```bash
python examples/demo_scatter.py
julia examples/demo_scatter.jl
```

The scripts save their figures to:

- `examples/output/rcor_scatter_python.png`
- `examples/output/rcor_scatter_julia.png`

Both figures are included in `assets/`. Re-running either example regenerates its plot in `examples/output/`.

A Jupyter notebook is also available at [`examples/demo_scatter.ipynb`](examples/demo_scatter.ipynb).

## Repository structure

```text
RCoR/
├── assets/
│   ├── intense-logo.png
│   ├── scatter_example_julia.png
│   └── scatter_example_python.png
├── examples/
│   ├── data/demo_variable_spectra_n12_d8/
│   ├── demo_scatter.jl
│   ├── demo_scatter.py
│   └── demo_scatter.ipynb
├── julia/RR_nuclear_optimized_final.jl
├── src/rcor_nuclear/
│   ├── __init__.py
│   └── core.py
├── pyproject.toml
└── README.md
```

## Acknowledgment

The present work was carried out as part of the “INTENSE” project, financially supported by the Ministry of Enterprise and Made in Italy (MIMI).
