# spectral-graph-core

**Spectral graph theory in pure Rust. Jacobi eigenvalues, conservation ratios, Fiedler vectors, Cheeger constants, and perturbation analysis — zero dependencies.**

## The Big Idea: The Laplacian's Eigenvalues Encode EVERYTHING

Every undirected graph has a Laplacian matrix L = D − A (degree matrix minus adjacency). This is a real symmetric matrix, so it has real eigenvalues and orthogonal eigenvectors. And those eigenvalues tell you *everything* about the graph's structure:

```
λ₁ = 0           → always (constant vector is the eigenvector)
λ₂               → how connected the graph is (algebraic connectivity)
λ_max            → total spectral energy
λ₂ / λ_max = CR  → the conservation ratio: overall coherence
```

- **λ₂ = 0?** The graph is disconnected. Period.
- **λ₂ small?** There's a bottleneck. One or two edges hold the graph together.
- **λ₂ large?** The graph is highly connected. Information flows freely.
- **CR = λ₂/λ_max ≈ 1?** Maximum coherence (complete graph).
- **CR ≈ 0?** The graph is barely holding together.

The Fiedler vector (eigenvector for λ₂) even tells you *where* the bottleneck is — its sign pattern gives you a spectral partition.

## Quick Start

```rust
use spectral_graph_core::{Graph, laplacian, eigenvalues, conservation_ratio};

// A path graph: 0 — 1 — 2 — 3 — 4
let g = Graph::path(5);

// Eigenvalues of the Laplacian
let eigs = eigenvalues(&laplacian(&g));
println!("λ₁ = {:.4}  (should be ≈ 0)", eigs[0]);
println!("λ₂ = {:.4}  (algebraic connectivity)", eigs[1]);
println!("λ_max = {:.4}", eigs[eigs.len()-1]);

// Conservation ratio
let cr = conservation_ratio(&g);
println!("CR = {:.4}", cr);  // path(5): ≈ 0.0955

// Compare: complete graph
let complete = Graph::complete(5);
let cr_complete = conservation_ratio(&complete);
println!("CR(K₅) = {:.4}", cr_complete);  // 1.0
```

## Graph Construction

The library provides a `Graph` type with convenient constructors:

```rust
use spectral_graph_core::Graph;

// Named constructors
let path     = Graph::path(5);       // P₅: line graph
let cycle    = Graph::cycle(5);      // C₅: ring graph
let complete = Graph::complete(5);   // K₅: fully connected
let star     = Graph::star(5);       // S₅: one center, 4 leaves

// Build your own
let mut g = Graph::new(4);
g.add_edge(0, 1, 1.0);  // unit weight
g.add_edge(1, 2, 2.5);  // weighted edge
g.add_edge(2, 3, 1.0);

// Inspect
println!("vertices: {}", g.order());       // 4
println!("edges: {}", g.size());           // 3
println!("degree(1): {}", g.degree(1));    // 3.5
println!("volume: {}", g.volume());        // 7.0
```

## The Jacobi Eigenvalue Algorithm

This library implements eigenvalue decomposition from scratch using **Jacobi rotations** — the simplest correct algorithm for symmetric matrices.

The idea is elegant: repeatedly find the largest off-diagonal element, and apply a Givens rotation to annihilate it. Each rotation is orthogonal (so eigenvalues are preserved), and each step reduces the total off-diagonal energy. Convergence is guaranteed.

```
Before rotation:           After rotation:

    ┌ a  b ┐       →          ┌ α  0 ┐
    └ b  c ┘                  └ 0  γ ┘

where α = a − τ·b, γ = c + τ·b, τ = tan(2θ)
```

Each sweep is O(n³). For the small-to-medium graphs this library targets (up to ~100 vertices), convergence is fast (typically <20 sweeps).

```rust
use spectral_graph_core::{Graph, laplacian, eigenvalues, eigenvectors};

let g = Graph::cycle(6);
let l = laplacian(&g);

// Just eigenvalues
let vals = eigenvalues(&l);
// [0.0, 0.586, 0.586, 2.0, 3.414, 3.414]

// Eigenvalues AND eigenvectors
let (vals, vecs) = eigenvectors(&l);
// vecs[i][j] = i-th component of the j-th eigenvector
// The first eigenvector is the constant vector (1/√n, ..., 1/√n)
```

## Three Laplacians

```rust
use spectral_graph_core::{Graph, laplacian, normalized_laplacian, signless_laplacian};

let g = Graph::star(5);

// Combinatorial Laplacian: L = D − A
let l = laplacian(&g);
// Star(5): diag(4, 1, 1, 1, 1) − adjacency

// Normalized Laplacian: ℒ = D^{-1/2} L D^{-1/2}
// Eigenvalues always in [0, 2] for connected graphs
let nl = normalized_laplacian(&g).unwrap();
// Returns None if any vertex has degree 0

// Signless Laplacian: Q = D + A
// Eigenvalues relate to graph bipartiteness
let q = signless_laplacian(&g);
```

The normalized Laplacian is what we use for CR — its eigenvalues are bounded in [0, 2], making the ratio well-defined and comparable across graphs of different sizes.

## The Conservation Ratio: CR = λ₂ / λ_max

The crown jewel of the library. A single number that captures graph coherence:

```rust
use spectral_graph_core::{Graph, conservation_ratio};

// Complete graph: CR = 1.0 (perfect)
assert!((conservation_ratio(&Graph::complete(5)) - 1.0).abs() < 1e-8);

// Path graph: CR ≈ 0.095 (fragile)
let cr_path = conservation_ratio(&Graph::path(5));

// Star graph: CR ≈ 0.25 (single point of failure)
let cr_star = conservation_ratio(&Graph::star(5));

// Add edges → CR goes up
let mut g = Graph::path(5);
let cr_before = conservation_ratio(&g);
g.add_edge(0, 4, 1.0);  // make it a cycle
let cr_after = conservation_ratio(&g);
assert!(cr_after > cr_before);
```

**What CR tells you:**
- CR = 1.0 → complete graph, maximum robustness
- CR ≈ 0.5 → reasonably well-connected
- CR < 0.1 → graph has serious bottlenecks
- CR = 0 → disconnected

CR is related to the Cheeger constant h by Cheeger's inequality:
```
CR/2 ≤ h ≤ √(2·CR)
```

## The Fiedler Vector: Where the Graph Wants to Split

The eigenvector for λ₂ (the Fiedler vector) encodes the graph's "natural partition":

```rust
use spectral_graph_core::{Graph, fiedler_vector};

// Bipartite graph: {0,1} vs {2,3}
let mut g = Graph::new(4);
g.add_edge(0, 2, 1.0);
g.add_edge(0, 3, 1.0);
g.add_edge(1, 2, 1.0);
g.add_edge(1, 3, 1.0);

let fv = fiedler_vector(&g);
// fv ≈ [0.5, 0.5, -0.5, -0.5]
// Sign pattern: {0,1} positive, {2,3} negative → the natural cut!
```

The sign pattern of the Fiedler vector gives you a spectral bisection of the graph. Vertices with the same sign tend to be in the same community. This is the foundation of spectral clustering.

## Spectral Perturbation: Which Edge Matters Most?

Not all edges are created equal. Removing some edges barely changes λ₂; removing others is catastrophic.

```rust
use spectral_graph_core::{Graph, edge_sensitivity, most_influential_edge, optimal_edge_to_add};

let mut g = Graph::path(6);
// Path: 0-1-2-3-4-5

// How sensitive is λ₂ to each edge?
for (i, j, _w) in g.edges() {
    let sens = edge_sensitivity(&g, i, j);
    println!("Edge ({},{}): sensitivity = {:.4}", i, j, sens);
}
// Center edges are most sensitive — they're the bottlenecks

// Which existing edge, if removed, hurts connectivity most?
let (i, j, delta) = most_influential_edge(&g);
println!("Most influential: edge ({},{}), Δλ₂ = {:.4}", i, j, delta);

// Which non-edge, if added, helps connectivity most?
let (i, j, delta) = optimal_edge_to_add(&g);
println!("Best edge to add: ({},{}), Δλ₂ = +{:.4}", i, j, delta);
// For path(6): adding (0,5) to close the cycle
```

This is incredibly useful for:
- **Network design** — where to add a link to maximize robustness
- **Vulnerability analysis** — which link's failure is most damaging
- **Community detection** — edges with low sensitivity are inter-community bridges

## Alignment Coefficient: Comparing Graphs

```rust
use spectral_graph_core::{Graph, alignment_coefficient};

let g1 = Graph::path(5);
let g2 = Graph::path(5);
let g3 = Graph::complete(5);

let ac_identical = alignment_coefficient(&g1, &g2);  // = 1.0
let ac_different = alignment_coefficient(&g1, &g3);   // much lower
```

The alignment coefficient is the cosine similarity of the eigenvalue spectra of two graphs' normalized Laplacians. Same structure → high alignment. Different structure → low alignment.

## Full API

### Graph Operations

| Method | Description |
|--------|-------------|
| `Graph::new(n)` | Empty graph on n vertices |
| `Graph::complete(n)` | K_n |
| `Graph::path(n)` | P_n |
| `Graph::cycle(n)` | C_n |
| `Graph::star(n)` | S_n (center = vertex 0) |
| `g.add_edge(i, j, w)` | Add weighted undirected edge |
| `g.remove_edge(i, j)` | Remove edge |
| `g.order()` | Number of vertices |
| `g.size()` | Number of edges |
| `g.degree(i)` | Weighted degree |
| `g.volume()` | Total volume (sum of degrees) |
| `g.edges()` | Iterator: (i, j, weight) |

### Laplacians

| Function | Formula | Use |
|----------|---------|-----|
| `laplacian` | L = D − A | General spectral analysis |
| `normalized_laplacian` | ℒ = D^{-1/2}LD^{-1/2} | Size-invariant comparison |
| `signless_laplacian` | Q = D + A | Bipartiteness detection |

### Eigenvalues

| Function | Returns | Complexity |
|----------|---------|------------|
| `eigenvalues` | Vec<f64> (sorted ascending) | O(n³) per sweep |
| `eigenvectors` | (eigenvalues, eigenvectors) | O(n³) per sweep |

### Spectral Measures

| Function | Returns | Meaning |
|----------|---------|---------|
| `conservation_ratio` | f64 ∈ [0, 1] | λ₂/λ_max of normalized Laplacian |
| `algebraic_connectivity` | f64 | λ₂ of combinatorial Laplacian |
| `spectral_gap` | f64 | Same as algebraic connectivity (λ₁ = 0) |
| `cheeger_constant` | f64 | Lower bound: λ₂/2 |
| `fiedler_vector` | Vec<f64> | Eigenvector for λ₂ |
| `alignment_coefficient` | f64 ∈ [0, 1] | Spectral cosine similarity |

### Perturbation Analysis

| Function | Returns |
|----------|---------|
| `edge_sensitivity(i, j)` | |Δλ₂| from toggling edge |
| `most_influential_edge()` | (i, j, Δλ₂) — removal hurts most |
| `optimal_edge_to_add()` | (i, j, Δλ₂) — addition helps most |

## Honest Limitations

- **Dense representation.** The adjacency matrix is stored as Vec<Vec<f64>>. Memory usage is O(n²). Fine for n ≤ ~1000. For sparse social/network graphs with millions of vertices, use a sparse library.
- **Jacobi eigenvalue convergence.** O(n³) per sweep, typically 10-20 sweeps. For n > 500, consider Lanczos or ARPACK-based alternatives. Jacobi wins on simplicity and numerical stability, not speed.
- **No directed graphs.** The Laplacian theory here assumes undirected graphs. Directed graphs require Hermitian Laplacians or random-walk Laplacians.
- **No generalized eigenvalue problems.** Can't solve Lx = λDx (normalized Laplacian as a generalized problem). Instead, we form D^{-1/2}LD^{-1/2} explicitly.
- **No graph generators.** No Barabási–Albert, Erdős–Rényi, stochastic block model, etc. You build graphs by hand or bring your own generator.

## The Foundation

This crate is the mathematical foundation for the SuperInstance ecosystem:

- **conservation-regime** uses `conservation_ratio` and `eigenvalues` for time-series regime detection
- **sheaf-cohomology** builds on the same Laplacian theory but generalizes vertex data to vector spaces
- **symplectic-geometry** shares the dense matrix infrastructure for Hamiltonian mechanics

Everything starts here: a graph, its Laplacian, and the eigenvalues that encode its soul.

## Installation

```toml
[dependencies]
spectral-graph-core = "0.1"
```

## License

MIT

## Ecosystem Integration

- Core foundation of the SuperInstance spectral ecosystem — provides graph Laplacians, spectral decompositions, and eigenspace utilities
- Consumed by 13+ other SuperInstance crates for spectral analysis, clustering, and control
- Integrates with `conservation-regime` for regime detection via spectral signatures
- Feeds `spectral-mechanics` and `graph-thermodynamics` for physics-inspired graph analysis
- Powers the `spectral-deadband` ↔ `neyman-pearson-gap` statistical testing pipeline

