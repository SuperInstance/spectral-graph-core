//! # spectral-graph-core
//!
//! A mathematically elegant spectral graph theory library in pure Rust.
//!
//! Zero external dependencies. Dense adjacency representation. Jacobi eigenvalue
//! decomposition. Conservation ratios, Fiedler vectors, Cheeger constants,
//! and spectral perturbation analysis — all from first principles.

// ─────────────────────────────────────────────────────────────────────────────
// Module 1 — Graph
// ─────────────────────────────────────────────────────────────────────────────

/// A simple, dense, undirected, weighted graph indexed by `usize`.
#[derive(Clone, Debug)]
pub struct Graph {
    /// Adjacency matrix. `adj[i][j] = w` for edge weight `w` (0.0 if no edge).
    adj: Vec<Vec<f64>>,
    /// Number of vertices.
    n: usize,
}

impl Graph {
    /// Create an empty graph on `n` vertices (no self-loops).
    pub fn new(n: usize) -> Self {
        Self {
            adj: vec![vec![0.0; n]; n],
            n,
        }
    }

    /// Add an undirected edge `(i, j)` with weight `w`.
    /// Adds to existing weight if the edge already exists.
    pub fn add_edge(&mut self, i: usize, j: usize, w: f64) {
        assert!(i < self.n && j < self.n, "vertex index out of bounds");
        assert!(i != j, "self-loops are not supported");
        self.adj[i][j] += w;
        self.adj[j][i] += w;
    }

    /// Remove the undirected edge `(i, j)` by setting its weight to zero.
    pub fn remove_edge(&mut self, i: usize, j: usize) {
        assert!(i < self.n && j < self.n, "vertex index out of bounds");
        self.adj[i][j] = 0.0;
        self.adj[j][i] = 0.0;
    }

    /// Return `true` if edge `(i, j)` exists (weight > 0).
    pub fn has_edge(&self, i: usize, j: usize) -> bool {
        self.adj[i][j] > 0.0
    }

    /// Weight of edge `(i, j)`, or 0.0 if no edge.
    pub fn edge_weight(&self, i: usize, j: usize) -> f64 {
        self.adj[i][j]
    }

    /// Number of vertices.
    pub fn order(&self) -> usize {
        self.n
    }

    /// Number of edges (counting each undirected edge once).
    pub fn size(&self) -> usize {
        let mut count = 0;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                if self.adj[i][j] > 0.0 {
                    count += 1;
                }
            }
        }
        count
    }

    /// Weighted degree of vertex `i`: sum of incident edge weights.
    pub fn degree(&self, i: usize) -> f64 {
        self.adj[i].iter().sum()
    }

    /// Total volume: sum of all weighted degrees (= 2 × sum of all edge weights).
    pub fn volume(&self) -> f64 {
        (0..self.n).map(|i| self.degree(i)).sum()
    }

    /// Reference to the raw adjacency matrix.
    pub fn adjacency(&self) -> &[Vec<f64>] {
        &self.adj
    }

    /// Iterator over edges as `(i, j, weight)`, with `i < j`.
    pub fn edges(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        let n = self.n;
        (0..n).flat_map(move |i| {
            (i + 1..n).filter_map(move |j| {
                let w = self.adj[i][j];
                if w > 0.0 { Some((i, j, w)) } else { None }
            })
        })
    }

    /// Construct a complete graph K_n with uniform weight 1.0.
    pub fn complete(n: usize) -> Self {
        let mut g = Self::new(n);
        for i in 0..n {
            for j in (i + 1)..n {
                g.add_edge(i, j, 1.0);
            }
        }
        g
    }

    /// Construct a path graph P_n.
    pub fn path(n: usize) -> Self {
        let mut g = Self::new(n);
        for i in 0..n.saturating_sub(1) {
            g.add_edge(i, i + 1, 1.0);
        }
        g
    }

    /// Construct a cycle graph C_n.
    pub fn cycle(n: usize) -> Self {
        let mut g = Self::path(n);
        if n > 2 {
            g.add_edge(n - 1, 0, 1.0);
        }
        g
    }

    /// Construct a star graph on `n` vertices (center = vertex 0).
    pub fn star(n: usize) -> Self {
        let mut g = Self::new(n);
        for i in 1..n {
            g.add_edge(0, i, 1.0);
        }
        g
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 2 — Laplacian
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the **combinatorial Laplacian** `L = D − A`.
pub fn laplacian(g: &Graph) -> Vec<Vec<f64>> {
    let n = g.order();
    let adj = g.adjacency();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        let deg = g.degree(i);
        l[i][i] = deg;
        for j in 0..n {
            if i != j {
                l[i][j] = -adj[i][j];
            }
        }
    }
    l
}

/// Compute the **normalized Laplacian** `ℒ = D^{-1/2} L D^{-1/2}`.
///
/// Returns `None` if any vertex has degree zero (isolated vertex).
pub fn normalized_laplacian(g: &Graph) -> Option<Vec<Vec<f64>>> {
    let n = g.order();
    let adj = g.adjacency();
    let mut inv_sqrt_deg = vec![0.0; n];
    for i in 0..n {
        let d = g.degree(i);
        if d <= 0.0 {
            return None;
        }
        inv_sqrt_deg[i] = 1.0 / d.sqrt();
    }

    let mut nl = vec![vec![0.0; n]; n];
    for i in 0..n {
        nl[i][i] = 1.0;
        for j in 0..n {
            if i != j && adj[i][j] > 0.0 {
                nl[i][j] = -inv_sqrt_deg[i] * inv_sqrt_deg[j] * adj[i][j];
            }
        }
    }
    Some(nl)
}

/// Compute the **signless Laplacian** `Q = D + A`.
pub fn signless_laplacian(g: &Graph) -> Vec<Vec<f64>> {
    let n = g.order();
    let adj = g.adjacency();
    let mut q = vec![vec![0.0; n]; n];
    for i in 0..n {
        q[i][i] = g.degree(i);
        for j in 0..n {
            if i != j {
                q[i][j] = adj[i][j];
            }
        }
    }
    q
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 3 — Eigen (Jacobi rotation method)
// ─────────────────────────────────────────────────────────────────────────────

/// Jacobi eigenvalue algorithm for real symmetric matrices.
///
/// Finds the largest off-diagonal element, applies a Givens rotation to
/// annihilate it, and repeats until convergence. O(n³) per sweep, numerically
/// stable, and beautifully simple.
struct JacobiSolver {
    a: Vec<Vec<f64>>,
    n: usize,
}

impl JacobiSolver {
    fn new(matrix: &[Vec<f64>]) -> Self {
        let n = matrix.len();
        let mut a = matrix.to_vec();
        // Symmetrize (handle floating-point asymmetry)
        for i in 0..n {
            for j in (i + 1)..n {
                let avg = (a[i][j] + a[j][i]) / 2.0;
                a[i][j] = avg;
                a[j][i] = avg;
            }
        }
        Self { a, n }
    }

    /// Find the index `(p, q)` of the largest off-diagonal element.
    fn find_pivot(&self) -> (usize, usize, f64) {
        let mut p = 0;
        let mut q = 1;
        let mut max_val = 0.0f64;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let v = self.a[i][j].abs();
                if v > max_val {
                    max_val = v;
                    p = i;
                    q = j;
                }
            }
        }
        (p, q, max_val)
    }

    /// Apply one Jacobi rotation to zero out `a[p][q]`.
    fn rotate(&mut self, v: &mut Vec<Vec<f64>>, p: usize, q: usize) {
        let (theta, tau);
        let apq = self.a[p][q];

        if apq.abs() < 1e-30 {
            theta = 0.0;
            tau = 0.0;
        } else {
            let app = self.a[p][p];
            let aqq = self.a[q][q];
            let diff = aqq - app;

            if diff.abs() < 1e-30 {
                theta = std::f64::consts::FRAC_PI_4;
            } else {
                theta = 0.5 * (2.0 * apq / diff).atan();
            }
            tau = theta.tan();
        }

        let cos = theta.cos();
        let sin = theta.sin();
        let _one_minus_cos_t_tau = 1.0 - cos - sin * tau;

        // Update matrix elements
        // a'[p][p] = a[p][p] - tan(θ) * a[p][q]
        // a'[q][q] = a[q][q] + tan(θ) * a[p][q]
        // a'[p][q] = 0
        self.a[p][p] -= tau * apq;
        self.a[q][q] += tau * apq;
        self.a[p][q] = 0.0;
        self.a[q][p] = 0.0;

        // Update remaining elements
        for r in 0..self.n {
            if r == p || r == q {
                continue;
            }
            let arp = self.a[r][p];
            let arq = self.a[r][q];
            self.a[r][p] = cos * arp - sin * arq;
            self.a[p][r] = self.a[r][p];
            self.a[r][q] = sin * arp + cos * arq;
            self.a[q][r] = self.a[r][q];
        }

        // Update eigenvectors
        for r in 0..self.n {
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = cos * vrp - sin * vrq;
            v[r][q] = sin * vrp + cos * vrq;
        }
    }

    /// Run the Jacobi iteration to convergence.
    fn solve(mut self) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = self.n;
        let mut v = vec![vec![0.0; n]; n];
        for i in 0..n {
            v[i][i] = 1.0;
        }

        let tol = 1e-12;
        let max_iter = 100 * n * n;

        for _ in 0..max_iter {
            let (p, q, max_off) = self.find_pivot();
            if max_off < tol {
                break;
            }
            self.rotate(&mut v, p, q);
        }

        // Extract eigenvalues from diagonal
        let eigenvalues: Vec<f64> = (0..n).map(|i| self.a[i][i]).collect();

        // Sort by eigenvalue (ascending), rearranging eigenvectors accordingly
        let mut indices: Vec<usize> = (0..n).collect();
        indices.sort_by(|&a, &b| eigenvalues[a].partial_cmp(&eigenvalues[b]).unwrap());

        let sorted_vals: Vec<f64> = indices.iter().map(|&i| eigenvalues[i]).collect();
        let sorted_vecs: Vec<Vec<f64>> = (0..n)
            .map(|r| indices.iter().map(|&c| v[r][c]).collect())
            .collect();

        (sorted_vals, sorted_vecs)
    }
}

/// Compute eigenvalues of a real symmetric matrix using Jacobi rotation.
pub fn eigenvalues(matrix: &[Vec<f64>]) -> Vec<f64> {
    JacobiSolver::new(matrix).solve().0
}

/// Compute eigenvalues and eigenvectors of a real symmetric matrix.
///
/// Returns `(eigenvalues, eigenvectors)` where `eigenvectors[i][j]` is the
/// `i`-th component of the `j`-th eigenvector (corresponding to `eigenvalues[j]`).
pub fn eigenvectors(matrix: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    JacobiSolver::new(matrix).solve()
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 4 — Spectral
// ─────────────────────────────────────────────────────────────────────────────

/// **Conservation ratio** λ₂ / λ_max of the normalized Laplacian.
///
/// - Complete graph → 1.0 (maximum connectivity)
/// - Disconnected graph → 0.0
pub fn conservation_ratio(g: &Graph) -> f64 {
    match normalized_laplacian(g) {
        Some(nl) => {
            let eigs = eigenvalues(&nl);
            if eigs.len() < 2 {
                return 0.0;
            }
            let lambda2 = eigs[1];
            let lambda_max = eigs[eigs.len() - 1];
            if lambda_max.abs() < 1e-14 {
                0.0
            } else {
                lambda2 / lambda_max
            }
        }
        None => 0.0, // isolated vertices
    }
}

/// **Fiedler vector**: the eigenvector corresponding to λ₂ of the Laplacian.
///
/// The sign pattern of this vector gives a spectral partition of the graph.
pub fn fiedler_vector(g: &Graph) -> Vec<f64> {
    let l = laplacian(g);
    let (_vals, vecs) = eigenvectors(&l);
    let n = g.order();

    // Find λ₂ (second smallest eigenvalue)
    // Eigenvalues are sorted ascending, so index 1 is λ₂
    // But verify: λ₁ should be ≈ 0
    if n < 2 {
        return vec![0.0; n];
    }

    // Extract the eigenvector for λ₂
    let fiedler_col: Vec<f64> = vecs.iter().map(|row| row[1]).collect();

    // Normalize so that the first nonzero entry is positive (canonical form)
    let first_nonzero = fiedler_col.iter().find(|&&x| x.abs() > 1e-14).copied().unwrap_or(0.0);
    if first_nonzero < 0.0 {
        return fiedler_col.iter().map(|x| -x).collect();
    }
    fiedler_col
}

/// **Spectral gap**: λ₂ − λ₁ of the Laplacian.
///
/// Since λ₁ = 0 for any graph, this equals λ₂ (= algebraic connectivity).
pub fn spectral_gap(g: &Graph) -> f64 {
    algebraic_connectivity(g)
}

/// **Algebraic connectivity**: λ₂ of the Laplacian.
///
/// Positive iff the graph is connected. Larger values mean better connectivity.
pub fn algebraic_connectivity(g: &Graph) -> f64 {
    let l = laplacian(g);
    let vals = eigenvalues(&l);
    if vals.len() < 2 {
        0.0
    } else {
        // λ₁ ≈ 0, so λ₂ = vals[1]
        vals[1].max(0.0) // clamp numerical noise
    }
}

/// **Cheeger constant** estimate: h ≥ λ₂ / 2.
///
/// This is the lower bound from Cheeger's inequality. The true Cheeger constant
/// satisfies λ₂/2 ≤ h ≤ √(2λ₂).
pub fn cheeger_constant(g: &Graph) -> f64 {
    let ac = algebraic_connectivity(g);
    ac / 2.0
}

/// **Alignment coefficient**: cosine similarity of the eigenvalue spectra
/// of two graphs' normalized Laplacians.
///
/// - Identical graphs → 1.0
/// - Completely different structure → ≈ 0.0
pub fn alignment_coefficient(g1: &Graph, g2: &Graph) -> f64 {
    let spec1 = match normalized_laplacian(g1) {
        Some(nl) => eigenvalues(&nl),
        None => vec![0.0; g1.order()],
    };
    let spec2 = match normalized_laplacian(g2) {
        Some(nl) => eigenvalues(&nl),
        None => vec![0.0; g2.order()],
    };

    // Pad shorter spectrum with zeros
    let max_len = spec1.len().max(spec2.len());
    let mut s1 = spec1;
    let mut s2 = spec2;
    s1.resize(max_len, 0.0);
    s2.resize(max_len, 0.0);

    // Cosine similarity
    let dot: f64 = s1.iter().zip(&s2).map(|(a, b)| a * b).sum();
    let norm1: f64 = s1.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm2: f64 = s2.iter().map(|x| x * x).sum::<f64>().sqrt();

    if norm1 < 1e-14 || norm2 < 1e-14 {
        0.0
    } else {
        (dot / (norm1 * norm2)).clamp(0.0, 1.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module 5 — Perturbation (spectral sensitivity)
// ─────────────────────────────────────────────────────────────────────────────

/// **Edge sensitivity**: change in λ₂ from toggling edge `(i, j)`.
///
/// Returns the magnitude of Δλ₂ from adding/removing the edge.
pub fn edge_sensitivity(g: &Graph, i: usize, j: usize) -> f64 {
    let lambda2_before = algebraic_connectivity(g);

    let mut g_mod = g.clone();
    if g.has_edge(i, j) {
        g_mod.remove_edge(i, j);
    } else {
        g_mod.add_edge(i, j, 1.0);
    }

    let lambda2_after = algebraic_connectivity(&g_mod);
    (lambda2_after - lambda2_before).abs()
}

/// **Most influential edge**: the existing edge whose removal most decreases λ₂.
///
/// Returns `(i, j, delta)` where `delta` is the decrease in algebraic connectivity.
pub fn most_influential_edge(g: &Graph) -> (usize, usize, f64) {
    let lambda2_orig = algebraic_connectivity(g);
    let mut best = (0, 1, 0.0f64);

    for (i, j, _w) in g.edges() {
        let mut g_mod = g.clone();
        g_mod.remove_edge(i, j);
        let lambda2_new = algebraic_connectivity(&g_mod);
        let delta = lambda2_orig - lambda2_new;
        if delta > best.2 {
            best = (i, j, delta);
        }
    }

    best
}

/// **Optimal edge to add**: the non-edge whose addition most increases λ₂.
///
/// Returns `(i, j, delta)` where `delta` is the increase in algebraic connectivity.
pub fn optimal_edge_to_add(g: &Graph) -> (usize, usize, f64) {
    let lambda2_orig = algebraic_connectivity(g);
    let n = g.order();
    let mut best = (0, 0, 0.0f64);
    let mut found = false;

    for i in 0..n {
        for j in (i + 1)..n {
            if !g.has_edge(i, j) {
                let mut g_mod = g.clone();
                g_mod.add_edge(i, j, 1.0);
                let lambda2_new = algebraic_connectivity(&g_mod);
                let delta = lambda2_new - lambda2_orig;
                if !found || delta > best.2 {
                    best = (i, j, delta);
                    found = true;
                }
            }
        }
    }

    best
}

// ─────────────────────────────────────────────────────────────────────────────
// Re-exports: a flat namespace for ergonomic use
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complete_graph_eigenvalues() {
        let g = Graph::complete(4);
        let l = laplacian(&g);
        let eigs = eigenvalues(&l);
        // K_4: eigenvalues are 0, 4, 4, 4
        assert!(eigs[0].abs() < 1e-8, "λ₁ should be 0, got {}", eigs[0]);
        assert!((eigs[1] - 4.0).abs() < 1e-8, "λ₂ should be 4, got {}", eigs[1]);
    }

    #[test]
    fn test_path_graph_eigenvalues() {
        let g = Graph::path(4);
        let l = laplacian(&g);
        let eigs = eigenvalues(&l);
        let expected_lambda2 = 2.0 * (1.0 - (std::f64::consts::PI / 4.0).cos());
        assert!(
            (eigs[1] - expected_lambda2).abs() < 1e-6,
            "λ₂ for P_4 should be {}, got {}",
            expected_lambda2,
            eigs[1]
        );
    }

    #[test]
    fn test_star_graph_lambda2() {
        let g = Graph::star(5);
        let ac = algebraic_connectivity(&g);
        // Star S_5: λ₂ = 1
        assert!(
            (ac - 1.0).abs() < 1e-6,
            "λ₂ for star(5) should be 1, got {}",
            ac
        );
    }

    #[test]
    fn test_conservation_ratio_complete() {
        let g = Graph::complete(5);
        let cr = conservation_ratio(&g);
        assert!(
            (cr - 1.0).abs() < 1e-8,
            "conservation ratio of K_5 should be 1.0, got {}",
            cr
        );
    }

    #[test]
    fn test_fiedler_bipartite() {
        // Create a bipartite graph: {0,1} vs {2,3}
        let mut g = Graph::new(4);
        g.add_edge(0, 2, 1.0);
        g.add_edge(0, 3, 1.0);
        g.add_edge(1, 2, 1.0);
        g.add_edge(1, 3, 1.0);

        let fv = fiedler_vector(&g);
        // Fiedler vector should split into positive and negative groups
        let pos: Vec<_> = fv.iter().filter(|&&x| x > 1e-8).collect();
        let neg: Vec<_> = fv.iter().filter(|&&x| x < -1e-8).collect();
        assert!(pos.len() > 0 && neg.len() > 0, "Fiedler vector should have both signs");
    }

    #[test]
    fn test_alignment_identical() {
        let g1 = Graph::path(5);
        let g2 = Graph::path(5);
        let ac = alignment_coefficient(&g1, &g2);
        assert!(
            (ac - 1.0).abs() < 1e-8,
            "alignment of identical graphs should be 1.0, got {}",
            ac
        );
    }

    #[test]
    fn test_perturbation_adding_edge_increases_connectivity() {
        let mut g = Graph::path(5);
        let ac_before = algebraic_connectivity(&g);
        // Add edge to make it more connected
        g.add_edge(0, 4, 1.0); // create a cycle
        let ac_after = algebraic_connectivity(&g);
        assert!(
            ac_after > ac_before - 1e-10,
            "adding an edge should not decrease algebraic connectivity"
        );
    }
}
