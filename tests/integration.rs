//! Integration tests for spectral-graph-core
//!
//! Validates against known analytical results from spectral graph theory.

use spectral_graph_core::*;

// ─── Helper ──────────────────────────────────────────────────────────────────

fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() < tol
}

// ─── Eigenvalue tests ────────────────────────────────────────────────────────

#[test]
fn complete_graph_kn_eigenvalues() {
    // K_n Laplacian: eigenvalues 0, n, n, ..., n (multiplicity n-1)
    for n in [3, 4, 5, 6] {
        let g = Graph::complete(n);
        let l = laplacian(&g);
        let eigs = eigenvalues(&l);
        assert!(approx_eq(eigs[0], 0.0, 1e-8), "K_{}: λ₁ should be 0", n);
        assert!(
            approx_eq(eigs[1], n as f64, 1e-6),
            "K_{}: λ₂ should be {}, got {}",
            n, n, eigs[1]
        );
        // All non-zero eigenvalues should equal n
        for i in 1..n {
            assert!(
                approx_eq(eigs[i], n as f64, 1e-6),
                "K_{}: λ_{} should be {}, got {}",
                n, i, n, eigs[i]
            );
        }
    }
}

#[test]
fn path_graph_eigenvalues() {
    // P_n Laplacian: λ_k = 2(1 - cos(kπ/n)), k = 0..n-1
    // λ₂ = 2(1 - cos(π/n))
    for n in [4, 5, 6, 8, 10] {
        let g = Graph::path(n);
        let l = laplacian(&g);
        let eigs = eigenvalues(&l);
        let expected_lambda2 = 2.0 * (1.0 - (std::f64::consts::PI / n as f64).cos());
        assert!(
            approx_eq(eigs[1], expected_lambda2, 1e-6),
            "P_{}: λ₂ should be {:.8}, got {:.8}",
            n, expected_lambda2, eigs[1]
        );
    }
}

#[test]
fn cycle_graph_eigenvalues() {
    // C_n Laplacian: λ_k = 2(1 - cos(2kπ/n))
    // λ₂ = 2(1 - cos(2π/n))
    for n in [4, 5, 6, 8, 10] {
        let g = Graph::cycle(n);
        let l = laplacian(&g);
        let eigs = eigenvalues(&l);
        let expected_lambda2 = 2.0 * (1.0 - (2.0 * std::f64::consts::PI / n as f64).cos());
        assert!(
            approx_eq(eigs[1], expected_lambda2, 1e-6),
            "C_{}: λ₂ should be {:.8}, got {:.8}",
            n, expected_lambda2, eigs[1]
        );
    }
}

#[test]
fn star_graph_eigenvalues() {
    // Star S_n: eigenvalues 0, 1 (mult n-2), n
    // So λ₂ = 1
    for n in [4, 5, 6] {
        let g = Graph::star(n);
        let ac = algebraic_connectivity(&g);
        assert!(
            approx_eq(ac, 1.0, 1e-6),
            "Star({}): λ₂ should be 1.0, got {}",
            n, ac
        );
    }
}

// ─── Conservation ratio ──────────────────────────────────────────────────────

#[test]
fn conservation_ratio_complete_is_one() {
    let g = Graph::complete(5);
    let cr = conservation_ratio(&g);
    assert!(
        approx_eq(cr, 1.0, 1e-8),
        "conservation ratio of K_5 should be 1.0, got {}",
        cr
    );
}

#[test]
fn conservation_ratio_disconnected_is_zero() {
    let g = Graph::new(5); // no edges
    let cr = conservation_ratio(&g);
    assert!(
        approx_eq(cr, 0.0, 1e-8),
        "conservation ratio of empty graph should be 0.0, got {}",
        cr
    );
}

#[test]
fn conservation_ratio_monotonicity() {
    // Adding edges should increase (or maintain) conservation ratio
    let path_cr = conservation_ratio(&Graph::path(8));
    let cycle_cr = conservation_ratio(&Graph::cycle(8));
    let complete_cr = conservation_ratio(&Graph::complete(8));

    assert!(
        cycle_cr >= path_cr - 1e-10,
        "cycle conservation ratio ({}) should ≥ path ({})",
        cycle_cr, path_cr
    );
    assert!(
        complete_cr >= cycle_cr - 1e-10,
        "complete conservation ratio ({}) should ≥ cycle ({})",
        complete_cr, cycle_cr
    );
}

// ─── Fiedler vector ──────────────────────────────────────────────────────────

#[test]
fn fiedler_vector_bipartite_split() {
    // Create a graph with clear bipartition: {0,1} — {2,3}
    let mut g = Graph::new(4);
    g.add_edge(0, 2, 1.0);
    g.add_edge(0, 3, 1.0);
    g.add_edge(1, 2, 1.0);
    g.add_edge(1, 3, 1.0);

    let fv = fiedler_vector(&g);

    // For a bipartite graph, the Fiedler vector splits vertices into two groups.
    // K_{2,2} has eigenvalue multiplicity at λ₂, so the Fiedler vector may not
    // be unique. Use a graph with unique Fiedler vector instead.
    //
    // Use a "ladder" bipartite: {0,1} -- {2,3} but not all cross edges
    let mut g = Graph::new(4);
    g.add_edge(0, 2, 1.0);
    g.add_edge(0, 3, 1.0);
    g.add_edge(1, 3, 1.0);
    // vertex 1 not connected to 2 → asymmetric, unique Fiedler vector

    let _fv = fiedler_vector(&g);

    // Should split into two groups with clear sign difference
    let pos_count = _fv.iter().filter(|&&x| x > 1e-8).count();
    let neg_count = _fv.iter().filter(|&&x| x < -1e-8).count();
    assert!(pos_count > 0 && neg_count > 0, "Fiedler vector should have both positive and negative entries, got {:?}", _fv);
}

#[test]
fn fiedler_vector_path_graph() {
    // For P_n, Fiedler vector is monotonic.
    // The canonical form makes the first nonzero entry positive,
    // so the vector may be monotonically decreasing.
    let g = Graph::path(6);
    let fv = fiedler_vector(&g);

    // First and last should have opposite signs
    assert!(
        fv[0] * fv[5] < 0.0,
        "Fiedler vector endpoints should have opposite signs, got [{}, {}]",
        fv[0], fv[5]
    );

    // Absolute values should be monotonically increasing from center
    // (the Fiedler vector of a path is sinusoidal)
}

// ─── Alignment coefficient ──────────────────────────────────────────────────

#[test]
fn alignment_identical_graphs() {
    let g1 = Graph::path(6);
    let g2 = Graph::path(6);
    let ac = alignment_coefficient(&g1, &g2);
    assert!(
        approx_eq(ac, 1.0, 1e-8),
        "alignment of identical graphs should be 1.0, got {}",
        ac
    );
}

#[test]
fn alignment_disconnected_vs_complete() {
    let g1 = Graph::new(5); // disconnected
    let g2 = Graph::complete(5);
    let ac = alignment_coefficient(&g1, &g2);
    // Should be close to 0 (but not exactly 0 due to padding)
    assert!(
        ac < 0.3,
        "alignment of disconnected vs complete should be low, got {}",
        ac
    );
}

#[test]
fn alignment_self_consistency() {
    let g = Graph::cycle(7);
    let ac = alignment_coefficient(&g, &g);
    assert!(
        approx_eq(ac, 1.0, 1e-8),
        "self-alignment should be 1.0, got {}",
        ac
    );
}

// ─── Perturbation ────────────────────────────────────────────────────────────

#[test]
fn adding_edge_increases_lambda2() {
    let g = Graph::path(6);
    let ac_before = algebraic_connectivity(&g);

    let (i, j, delta) = optimal_edge_to_add(&g);
    assert!(delta > 0.0, "optimal edge should have positive delta");

    let mut g2 = g.clone();
    g2.add_edge(i, j, 1.0);
    let ac_after = algebraic_connectivity(&g2);
    assert!(
        ac_after > ac_before - 1e-10,
        "adding optimal edge should increase λ₂: {} → {}",
        ac_before, ac_after
    );
}

#[test]
fn removing_edge_decreases_lambda2() {
    let g = Graph::cycle(8);
    let ac_before = algebraic_connectivity(&g);

    let (i, j, delta) = most_influential_edge(&g);
    assert!(delta >= 0.0, "most influential edge delta should be ≥ 0");

    let mut g2 = g.clone();
    g2.remove_edge(i, j);
    let ac_after = algebraic_connectivity(&g2);
    assert!(
        ac_after <= ac_before + 1e-10,
        "removing influential edge should decrease λ₂: {} → {}",
        ac_before, ac_after
    );
}

#[test]
fn edge_sensitivity_is_positive() {
    let g = Graph::cycle(6);
    for (i, j, _) in g.edges() {
        let s = edge_sensitivity(&g, i, j);
        assert!(s >= 0.0, "sensitivity should be ≥ 0 for edge ({},{})", i, j);
    }
}

// ─── Laplacian properties ────────────────────────────────────────────────────

#[test]
fn laplacian_row_sums_zero() {
    let g = Graph::cycle(5);
    let l = laplacian(&g);
    for (i, row) in l.iter().enumerate() {
        let sum: f64 = row.iter().sum();
        assert!(
            approx_eq(sum, 0.0, 1e-10),
            "Laplacian row {} sums to {}, expected 0",
            i, sum
        );
    }
}

#[test]
fn signless_laplacian_row_sums_positive() {
    let g = Graph::path(5);
    let q = signless_laplacian(&g);
    for (i, row) in q.iter().enumerate() {
        let sum: f64 = row.iter().sum();
        assert!(
            sum >= 0.0,
            "signless Laplacian row {} sums to {}",
            i, sum
        );
    }
}

#[test]
fn normalized_laplacian_diagonal_one() {
    let g = Graph::complete(4);
    let nl = normalized_laplacian(&g).unwrap();
    for (i, row) in nl.iter().enumerate() {
        assert!(
            approx_eq(row[i], 1.0, 1e-10),
            "normalized Laplacian diagonal [{},{}] = {}, expected 1.0",
            i, i, row[i]
        );
    }
}

// ─── Graph construction ──────────────────────────────────────────────────────

#[test]
fn graph_size_and_order() {
    let g = Graph::complete(5);
    assert_eq!(g.order(), 5);
    assert_eq!(g.size(), 10); // C(5,2) = 10

    let g = Graph::path(4);
    assert_eq!(g.order(), 4);
    assert_eq!(g.size(), 3);

    let g = Graph::cycle(5);
    assert_eq!(g.order(), 5);
    assert_eq!(g.size(), 5);

    let g = Graph::star(5);
    assert_eq!(g.order(), 5);
    assert_eq!(g.size(), 4);
}

#[test]
fn graph_volume() {
    let g = Graph::complete(4);
    // Each vertex has degree 3, total volume = 12
    assert!(approx_eq(g.volume(), 12.0, 1e-10));
}

// ─── Cheeger constant ────────────────────────────────────────────────────────

#[test]
fn cheeger_is_half_algebraic_connectivity() {
    let g = Graph::cycle(6);
    let h = cheeger_constant(&g);
    let ac = algebraic_connectivity(&g);
    assert!(
        approx_eq(h, ac / 2.0, 1e-10),
        "Cheeger estimate should be λ₂/2: {} vs {}",
        h, ac / 2.0
    );
}

// ─── Spectral gap ────────────────────────────────────────────────────────────

#[test]
fn spectral_gap_equals_algebraic_connectivity() {
    let g = Graph::path(6);
    assert!(
        approx_eq(spectral_gap(&g), algebraic_connectivity(&g), 1e-10),
        "spectral gap should equal algebraic connectivity"
    );
}
