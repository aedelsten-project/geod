// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Complementary charges that cancel the multipole moments of a unit cell.
//!
//! Implements the method of P. V. Sushko and I. V. Abarenkov, "General Purpose Electrostatic
//! Embedding Potential", J. Chem. Theory Comput. 6, 1323 (2010). Equation numbers in the
//! comments refer to that paper.
//!
//! All components of all multipole moments up to order `M` are labelled by the integer points
//! `n = (n1, n2, n3)` with `n1 + n2 + n3 <= M`, i.e. the points of the tetrahedron `T_M` with
//! vertices (0,0,0), (M,0,0), (0,M,0) and (0,0,M). Placing one charge `e(n)` on each of these
//! points (in units of the lattice vectors) gives exactly as many unknowns as there are moment
//! components, so the charges can be chosen to cancel all moments up to order `M` (eq 5).
//!
//! To keep the symmetry of the lattice, one such tetrahedron `T_M^α` is attached to each of the
//! eight corners of the unit cell, each cancelling 1/8 of the moments. The charges of all
//! tetrahedra are then merged into a single set of auxiliary charges.
//!
//! All positions are fractional co-ordinates.

use crate::crystal::WyckoffBlock;

/// Highest multipole order we allow. Higher orders work in principle, but have no realistic use
/// case and the integer tables grow quickly (`G_k(n)` contains `n!/(n-k)!`).
const MAX_ORDER: usize = 4;

/// Orientations (signs along a1, a2, a3) of the eight tetrahedra `T_M^α`, one per unit-cell
/// corner. Each tetrahedron points from its corner into the cell.
const ORIENTATIONS: [[i32; 3]; 8] = [
    [ 1,  1,  1],
    [-1,  1,  1],
    [ 1, -1,  1],
    [ 1,  1, -1],
    [-1, -1,  1],
    [-1,  1, -1],
    [ 1, -1, -1],
    [-1, -1, -1],
];

/// A point charge at a fractional position.
#[derive(Clone, Debug, Default)]
pub struct PointCharge {
    pub charge: f64,
    pub position: [f64; 3],
}

/// One of the tetrahedra `T_M^α` of charges.
struct Tetrahedron {
    /// Direction (±1) of the tetrahedron edges along each lattice vector.
    orientation: [i32; 3],
    /// Main vertex `r_0^α` (eq 6), in the frame where the unit cell spans [-0.5, 0.5].
    main_vertex: [f64; 3],
}

impl Tetrahedron {
    /// Places the main vertex at the unit-cell corner opposite to `orientation`, moved a further
    /// `shift` (the paper's `V`) away from the cell (eq 6).
    fn new(orientation: [i32; 3], shift: f64) -> Self {
        let main_vertex = orientation.map(|o| -(0.5 + shift) * (o as f64));
        Self { orientation, main_vertex }
    }

    /// Position of the charge `e(n)` (eq 7).
    fn charge_position(&self, n: [usize; 3]) -> [f64; 3] {
        std::array::from_fn(|a| self.main_vertex[a] + (self.orientation[a] * n[a] as i32) as f64)
    }
}

/// Everything that only depends on the multipole order `M`, shared between all tetrahedra.
///
/// Matrices are stored column-major, see [`col_idx`].
struct TetrahedronBasis {
    /// Highest multipole order `M` to cancel.
    max_order: usize,
    /// Points `n` of `T_M`, in the order given by [`for_each_point`]. The same list labels the
    /// moment components `m`, the charges `e(n)` and the equations `k` (eq 18).
    points: Vec<[usize; 3]>,
    /// Coefficients `g(k, m)` of the falling factorial `G_k(x) = x(x-1)...(x-k+1)` in powers of
    /// `x` (eqs 9, 11, 12); `(M+1) x (M+1)`, row `k`, column `m`. These are the signed Stirling
    /// numbers of the first kind.
    falling_factorial_coeffs: Vec<i64>,
    /// System matrix `G̃_k(n) = G_k1(n1) G_k2(n2) G_k3(n3)` (eq 19); `N_M x N_M`, row `k`,
    /// column `n`. Upper triangular in the point ordering, so the system is solved by back
    /// substitution.
    system_matrix: Vec<i64>,
}

impl TetrahedronBasis {
    fn new(max_order: usize) -> Self {
        let points = tetrahedron_points(max_order);
        let falling_factorial_coeffs = falling_factorial_coeffs(max_order);
        let falling_factorial_values = falling_factorial_values(&falling_factorial_coeffs, max_order);
        let system_matrix = system_matrix(&falling_factorial_values, &points, max_order);
        Self { max_order, points, falling_factorial_coeffs, system_matrix }
    }

    /// Number of points `N_M` in `T_M`, i.e. number of moment components and of charges.
    fn n_points(&self) -> usize {
        self.points.len()
    }

    /// `g(k, m)`.
    fn g(&self, k: usize, m: usize) -> i64 {
        self.falling_factorial_coeffs[col_idx(k, m, self.max_order + 1)]
    }

    /// `G̃_k(n)`, with `k` and `n` given as point indices.
    fn system_element(&self, k: usize, n: usize) -> f64 {
        self.system_matrix[col_idx(k, n, self.n_points())] as f64
    }
}

/// Calls `f(i, j, k)` for every integer point of `T_M`, i.e. all `i + j + k <= max_order`.
///
/// This fixes the ordering of points, moment components and charges throughout. The first point
/// is (0,0,0) and the last is (M,0,0), as in the paper's numbering (eq 18).
pub fn for_each_point(max_order: usize, mut f: impl FnMut(usize, usize, usize)) {
    for i in 0..=max_order {
        for j in 0..=(max_order-i) {
            for k in 0..=(max_order-i-j) {
                f(i, j, k);
            }
        }
    }
}

/// Index of element (`i`, `j`) in a column-major matrix with `rows` rows.
fn col_idx(i: usize, j: usize, rows: usize) -> usize {
    j * rows + i
}

/// Number of points `N_M = (M+1)(M+2)(M+3)/6` in `T_M`.
fn n_points(max_order: usize) -> usize {
    (max_order+1) * (max_order+2) * (max_order+3) / 6
}

/// All points of `T_M` in [`for_each_point`] order.
fn tetrahedron_points(max_order: usize) -> Vec<[usize; 3]> {
    let mut points = Vec::with_capacity(n_points(max_order));
    for_each_point(max_order, |i, j, k| points.push([i, j, k]));
    points
}

/// Coefficients `g(k, m)` from the recurrence (eq 12):
/// `g(k, k) = 1`, `g(k, 0) = 0` for `k > 0`, and `g(k, m) = g(k-1, m-1) - (k-1) g(k-1, m)`.
fn falling_factorial_coeffs(max_order: usize) -> Vec<i64> {
    let rows = max_order + 1;
    let mut g = vec![0; rows * rows];
    g[col_idx(0, 0, rows)] = 1;
    // Filled column by column; g(k, 0) = 0 for k > 0 is already set by the zero init.
    for m in 1..=max_order {
        g[col_idx(m, m, rows)] = 1;
        for k in m+1..=max_order {
            g[col_idx(k, m, rows)] = g[col_idx(k-1, m-1, rows)] - ((k as i64) - 1) * g[col_idx(k-1, m, rows)];
        }
    }
    g
}

/// `G_k(n) = sum_m g(k, m) n^m` for integer `n = 0..=M` (eqs 10, 11), which equals `n!/(n-k)!`
/// for `n >= k` and zero otherwise; `(M+1) x (M+1)`, row `k`, column `n`.
fn falling_factorial_values(g: &[i64], max_order: usize) -> Vec<i64> {
    let rows = max_order + 1;
    let mut values = vec![0; rows * rows];
    for n in 0..=max_order {
        for k in 0..=max_order {
            values[col_idx(k, n, rows)] = (0..=k)
                .map(|m| g[col_idx(k, m, rows)] * (n as i64).pow(m as u32))
                .sum();
        }
    }
    values
}

/// System matrix `G̃_k(n)` (eq 19).
fn system_matrix(values: &[i64], points: &[[usize; 3]], max_order: usize) -> Vec<i64> {
    let rows = max_order + 1;
    let n_points = points.len();
    let mut matrix = vec![0; n_points * n_points];
    for (n_idx, n) in points.iter().enumerate() {
        for (k_idx, k) in points.iter().enumerate() {
            matrix[col_idx(k_idx, n_idx, n_points)] =
                values[col_idx(k[0], n[0], rows)] *
                values[col_idx(k[1], n[1], rows)] *
                values[col_idx(k[2], n[2], rows)];
        }
    }
    matrix
}

/// All atoms of the unit cell as (charge, position) pairs.
fn cell_charges(blocks: &[WyckoffBlock]) -> impl Iterator<Item = (f64, &[f64; 3])> {
    blocks.iter().flat_map(|block| {
        let charge = block.site.oxidation.expect("Need oxidation state to run the code...");
        block.positions.iter().map(move |position| (charge, position))
    })
}

/// `x^m1 y^m2 z^m3` of `position - origin`.
fn monomial(position: &[f64; 3], origin: &[f64; 3], m: [usize; 3]) -> f64 {
    (0..3).map(|a| (position[a] - origin[a]).powi(m[a] as i32)).product()
}

/// Multipole moment components `P(m) = sum_j e_j x_j^m1 y_j^m2 z_j^m3` of the unit cell, taken
/// about `origin` (eq 3), for every `m` in `basis.points`.
fn compute_multipole_moments(moments: &mut [f64], blocks: &[WyckoffBlock], origin: &[f64; 3], basis: &TetrahedronBasis) {
    for (moment, m) in moments.iter_mut().zip(basis.points.iter()) {
        *moment = cell_charges(blocks)
            .map(|(charge, position)| charge * monomial(position, origin, *m))
            .sum();
    }
}

/// Right-hand side `f̃(k)` of the system for one tetrahedron (eq 20):
///
/// `f̃(k) = - sum_m P(m) prod_a (±1)^m_a g(k_a, m_a)`,
///
/// where the sign along axis `a` is the tetrahedron's orientation, and `P(m)` is taken about its
/// main vertex. The minus sign makes the charges cancel, rather than reproduce, the moments.
fn compute_right_hand_side(rhs: &mut [f64], moments: &[f64], tetrahedron: &Tetrahedron, basis: &TetrahedronBasis) {
    // (±1)^m g(k, m) along one axis.
    let signed_g = |axis: usize, k: usize, m: usize| tetrahedron.orientation[axis].pow(m as u32) as i64 * basis.g(k, m);
    for (f, k) in rhs.iter_mut().zip(basis.points.iter()) {
        let sum: f64 = moments.iter().zip(basis.points.iter())
            .map(|(moment, m)| {
                let coeff = signed_g(0, k[0], m[0]) * signed_g(1, k[1], m[1]) * signed_g(2, k[2], m[2]);
                (coeff as f64) * moment
            })
            .sum();
        *f = -sum;
    }
}

/// Solves `sum_n G̃_k(n) e(n) = f̃(k)` for the charges `e(n)` by back substitution, starting from
/// the last point (M,0,0) (eq 21; the 3D analogue of eqs 16-17).
fn solve_charges(charges: &mut [PointCharge], rhs: &[f64], basis: &TetrahedronBasis) {
    let n = basis.n_points();
    assert!(charges.len() == n);
    assert!(rhs.len() == n);
    for k in (0..n).rev() {
        let sum: f64 = (k+1..n).map(|j| basis.system_element(k, j) * charges[j].charge).sum();
        charges[k].charge = (rhs[k] - sum) / basis.system_element(k, k);
    }
}

/// Sets the position of every charge `e(n)` of the tetrahedron (eq 7).
fn place_charges(charges: &mut [PointCharge], tetrahedron: &Tetrahedron, basis: &TetrahedronBasis) {
    for (c, n) in charges.iter_mut().zip(basis.points.iter()) {
        c.position = tetrahedron.charge_position(*n);
    }
}

/// Merges charges of different tetrahedra that sit on the same position, divides by the number
/// of tetrahedra (each tetrahedron cancels the full moments, but we want 1/8 each), and drops
/// charges that cancel out.
fn combine_charges(charges: &[PointCharge], n_tetrahedra: usize) -> Vec<PointCharge> {
    const DIST_CUTOFF: f64 = 1.0e-6;
    const CHARGE_CUTOFF: f64 = 1.0e-15;
    let mut combined: Vec<PointCharge> = Vec::new();
    for c in charges.iter().filter(|c| c.charge != 0.0) {
        let existing = combined.iter_mut().find(|cc| {
            let d2: f64 = cc.position.iter().zip(c.position.iter()).map(|(a, b)| (a - b) * (a - b)).sum();
            d2.sqrt() < DIST_CUTOFF
        });
        match existing {
            Some(cc) => cc.charge += c.charge,
            None => combined.push(c.clone()),
        }
    }
    for c in combined.iter_mut() {
        c.charge /= n_tetrahedra as f64;
    }
    combined.retain(|c| c.charge.abs() >= CHARGE_CUTOFF);
    combined
}

/// Prints the summed components of each multipole order `1..=M+2` of a set of charges about
/// `origin`. If the set is a unit cell plus its auxiliary charges (or a cluster built from such
/// cells), orders up to `M` should vanish.
///
/// `charges` is called once per moment component and must yield the same (charge, position)
/// pairs every time; this avoids collecting them into a buffer.
pub fn print_moment_deviations<I>(charges: impl Fn() -> I, origin: &[f64; 3], max_order: usize)
where
    I: Iterator<Item = (f64, [f64; 3])>,
{
    println!("  Deviation of multipoles from zero.");
    println!("  Cancellation up to order {max_order} was requested:");
    for order in 1..=max_order+2 {
        let mut sum = 0.0;
        for m1 in 0..=order {
            for m2 in 0..=(order-m1) {
                let m = [m1, m2, order - m1 - m2];
                sum += charges()
                    .map(|(charge, position)| charge * monomial(&position, origin, m))
                    .sum::<f64>();
            }
        }
        println!("    P({order}) = {sum:15.5e}");
    }
}

/// Computes the auxiliary charges that, added to the unit cell, cancel all of its multipole
/// moments up to and including order `max_order` (the paper's `M`).
///
/// `shift` is the paper's `V`: how far the tetrahedra are moved out from the unit-cell corners.
///
/// The cell must be charge neutral. The returned charges are in the same fractional frame as the
/// positions in `blocks`.
pub fn compute_auxiliary_charges(blocks: &[WyckoffBlock], max_order: usize, shift: f64) -> Result<Vec<PointCharge>, Box<dyn std::error::Error>> {
    if max_order > MAX_ORDER {
        return Err(format!("Requested a multipole order of {max_order}, maximum allowed is {MAX_ORDER}").into());
    }

    // The paper puts the unit-cell corners at (±1/2, ±1/2, ±1/2) (eq 6), so the tetrahedron
    // geometry is set up in a frame where the cell spans [-0.5, 0.5]; `cell_centre` is the origin
    // of that frame in the frame of `blocks`.
    let cell_centre = [0.5; 3];
    let to_cell_frame = |p: [f64; 3]| std::array::from_fn(|a| p[a] + cell_centre[a]);

    let basis = TetrahedronBasis::new(max_order);
    let tetrahedra: Vec<Tetrahedron> = ORIENTATIONS.iter().map(|&o| Tetrahedron::new(o, shift)).collect();

    let n = basis.n_points();
    let mut charges = vec![PointCharge::default(); n * tetrahedra.len()];
    let mut moments = vec![0.0; n];
    let mut rhs = vec![0.0; n];

    for (tetrahedron, tetrahedron_charges) in tetrahedra.iter().zip(charges.chunks_mut(n)) {
        // The system is set up with the main vertex as origin (r_0^α = 0).
        let origin = to_cell_frame(tetrahedron.main_vertex);
        compute_multipole_moments(&mut moments, blocks, &origin, &basis);
        compute_right_hand_side(&mut rhs, &moments, tetrahedron, &basis);
        solve_charges(tetrahedron_charges, &rhs, &basis);
        place_charges(tetrahedron_charges, tetrahedron, &basis);
    }

    let mut aux_charges = combine_charges(&charges, tetrahedra.len());
    for c in aux_charges.iter_mut() {
        c.position = to_cell_frame(c.position);
    }

    let all_charges = || cell_charges(blocks)
        .map(|(charge, position)| (charge, *position))
        .chain(aux_charges.iter().map(|c| (c.charge, c.position)));
    print_moment_deviations(all_charges, &cell_centre, max_order);

    Ok(aux_charges)
}
