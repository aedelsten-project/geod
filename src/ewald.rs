// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Electrostatic potential and field of the infinite crystal, by Ewald summation.
//!
//! Used as the reference the finite embedding cluster should reproduce. Atomic units throughout:
//! positions in Bohr, potential in Eh/e, field in Eh/(e Bohr). The unit cell must be neutral; the
//! average potential over the cell is zero (the usual "tin-foil" convention).

use std::f64::consts::PI;

use crate::crystal::{CellParameters, WyckoffBlock, ANGSTROM_TO_BOHR};

/// Target accuracy: the real- and reciprocal-space sums are truncated where their terms drop
/// below roughly `exp(-EWALD_CUTOFF^2)`.
const EWALD_CUTOFF: f64 = 5.5;
/// Terms closer than this (Bohr) are the evaluation point's own charge.
const SELF_TOLERANCE: f64 = 1.0e-6;

/// Complementary error function, accurate to ~1e-15 for `x >= 0`.
pub fn erfc(x: f64) -> f64 {
    assert!(x >= 0.0);
    if x < 2.0 {
        // erf(x) = 2/sqrt(pi) sum_n (-1)^n x^(2n+1) / (n! (2n+1))
        let mut term = x;
        let mut sum = x;
        let mut n = 0.0;
        while term.abs() > 1.0e-17 * sum.abs() {
            n += 1.0;
            term *= -x * x / n;
            sum += term / (2.0 * n + 1.0);
        }
        1.0 - 2.0 / PI.sqrt() * sum
    } else {
        // erfc(x) = exp(-x^2)/sqrt(pi) / (x + (1/2)/(x + 1/(x + (3/2)/(x + 2/(x + ...)))))
        let mut t = x;
        for n in (1..=60).rev() {
            t = x + (n as f64 / 2.0) / t;
        }
        (-x * x).exp() / (PI.sqrt() * t)
    }
}

/// Precomputed Ewald sum for the atoms of one unit cell.
pub struct Ewald {
    /// Unit-cell charges and Cartesian positions (Bohr).
    charges: Vec<(f64, [f64; 3])>,
    alpha: f64,
    /// Lattice translations within the real-space cutoff.
    translations: Vec<[f64; 3]>,
    /// Reciprocal vectors of one half-space, with `2 * 4 pi/V exp(-k^2/4 alpha^2)/k^2` and the
    /// structure factor sums `sum_j q_j cos(k.r_j)` and `sum_j q_j sin(k.r_j)`.
    reciprocal: Vec<([f64; 3], f64, f64, f64)>,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

impl Ewald {
    pub fn new(blocks: &[WyckoffBlock], cell: &CellParameters) -> Self {
        let charges: Vec<(f64, [f64; 3])> = blocks.iter().flat_map(|block| {
            let q = block.site.oxidation.expect("Need oxidation state to run the code...");
            block.positions.iter().map(move |p| (q, cell.to_cartesian(*p).map(|x| x * ANGSTROM_TO_BOHR)))
        }).collect();
        let a = cell.vectors.map(|v| v.map(|x| x * ANGSTROM_TO_BOHR));
        let volume = dot(a[0], cross(a[1], a[2])).abs();
        // Reciprocal vectors b_i with a_i . b_j = delta_ij; |b_i| = 1 / (spacing of planes i).
        let b = [cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])].map(|v| v.map(|x| x / volume));
        let spacing = b.map(|v| 1.0 / dot(v, v).sqrt());

        // Real-space cutoff of about one cell; alpha follows from it.
        let r_cut = spacing.iter().cloned().fold(0.0, f64::max);
        let alpha = EWALD_CUTOFF / r_cut;
        let k_cut = 2.0 * alpha * EWALD_CUTOFF;

        let n = spacing.map(|d| (r_cut / d).ceil() as i32 + 1);
        let mut translations = Vec::new();
        for i in -n[0]..=n[0] {
            for j in -n[1]..=n[1] {
                for k in -n[2]..=n[2] {
                    translations.push(std::array::from_fn(|c| i as f64 * a[0][c] + j as f64 * a[1][c] + k as f64 * a[2][c]));
                }
            }
        }

        // k = 2 pi (m0 b0 + m1 b1 + m2 b2); |k . a_i| = 2 pi |m_i| <= k_cut |a_i|.
        let m = a.map(|v| (k_cut * dot(v, v).sqrt() / (2.0 * PI)).ceil() as i32);
        let mut reciprocal = Vec::new();
        for m0 in 0..=m[0] {
            for m1 in -m[1]..=m[1] {
                for m2 in -m[2]..=m[2] {
                    // One half-space: k and -k give the same contribution.
                    if m0 == 0 && (m1 < 0 || (m1 == 0 && m2 <= 0)) {
                        continue;
                    }
                    let k: [f64; 3] = std::array::from_fn(|c| 2.0 * PI * (m0 as f64 * b[0][c] + m1 as f64 * b[1][c] + m2 as f64 * b[2][c]));
                    let k2 = dot(k, k);
                    if k2 > k_cut * k_cut {
                        continue;
                    }
                    let coeff = 2.0 * 4.0 * PI / volume * (-k2 / (4.0 * alpha * alpha)).exp() / k2;
                    let (c, s) = charges.iter().fold((0.0, 0.0), |(c, s), (q, r)| {
                        let kr = dot(k, *r);
                        (c + q * kr.cos(), s + q * kr.sin())
                    });
                    reciprocal.push((k, coeff, c, s));
                }
            }
        }
        Self { charges, alpha, translations, reciprocal }
    }

    /// Potential and field at `r` (Cartesian, Bohr). If `r` is an atom site, that atom's own
    /// charge is left out, as for an on-site potential.
    pub fn potential_and_field(&self, r: [f64; 3]) -> (f64, [f64; 3]) {
        let alpha = self.alpha;
        let mut v = 0.0;
        let mut e = [0.0; 3];
        for (q, rj) in self.charges.iter() {
            for t in self.translations.iter() {
                let d: [f64; 3] = std::array::from_fn(|c| r[c] - rj[c] - t[c]);
                let dist = dot(d, d).sqrt();
                if dist < SELF_TOLERANCE {
                    // Own charge: erfc(a d)/d - 1/d -> -2a/sqrt(pi); no field by symmetry.
                    v -= q * 2.0 * alpha / PI.sqrt();
                    continue;
                }
                let erfc_term = erfc(alpha * dist) / dist;
                v += q * erfc_term;
                let radial = q * (erfc_term + 2.0 * alpha / PI.sqrt() * (-alpha * alpha * dist * dist).exp()) / (dist * dist);
                for c in 0..3 {
                    e[c] += radial * d[c];
                }
            }
        }
        for (k, coeff, c, s) in self.reciprocal.iter() {
            let kr = dot(*k, r);
            let (sin_kr, cos_kr) = kr.sin_cos();
            v += coeff * (cos_kr * c + sin_kr * s);
            let ek = coeff * (sin_kr * c - cos_kr * s);
            for i in 0..3 {
                e[i] += ek * k[i];
            }
        }
        (v, e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::WyckoffSite;

    #[test]
    fn erfc_values() {
        for (x, expected) in [
            (0.0, 1.0),
            (0.5, 0.4795001221869535),
            (1.0, 0.15729920705028513),
            (1.9999999, 0.0046777370477462145),
            (2.0, 0.004677734981047266),
            (3.0, 2.209049699858544e-05),
            (4.0, 1.541725790028002e-08),
            (5.0, 1.5374597944280347e-12),
        ] {
            assert!((erfc(x) - expected).abs() < 1.0e-14 + 1.0e-12 * expected, "erfc({x}) = {}", erfc(x));
        }
    }

    fn block(element: &str, q: f64, positions: Vec<[f64; 3]>) -> WyckoffBlock {
        WyckoffBlock {
            site: WyckoffSite {
                label: element.to_string(),
                type_symbol: element.to_string(),
                element: element.to_string(),
                oxidation: Some(q),
                occupancy: 1.0,
                multiplicity: None,
                wyckoff_letter: None,
                position: positions[0],
            },
            positions,
        }
    }

    /// Rock salt: the on-site potential is -M / r0 with the Madelung constant M = 1.747564594633.
    #[test]
    fn rock_salt_madelung_constant() {
        let fcc = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]];
        let shifted = fcc.map(|p| [p[0] + 0.5, p[1], p[2]].map(crate::crystal::wrap));
        let blocks = [block("Na", 1.0, fcc.to_vec()), block("Cl", -1.0, shifted.to_vec())];
        let a = 5.64; // Angstrom
        let cell = CellParameters::new(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let ewald = Ewald::new(&blocks, &cell);
        let r0 = a / 2.0 * ANGSTROM_TO_BOHR;
        let na = cell.to_cartesian(fcc[1]).map(|x| x * ANGSTROM_TO_BOHR);
        let (v, e) = ewald.potential_and_field(na);
        assert!((v + 1.747564594633 / r0).abs() < 1.0e-10, "v = {v}, expected {}", -1.747564594633 / r0);
        assert!(e.iter().all(|x| x.abs() < 1.0e-10));
    }

    /// The field is minus the gradient of the potential.
    #[test]
    fn field_is_gradient() {
        let blocks = [
            block("A", 2.0, vec![[0.1, 0.2, 0.3], [0.6, 0.1, 0.8]]),
            block("B", -1.0, vec![[0.4, 0.7, 0.2], [0.9, 0.5, 0.5], [0.2, 0.9, 0.6], [0.7, 0.3, 0.1]]),
        ];
        let cell = CellParameters::new(4.0, 5.0, 6.0, 80.0, 95.0, 105.0).unwrap();
        let ewald = Ewald::new(&blocks, &cell);
        let r = [1.3, 2.1, 0.7];
        let (_, e) = ewald.potential_and_field(r);
        let h = 1.0e-5;
        for c in 0..3 {
            let mut rp = r;
            let mut rm = r;
            rp[c] += h;
            rm[c] -= h;
            let grad = (ewald.potential_and_field(rp).0 - ewald.potential_and_field(rm).0) / (2.0 * h);
            assert!((e[c] + grad).abs() < 1.0e-7, "component {c}: E = {}, -dV = {}", e[c], -grad);
        }
    }
}
