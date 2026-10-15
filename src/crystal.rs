// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Crystallography on top of the raw CIF data: cell, symmetry operations,
//! asymmetric-unit sites and their unfolded "Wyckoff blocks".

use std::collections::HashMap;
use std::path::Path;

use crate::cif::{self, err, parse_number, CifError, DataBlock, Result};

/// Bohr per Angstrom.
pub const ANGSTROM_TO_BOHR: f64 = 1.0 / 0.529177210903;

#[derive(Debug, Clone)]
pub struct CellParameters {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    /// Angles in degrees, as in the CIF.
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    /// Lattice vectors as rows: vectors[0] = a, vectors[1] = b, vectors[2] = c.
    /// Standard orientation: a along x, b in the xy-plane.
    pub vectors: [[f64; 3]; 3],
}

impl CellParameters {
    pub fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Result<Self> {
        // Purge round-off like cos(90 deg) = 6e-17.
        let clean = |x: f64| if x.abs() < 1.0e-8 { 0.0 } else { x };
        let (cos_a, cos_b) = (alpha.to_radians().cos(), beta.to_radians().cos());
        let (cos_g, sin_g) = (gamma.to_radians().cos(), gamma.to_radians().sin());
        // cy / c and (cz / c)^2
        let cy_frac = (cos_a - cos_b * cos_g) / sin_g;
        let cz_frac2 = 1.0 - cos_b * cos_b - cy_frac * cy_frac;
        if !(a > 0.0 && b > 0.0 && c > 0.0) || sin_g.abs() < 1.0e-8 || !(cz_frac2 > 0.0) {
            return err(format!(
                "impossible cell: a={a} b={b} c={c} alpha={alpha} beta={beta} gamma={gamma}"
            ));
        }
        let vectors = [
            [a, 0.0, 0.0],
            [clean(b * cos_g), clean(b * sin_g), 0.0],
            [clean(c * cos_b), clean(c * cy_frac), clean(c * cz_frac2.sqrt())],
        ];
        Ok(CellParameters { a, b, c, alpha, beta, gamma, vectors })
    }

    pub fn to_cartesian(&self, frac: [f64; 3]) -> [f64; 3] {
        let v = &self.vectors;
        [0, 1, 2].map(|k| frac[0] * v[0][k] + frac[1] * v[1][k] + frac[2] * v[2][k])
    }

    /// Inverse of [`CellParameters::to_cartesian`]. Relies on the standard orientation
    /// (a along x, b in the xy-plane), which makes the vector matrix triangular.
    pub fn to_fractional(&self, cart: [f64; 3]) -> [f64; 3] {
        let v = &self.vectors;
        let f2 = cart[2] / v[2][2];
        let f1 = (cart[1] - f2 * v[2][1]) / v[1][1];
        let f0 = (cart[0] - f1 * v[1][0] - f2 * v[2][0]) / v[0][0];
        [f0, f1, f2]
    }

    /// Spacing (Angstrom) between neighbouring lattice planes spanned by (b, c), (c, a) and (a, b).
    /// A sphere of radius r spans at most r / spacing[i] cells along lattice vector i.
    pub fn plane_spacings(&self) -> [f64; 3] {
        let v = &self.vectors;
        let cross = |p: [f64; 3], q: [f64; 3]| [p[1]*q[2] - p[2]*q[1], p[2]*q[0] - p[0]*q[2], p[0]*q[1] - p[1]*q[0]];
        let norm = |p: [f64; 3]| (p[0]*p[0] + p[1]*p[1] + p[2]*p[2]).sqrt();
        let bc = cross(v[1], v[2]);
        let volume = (v[0][0]*bc[0] + v[0][1]*bc[1] + v[0][2]*bc[2]).abs();
        [volume / norm(bc), volume / norm(cross(v[2], v[0])), volume / norm(cross(v[0], v[1]))]
    }

    /// Distance in Angstrom between two fractional positions, using the
    /// nearest periodic image.
    pub fn periodic_distance(&self, p: [f64; 3], q: [f64; 3]) -> f64 {
        let d = [0, 1, 2].map(|k| {
            let x = p[k] - q[k];
            x - x.round()
        });
        let r = self.to_cartesian(d);
        (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt()
    }
}

/// x' = rotation * x + translation, in fractional coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct SymOp {
    pub rotation: [[f64; 3]; 3],
    pub translation: [f64; 3],
}

impl SymOp {
    /// Parses e.g. "-x+1/2, y, z", "x-y,x,z+1/3", "1/2+x, -y, 0.25-z".
    pub fn parse(s: &str) -> Result<Self> {
        let parts: Vec<&str> = s.split(',').collect();
        if parts.len() != 3 {
            return err(format!("symmetry operation '{s}' does not have 3 components"));
        }
        let mut rotation = [[0.0; 3]; 3];
        let mut translation = [0.0; 3];
        for (k, part) in parts.iter().enumerate() {
            (rotation[k], translation[k]) = parse_component(part)
                .map_err(|e| CifError(format!("symmetry operation '{s}': {e}")))?;
        }
        Ok(SymOp { rotation, translation })
    }

    pub fn identity() -> Self {
        SymOp {
            rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            translation: [0.0; 3],
        }
    }

    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        [0, 1, 2].map(|k| {
            let r = &self.rotation[k];
            r[0] * p[0] + r[1] * p[1] + r[2] * p[2] + self.translation[k]
        })
    }
}

/// One component like "-x+y+1/3" -> ([-1, 1, 0], 1/3).
fn parse_component(expr: &str) -> Result<([f64; 3], f64)> {
    let chars: Vec<char> = expr
        .chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| c.to_ascii_lowercase())
        .collect();
    if chars.is_empty() {
        return err("empty component");
    }
    let mut row = [0.0; 3];
    let mut shift = 0.0;
    let mut i = 0;
    while i < chars.len() {
        let mut sign = 1.0;
        match chars[i] {
            '+' => i += 1,
            '-' => {
                sign = -1.0;
                i += 1;
            }
            _ => {}
        }
        // Optional number: "1/2", "0.25", "2" (as in "2x").
        let start = i;
        while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.' || chars[i] == '/') {
            i += 1;
        }
        let number: Option<f64> = if i > start {
            let text: String = chars[start..i].iter().collect();
            Some(parse_fraction(&text).ok_or_else(|| CifError(format!("bad number '{text}'")))?)
        } else {
            None
        };
        if i < chars.len() && chars[i] == '*' {
            i += 1;
        }
        match chars.get(i) {
            Some('x') | Some('y') | Some('z') => {
                let axis = (chars[i] as u8 - b'x') as usize;
                row[axis] += sign * number.unwrap_or(1.0);
                i += 1;
            }
            _ => match number {
                Some(v) => shift += sign * v,
                None => return err(format!("cannot read '{expr}'")),
            },
        }
    }
    Ok((row, shift))
}

fn parse_fraction(s: &str) -> Option<f64> {
    match s.split_once('/') {
        Some((num, den)) => {
            let (num, den): (f64, f64) = (num.parse().ok()?, den.parse().ok()?);
            if den == 0.0 { None } else { Some(num / den) }
        }
        None => s.parse().ok(),
    }
}

/// Wraps into [0, 1), making sure round-off never leaves exactly 1.0.
pub fn wrap(x: f64) -> f64 {
    let r = x.rem_euclid(1.0);
    if r > 1.0 - 1.0e-10 { 0.0 } else { r }
}

/// One row of the `_atom_site_` loop (one site of the asymmetric unit).
#[derive(Debug, Clone)]
pub struct WyckoffSite {
    /// `_atom_site_label`, e.g. "Y0". Falls back to the type symbol.
    pub label: String,
    /// `_atom_site_type_symbol`, e.g. "Y3+". Falls back to the label.
    pub type_symbol: String,
    /// Chemical element guessed from the type symbol/label, e.g. "Y", "Pr".
    pub element: String,
    pub oxidation: Option<f64>,
    /// 1.0 if the CIF does not give it.
    pub occupancy: f64,
    /// Multiplicity stated in the CIF, useful to check the unfolding.
    pub multiplicity: Option<usize>,
    pub wyckoff_letter: Option<String>,
    pub position: [f64; 3],
}

/// A site together with all its symmetry-equivalent positions in the cell.
#[derive(Debug, Clone)]
pub struct WyckoffBlock {
    pub site: WyckoffSite,
    pub positions: Vec<[f64; 3]>,
}

#[derive(Debug, Clone)]
pub struct Crystal {
    /// Name of the CIF data block, without the `data_` prefix.
    pub name: String,
    pub cell: CellParameters,
    pub symmetry_operations: Vec<SymOp>,
    pub sites: Vec<WyckoffSite>,
}

impl Crystal {
    /// Reads the first data block of a CIF file.
    pub fn read(path: impl AsRef<Path>) -> Result<Self> {
        let blocks = cif::read_file(path)?;
        match blocks.first() {
            Some(block) => Crystal::from_block(block),
            None => err("no data_ block in file"),
        }
    }

    pub fn from_block(block: &DataBlock) -> Result<Self> {
        Ok(Crystal {
            name: block.name.clone(),
            cell: read_cell(block)?,
            symmetry_operations: read_symmetry_operations(block)?,
            sites: read_sites(block)?,
        })
    }

    /// Applies all symmetry operations to every site, wraps into the cell and
    /// drops duplicates closer than `tolerance` Angstrom.
    pub fn unfold(&self, tolerance: f64) -> Vec<WyckoffBlock> {
        self.sites
            .iter()
            .map(|site| {
                let mut positions: Vec<[f64; 3]> = Vec::new();
                for op in &self.symmetry_operations {
                    let p = op.apply(site.position).map(wrap);
                    if !positions.iter().any(|q| self.cell.periodic_distance(p, *q) < tolerance) {
                        positions.push(p);
                    }
                }
                WyckoffBlock { site: site.clone(), positions }
            })
            .collect()
    }
}

fn required_number(block: &DataBlock, tag: &str) -> Result<f64> {
    match block.get(tag) {
        Some(v) => parse_number(v).ok_or_else(|| CifError(format!("{tag}: '{v}' is not a number"))),
        None => err(format!("{tag} is missing")),
    }
}

fn read_cell(block: &DataBlock) -> Result<CellParameters> {
    CellParameters::new(
        required_number(block, "_cell_length_a")?,
        required_number(block, "_cell_length_b")?,
        required_number(block, "_cell_length_c")?,
        required_number(block, "_cell_angle_alpha")?,
        required_number(block, "_cell_angle_beta")?,
        required_number(block, "_cell_angle_gamma")?,
    )
}

fn read_symmetry_operations(block: &DataBlock) -> Result<Vec<SymOp>> {
    // Modern name first, then the older one pymatgen writes.
    for tag in ["_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"] {
        if let Some(lp) = block.loop_with(tag) {
            let col = lp.column(tag).unwrap();
            return lp.rows.iter().map(|row| SymOp::parse(&row[col])).collect();
        }
        // A single operation can be written without loop_.
        if let Some(op) = block.get(tag) {
            return Ok(vec![SymOp::parse(op)?]);
        }
    }
    // No operations listed: only safe if the group is P1.
    let group = ["_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m"]
        .iter()
        .find_map(|t| block.get(t))
        .map(|g| g.replace(' ', ""));
    match group.as_deref() {
        Some("P1") => Ok(vec![SymOp::identity()]),
        Some(g) => err(format!(
            "no symmetry operations listed for space group {g}; generating them from the name is not supported"
        )),
        None => err("no symmetry operations and no space group in file"),
    }
}

fn read_sites(block: &DataBlock) -> Result<Vec<WyckoffSite>> {
    let lp = match block.loop_with("_atom_site_fract_x") {
        Some(lp) => lp,
        None => return err("no _atom_site_fract_x loop (Cartesian-only CIFs are not supported)"),
    };
    let required = |tag: &str| lp.column(tag).ok_or_else(|| CifError(format!("{tag} is missing")));
    let (ix, iy, iz) = (
        required("_atom_site_fract_x")?,
        required("_atom_site_fract_y")?,
        required("_atom_site_fract_z")?,
    );
    let i_label = lp.column("_atom_site_label");
    let i_type = lp.column("_atom_site_type_symbol");
    let i_occ = lp.column("_atom_site_occupancy");
    let i_mult = lp
        .column("_atom_site_symmetry_multiplicity")
        .or(lp.column("_atom_site_site_symmetry_multiplicity"));
    let i_wyck = lp.column("_atom_site_wyckoff_symbol");
    if i_label.is_none() && i_type.is_none() {
        return err("atom sites have neither _atom_site_label nor _atom_site_type_symbol");
    }
    let oxidation_table = read_oxidation_table(block);

    let mut sites = Vec::new();
    for row in &lp.rows {
        // "?" and "." mean unknown / not applicable.
        let value = |i: Option<usize>| i.map(|i| row[i].as_str()).filter(|v| *v != "?" && *v != ".");
        let coordinate = |i: usize| {
            parse_number(&row[i]).ok_or_else(|| CifError(format!("bad coordinate '{}' in row {row:?}", row[i])))
        };
        let label = value(i_label).or(value(i_type)).unwrap_or("").to_string();
        let type_symbol = value(i_type).unwrap_or(&label).to_string();
        let element = element_from(&type_symbol)
            .ok_or_else(|| CifError(format!("cannot tell the element of site '{label}'")))?;
        let oxidation = oxidation_table
            .get(&type_symbol)
            .copied()
            .or_else(|| oxidation_from_symbol(&type_symbol));
        let occupancy = match value(i_occ) {
            Some(v) => parse_number(v).ok_or_else(|| CifError(format!("bad occupancy '{v}'")))?,
            None => 1.0,
        };
        sites.push(WyckoffSite {
            label,
            type_symbol,
            element,
            oxidation,
            occupancy,
            multiplicity: value(i_mult).and_then(|v| v.parse().ok()),
            wyckoff_letter: value(i_wyck).map(str::to_string),
            position: [coordinate(ix)?, coordinate(iy)?, coordinate(iz)?],
        });
    }
    Ok(sites)
}

fn read_oxidation_table(block: &DataBlock) -> HashMap<String, f64> {
    let mut table = HashMap::new();
    if let Some(lp) = block.loop_with("_atom_type_oxidation_number") {
        if let (Some(is), Some(io)) = (lp.column("_atom_type_symbol"), lp.column("_atom_type_oxidation_number")) {
            for row in &lp.rows {
                if let Some(ox) = parse_number(&row[io]) {
                    table.insert(row[is].clone(), ox);
                }
            }
        }
    }
    table
}

/// "Y3+" -> "Y", "Pr1" -> "Pr", "O2-" -> "O". Not checked against the periodic table.
fn element_from(s: &str) -> Option<String> {
    let mut chars = s.chars();
    let first = chars.next().filter(|c| c.is_ascii_alphabetic())?;
    let mut element = first.to_ascii_uppercase().to_string();
    element.extend(chars.take_while(|c| c.is_ascii_lowercase()));
    Some(element)
}

/// "Y3+" -> 3, "O2-" -> -2, "Pr" -> None.
fn oxidation_from_symbol(s: &str) -> Option<f64> {
    let tail = s.trim_start_matches(|c: char| c.is_ascii_alphabetic());
    let sign = match tail.chars().last()? {
        '+' => 1.0,
        '-' => -1.0,
        _ => return None,
    };
    let digits = &tail[..tail.len() - 1];
    let n: f64 = if digits.is_empty() { 1.0 } else { digits.parse().ok()? };
    Some(sign * n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symop_parsing() {
        let op = SymOp::parse("-x+1/2, y, z+1/2").unwrap();
        assert_eq!(op.rotation, [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        assert_eq!(op.translation, [0.5, 0.0, 0.5]);

        let op = SymOp::parse("x-y,x,z+1/3").unwrap();
        assert_eq!(op.rotation[0], [1.0, -1.0, 0.0]);
        assert!((op.translation[2] - 1.0 / 3.0).abs() < 1e-15);

        // The Julia version read "x-1/2" as -x.
        let op = SymOp::parse("x-1/2, 1/2+y, 0.25-z").unwrap();
        assert_eq!(op.rotation, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]]);
        assert_eq!(op.translation, [-0.5, 0.5, 0.25]);

        assert!(SymOp::parse("x, y").is_err());
        assert!(SymOp::parse("x, y, w").is_err());
    }

    #[test]
    fn wrap_never_gives_one() {
        assert_eq!(wrap(1.0), 0.0);
        assert_eq!(wrap(-1.0e-17), 0.0);
        assert_eq!(wrap(-0.25), 0.75);
    }

    #[test]
    fn symbols() {
        assert_eq!(element_from("Pr3+").as_deref(), Some("Pr"));
        assert_eq!(element_from("Y0").as_deref(), Some("Y"));
        assert_eq!(oxidation_from_symbol("O2-"), Some(-2.0));
        assert_eq!(oxidation_from_symbol("Y"), None);
    }

    #[test]
    fn hexagonal_cell() {
        let cell = CellParameters::new(3.0, 3.0, 5.0, 90.0, 90.0, 120.0).unwrap();
        assert!((cell.vectors[1][0] + 1.5).abs() < 1e-12);
        assert!((cell.vectors[2][2] - 5.0).abs() < 1e-12);
        assert!(CellParameters::new(1.0, 1.0, 1.0, 0.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn fractional_round_trip() {
        let cell = CellParameters::new(3.0, 4.0, 5.0, 80.0, 95.0, 110.0).unwrap();
        let frac = [0.3, -1.2, 2.7];
        let back = cell.to_fractional(cell.to_cartesian(frac));
        for k in 0..3 {
            assert!((back[k] - frac[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn y2o3() {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        assert_eq!(crystal.symmetry_operations.len(), 48);
        let unfolded = crystal.unfold(1.0e-2);
        let counts: Vec<usize> = unfolded.iter().map(|b| b.positions.len()).collect();
        assert_eq!(counts, vec![24, 8, 48]);
        for b in &unfolded {
            assert_eq!(Some(b.positions.len()), b.site.multiplicity);
        }
        assert_eq!(unfolded[0].site.oxidation, Some(3.0));
        assert_eq!(unfolded[2].site.element, "O");
    }
}
