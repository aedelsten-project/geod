// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Spherical cuts of the pristine crystal, e.g. as a starting point for a QM region.
//!
//! All atoms within a radius of a centre (a Wyckoff site or an arbitrary point) are taken from a
//! supercell just large enough to contain the sphere, sorted by distance from the centre, and
//! written as .xyz in Angstrom, in the Cartesian frame of the CIF cell. That is the frame the
//! production mode expects for its `qm` file, so a cut can be used as is or edited first.

use std::error::Error;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::crystal::{CellParameters, WyckoffBlock};
use crate::input::CutCentre;

/// Atoms within this distance (Angstrom) of each other's distance from the centre count as one
/// shell; within a shell, atoms are ordered by Wyckoff site.
const SHELL_TOLERANCE: f64 = 1.0e-4;
/// Atoms this close (Angstrom) beyond the radius are still included, so that a shell exactly at
/// the radius is not split by rounding.
const RADIUS_TOLERANCE: f64 = 1.0e-6;

/// One atom of a cut.
pub struct CutAtom {
    /// Index of the Wyckoff block (site).
    pub block: usize,
    /// Fractional position.
    pub position: [f64; 3],
    /// Distance from the centre, Angstrom.
    pub distance: f64,
}

/// Resolves the centre of a cut to a fractional position.
pub fn centre_position(centre: &CutCentre, blocks: &[WyckoffBlock], cell: &CellParameters) -> Result<[f64; 3], Box<dyn Error>> {
    match centre {
        CutCentre::Site(label) => blocks.iter()
            .find(|b| b.site.label == *label)
            .map(|b| b.positions[0])
            .ok_or_else(|| {
                let labels: Vec<&str> = blocks.iter().map(|b| b.site.label.as_str()).collect();
                format!("no Wyckoff site '{label}'; the sites are {}", labels.join(", ")).into()
            }),
        CutCentre::Fractional(p) => Ok(*p),
        CutCentre::Cartesian(r) => Ok(cell.to_fractional(*r)),
    }
}

/// All atoms within `radius` (Angstrom) of `centre` (fractional), sorted by distance, then by
/// Wyckoff site within each shell, then by position.
pub fn cut_sphere(blocks: &[WyckoffBlock], cell: &CellParameters, centre: [f64; 3], radius: f64) -> Vec<CutAtom> {
    let span = cell.plane_spacings().map(|d| (radius / d).ceil() as i32 + 1);
    let base = centre.map(|x| x.floor() as i32);
    let mut atoms = Vec::new();
    for i in base[0] - span[0]..=base[0] + span[0] {
        for j in base[1] - span[1]..=base[1] + span[1] {
            for k in base[2] - span[2]..=base[2] + span[2] {
                let d = [i, j, k];
                for (b, block) in blocks.iter().enumerate() {
                    for p in block.positions.iter() {
                        let position: [f64; 3] = std::array::from_fn(|a| p[a] + d[a] as f64);
                        let r = cell.to_cartesian(std::array::from_fn(|a| position[a] - centre[a]));
                        let distance = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
                        if distance <= radius + RADIUS_TOLERANCE {
                            atoms.push(CutAtom { block: b, position, distance });
                        }
                    }
                }
            }
        }
    }
    let shell = |a: &CutAtom| (a.distance / SHELL_TOLERANCE).round() as i64;
    atoms.sort_by(|a, b| {
        shell(a).cmp(&shell(b))
            .then(a.block.cmp(&b.block))
            .then(a.position.partial_cmp(&b.position).unwrap())
    });
    atoms
}

/// Prints the cut shell by shell: distance, and how many atoms of which site.
pub fn print_shells(atoms: &[CutAtom], blocks: &[WyckoffBlock]) {
    println!("  {:>12}  Atoms", "Distance (A)");
    let mut i = 0;
    while i < atoms.len() {
        let d = atoms[i].distance;
        let mut parts: Vec<String> = Vec::new();
        while i < atoms.len() && (atoms[i].distance - d).abs() < SHELL_TOLERANCE {
            let b = atoms[i].block;
            let n = atoms[i..].iter().take_while(|a| a.block == b && (a.distance - d).abs() < SHELL_TOLERANCE).count();
            parts.push(format!("{n} x {} ({})", blocks[b].site.element, blocks[b].site.label));
            i += n;
        }
        println!("  {d:>12.4}  {}", parts.join(", "));
    }
}

/// Net formal charge of the cut, if all its sites have an oxidation state.
pub fn formal_charge(atoms: &[CutAtom], blocks: &[WyckoffBlock]) -> Option<f64> {
    atoms.iter().map(|a| blocks[a.block].site.oxidation).sum()
}

/// Writes the cut as .xyz: element and Cartesian position in Angstrom, followed by the Wyckoff
/// site and the distance from the centre as extra columns (ignored by the .xyz reader).
pub fn write_xyz(path: impl AsRef<Path>, atoms: &[CutAtom], blocks: &[WyckoffBlock], cell: &CellParameters, comment: &str) -> io::Result<()> {
    let mut file = BufWriter::new(File::create(path)?);
    writeln!(file, "{}", atoms.len())?;
    writeln!(file, "{comment}")?;
    for a in atoms.iter() {
        let r = cell.to_cartesian(a.position);
        let site = &blocks[a.block].site;
        writeln!(file, "{:<3}{:16.8}{:16.8}{:16.8}   {:<6}{:12.6}", site.element, r[0], r[1], r[2], site.label, a.distance)?;
    }
    file.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cif;
    use crate::cluster::{build_cluster, select_cells, ClusterShape};
    use crate::crystal::Crystal;
    use crate::tetrahedron::compute_auxiliary_charges;

    fn y2o3() -> (Vec<WyckoffBlock>, CellParameters) {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        (crystal.unfold(1.0e-2), crystal.cell)
    }

    #[test]
    fn cut_around_site() {
        let (blocks, cell) = y2o3();
        let centre = centre_position(&CutCentre::Site("Y0".into()), &blocks, &cell).unwrap();
        let atoms = cut_sphere(&blocks, &cell, centre, 6.0);
        assert_eq!(atoms[0].block, 0);
        assert!(atoms[0].distance < 1e-12);
        assert!(atoms.windows(2).all(|w| w[0].distance <= w[1].distance + SHELL_TOLERANCE));
        assert!(atoms.iter().all(|a| a.distance <= 6.0 + RADIUS_TOLERANCE));
        // Brute force over a generous supercell gives the same number of atoms.
        let mut n = 0;
        for d in (-3..=3).flat_map(|i| (-3..=3).flat_map(move |j| (-3..=3).map(move |k| [i, j, k]))) {
            for block in blocks.iter() {
                for p in block.positions.iter() {
                    let r = cell.to_cartesian(std::array::from_fn(|a| p[a] + d[a] as f64 - centre[a]));
                    if (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt() <= 6.0 + RADIUS_TOLERANCE {
                        n += 1;
                    }
                }
            }
        }
        assert_eq!(atoms.len(), n);
        assert!(blocks.iter().any(|b| b.site.label == "Y0"));
        assert!(centre_position(&CutCentre::Site("X9".into()), &blocks, &cell).is_err());
    }

    /// A cut written to .xyz can be read back as the QM region of an embedding.
    #[test]
    fn cut_is_a_valid_qm_region() {
        let (blocks, cell) = y2o3();
        let atoms = cut_sphere(&blocks, &cell, [0.4, 0.6, 0.3], 4.0);
        let path = std::env::temp_dir().join(format!("qm_cut_test_{}.xyz", std::process::id()));
        write_xyz(&path, &atoms, &blocks, &cell, "test").unwrap();
        let qm = crate::xyz::read_file(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(qm.len(), atoms.len());

        let aux = compute_auxiliary_charges(&blocks, 2, 0.0).unwrap();
        let cells = select_cells(&ClusterShape::Cuboid([[1, 1]; 3]), &cell, &blocks);
        let mut cluster = build_cluster(&blocks, &aux, cells, &cell);
        cluster.mark_qm_region(&qm, &blocks, &cell).unwrap();
        for (k, &i) in cluster.qm.iter().enumerate() {
            assert_eq!(cluster.atoms[i].block, atoms[k].block);
        }
    }
}
