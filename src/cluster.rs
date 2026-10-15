// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Builds the finite embedding cluster from the unit cell and its auxiliary charges.
//!
//! The cluster is described by the list of integer cell displacements that make it up; the
//! shape (cuboid, sphere, ...) only decides which displacements are in the list. From that list:
//!
//! - The *main field* is every unit-cell atom copied into every cell. Its size is known up
//!   front, and each entry remembers its Wyckoff block, so regions (e.g. AIMPs) can later be
//!   assigned per Wyckoff site.
//! - The *auxiliary field* is the sum of the auxiliary charges of all cells. Inside the cluster
//!   these cancel exactly (Sushko & Abarenkov, end of section 2), so only a shell at the surface
//!   remains.
//!
//! All auxiliary charges lie on a few fixed sublattices: along each axis at `Z - V` or `Z + V`
//! (`V` being the tetrahedron shift; both coincide for `V = 0` and `V = 1/2`). The auxiliary
//! field is therefore accumulated on a dense integer grid over the cluster, with one layer per
//! sublattice. No position matching is needed: overlapping charges simply land on the same grid
//! point. Atoms that coincide with an auxiliary sublattice are found once for the unit cell and
//! take over the auxiliary charge on their grid point.
//!
//! All positions are fractional co-ordinates.

use crate::crystal::{wrap, CellParameters, WyckoffBlock};
use crate::tetrahedron::PointCharge;
use crate::xyz::XyzAtom;

/// Two fractional offsets closer than this along an axis belong to the same sublattice.
const SUBLATTICE_TOLERANCE: f64 = 1.0e-6;
/// An atom closer than this (Angstrom) to an auxiliary sublattice point sits on it, and takes over
/// its charge.
const COINCIDENCE_TOLERANCE: f64 = 1.0e-3;
/// An atom closer than this (Angstrom) to an auxiliary sublattice point, without coinciding with
/// it, would produce a near-singular potential. We warn about these.
const NEAR_MISS_TOLERANCE: f64 = 0.5;
/// Auxiliary charges smaller than this, relative to the largest auxiliary charge of the unit
/// cell, are considered cancelled. The interior points are sums of many O(1) charges that cancel
/// analytically, so the rounding noise is well above machine epsilon.
const PRUNE_RELATIVE_TOLERANCE: f64 = 1.0e-10;
/// A QM atom from the input matches a cluster atom if they are closer than this (Angstrom).
const QM_MATCH_TOLERANCE: f64 = 1.0e-2;

/// Which cells make up the cluster.
pub enum ClusterShape {
    /// Cells from `-n[a][0]` to `+n[a][1]` along each lattice vector `a`.
    Cuboid([[u32; 2]; 3]),
    /// All cells whose atom centroid is within `radius` (Angstrom) of `reference` (fractional;
    /// defaults to the atom centroid of the unit cell at the origin).
    Sphere { radius: f64, reference: Option<[f64; 3]> },
}

/// Which atoms around the QM region are represented by AIMPs.
pub enum AimpRegion {
    /// All atoms in the cells from `-n[a][0]` to `+n[a][1]` along each lattice vector `a`.
    Cells([[u32; 2]; 3]),
    /// All atoms within `radius` (Angstrom) of `reference` (fractional; defaults to the centroid
    /// of the QM region, or the centre of the cell at the origin if there is none).
    Radius { radius: f64, reference: Option<[f64; 3]> },
}

/// Which part of the embedding an atom belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    /// Treated quantum mechanically.
    Qm,
    /// Ab initio model potential.
    Aimp,
    /// Plain point charge.
    PointCharge,
}

/// One atom of the main field.
#[derive(Clone, Debug)]
pub struct ClusterAtom {
    /// Index of the Wyckoff block (site) in the unit cell.
    pub block: usize,
    /// Index of the position within the Wyckoff block.
    pub atom: usize,
    /// Displacement of the cell the atom sits in.
    pub cell: [i32; 3],
    pub position: [f64; 3],
    /// Formal charge, plus any auxiliary charge sitting on the atom.
    pub charge: f64,
    pub region: Region,
}

/// The finite cluster: main field and auxiliary field.
pub struct EmbeddingCluster {
    /// Displacements of the cells in the cluster.
    pub cells: Vec<[i32; 3]>,
    /// Main field, cell by cell, with the atoms of each cell in unit-cell order.
    pub atoms: Vec<ClusterAtom>,
    /// Auxiliary charges left after cancellation and after merging into atoms.
    pub aux: Vec<PointCharge>,
    /// Indices into `atoms` of the QM region, in the order of the QM input.
    pub qm: Vec<usize>,
    /// Indices into `atoms` of the AIMP region, grouped by Wyckoff block (in unit-cell order) and
    /// by distance from the AIMP reference point within each block.
    pub aimp: Vec<usize>,
}

impl EmbeddingCluster {
    /// All charges of the cluster as (charge, position) pairs.
    pub fn charges(&self) -> impl Iterator<Item = (f64, [f64; 3])> + '_ {
        self.atoms.iter().map(|a| (a.charge, a.position))
            .chain(self.aux.iter().map(|c| (c.charge, c.position)))
    }
}

/// An atom of the unit cell, flattened out of the Wyckoff blocks.
struct UnitCellAtom {
    block: usize,
    atom: usize,
    position: [f64; 3],
    charge: f64,
}

fn unit_cell_atoms(blocks: &[WyckoffBlock]) -> Vec<UnitCellAtom> {
    blocks.iter().enumerate().flat_map(|(b, block)| {
        let charge = block.site.oxidation.expect("Need oxidation state to run the code...");
        block.positions.iter().enumerate().map(move |(i, p)| UnitCellAtom { block: b, atom: i, position: *p, charge })
    }).collect()
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Distance between `x` and the nearest point of `offset + Z`.
fn periodic_offset(x: f64, offset: f64) -> f64 {
    let d = x - offset;
    d - d.round()
}

/// Cell displacements making up a cluster of the given shape.
pub fn select_cells(shape: &ClusterShape, cell: &CellParameters, blocks: &[WyckoffBlock]) -> Vec<[i32; 3]> {
    match shape {
        ClusterShape::Cuboid(n) => {
            let range = |a: usize| -(n[a][0] as i32)..=(n[a][1] as i32);
            let mut cells = Vec::with_capacity(n.iter().map(|r| (r[0] + r[1] + 1) as usize).product());
            for i in range(0) {
                for j in range(1) {
                    for k in range(2) {
                        cells.push([i, j, k]);
                    }
                }
            }
            cells
        }
        ClusterShape::Sphere { radius, reference } => {
            let atoms = unit_cell_atoms(blocks);
            let mut centroid = [0.0; 3];
            for atom in atoms.iter() {
                for a in 0..3 {
                    centroid[a] += atom.position[a] / atoms.len() as f64;
                }
            }
            let reference = reference.unwrap_or(centroid);
            // The sphere spans at most radius / spacing cells along each lattice vector.
            let n = cell.plane_spacings().map(|d| (radius / d).ceil() as i32 + 1);
            let mut cells = Vec::new();
            for i in -n[0]..=n[0] {
                for j in -n[1]..=n[1] {
                    for k in -n[2]..=n[2] {
                        let d = [i, j, k];
                        let r = cell.to_cartesian(std::array::from_fn(|a| centroid[a] + d[a] as f64 - reference[a]));
                        if norm(r) <= radius + 1.0e-8 {
                            cells.push(d);
                        }
                    }
                }
            }
            cells
        }
    }
}

/// The sublattices `offset + Z^3` the auxiliary charges sit on.
struct AuxLattice {
    /// Distinct fractional offsets in [0, 1) along each axis; one or two each.
    offsets: [Vec<f64>; 3],
}

impl AuxLattice {
    fn from_charges(aux: &[PointCharge]) -> Self {
        let mut offsets: [Vec<f64>; 3] = Default::default();
        for c in aux.iter() {
            for a in 0..3 {
                let f = wrap(c.position[a]);
                if !offsets[a].iter().any(|&o| periodic_offset(f, o).abs() < SUBLATTICE_TOLERANCE) {
                    offsets[a].push(f);
                }
            }
        }
        Self { offsets }
    }

    fn n_sublattices(&self) -> usize {
        self.offsets.iter().map(|o| o.len()).product()
    }

    /// Per-axis offset indices of sublattice `s`.
    fn split(&self, s: usize) -> [usize; 3] {
        let n0 = self.offsets[0].len();
        let n1 = self.offsets[1].len();
        [s % n0, (s / n0) % n1, s / (n0 * n1)]
    }

    fn offset(&self, s: usize) -> [f64; 3] {
        let j = self.split(s);
        std::array::from_fn(|a| self.offsets[a][j[a]])
    }

    /// Sublattice and integer co-ordinates of a position known to be on the lattice.
    fn locate(&self, p: [f64; 3]) -> (usize, [i32; 3]) {
        let j: [usize; 3] = std::array::from_fn(|a| {
            self.offsets[a].iter()
                .position(|&o| periodic_offset(p[a], o).abs() < SUBLATTICE_TOLERANCE)
                .expect("auxiliary charge off its own lattice")
        });
        let n0 = self.offsets[0].len();
        let n1 = self.offsets[1].len();
        let s = j[0] + n0 * (j[1] + n1 * j[2]);
        let off = self.offset(s);
        (s, std::array::from_fn(|a| (p[a] - off[a]).round() as i32))
    }
}

/// Dense accumulation grid for the auxiliary field: one layer per sublattice over an integer
/// box `lo..lo+dims`.
struct AuxGrid {
    lo: [i32; 3],
    dims: [usize; 3],
    n_sublattices: usize,
    charges: Vec<f64>,
}

impl AuxGrid {
    fn new(lo: [i32; 3], hi: [i32; 3], n_sublattices: usize) -> Self {
        let dims: [usize; 3] = std::array::from_fn(|a| (hi[a] - lo[a] + 1) as usize);
        let charges = vec![0.0; n_sublattices * dims.iter().product::<usize>()];
        Self { lo, dims, n_sublattices, charges }
    }

    fn index(&self, s: usize, n: [i32; 3]) -> usize {
        let i: [usize; 3] = std::array::from_fn(|a| (n[a] - self.lo[a]) as usize);
        s + self.n_sublattices * (i[0] + self.dims[0] * (i[1] + self.dims[1] * i[2]))
    }

    /// Inverse of [`AuxGrid::index`].
    fn point(&self, idx: usize) -> (usize, [i32; 3]) {
        let s = idx % self.n_sublattices;
        let rest = idx / self.n_sublattices;
        let i = [rest % self.dims[0], (rest / self.dims[0]) % self.dims[1], rest / (self.dims[0] * self.dims[1])];
        (s, std::array::from_fn(|a| self.lo[a] + i[a] as i32))
    }
}

/// A unit-cell atom sitting on an auxiliary sublattice point.
struct Coincidence {
    /// Index into the flattened unit-cell atoms.
    atom: usize,
    sublattice: usize,
    /// Integer co-ordinates of the sublattice point, in the cell at the origin.
    n: [i32; 3],
}

/// Finds the unit-cell atoms that sit on an auxiliary sublattice point, and warns about atoms
/// that come close to one. Both atoms and sublattices are periodic, so this holds in every cell.
fn find_coincidences(atoms: &[UnitCellAtom], lattice: &AuxLattice, cell: &CellParameters, blocks: &[WyckoffBlock]) -> Vec<Coincidence> {
    let mut coincidences = Vec::new();
    for (i, atom) in atoms.iter().enumerate() {
        for s in 0..lattice.n_sublattices() {
            let off = lattice.offset(s);
            let delta: [f64; 3] = std::array::from_fn(|a| periodic_offset(atom.position[a], off[a]));
            let dist = norm(cell.to_cartesian(delta));
            if dist < COINCIDENCE_TOLERANCE {
                let n = std::array::from_fn(|a| (atom.position[a] - off[a]).round() as i32);
                coincidences.push(Coincidence { atom: i, sublattice: s, n });
            } else if dist < NEAR_MISS_TOLERANCE {
                eprintln!(
                    "warning: {} at {:?} is {dist:.3} Angstrom from an auxiliary charge site; consider another shift",
                    blocks[atom.block].site.label, atom.position
                );
            }
        }
    }
    coincidences
}

/// Builds the cluster made of `cells` copies of the unit cell (`blocks`), complemented with the
/// auxiliary charges `aux_cell` of each copy.
pub fn build_cluster(blocks: &[WyckoffBlock], aux_cell: &[PointCharge], cells: Vec<[i32; 3]>, cell: &CellParameters) -> EmbeddingCluster {
    let atoms_cell = unit_cell_atoms(blocks);
    let lattice = AuxLattice::from_charges(aux_cell);
    let aux_located: Vec<(usize, [i32; 3], f64)> = aux_cell.iter()
        .map(|c| {
            let (s, n) = lattice.locate(c.position);
            (s, n, c.charge)
        })
        .collect();
    let coincidences = find_coincidences(&atoms_cell, &lattice, cell, blocks);

    // Grid bounds: cell displacements plus the extent of the auxiliary charges of one cell.
    let mut lo = [i32::MAX; 3];
    let mut hi = [i32::MIN; 3];
    for d in cells.iter() {
        for (_, n, _) in aux_located.iter() {
            for a in 0..3 {
                lo[a] = lo[a].min(n[a] + d[a]);
                hi[a] = hi[a].max(n[a] + d[a]);
            }
        }
    }

    // Pass 1: accumulate the auxiliary field.
    let mut grid = AuxGrid::new(lo, hi, lattice.n_sublattices());
    for d in cells.iter() {
        for (s, n, q) in aux_located.iter() {
            let idx = grid.index(*s, std::array::from_fn(|a| n[a] + d[a]));
            grid.charges[idx] += q;
        }
    }

    // Pass 2: the main field, with auxiliary charges sitting on atoms moved onto them.
    let mut atoms = Vec::with_capacity(cells.len() * atoms_cell.len());
    for d in cells.iter() {
        let first = atoms.len();
        atoms.extend(atoms_cell.iter().map(|a| ClusterAtom {
            block: a.block,
            atom: a.atom,
            cell: *d,
            position: std::array::from_fn(|k| a.position[k] + d[k] as f64),
            charge: a.charge,
            region: Region::PointCharge,
        }));
        for c in coincidences.iter() {
            let idx = grid.index(c.sublattice, std::array::from_fn(|a| c.n[a] + d[a]));
            atoms[first + c.atom].charge += grid.charges[idx];
            grid.charges[idx] = 0.0;
        }
    }

    // Prune cancelled auxiliary charges.
    let max_aux = aux_cell.iter().map(|c| c.charge.abs()).fold(0.0, f64::max);
    let tolerance = PRUNE_RELATIVE_TOLERANCE * max_aux;
    let n_aux = grid.charges.iter().filter(|q| q.abs() > tolerance).count();
    let mut aux = Vec::with_capacity(n_aux);
    for (idx, q) in grid.charges.iter().enumerate() {
        if q.abs() > tolerance {
            let (s, n) = grid.point(idx);
            let off = lattice.offset(s);
            aux.push(PointCharge { charge: *q, position: std::array::from_fn(|a| off[a] + n[a] as f64) });
        }
    }

    EmbeddingCluster { cells, atoms, aux, qm: Vec::new(), aimp: Vec::new() }
}

impl EmbeddingCluster {
    /// Marks the atoms matching `qm_atoms` (Cartesian, Angstrom, in the frame of `cell`) as the
    /// QM region, keeping the order of `qm_atoms`.
    ///
    /// Every QM atom must match exactly one cluster atom of the same element. Which atoms make
    /// up a sensible QM region is left to the user; no completion by Wyckoff site is done.
    pub fn mark_qm_region(&mut self, qm_atoms: &[XyzAtom], blocks: &[WyckoffBlock], cell: &CellParameters) -> Result<(), Box<dyn std::error::Error>> {
        let atoms_cell = unit_cell_atoms(blocks);
        let mut qm: Vec<usize> = Vec::with_capacity(qm_atoms.len());
        for qm_atom in qm_atoms.iter() {
            let frac = cell.to_fractional(qm_atom.position);
            // Find the unit-cell atom and the cell displacement that put an atom on `frac`.
            let found = atoms_cell.iter().enumerate().find_map(|(i, a)| {
                let delta: [f64; 3] = std::array::from_fn(|k| frac[k] - a.position[k]);
                let d = delta.map(|x| x.round());
                let residual = std::array::from_fn(|k| delta[k] - d[k]);
                (norm(cell.to_cartesian(residual)) < QM_MATCH_TOLERANCE).then(|| (i, d.map(|x| x as i32)))
            });
            let (i, d) = found.ok_or_else(|| format!(
                "QM atom {} at {:?} Angstrom does not match any atom of the crystal", qm_atom.label, qm_atom.position
            ))?;
            let cell_idx = self.cells.iter().position(|c| *c == d).ok_or_else(|| format!(
                "QM atom {} at {:?} Angstrom lies in cell {d:?}, outside the cluster", qm_atom.label, qm_atom.position
            ))?;
            let idx = cell_idx * atoms_cell.len() + i;
            let element = &blocks[atoms_cell[i].block].site.element;
            if *element != qm_atom.element {
                return Err(format!(
                    "QM atom {} at {:?} Angstrom sits on a {element} site ({})",
                    qm_atom.label, qm_atom.position, blocks[atoms_cell[i].block].site.label
                ).into());
            }
            if qm.contains(&idx) {
                return Err(format!("QM atom {} at {:?} Angstrom is given twice", qm_atom.label, qm_atom.position).into());
            }
            qm.push(idx);
        }
        self.mark_qm_atoms(qm);
        Ok(())
    }

    /// Makes `indices` (into `atoms`, in output order) the QM region, replacing any earlier one.
    pub fn mark_qm_atoms(&mut self, indices: Vec<usize>) {
        for &i in self.qm.iter() {
            self.atoms[i].region = Region::PointCharge;
        }
        for &i in indices.iter() {
            self.atoms[i].region = Region::Qm;
        }
        self.qm = indices;
    }

    /// Index into `atoms` of position `atom` of Wyckoff block `block` in cell `cell`, if that cell
    /// is part of the cluster.
    pub fn atom_index(&self, cell: [i32; 3], block: usize, atom: usize, blocks: &[WyckoffBlock]) -> Option<usize> {
        let n_atoms: usize = blocks.iter().map(|b| b.positions.len()).sum();
        let offset: usize = blocks[..block].iter().map(|b| b.positions.len()).sum();
        self.cells.iter().position(|c| *c == cell).map(|c| c * n_atoms + offset + atom)
    }

    /// Marks the atoms of `region` as AIMPs, except those already in the QM region. Call after
    /// [`EmbeddingCluster::mark_qm_region`].
    ///
    /// The AIMPs are grouped by Wyckoff block, since each crystallographic site gets its own
    /// model potential.
    pub fn mark_aimp_region(&mut self, region: &AimpRegion, cell: &CellParameters) {
        let reference = match region {
            AimpRegion::Radius { reference: Some(r), .. } => *r,
            _ if !self.qm.is_empty() => {
                let mut centroid = [0.0; 3];
                for &i in self.qm.iter() {
                    for a in 0..3 {
                        centroid[a] += self.atoms[i].position[a] / self.qm.len() as f64;
                    }
                }
                centroid
            }
            _ => [0.5; 3],
        };
        let distance = |p: &[f64; 3]| norm(cell.to_cartesian(std::array::from_fn(|a| p[a] - reference[a])));

        if let AimpRegion::Cells(n) = region {
            let outside = (-(n[0][0] as i32)..=n[0][1] as i32)
                .flat_map(|i| (-(n[1][0] as i32)..=n[1][1] as i32).flat_map(move |j| (-(n[2][0] as i32)..=n[2][1] as i32).map(move |k| [i, j, k])))
                .filter(|d| !self.cells.contains(d))
                .count();
            if outside > 0 {
                eprintln!("warning: {outside} cells of the AIMP region lie outside the cluster");
            }
        }
        let in_region = |atom: &ClusterAtom| match region {
            AimpRegion::Cells(n) => (0..3).all(|a| -(n[a][0] as i32) <= atom.cell[a] && atom.cell[a] <= n[a][1] as i32),
            AimpRegion::Radius { radius, .. } => distance(&atom.position) <= radius + 1.0e-8,
        };

        let mut aimp: Vec<usize> = Vec::new();
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            if atom.region == Region::Qm {
                continue;
            }
            if in_region(atom) {
                atom.region = Region::Aimp;
                aimp.push(i);
            } else if atom.region == Region::Aimp {
                atom.region = Region::PointCharge;
            }
        }
        let atoms = &self.atoms;
        aimp.sort_by(|&i, &j| {
            let (a, b) = (&atoms[i], &atoms[j]);
            a.block.cmp(&b.block)
                .then(distance(&a.position).total_cmp(&distance(&b.position)))
                .then(i.cmp(&j))
        });
        self.aimp = aimp;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cif;
    use crate::crystal::Crystal;
    use crate::tetrahedron::compute_auxiliary_charges;

    /// The grid construction must give the same result as copying every charge and merging
    /// overlapping positions by brute force.
    fn compare_with_brute_force(shift: f64) {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        let blocks = crystal.unfold(1.0e-2);
        let aux_cell = compute_auxiliary_charges(&blocks, 4, shift).unwrap();
        let cells = select_cells(&ClusterShape::Cuboid([[1, 1]; 3]), &crystal.cell, &blocks);
        let cluster = build_cluster(&blocks, &aux_cell, cells.clone(), &crystal.cell);

        let same = |p: &[f64; 3], q: &[f64; 3]| (0..3).all(|a| (p[a] - q[a]).abs() < 1.0e-6);
        let mut merged: Vec<PointCharge> = Vec::new();
        for d in cells.iter() {
            for c in aux_cell.iter() {
                let p = std::array::from_fn(|a| c.position[a] + d[a] as f64);
                match merged.iter_mut().find(|m| same(&m.position, &p)) {
                    Some(m) => m.charge += c.charge,
                    None => merged.push(PointCharge { charge: c.charge, position: p }),
                }
            }
        }
        let mut atom_charges: Vec<f64> = cluster.atoms.iter().map(|a| blocks[a.block].site.oxidation.unwrap()).collect();
        for m in merged.iter_mut() {
            if let Some(i) = cluster.atoms.iter().position(|a| same(&a.position, &m.position)) {
                atom_charges[i] += m.charge;
                m.charge = 0.0;
            }
        }
        merged.retain(|m| m.charge.abs() > 1.0e-10);

        assert_eq!(merged.len(), cluster.aux.len());
        for m in merged.iter() {
            let c = cluster.aux.iter().find(|c| same(&c.position, &m.position)).expect("missing aux charge");
            assert!((c.charge - m.charge).abs() < 1.0e-10);
        }
        for (a, q) in cluster.atoms.iter().zip(atom_charges.iter()) {
            assert!((a.charge - q).abs() < 1.0e-10);
        }
        let total: f64 = cluster.charges().map(|(q, _)| q).sum();
        assert!(total.abs() < 1.0e-9);
    }

    #[test]
    fn grid_matches_brute_force_unshifted() {
        compare_with_brute_force(0.0);
    }

    #[test]
    fn sphere_selects_cells_within_radius() {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        let blocks = crystal.unfold(1.0e-2);
        // a = 10.6 Angstrom (cubic): radius 25 includes the cells with |d| <= 25/a = 2.36.
        let a = crystal.cell.a;
        let cells = select_cells(&ClusterShape::Sphere { radius: 25.0, reference: None }, &crystal.cell, &blocks);
        let expected = (-3..=3i32).flat_map(|i| (-3..=3i32).flat_map(move |j| (-3..=3i32).map(move |k| [i, j, k])))
            .filter(|d| a * ((d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) as f64).sqrt() <= 25.0)
            .count();
        assert_eq!(cells.len(), expected);
        assert!(cells.contains(&[0, 0, 0]) && cells.contains(&[2, 0, 0]) && !cells.contains(&[2, 2, 0]));
    }

    #[test]
    fn qm_region_keeps_input_order() {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        let blocks = crystal.unfold(1.0e-2);
        let aux_cell = compute_auxiliary_charges(&blocks, 2, 0.0).unwrap();
        let cells = select_cells(&ClusterShape::Cuboid([[1, 1]; 3]), &crystal.cell, &blocks);
        let mut cluster = build_cluster(&blocks, &aux_cell, cells, &crystal.cell);

        // An O from the cell at (-1, 0, 0) first, then a Y from the origin cell.
        let o = blocks[2].positions[3];
        let y = blocks[0].positions[5];
        let xyz = |label: &str, frac: [f64; 3]| XyzAtom {
            label: label.to_string(),
            element: label.to_string(),
            position: crystal.cell.to_cartesian(frac),
        };
        let qm = [xyz("O", [o[0] - 1.0, o[1], o[2]]), xyz("Y", y)];
        cluster.mark_qm_region(&qm, &blocks, &crystal.cell).unwrap();
        assert_eq!(cluster.qm.len(), 2);
        let first = &cluster.atoms[cluster.qm[0]];
        assert_eq!((first.block, first.atom, first.cell), (2, 3, [-1, 0, 0]));
        let second = &cluster.atoms[cluster.qm[1]];
        assert_eq!((second.block, second.atom, second.cell), (0, 5, [0, 0, 0]));
        assert_eq!(cluster.atoms.iter().filter(|a| a.region == Region::Qm).count(), 2);

        // Wrong element, outside the cluster, off any site, given twice.
        assert!(cluster.mark_qm_region(&[xyz("Y", o)], &blocks, &crystal.cell).is_err());
        assert!(cluster.mark_qm_region(&[xyz("Y", [y[0] + 2.0, y[1], y[2]])], &blocks, &crystal.cell).is_err());
        assert!(cluster.mark_qm_region(&[xyz("Y", [0.123, 0.456, 0.789])], &blocks, &crystal.cell).is_err());
        let mut fresh = build_cluster(&blocks, &aux_cell, cluster.cells.clone(), &crystal.cell);
        assert!(fresh.mark_qm_region(&[xyz("Y", y), xyz("Y", y)], &blocks, &crystal.cell).is_err());
    }

    #[test]
    fn aimp_region_grouped_by_wyckoff_block() {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        let blocks = crystal.unfold(1.0e-2);
        let aux_cell = compute_auxiliary_charges(&blocks, 2, 0.0).unwrap();
        let cells = select_cells(&ClusterShape::Cuboid([[2, 2]; 3]), &crystal.cell, &blocks);
        let mut cluster = build_cluster(&blocks, &aux_cell, cells, &crystal.cell);
        let y = blocks[0].positions[0];
        let qm = [XyzAtom { label: "Y".into(), element: "Y".into(), position: crystal.cell.to_cartesian(y) }];
        cluster.mark_qm_region(&qm, &blocks, &crystal.cell).unwrap();

        // Whole cells: the origin cell and its +a neighbour, minus the QM atom.
        let n_cell: usize = blocks.iter().map(|b| b.positions.len()).sum();
        cluster.mark_aimp_region(&AimpRegion::Cells([[0, 1], [0, 0], [0, 0]]), &crystal.cell);
        assert_eq!(cluster.aimp.len(), 2 * n_cell - 1);
        assert!(!cluster.aimp.contains(&cluster.qm[0]));

        // Radius around the QM atom: every atom within 6 Angstrom, grouped by block.
        cluster.mark_aimp_region(&AimpRegion::Radius { radius: 6.0, reference: None }, &crystal.cell);
        let dist = |p: [f64; 3]| norm(crystal.cell.to_cartesian(std::array::from_fn(|a| p[a] - y[a])));
        let expected = cluster.atoms.iter().filter(|a| a.region != Region::Qm && dist(a.position) <= 6.0).count();
        assert_eq!(cluster.aimp.len(), expected);
        assert_eq!(cluster.atoms.iter().filter(|a| a.region == Region::Aimp).count(), expected);
        let block_order: Vec<usize> = cluster.aimp.iter().map(|&i| cluster.atoms[i].block).collect();
        assert!(block_order.windows(2).all(|w| w[0] <= w[1]));
        for w in cluster.aimp.windows(2) {
            let (a, b) = (&cluster.atoms[w[0]], &cluster.atoms[w[1]]);
            assert!(a.block < b.block || dist(a.position) <= dist(b.position));
        }
    }

    #[test]
    fn grid_matches_brute_force_shifted() {
        // Y (8b) at (1/4, 1/4, 1/4) sits on an auxiliary sublattice for V = 1/4.
        compare_with_brute_force(0.25);
    }
}
