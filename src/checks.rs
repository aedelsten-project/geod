// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Consistency checks of the embedding cluster.
//!
//! Two kinds:
//!
//! - Structural (errors): the QM and AIMP atoms must carry exactly their formal charges. Interior
//!   auxiliary charges cancel exactly, so an auxiliary charge folded onto a QM or AIMP atom means
//!   that region reaches the boundary shell of the cluster. The QM calculation replaces those
//!   atoms by their formal charges, so the folded charge would be lost silently.
//! - Accuracy (report): the potential and field of the cluster at the QM (and AIMP) sites,
//!   compared with the Ewald potential of the infinite crystal, as in Sushko & Abarenkov,
//!   section 3. The two may differ by a constant (which does not affect the QM calculation), so
//!   the mean offset is reported separately from the spread.

use std::error::Error;

use crate::cluster::{EmbeddingCluster, Region};
use crate::crystal::{CellParameters, WyckoffBlock, ANGSTROM_TO_BOHR};
use crate::ewald::Ewald;

/// Relative tolerance for "carries exactly its formal charge"; folded auxiliary charges in the
/// interior cancel to rounding noise.
const FORMAL_CHARGE_TOLERANCE: f64 = 1.0e-8;

/// Fails if any QM or AIMP atom carries auxiliary charge on top of its formal charge.
pub fn check_regions_free_of_aux(cluster: &EmbeddingCluster, blocks: &[WyckoffBlock]) -> Result<(), Box<dyn Error>> {
    for (name, region, indices) in [("QM", Region::Qm, &cluster.qm), ("AIMP", Region::Aimp, &cluster.aimp)] {
        let offending: Vec<usize> = indices.iter().copied()
            .filter(|&i| {
                let atom = &cluster.atoms[i];
                debug_assert!(atom.region == region);
                let formal = blocks[atom.block].site.oxidation.expect("Need oxidation state to run the code...");
                (atom.charge - formal).abs() > FORMAL_CHARGE_TOLERANCE * formal.abs().max(1.0)
            })
            .collect();
        if let Some(&first) = offending.first() {
            let atom = &cluster.atoms[first];
            return Err(format!(
                "{} {name} atoms carry auxiliary charge (e.g. {} position {} in cell {:?}, charge {:.6}): the {name} region reaches the boundary of the cluster; enlarge the cluster or shrink the region",
                offending.len(), blocks[atom.block].site.label, atom.atom, atom.cell, atom.charge
            ).into());
        }
    }
    Ok(())
}

/// Potential and field (atomic units) of all cluster charges at atom `site`, leaving out the atom
/// itself.
fn cluster_potential_and_field(cluster: &EmbeddingCluster, site: usize, cell: &CellParameters) -> (f64, [f64; 3]) {
    let to_bohr = |p: [f64; 3]| cell.to_cartesian(p).map(|x| x * ANGSTROM_TO_BOHR);
    let r = to_bohr(cluster.atoms[site].position);
    let mut v = 0.0;
    let mut e = [0.0; 3];
    let others = cluster.atoms.iter().enumerate()
        .filter(|(i, _)| *i != site)
        .map(|(_, a)| (a.charge, a.position))
        .chain(cluster.aux.iter().map(|c| (c.charge, c.position)));
    for (q, p) in others {
        let rj = to_bohr(p);
        let d: [f64; 3] = std::array::from_fn(|c| r[c] - rj[c]);
        let dist = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        v += q / dist;
        for c in 0..3 {
            e[c] += q * d[c] / (dist * dist * dist);
        }
    }
    (v, e)
}

/// Differences between the cluster and the infinite crystal at a set of sites.
struct Deviations {
    /// Cluster (on-site) potential at each site.
    v: Vec<f64>,
    /// Cluster minus Ewald potential at each site.
    dv: Vec<f64>,
    /// |cluster minus Ewald field| at each site.
    de: Vec<f64>,
}

impl Deviations {
    fn compute(cluster: &EmbeddingCluster, sites: &[usize], ewald: &Ewald, cell: &CellParameters) -> Self {
        let mut v = Vec::with_capacity(sites.len());
        let mut dv = Vec::with_capacity(sites.len());
        let mut de = Vec::with_capacity(sites.len());
        for &i in sites.iter() {
            let (v_cluster, e_cluster) = cluster_potential_and_field(cluster, i, cell);
            let r = cell.to_cartesian(cluster.atoms[i].position).map(|x| x * ANGSTROM_TO_BOHR);
            let (v_ewald, e_ewald) = ewald.potential_and_field(r);
            v.push(v_cluster);
            dv.push(v_cluster - v_ewald);
            let d: [f64; 3] = std::array::from_fn(|c| e_cluster[c] - e_ewald[c]);
            de.push((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
        }
        Self { v, dv, de }
    }

    fn offset(&self) -> f64 {
        self.dv.iter().sum::<f64>() / self.dv.len() as f64
    }

    /// RMS of the potential deviation after removing the mean offset.
    fn rms_spread(&self) -> f64 {
        let offset = self.offset();
        (self.dv.iter().map(|x| (x - offset).powi(2)).sum::<f64>() / self.dv.len() as f64).sqrt()
    }

    fn max_field(&self) -> f64 {
        self.de.iter().cloned().fold(0.0, f64::max)
    }
}

/// Reports how well the cluster reproduces the potential and field of the infinite crystal at the
/// QM and AIMP sites.
pub fn report_embedding_accuracy(cluster: &EmbeddingCluster, blocks: &[WyckoffBlock], cell: &CellParameters) {
    let ewald = Ewald::new(blocks, cell);
    println!("Embedding accuracy (cluster minus infinite crystal, atomic units):");
    if !cluster.qm.is_empty() {
        let qm = Deviations::compute(cluster, &cluster.qm, &ewald, cell);
        println!("  QM sites:");
        println!("    {:<6} {:>16} {:>16} {:>16}", "Site", "V (Eh/e)", "dV (Eh/e)", "|dE| (Eh/e/a0)");
        for (k, &i) in cluster.qm.iter().enumerate() {
            let label = &blocks[cluster.atoms[i].block].site.label;
            println!("    {label:<6} {:>16.10} {:>16.3e} {:>16.3e}", qm.v[k], qm.dv[k], qm.de[k]);
        }
        println!("    Mean offset        = {:12.3e} Eh/e", qm.offset());
        println!("    RMS spread of dV   = {:12.3e} Eh/e", qm.rms_spread());
        println!("    Max |dE|           = {:12.3e} Eh/(e a0)", qm.max_field());
    }
    if !cluster.aimp.is_empty() {
        let aimp = Deviations::compute(cluster, &cluster.aimp, &ewald, cell);
        println!("  AIMP sites ({}):", cluster.aimp.len());
        println!("    Mean offset        = {:12.3e} Eh/e", aimp.offset());
        println!("    RMS spread of dV   = {:12.3e} Eh/e", aimp.rms_spread());
        println!("    Max |dE|           = {:12.3e} Eh/(e a0)", aimp.max_field());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cif;
    use crate::cluster::{build_cluster, select_cells, AimpRegion, ClusterShape};
    use crate::crystal::Crystal;
    use crate::tetrahedron::compute_auxiliary_charges;
    use crate::xyz::XyzAtom;

    fn y2o3_cluster(radius: f64, aimp_radius: f64) -> (EmbeddingCluster, Vec<WyckoffBlock>, CellParameters) {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        let blocks = crystal.unfold(1.0e-2);
        let aux_cell = compute_auxiliary_charges(&blocks, 4, 0.0).unwrap();
        let cells = select_cells(&ClusterShape::Sphere { radius, reference: None }, &crystal.cell, &blocks);
        let mut cluster = build_cluster(&blocks, &aux_cell, cells, &crystal.cell);
        let y = blocks[0].positions[0];
        let qm = [XyzAtom { label: "Y".into(), element: "Y".into(), position: crystal.cell.to_cartesian(y) }];
        cluster.mark_qm_region(&qm, &blocks, &crystal.cell).unwrap();
        cluster.mark_aimp_region(&AimpRegion::Radius { radius: aimp_radius, reference: None }, &crystal.cell);
        (cluster, blocks, crystal.cell)
    }

    #[test]
    fn aimp_region_reaching_the_boundary_shell_is_an_error() {
        // M = 4: the shell of uncancelled auxiliary charge is about four cells thick.
        let (cluster, blocks, _) = y2o3_cluster(50.0, 15.0);
        assert!(check_regions_free_of_aux(&cluster, &blocks).is_err());
        let (cluster, blocks, _) = y2o3_cluster(60.0, 15.0);
        assert!(check_regions_free_of_aux(&cluster, &blocks).is_ok());
    }

    #[test]
    fn cluster_reproduces_infinite_crystal() {
        let (cluster, blocks, cell) = y2o3_cluster(60.0, 5.0);
        let ewald = Ewald::new(&blocks, &cell);
        let qm = Deviations::compute(&cluster, &cluster.qm, &ewald, &cell);
        assert!(qm.dv[0].abs() < 1.0e-3 && qm.de[0] < 1.0e-4);
        let aimp = Deviations::compute(&cluster, &cluster.aimp, &ewald, &cell);
        assert!(aimp.rms_spread() < 1.0e-3 && aimp.max_field() < 1.0e-3);
    }
}
