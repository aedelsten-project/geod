// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Input files for the quantum-chemistry code: the QM and AIMP regions as basis-set blocks
//! (`.inp`) and the point-charge field (`.xfield`). Both follow the OpenMolcas conventions, so the
//! files work with OpenMolcas as well.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use crate::cluster::{EmbeddingCluster, Region};
use crate::crystal::{CellParameters, WyckoffBlock, ANGSTROM_TO_BOHR};
use crate::input::BasisInput;

/// One `label x y z` line, position in Bohr (the default unit, so none is given). Six
/// decimals is already far beyond the precision of the crystal structure.
fn write_centre(file: &mut impl Write, label: &str, position: [f64; 3], cell: &CellParameters) -> io::Result<()> {
    let r = cell.to_cartesian(position).map(|x| x * ANGSTROM_TO_BOHR);
    writeln!(file, "    {label:<8}{:16.6}{:16.6}{:16.6}", r[0], r[1], r[2])
}

/// Writes the QM region followed by the AIMP region as `Basis set` blocks.
///
/// QM atoms are written in the order of the QM input. A new basis block starts whenever the
/// element changes from the previous atom; atoms are numbered per element (Y1, Y2, O1, ...) across
/// blocks.
///
/// AIMPs get one basis block per Wyckoff site, with the site's own embedding potential
/// (by default `{element}.EMB-AIMP.author.0s.0s.ECP.{site label}`), and are numbered A1, A2, ...
/// throughout.
pub fn write_inp(path: impl AsRef<Path>, cluster: &EmbeddingCluster, blocks: &[WyckoffBlock], cell: &CellParameters, basis: &BasisInput) -> io::Result<()> {
    let mut file = BufWriter::new(File::create(path)?);

    let mut counts: HashMap<&str, usize> = HashMap::new();
    let mut current: Option<&str> = None;
    for &i in cluster.qm.iter() {
        let atom = &cluster.atoms[i];
        let site = &blocks[atom.block].site;
        let element = site.element.as_str();
        if current != Some(element) {
            if current.is_some() {
                writeln!(file, "  End of Basis")?;
            }
            writeln!(file, "  Basis set")?;
            writeln!(file, "    {}", basis.qm_label(element, &site.label))?;
            current = Some(element);
        }
        let n = counts.entry(element).or_insert(0);
        *n += 1;
        write_centre(&mut file, &format!("{element}{n}"), atom.position, cell)?;
    }
    if current.is_some() {
        writeln!(file, "  End of Basis")?;
    }

    // `cluster.aimp` is grouped by Wyckoff block, so each block is opened once.
    let mut current_block: Option<usize> = None;
    for (n, &i) in cluster.aimp.iter().enumerate() {
        let atom = &cluster.atoms[i];
        if current_block != Some(atom.block) {
            if current_block.is_some() {
                writeln!(file, "  End of Basis")?;
            }
            let site = &blocks[atom.block].site;
            writeln!(file, "  Basis set")?;
            writeln!(file, "    {}", basis.aimp_label(&site.element, &site.label))?;
            writeln!(file, "    pseudocharge")?;
            current_block = Some(atom.block);
        }
        write_centre(&mut file, &format!("A{}", n + 1), atom.position, cell)?;
    }
    if current_block.is_some() {
        writeln!(file, "  End of Basis")?;
    }

    file.flush()
}

/// Writes the point-charge field (OpenMolcas-style `XField`): all atoms outside the QM and AIMP regions,
/// followed by the auxiliary charges.
///
/// Header `N 0` (N charges, multipole order 0), then `x y z q` per line, positions in Bohr. The
/// charges get 12 decimals: the auxiliary charges are large and only cancel the cell moments in
/// sums over many of them.
pub fn write_xfield(path: impl AsRef<Path>, cluster: &EmbeddingCluster, cell: &CellParameters) -> io::Result<()> {
    let mut file = BufWriter::new(File::create(path)?);
    let point_charges = || cluster.atoms.iter()
        .filter(|a| a.region == Region::PointCharge)
        .map(|a| (a.charge, a.position))
        .chain(cluster.aux.iter().map(|c| (c.charge, c.position)));
    writeln!(file, "{} 0", point_charges().count())?;
    for (q, p) in point_charges() {
        let r = cell.to_cartesian(p).map(|x| x * ANGSTROM_TO_BOHR);
        writeln!(file, "{:16.6}{:16.6}{:16.6}{:20.12}", r[0], r[1], r[2], q)?;
    }
    file.flush()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::cif;
    use crate::cluster::{build_cluster, select_cells, AimpRegion, ClusterShape};
    use crate::crystal::Crystal;
    use crate::tetrahedron::compute_auxiliary_charges;
    use crate::xyz::XyzAtom;

    #[test]
    fn blocks_follow_input_order() {
        let blocks = cif::parse(include_str!("../tests/data/Y2O3.cif")).unwrap();
        let crystal = Crystal::from_block(&blocks[0]).unwrap();
        let blocks = crystal.unfold(1.0e-2);
        let aux_cell = compute_auxiliary_charges(&blocks, 2, 0.0).unwrap();
        let cells = select_cells(&ClusterShape::Cuboid([[0, 0]; 3]), &crystal.cell, &blocks);
        let mut cluster = build_cluster(&blocks, &aux_cell, cells, &crystal.cell);
        let xyz = |label: &str, frac: [f64; 3]| XyzAtom {
            label: label.to_string(),
            element: label.to_string(),
            position: crystal.cell.to_cartesian(frac),
        };
        let qm = [
            xyz("Y", blocks[0].positions[0]),
            xyz("O", blocks[2].positions[0]),
            xyz("O", blocks[2].positions[1]),
            xyz("Y", blocks[1].positions[0]),
        ];
        cluster.mark_qm_region(&qm, &blocks, &crystal.cell).unwrap();
        cluster.mark_aimp_region(&AimpRegion::Cells([[0, 0]; 3]), &crystal.cell);

        let path = std::env::temp_dir().join(format!("output_test_{}.inp", std::process::id()));
        write_inp(&path, &cluster, &blocks, &crystal.cell, &BasisInput::default()).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();

        let labels: Vec<&str> = text.lines()
            .filter(|l| l.starts_with("    ") && !l.contains("ANO") && !l.contains("AIMP") && !l.contains("pseudocharge"))
            .map(|l| l.split_whitespace().next().unwrap())
            .collect();
        let n_cell: usize = blocks.iter().map(|b| b.positions.len()).sum();
        let aimp_labels: Vec<String> = (1..=n_cell - 4).map(|n| format!("A{n}")).collect();
        assert_eq!(labels[..4], ["Y1", "O1", "O2", "Y2"]);
        assert_eq!(labels[4..], aimp_labels);
        // Three QM blocks, then one AIMP block per Wyckoff site.
        assert_eq!(text.matches("Basis set").count(), 3 + blocks.len());
        assert_eq!(text.matches("End of Basis").count(), 3 + blocks.len());
        for block in blocks.iter() {
            let header = format!("    {}.EMB-AIMP.author.0s.0s.ECP.{}\n    pseudocharge\n", block.site.element, block.site.label);
            assert_eq!(text.matches(&header).count(), 1);
        }
        assert!(!text.contains("Angstrom") && !text.contains("AIMPLIB"));

        // The field holds everything else: one line per point charge, x y z q.
        write_xfield(&path, &cluster, &crystal.cell).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(&path).unwrap();
        let mut lines = text.lines();
        let n_field = cluster.atoms.len() - cluster.qm.len() - cluster.aimp.len() + cluster.aux.len();
        assert_eq!(lines.next().unwrap(), format!("{n_field} 0"));
        let rows: Vec<Vec<f64>> = lines.map(|l| l.split_whitespace().map(|x| x.parse().unwrap()).collect()).collect();
        assert_eq!(rows.len(), n_field);
        assert!(rows.iter().all(|r| r.len() == 4));
        // QM, AIMP and field together are neutral.
        let field: f64 = rows.iter().map(|r| r[3]).sum();
        let regions: f64 = cluster.qm.iter().chain(cluster.aimp.iter()).map(|&i| cluster.atoms[i].charge).sum();
        assert!((field + regions).abs() < 1.0e-8);
    }
}
