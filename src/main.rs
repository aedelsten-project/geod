// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Geod: electrostatic embedding of ionic crystals from CIF files: moment-cancelled point-charge fields,
//! QM and AIMP regions, and the inputs for the quantum-chemistry code.
//!
//! The method is that of P. V. Sushko and I. V. Abarenkov, J. Chem. Theory Comput. 6, 1323 (2010),
//! originally implemented by P. V. Sushko in Fortran (EMBQ), then by E. D. Larsson in Julia as
//! part of SCEPIC (with Sushko's permission). This Rust version, part of the Ädelsten project, is
//! by E. D. Larsson, written together with Claude (Anthropic). See README.md and docs/MANUAL.md.

mod checks;
mod cif;
mod cluster;
mod crystal;
mod ewald;
mod input;
mod output;
mod qm_cut;
mod tetrahedron;
mod xyz;

use std::fs::File;
use std::io::{BufWriter, Write};

use std::error::Error;

use crystal::{CellParameters, Crystal, WyckoffBlock};

use crate::cluster::{build_cluster, select_cells, AimpRegion, ClusterShape, EmbeddingCluster};
use crate::input::{Input, Mode, RegionSpec};
use crate::tetrahedron::{compute_auxiliary_charges, print_moment_deviations, PointCharge};

/// Two positions closer than this (Angstrom) count as the same site.
const DUPLICATE_TOLERANCE: f64 = 1.0e-2;
/// Papers to cite when using results of this program.
const CITATIONS: [&str; 2] = [
    "E. D. Larsson, M. Krosnicki, V. Veryazov, Chem. Phys. 562, 111549 (2022), doi:10.1016/j.chemphys.2022.111549",
    "P. V. Sushko, I. V. Abarenkov, J. Chem. Theory Comput. 6, 1323 (2010), doi:10.1021/ct900480p",
];

/// Largest net charge of the unit cell still considered neutral.
const NEUTRALITY_TOLERANCE: f64 = 1.0e-8;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: {} <input.toml>", args[0]);
        std::process::exit(1);
    }
    if let Err(e) = run(&args[1]) {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
    println!("\nPlease cite the following papers:");
    for citation in CITATIONS {
        println!("  {citation}");
    }
}

fn run(input_path: &str) -> Result<(), Box<dyn Error>> {
    let input = Input::read(input_path)?;

    // Generate crystallographic information from the cif-file.
    let (crystal, mut blocks) = generate_unit_cell_data(&input.cif)?;
    input.apply_oxidation(&mut blocks)?;

    if input.mode == Mode::QmCut {
        return write_qm_cut(&input, &blocks, &crystal.cell);
    }

    // Ensure we have a charge neutral system, anything else would be nonsense.
    check_charge_neutrality(&blocks)?;

    // Point charges that cancel the unit-cell multipole moments up to order max_order.
    let aux_charges = compute_auxiliary_charges(&blocks, input.max_order, input.shift)?;
    println!("{} auxiliary charges", aux_charges.len());

    match input.mode {
        Mode::Production => {
            // QM region, given as Cartesian co-ordinates (Angstrom) in the frame of the CIF cell.
            let qm_file = input.qm.as_deref().expect("checked when reading the input");
            let qm_atoms = xyz::read_file(qm_file)?;
            // Centre the cluster on the QM region, where the margin to the boundary matters.
            let mut centre = [0.0; 3];
            for atom in qm_atoms.iter() {
                let f = crystal.cell.to_fractional(atom.position);
                for a in 0..3 {
                    centre[a] += f[a] / qm_atoms.len() as f64;
                }
            }
            let mut cluster = build_embedding(&input, &blocks, &aux_charges, &crystal.cell, centre);
            cluster.mark_qm_region(&qm_atoms, &blocks, &crystal.cell)?;
            finish_embedding(&mut cluster, &input, &blocks, &crystal.cell, &input.project)?;
        }
        Mode::AimpOptimization => {
            // One embedding per Wyckoff site, with one atom of the site in the origin cell as the
            // QM region and the cluster centred on it.
            for (b, block) in blocks.iter().enumerate() {
                println!("\nWyckoff site {}:", block.site.label);
                let mut cluster = build_embedding(&input, &blocks, &aux_charges, &crystal.cell, block.positions[0]);
                let qm = cluster.atom_index([0, 0, 0], b, 0, &blocks)
                    .ok_or("the origin cell is not part of the cluster")?;
                cluster.mark_qm_atoms(vec![qm]);
                let name = format!("{}_{}", input.project, block.site.label);
                finish_embedding(&mut cluster, &input, &blocks, &crystal.cell, &name)?;
            }
        }
        Mode::QmCut => unreachable!("handled before the auxiliary charges"),
    }
    Ok(())
}

/// Builds the cluster (main and auxiliary field) of the requested shape. A spherical cluster is
/// centred on `centre` (fractional); a cuboid one spans the given cells around the origin cell.
fn build_embedding(input: &Input, blocks: &[WyckoffBlock], aux_charges: &[PointCharge], cell: &CellParameters, centre: [f64; 3]) -> EmbeddingCluster {
    let shape = match *input.cluster.as_ref().expect("checked when reading the input") {
        RegionSpec::Radius(radius) => ClusterShape::Sphere { radius, reference: Some(centre) },
        RegionSpec::Cells(n) => ClusterShape::Cuboid(n),
    };
    let cells = select_cells(&shape, cell, blocks);
    let cluster = build_cluster(blocks, aux_charges, cells, cell);
    println!(
        "Cluster: {} cells, {} atoms, {} auxiliary charges",
        cluster.cells.len(), cluster.atoms.len(), cluster.aux.len()
    );
    let total_charge: f64 = cluster.charges().map(|(q, _)| q).sum();
    println!("  Total charge = {total_charge:.3e}");
    print_moment_deviations(|| cluster.charges(), &centre, input.max_order);
    cluster
}

/// Marks the AIMP region around the (already marked) QM region, checks the embedding and writes
/// `{name}.inp` and `{name}.xfield`.
fn finish_embedding(cluster: &mut EmbeddingCluster, input: &Input, blocks: &[WyckoffBlock], cell: &CellParameters, name: &str) -> Result<(), Box<dyn Error>> {
    // AIMPs around the QM region, blocked per Wyckoff site in the output.
    if let Some(aimp) = &input.aimp {
        let region = match *aimp {
            RegionSpec::Radius(radius) => AimpRegion::Radius { radius, reference: None },
            RegionSpec::Cells(n) => AimpRegion::Cells(n),
        };
        cluster.mark_aimp_region(&region, cell);
    }
    println!("QM region: {} atoms, AIMP region: {} atoms", cluster.qm.len(), cluster.aimp.len());

    checks::check_regions_free_of_aux(cluster, blocks)?;
    checks::report_embedding_accuracy(cluster, blocks, cell);

    let inp = format!("{name}.inp");
    output::write_inp(&inp, cluster, blocks, cell, &input.basis).map_err(|e| format!("{inp}: {e}"))?;
    let xfield = format!("{name}.xfield");
    output::write_xfield(&xfield, cluster, cell).map_err(|e| format!("{xfield}: {e}"))?;
    println!("Written {inp} and {xfield}");
    Ok(())
}

/// Writes a spherical cut of the crystal to .xyz and prints its shells.
fn write_qm_cut(input: &Input, blocks: &[WyckoffBlock], cell: &CellParameters) -> Result<(), Box<dyn Error>> {
    let cut = input.cut.as_ref().expect("checked when reading the input");
    let centre = qm_cut::centre_position(&cut.centre, blocks, cell)?;
    let atoms = qm_cut::cut_sphere(blocks, cell, centre, cut.radius);
    let r = cell.to_cartesian(centre);
    println!(
        "Cut: {} atoms within {:.4} Angstrom of fractional ({:.6}, {:.6}, {:.6}) = Cartesian ({:.6}, {:.6}, {:.6}) Angstrom",
        atoms.len(), cut.radius, centre[0], centre[1], centre[2], r[0], r[1], r[2]
    );
    qm_cut::print_shells(&atoms, blocks);
    let charge = qm_cut::formal_charge(&atoms, blocks);
    match charge {
        Some(q) => println!("  Net formal charge = {q:+}"),
        None => println!("  Net formal charge unknown (missing oxidation states)"),
    }
    let comment = format!(
        "Cut of {}: {} atoms within {:.4} Angstrom of ({:.6}, {:.6}, {:.6}) Angstrom{}; columns: element x y z (Angstrom) site distance",
        input.cif, atoms.len(), cut.radius, r[0], r[1], r[2],
        charge.map(|q| format!(", net formal charge {q:+}")).unwrap_or_default()
    );
    qm_cut::write_xyz(&cut.output, &atoms, blocks, cell, &comment).map_err(|e| format!("{}: {e}", cut.output))?;
    println!("Written {}", cut.output);
    Ok(())
}

fn check_charge_neutrality(blocks: &[WyckoffBlock]) -> Result<(), Box<dyn Error>> {
    let mut total_charge = 0.0;
    for block in blocks.iter() {
        let ox = block.site.oxidation.ok_or_else(|| format!(
            "no oxidation state for site {}; give it in the CIF or under [oxidation]", block.site.label
        ))?;
        total_charge += ox * block.positions.len() as f64;
    }
    if total_charge.abs() > NEUTRALITY_TOLERANCE {
        return Err(format!("the unit cell has a net charge of {total_charge}; it must be neutral").into());
    }
    Ok(())
}

fn generate_unit_cell_data(cif_path: &str) -> Result<(Crystal, Vec<WyckoffBlock>), Box<dyn Error>> {
    let crystal = Crystal::read(cif_path)?;
    let cell = &crystal.cell;
    println!("Crystal {} from {cif_path}:", crystal.name);
    println!(
        "  a = {:.6}, b = {:.6}, c = {:.6} Angstrom; alpha = {:.4}, beta = {:.4}, gamma = {:.4} degrees",
        cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma
    );
    println!("  {} symmetry operations", crystal.symmetry_operations.len());
    let blocks = crystal.unfold(DUPLICATE_TOLERANCE);
    for block in &blocks {
        let site = &block.site;
        let wyckoff = site.wyckoff_letter.as_deref().map(|w| format!(" ({w})")).unwrap_or_default();
        println!("  {:<6} {:<3} {:>4} positions{wyckoff}", site.label, site.element, block.positions.len());
        if let Some(m) = site.multiplicity {
            if m != block.positions.len() {
                eprintln!("warning: {} has multiplicity {m} in the CIF", site.label);
            }
        }
        if site.occupancy < 1.0 {
            eprintln!("warning: {} has occupancy {}", site.label, site.occupancy);
        }
    }
    Ok((crystal, blocks))
}

/// Only for checking the read by eye / in VESTA.
#[allow(dead_code)]
fn write_poscar(path: &str, cell: &CellParameters, blocks: &[WyckoffBlock]) -> std::io::Result<()> {
    let mut file = BufWriter::new(File::create(path)?);
    writeln!(file, "comment")?;
    writeln!(file, "1.0")?;
    for v in &cell.vectors {
        writeln!(file, "{:.6} {:.6} {:.6}", v[0], v[1], v[2])?;
    }
    let elements: Vec<&str> = blocks.iter().map(|b| b.site.element.as_str()).collect();
    writeln!(file, "{}", elements.join(" "))?;
    let counts: Vec<String> = blocks.iter().map(|b| b.positions.len().to_string()).collect();
    writeln!(file, "{}", counts.join(" "))?;
    writeln!(file, "Direct")?;
    for block in blocks {
        for p in &block.positions {
            writeln!(file, "{:.6} {:.6} {:.6}", p[0], p[1], p[2])?;
        }
    }
    file.flush()
}
