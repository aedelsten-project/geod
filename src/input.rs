// SPDX-FileCopyrightText: 2026 Ernst D. Larsson
// SPDX-License-Identifier: MIT

//! Program input, read from a TOML file.
//!
//! ```toml
//! mode      = "production"        # or "aimp-optimization"
//! cif       = "Y2O3.cif"
//! project   = "project"           # base name of the output files
//! max_order = 4                   # highest multipole order cancelled (M)
//! shift     = 0.0                 # tetrahedron shift V
//! qm        = "qm.xyz"            # production only: QM region, Cartesian Angstrom
//!
//! [cluster]
//! radius = 60.0                   # Angstrom; or: cells = [4, 4, 4] / [[4, 4], [4, 4], [4, 4]]
//!
//! [aimp]                          # optional
//! radius = 15.0                   # Angstrom; or: cells = [1, 1, 1]
//!
//! [basis]                         # optional; {element} and {site} are substituted
//! qm   = "{element}.ANO-S-MB"
//! aimp = "{element}.EMB-AIMP.author.0s.0s.ECP.{site}"
//!
//! [oxidation]                     # optional overrides, by site label, type symbol or element
//! Y = 3
//! O = -2
//! ```
//!
//! The `qm-cut` mode only reads the CIF and writes a spherical cut of the crystal as an .xyz
//! file (Angstrom), e.g. as a starting point for the `qm` file of a production run:
//!
//! ```toml
//! mode = "qm-cut"
//! cif  = "Y2O3.cif"
//!
//! [cut]
//! site   = "Y0"                   # or: point = [0.5, 0.5, 0.5] (fractional)
//!                                 # or: cartesian = [x, y, z] (in `unit`)
//! radius = 6.0
//! unit   = "angstrom"             # optional, or "bohr"; applies to radius and cartesian
//! output = "qm.xyz"               # optional, default "{project}_cut.xyz"
//! ```
//!
//! Paths are relative to the working directory, and output files are written there.

use std::collections::HashMap;
use std::error::Error;
use std::fs;
use std::path::Path;

use serde::Deserialize;

use crate::crystal::{WyckoffBlock, ANGSTROM_TO_BOHR};

#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "kebab-case")]
pub enum Mode {
    /// One embedding around the QM region given in `qm`.
    Production,
    /// One embedding per Wyckoff site, each with a single atom of that site as QM region, for
    /// optimizing the embedding AIMPs.
    AimpOptimization,
    /// A spherical cut of the crystal, written as .xyz.
    QmCut,
}

#[derive(Debug, Deserialize, Clone, Copy, PartialEq, Eq, Default)]
#[serde(rename_all = "lowercase")]
pub enum LengthUnit {
    #[default]
    Angstrom,
    Bohr,
}

/// The `[cut]` section as written.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct CutRaw {
    site: Option<String>,
    point: Option<[f64; 3]>,
    cartesian: Option<[f64; 3]>,
    radius: f64,
    #[serde(default)]
    unit: LengthUnit,
    output: Option<String>,
}

/// Centre of a cut.
#[derive(Debug, Clone, PartialEq)]
pub enum CutCentre {
    /// The first position of the Wyckoff site with this label, in the origin cell.
    Site(String),
    Fractional([f64; 3]),
    /// Cartesian, Angstrom.
    Cartesian([f64; 3]),
}

/// Validated `[cut]` section; lengths in Angstrom.
#[derive(Debug, Clone)]
pub struct CutInput {
    pub centre: CutCentre,
    pub radius: f64,
    pub output: String,
}

/// Cells to repeat along each lattice vector: the same number in both directions, or
/// `[negative, positive]` per axis.
#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(untagged)]
pub enum CellRange {
    Symmetric([u32; 3]),
    PerDirection([[u32; 2]; 3]),
}

impl CellRange {
    pub fn per_direction(self) -> [[u32; 2]; 3] {
        match self {
            CellRange::Symmetric(n) => n.map(|n| [n, n]),
            CellRange::PerDirection(n) => n,
        }
    }
}

/// A region given either by a radius (Angstrom) or by a range of cells.
#[derive(Debug, Deserialize, Clone, Copy)]
#[serde(deny_unknown_fields)]
pub struct RegionInput {
    pub radius: Option<f64>,
    pub cells: Option<CellRange>,
}

/// A region resolved from [`RegionInput`].
pub enum RegionSpec {
    Radius(f64),
    Cells([[u32; 2]; 3]),
}

impl RegionInput {
    fn resolve(self, name: &str) -> Result<RegionSpec, Box<dyn Error>> {
        match (self.radius, self.cells) {
            (Some(r), None) if r > 0.0 => Ok(RegionSpec::Radius(r)),
            (Some(r), None) => Err(format!("[{name}] radius must be positive, got {r}").into()),
            (None, Some(c)) => Ok(RegionSpec::Cells(c.per_direction())),
            _ => Err(format!("[{name}] needs exactly one of 'radius' or 'cells'").into()),
        }
    }
}

/// Basis-set labels; `{element}` and `{site}` (the CIF site label) are substituted.
#[derive(Debug, Deserialize, Clone)]
#[serde(deny_unknown_fields, default)]
pub struct BasisInput {
    pub qm: String,
    pub aimp: String,
}

impl Default for BasisInput {
    fn default() -> Self {
        Self {
            qm: "{element}.ANO-S-MB".to_string(),
            aimp: "{element}.EMB-AIMP.author.0s.0s.ECP.{site}".to_string(),
        }
    }
}

impl BasisInput {
    pub fn qm_label(&self, element: &str, site: &str) -> String {
        substitute(&self.qm, element, site)
    }

    pub fn aimp_label(&self, element: &str, site: &str) -> String {
        substitute(&self.aimp, element, site)
    }
}

fn substitute(template: &str, element: &str, site: &str) -> String {
    template.replace("{element}", element).replace("{site}", site)
}

/// The input file as written.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct RawInput {
    mode: Mode,
    cif: String,
    #[serde(default = "default_project")]
    project: String,
    #[serde(default = "default_max_order")]
    max_order: usize,
    #[serde(default)]
    shift: f64,
    qm: Option<String>,
    cluster: Option<RegionInput>,
    aimp: Option<RegionInput>,
    cut: Option<CutRaw>,
    #[serde(default)]
    basis: BasisInput,
    #[serde(default)]
    oxidation: HashMap<String, f64>,
}

fn default_project() -> String {
    "project".to_string()
}

fn default_max_order() -> usize {
    4
}

fn resolve_cut(raw: CutRaw, project: &str) -> Result<CutInput, Box<dyn Error>> {
    let scale = match raw.unit {
        LengthUnit::Angstrom => 1.0,
        LengthUnit::Bohr => 1.0 / ANGSTROM_TO_BOHR,
    };
    let centre = match (raw.site, raw.point, raw.cartesian) {
        (Some(site), None, None) => CutCentre::Site(site),
        (None, Some(p), None) => CutCentre::Fractional(p),
        (None, None, Some(r)) => CutCentre::Cartesian(r.map(|x| x * scale)),
        _ => return Err("[cut] needs exactly one of 'site', 'point' or 'cartesian'".into()),
    };
    if !(raw.radius > 0.0) {
        return Err(format!("[cut] radius must be positive, got {}", raw.radius).into());
    }
    Ok(CutInput {
        centre,
        radius: raw.radius * scale,
        output: raw.output.unwrap_or_else(|| format!("{project}_cut.xyz")),
    })
}

/// Validated program input.
pub struct Input {
    pub mode: Mode,
    pub cif: String,
    pub project: String,
    pub max_order: usize,
    pub shift: f64,
    /// QM region file; set in production mode only.
    pub qm: Option<String>,
    /// Set in the embedding modes.
    pub cluster: Option<RegionSpec>,
    pub aimp: Option<RegionSpec>,
    /// Set in qm-cut mode only.
    pub cut: Option<CutInput>,
    pub basis: BasisInput,
    pub oxidation: HashMap<String, f64>,
}

impl Input {
    pub fn read(path: impl AsRef<Path>) -> Result<Self, Box<dyn Error>> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::parse(&text).map_err(|e| format!("{}: {e}", path.display()).into())
    }

    pub fn parse(text: &str) -> Result<Self, Box<dyn Error>> {
        let raw: RawInput = toml::from_str(text)?;
        match (raw.mode, &raw.qm) {
            (Mode::Production, None) => return Err("production mode needs a QM region ('qm')".into()),
            (Mode::AimpOptimization, Some(_)) => {
                return Err("aimp-optimization mode places the QM atom itself; remove 'qm'".into());
            }
            _ => {}
        }
        let embedding = raw.mode != Mode::QmCut;
        if embedding && raw.cluster.is_none() {
            return Err("the embedding modes need a [cluster] section".into());
        }
        if embedding && raw.cut.is_some() {
            return Err("[cut] is only used in qm-cut mode".into());
        }
        if !embedding && (raw.qm.is_some() || raw.cluster.is_some() || raw.aimp.is_some()) {
            return Err("qm-cut mode takes only a [cut] section; remove 'qm', [cluster] and [aimp]".into());
        }
        let cut = match raw.cut {
            Some(cut) => Some(resolve_cut(cut, &raw.project)?),
            None if !embedding => return Err("qm-cut mode needs a [cut] section".into()),
            None => None,
        };
        Ok(Self {
            mode: raw.mode,
            cif: raw.cif,
            project: raw.project,
            max_order: raw.max_order,
            shift: raw.shift,
            qm: raw.qm,
            cluster: raw.cluster.map(|c| c.resolve("cluster")).transpose()?,
            aimp: raw.aimp.map(|a| a.resolve("aimp")).transpose()?,
            cut,
            basis: raw.basis,
            oxidation: raw.oxidation,
        })
    }

    /// Applies the oxidation-state overrides: a key matches a site by label first, then by type
    /// symbol, then by element. Fails if a key matches no site.
    pub fn apply_oxidation(&self, blocks: &mut [WyckoffBlock]) -> Result<(), Box<dyn Error>> {
        for key in self.oxidation.keys() {
            let used = blocks.iter().any(|b| *key == b.site.label || *key == b.site.type_symbol || *key == b.site.element);
            if !used {
                return Err(format!("[oxidation] '{key}' matches no site of the crystal").into());
            }
        }
        for block in blocks.iter_mut() {
            let site = &block.site;
            let q = [&site.label, &site.type_symbol, &site.element]
                .iter()
                .find_map(|k| self.oxidation.get(*k).copied());
            if let Some(q) = q {
                block.site.oxidation = Some(q);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn production_input() {
        let input = Input::parse(r#"
            mode = "production"
            cif = "Y2O3.cif"
            qm = "qm.xyz"
            [cluster]
            radius = 60.0
            [aimp]
            cells = [[0, 1], [1, 1], [2, 0]]
            [oxidation]
            Y = 3
        "#).unwrap();
        assert_eq!(input.mode, Mode::Production);
        assert_eq!(input.project, "project");
        assert_eq!(input.max_order, 4);
        assert!(matches!(input.cluster, Some(RegionSpec::Radius(r)) if r == 60.0));
        assert!(matches!(input.aimp, Some(RegionSpec::Cells([[0, 1], [1, 1], [2, 0]]))));
        assert_eq!(input.basis.aimp_label("Y", "Y0"), "Y.EMB-AIMP.author.0s.0s.ECP.Y0");
        assert_eq!(input.oxidation["Y"], 3.0);
    }

    #[test]
    fn aimp_optimization_input() {
        let input = Input::parse(r#"
            mode = "aimp-optimization"
            cif = "Y2O3.cif"
            max_order = 3
            [cluster]
            cells = [4, 4, 4]
            [basis]
            qm = "{element}.ANO-RCC-VTZP"
        "#).unwrap();
        assert_eq!(input.mode, Mode::AimpOptimization);
        assert!(matches!(input.cluster, Some(RegionSpec::Cells([[4, 4], [4, 4], [4, 4]]))));
        assert!(input.aimp.is_none());
        assert_eq!(input.basis.qm_label("O", "O2"), "O.ANO-RCC-VTZP");
        assert_eq!(input.basis.aimp_label("O", "O2"), "O.EMB-AIMP.author.0s.0s.ECP.O2");
    }

    #[test]
    fn qm_cut_input() {
        let input = Input::parse(r#"
            mode = "qm-cut"
            cif = "Y2O3.cif"
            project = "y2o3"
            [cut]
            cartesian = [1.0, 0.0, 2.0]
            radius = 10.0
            unit = "bohr"
        "#).unwrap();
        let cut = input.cut.unwrap();
        let a = 1.0 / ANGSTROM_TO_BOHR;
        assert_eq!(cut.centre, CutCentre::Cartesian([a, 0.0, 2.0 * a]));
        assert!((cut.radius - 10.0 * a).abs() < 1e-12);
        assert_eq!(cut.output, "y2o3_cut.xyz");
        assert!(input.cluster.is_none());

        let base = "mode = \"qm-cut\"\ncif = \"x.cif\"\n";
        assert!(Input::parse(&format!("{base}[cut]\nsite = \"Y0\"\nradius = 5.0\n")).is_ok());
        // No centre, two centres, a [cluster] section, no [cut] at all.
        assert!(Input::parse(&format!("{base}[cut]\nradius = 5.0\n")).is_err());
        assert!(Input::parse(&format!("{base}[cut]\nsite = \"Y0\"\npoint = [0, 0, 0]\nradius = 5.0\n")).is_err());
        assert!(Input::parse(&format!("{base}[cut]\nsite = \"Y0\"\nradius = 5.0\n[cluster]\nradius = 50.0\n")).is_err());
        assert!(Input::parse(base).is_err());
    }

    #[test]
    fn invalid_inputs() {
        let base = "cif = \"x.cif\"\n[cluster]\nradius = 50.0\n";
        // Production without QM, optimization with QM.
        assert!(Input::parse(&format!("mode = \"production\"\n{base}")).is_err());
        assert!(Input::parse(&format!("mode = \"aimp-optimization\"\nqm = \"q.xyz\"\n{base}")).is_err());
        // Typos, both or neither of radius/cells, unknown mode.
        assert!(Input::parse(&format!("mode = \"production\"\nqm = \"q.xyz\"\nmax_ordr = 3\n{base}")).is_err());
        assert!(Input::parse("mode = \"aimp-optimization\"\ncif = \"x\"\n[cluster]\nradius = 5.0\ncells = [1, 1, 1]\n").is_err());
        assert!(Input::parse("mode = \"aimp-optimization\"\ncif = \"x\"\n[cluster]\n").is_err());
        assert!(Input::parse(&format!("mode = \"testing\"\n{base}")).is_err());
    }
}
