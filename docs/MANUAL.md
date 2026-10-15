# Geod manual

## Contents

1. [Overview](#1-overview)
2. [Running the program](#2-running-the-program)
3. [Input reference](#3-input-reference)
4. [Modes and workflows](#4-modes-and-workflows)
5. [Output files](#5-output-files)
6. [Diagnostics and checks](#6-diagnostics-and-checks)
7. [Conventions](#7-conventions)
8. [Limitations](#8-limitations)
9. [Citing, credits and licence](#9-citing-credits-and-licence)

## 1. Overview

Geod builds electrostatic embeddings for ionic crystals from a CIF file. It is part of the
Ädelsten project, where it covers the embedding side of SCEPIC; the self-consistent optimization of
the embedding AIMPs belongs to the Ädelsten SCF code. It is meant for
embedded-cluster calculations, in which a QM region is surrounded by embedding AIMPs (ab initio
model potentials) and a field of point charges that reproduces the crystal's Madelung potential.

### Method

A finite cluster of unit cells does not, in general, reproduce the potential of the infinite
crystal: if the cell has a dipole or quadrupole moment, the potential in the cluster depends on its
size and shape. Sushko and Abarenkov (2010) solve this by complementing the unit cell with a set of
point charges, chosen analytically, that cancel all components of all multipole moments of the
cell up to an order M. The lattice sum of the modified cells then converges absolutely.

The complementary (auxiliary) charges sit on the integer points of eight tetrahedra, one at each
corner of the unit cell, each cancelling one eighth of the moments. When the modified cells are
stacked into a cluster, the auxiliary charges of neighbouring cells cancel exactly in the interior.
Only a shell at the boundary of the cluster, about M cells thick, keeps a net auxiliary charge.

The program:

1. reads the CIF, unfolds the asymmetric unit with the symmetry operations, and assigns charges
   (oxidation states);
2. computes the auxiliary charges of the unit cell for the chosen M and tetrahedron shift V;
3. builds the cluster from a list of cells, accumulates the auxiliary field on a grid (so
   overlapping charges merge exactly), moves auxiliary charges that coincide with atoms onto those
   atoms, and drops the ones that cancel;
4. marks the QM region and the AIMP region;
5. checks that neither region reaches the boundary shell, and compares the potential and field of
   the cluster at the QM and AIMP sites with the Ewald potential of the infinite crystal;
6. writes the QM and AIMP regions as basis-set blocks and the remaining charges as a point-charge
   field.

## 2. Running the program

```sh
cargo build --release
target/release/geod input.toml
```

The only argument is a TOML input file (section 3). File paths in the input are relative to the
working directory, and output files are written there.

Run the test suite with `cargo test --release`.

**Name clash with PROJ.** The cartographic projection library PROJ also installs a command called
`geod` (a geodesic calculator), and it comes along with PROJ from Linux distributions, Homebrew,
conda and GDAL environment modules. If both are on your `PATH`, whichever comes first is run, so
check with `which -a geod`. To avoid the clash, call Geod by its full path
(`target/release/geod`), install it under another name (e.g. a symlink `aedelsten-geod`), or make
sure its directory comes before PROJ's in `PATH`.

## 3. Input reference

### Top level

| Key | Type | Default | Description |
|---|---|---|---|
| `mode` | string | required | `"production"`, `"aimp-optimization"` or `"qm-cut"` (section 4). |
| `cif` | string | required | The CIF file. Only the first data block is read. |
| `project` | string | `"project"` | Base name of the output files. |
| `max_order` | integer | `4` | Highest multipole order M cancelled by the auxiliary charges. At most 4. |
| `shift` | number | `0.0` | Tetrahedron shift V, in fractional units (section 6.4). |
| `qm` | string | — | Production mode only, and required there: the QM region as an .xyz file. |

### `[cluster]` (embedding modes, required)

Exactly one of:

| Key | Type | Description |
|---|---|---|
| `radius` | number | Spherical cluster: all cells whose atom centroid lies within this radius (Å) of the centre. The centre is the centroid of the QM region in production mode, and the QM atom in AIMP-optimization mode. |
| `cells` | array | Cuboid cluster around the origin cell: `[n1, n2, n3]` cells in both directions along each lattice vector, or `[[n1-, n1+], [n2-, n2+], [n3-, n3+]]`. |

### `[aimp]` (embedding modes, optional)

Atoms represented by embedding AIMPs. QM atoms are never AIMPs. Without this section there are no
AIMPs. Exactly one of:

| Key | Type | Description |
|---|---|---|
| `radius` | number | All atoms within this radius (Å) of the centroid of the QM region. |
| `cells` | array | All atoms in these cells, counted from the origin cell, same form as in `[cluster]`. |

### `[basis]` (embedding modes, optional)

Labels written in the basis-set blocks. `{element}` and `{site}` (the CIF site label) are
substituted.

| Key | Default |
|---|---|
| `qm` | `"{element}.ANO-S-MB"` |
| `aimp` | `"{element}.EMB-AIMP.author.0s.0s.ECP.{site}"` |

### `[oxidation]` (optional)

Charges overriding those from the CIF, as `key = charge`. A key matches a site by its label
(`Y0`), then by its type symbol (`Y3+`), then by its element (`Y`). A key that matches no site is
an error.

Without an override, the charge of a site comes from the CIF's `_atom_type_oxidation_number` loop,
or failing that from the type symbol (`Y3+` → +3, `O2-` → −2). The embedding modes require a
charge for every site and a neutral unit cell (net charge below 1e-8).

### `[cut]` (qm-cut mode, required)

| Key | Type | Default | Description |
|---|---|---|---|
| `site` | string | — | Centre on the first position of this Wyckoff site (CIF site label) in the origin cell. |
| `point` | array | — | Centre on this fractional position. |
| `cartesian` | array | — | Centre on this Cartesian position, in `unit`. |
| `radius` | number | required | Radius of the cut, in `unit`. |
| `unit` | string | `"angstrom"` | `"angstrom"` or `"bohr"`; applies to `radius` and `cartesian`. |
| `output` | string | `"{project}_cut.xyz"` | Output file. |

Exactly one of `site`, `point` and `cartesian` must be given.

Unknown keys, and sections that don't belong to the chosen mode, are errors. A misspelled key is
reported rather than silently ignored.

## 4. Modes and workflows

### Production

One embedding around a QM region you provide:

```toml
mode = "production"
cif  = "Y2O3.cif"
qm   = "qm.xyz"

[cluster]
radius = 60.0

[aimp]
radius = 15.0
```

The QM file is an .xyz file in Å (section 5.3) in the Cartesian frame of the CIF (section 7).
Each QM atom must coincide, within 0.01 Å, with an atom of the crystal of the same element, in a
cell that is part of the cluster. Each crystal atom may be given only once. The QM atoms are
written in the order of the file; no sorting by site is done. Which atoms make a sensible QM region
is left to the user.

Output: `{project}.inp` and `{project}.xfield`.

### AIMP optimization

One embedding per Wyckoff site, as the starting point for the self-consistent optimization of the
embedding AIMPs; the optimization itself is done by the Ädelsten SCF code. No QM file is given. For each site, the first position of the site in the origin cell becomes the single
QM atom, the cluster is centred on it, and the AIMP shell around it uses the AIMPs of all sites
(the ones being optimized):

```toml
mode = "aimp-optimization"
cif  = "Y2O3.cif"

[cluster]
radius = 60.0

[aimp]
radius = 10.0
```

Output: `{project}_{site}.inp` and `{project}_{site}.xfield` for every site, for example
`project_Y0.inp`, `project_Y1.inp` and `project_O2.inp` for Y₂O₃. All sites go through the same
checks; the run stops at the first site that fails one.

### QM cut

A spherical cut of the pristine crystal, as a starting point for a QM region:

```toml
mode = "qm-cut"
cif  = "Y2O3.cif"

[cut]
site   = "Y0"
radius = 4.0
output = "qm.xyz"
```

The program takes all atoms within the radius from a supercell just large enough to hold the
sphere, sorts them by distance from the centre (atoms within 1e-4 Å of the same distance form one
shell, ordered by Wyckoff site), prints the shells and the net formal charge, and writes the cut as
.xyz:

```
  Distance (A)  Atoms
        0.0000  1 x Y (Y0)
        2.2471  2 x O (O2)
        2.2697  2 x O (O2)
        2.3347  2 x O (O2)
        3.5148  2 x Y (Y1)
        3.5317  4 x Y (Y0)
  Net formal charge = +9
```

The file can be used directly as the `qm` file of a production run, or edited first, for example
to complete coordination spheres or adjust the charge.

## 5. Output files

### 5.1 Basis-set input (`.inp`)

The QM region followed by the AIMP region, as `Basis set` blocks. Positions are in Bohr with
6 decimals; no unit keyword is written, Bohr being the default.

QM atoms appear in the order of the QM input. A new block starts whenever the element changes
from the previous atom, and atoms are numbered per element across blocks (`Y1`, `O1`, `O2`, `Y2`):

```
  Basis set
    Y.ANO-S-MB
    Y1              0.000000        5.011471        5.663835
  End of Basis
```

AIMPs get one block per Wyckoff site, each with that site's own potential, and are numbered `A1`,
`A2`, ... throughout. Within a block they are sorted by distance from the AIMP reference point:

```
  Basis set
    Y.EMB-AIMP.author.0s.0s.ECP.Y0
    pseudocharge
    A1             -4.359106        0.000000        5.011471
    ...
  End of Basis
```

### 5.2 Point-charge field (`.xfield`)

All atoms outside the QM and AIMP regions, followed by the remaining auxiliary charges. The header
gives the number of charges and the multipole order (0, charges only). Each line is `x y z q`,
positions in Bohr with 6 decimals and charges with 12 decimals:

```
62057 0
     -120.275297      -15.034412        5.663835      3.000000000000
     ...
```

Charges carry more digits than positions because the auxiliary charges are large (up to ~20 e at
the boundary) and cancel the cell moments only in sums over many of them.

The QM region, AIMP region and field together are neutral.

### 5.3 .xyz files

The QM input and the qm-cut output: the number of atoms, a comment line, then one atom per line as
`label x y z` in Å. The element is read from the label (`Y1` → Y). Further columns are ignored;
qm-cut writes the site label and the distance from the centre there.

## 6. Diagnostics and checks

### 6.1 Multipole deviations

For the unit cell with its auxiliary charges, and again for the whole cluster, the program prints
the summed moment components of each order up to M + 2:

```
  Deviation of multipoles from zero.
  Cancellation up to order 4 was requested:
    P(1) =     1.91547e-12
    ...
    P(4) =     -2.10231e-8
    P(5) =       1.91723e4
```

Orders up to M should vanish to rounding error; the higher orders are not cancelled and are
printed for reference.

### 6.2 Region check (error)

Every QM and AIMP atom must carry exactly its formal charge. In the interior of the cluster the
auxiliary charges cancel, but an auxiliary charge that coincides with an atom is merged into that
atom's charge, and in the boundary shell that charge does not cancel. A QM or AIMP atom carrying
extra charge therefore means that the region reaches the boundary shell. The QM calculation would
replace the atom by its formal charge and silently lose the rest. The run stops with an error
naming the first such atom; enlarge the cluster or shrink the region.

The boundary shell is about M cells thick, so the cluster must extend roughly M cells beyond the
AIMP region in every direction. For Y₂O₃ (a = 10.6 Å) with M = 4 and a 15 Å AIMP radius, a 50 Å
cluster fails this check and a 60 Å cluster passes.

### 6.3 Embedding accuracy (report)

The potential and field of the cluster at the QM sites, and in summary at the AIMP sites, are
compared with those of the infinite crystal, computed by Ewald summation (atomic units):

```
Embedding accuracy (cluster minus infinite crystal, atomic units):
  QM sites:
    Site           V (Eh/e)        dV (Eh/e)   |dE| (Eh/e/a0)
    Y0        -1.1342556878        -1.398e-4         1.278e-6
    Mean offset        =    -1.398e-4 Eh/e
    RMS spread of dV   =      0.000e0 Eh/e
    Max |dE|           =     1.278e-6 Eh/(e a0)
  AIMP sites (928):
    Mean offset        =    -1.318e-4 Eh/e
    RMS spread of dV   =     1.620e-4 Eh/e
    Max |dE|           =     2.742e-4 Eh/(e a0)
```

`V` is the on-site potential at a QM site from all other charges. A constant offset between cluster
and crystal does not affect the QM calculation; the spread of ΔV across the region and the field
difference do. These numbers decrease as the cluster grows. They are the basis for choosing the
cluster size, rather than a fixed rule.

### 6.4 Warnings

- **Auxiliary charge near an atom.** With a nonzero shift V, the auxiliary charges sit on
  sublattices at ±V (fractional) from the cell corners. An atom that lies on one of them takes over
  its charge, as above. An atom closer than 0.5 Å to one, without lying on it, is reported, because
  the two would give a near-singular potential. Choose another shift.
- **CIF multiplicity or occupancy.** The number of unfolded positions of a site is compared with
  the multiplicity in the CIF, and partial occupancies are reported. Partially occupied sites are
  treated as fully occupied.
- **AIMP cells outside the cluster** (with `[aimp] cells`).

## 7. Conventions

- **Lattice orientation:** a along x, b in the xy-plane, from the CIF cell parameters. All
  Cartesian input (the QM file, `[cut] cartesian`) must be in this frame, with the origin at the
  cell origin. For non-orthogonal cells, an .xyz file from another program may be rotated relative
  to this frame and will then not match.
- **Units:** input lengths in Å (except `[cut]` with `unit = "bohr"`); `.inp` and `.xfield` in
  Bohr; Ewald report in atomic units.
- **Positions:** unit-cell positions are wrapped into [0, 1). Symmetry copies closer than 0.01 Å
  count as one position.
- **Sites:** a Wyckoff site is one row of the CIF's `_atom_site_` loop, identified by its label.
  AIMPs, the AIMP-optimization runs and the `{site}` substitution all use this label.

## 8. Limitations

- Bulk crystals only. Surfaces are not supported: the boundary shell of the auxiliary charges
  would lie at the surface itself.
- Partially occupied and disordered sites are treated as fully occupied.
- The multipole order is limited to M ≤ 4.
- The Cartesian frame is fixed by the CIF (section 7).
- The point-charge field is not symmetrized around the QM region.

## 9. Citing, credits and licence

If you use results from this program, please cite:

- E. D. Larsson, M. Krośnicki, V. Veryazov, "A program system for self-consistent embedded
  potentials for ionic crystals", Chem. Phys. **562**, 111549 (2022),
  doi:10.1016/j.chemphys.2022.111549
- P. V. Sushko, I. V. Abarenkov, "General Purpose Electrostatic Embedding Potential",
  J. Chem. Theory Comput. **6**, 1323 (2010), doi:10.1021/ct900480p

BibTeX entries are in `CITATION.bib`.

The method and its original Fortran implementation (EMBQ) are by Peter V. Sushko, with Igor V.
Abarenkov. The Julia implementation in SCEPIC is by Ernst D. Larsson, building on EMBQ with the
permission of P. V. Sushko. This Rust implementation is
by Ernst D. Larsson, written together with Claude (Anthropic).

Licensed under the MIT licence; see `LICENSE`.
