# Geod

Electrostatic embedding of ionic crystals, starting from a CIF file: the crystalline shell the QM
cluster sits in, as a geode's crystals line its cavity. Part of the Ädelsten
project: this program takes over the embedding side of SCEPIC (cluster construction, QM and AIMP
regions, point-charge fields), while the self-consistent optimization of the embedding AIMPs
belongs to the Ädelsten SCF code.

From the crystal structure, the program builds a finite cluster whose electrostatic potential in
its centre reproduces that of the infinite crystal:

- the unit cell is complemented with point charges that cancel all of its multipole moments up to
  a chosen order (Sushko and Abarenkov, 2010), so the lattice sum converges absolutely;
- copies of this modified cell form a large cluster, in which the complementary charges cancel
  everywhere except in a thin shell at the boundary;
- a QM region and a shell of embedding AIMPs are marked around the centre, checked against the
  boundary shell, and the result is compared with the Ewald potential of the infinite crystal;
- the QM and AIMP regions are written as basis-set blocks (`.inp`) and the rest as a point-charge
  field (`.xfield`), in OpenMolcas conventions.

There are three modes: a production embedding around a given QM region, one embedding per
Wyckoff site for optimizing the embedding AIMPs, and spherical cuts of the crystal as a starting
point for QM regions.

## Building and running

```sh
cargo build --release
target/release/geod input.toml
cargo test --release      # the test suite
```

A minimal production input:

```toml
mode = "production"
cif  = "Y2O3.cif"
qm   = "qm.xyz"

[cluster]
radius = 60.0

[aimp]
radius = 15.0
```

All modes, input keys, output formats and diagnostics are described in
[docs/MANUAL.md](docs/MANUAL.md).

## Citing

If you use results from this program, please cite:

- E. D. Larsson, M. Krośnicki, V. Veryazov, "A program system for self-consistent embedded
  potentials for ionic crystals", Chem. Phys. **562**, 111549 (2022),
  doi:[10.1016/j.chemphys.2022.111549](https://doi.org/10.1016/j.chemphys.2022.111549)
- P. V. Sushko, I. V. Abarenkov, "General Purpose Electrostatic Embedding Potential",
  J. Chem. Theory Comput. **6**, 1323 (2010),
  doi:[10.1021/ct900480p](https://doi.org/10.1021/ct900480p)

The program prints these at the end of every run. BibTeX entries are in
[CITATION.bib](CITATION.bib).

## Credits

- **Method and original implementation:** Peter V. Sushko (with Igor V. Abarenkov), in the Fortran
  code EMBQ.
- **Julia implementation:** Ernst D. Larsson, as part of SCEPIC, building on EMBQ with the
  permission of P. V. Sushko.
- **This Rust implementation:** Ernst D. Larsson, written together with Claude (Anthropic).

## Licence

MIT, see [LICENSE](LICENSE). The test structure `tests/data/Y2O3.cif` is from the Materials
Project, under CC BY 4.0. Every file's licence is declared in REUSE format (SPDX headers and
[REUSE.toml](REUSE.toml)); check with `reuse lint`.
