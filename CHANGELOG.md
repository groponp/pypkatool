# Changelog

All notable changes to `pypkatool` are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- `--pdbid` flag on `fixstructure`: fetches the deposited `SEQRES` from
  RCSB as a reference sequence for internal-gap detection, for PDBs that
  have no `SEQRES` records of their own (confirmed: without this,
  PDBFixer's `findMissingResidues()` has nothing to compare the chain
  against and silently repairs 0 internal gaps, even real ones).

### Changed
- **Breaking**: the `--ph` flag on `run` and `reprocess` is renamed to
  `--pH` (matching standard pKa/pH notation). Scripts calling `pypkatool`
  need updating; the old lowercase `--ph` is no longer accepted.

### Fixed
- `fixstructure` now prints an explicit `WARNING` when it cannot detect
  internal gaps at all (no `SEQRES` and no `--pdbid`), instead of silently
  reporting `internal_residues_added: 0` with no indication that this
  might mean "couldn't check" rather than "nothing missing".

## [1.0.0] - 2026-07-09

### Added
- `pypkatool run` - full PyPKA (Poisson-Boltzmann + Monte Carlo) pipeline
  with pKAI+ cross-validation and CHARMM36 RESI/PRES protonation-state
  mapping.
- `pypkatool reprocess` - regenerate reports from a previous run's raw
  PyPKA output without rerunning the PB+MC step.
- `pypkatool fixstructure` - repair a structurally incomplete PDB with
  PDBFixer before running it: fills missing heavy atoms everywhere
  (including chain termini), rebuilds missing residues only for internal
  chain gaps (never extends a chain), adds no hydrogens. Runs in a separate
  `pdbfixer` conda environment (PDBFixer/OpenMM need numpy>=2, incompatible
  with delphi4py's numpy<2 pin).
- `--charmm-input` flag on `run` - sets PyPKA's `ffinput=CHARMM` for PDBs
  that already carry CHARMM protonation-state residue names (e.g. a
  CHARMM-GUI PDB Reader round-trip). Confirmed empirically to have no
  effect on standard PDB input across all 12 benchmark proteins.
- Actionable error message when PyPKA's preprocessing step fails on a
  structurally incomplete PDB, instead of a raw Python traceback.
- `CITATION.cff` for GitHub's "Cite this repository" button.
- Three-layer test suite (`tests/test_pypkatool.py`): unit tests on pure
  mapping/validation functions, regression tests against frozen fixtures,
  and adversarial tests for the RTF/PDB cross-validation machinery.

### Fixed
- `DEFAULT_PARAMS` never set `ffID`, so PyPKA silently ran its own default
  (`G54A7`/GROMOS) instead of `CHARMM36m` - the force field the bundled RTF
  and every CHARMM RESI/PRES label this tool produces assume. `ffID` is now
  pinned to `"CHARMM36m"` explicitly. This also disables PyPKA's `SER`/`THR`
  titration (correct under CHARMM36m, which has no `PRES` for either
  deprotonated state) and shifts pKa values across the board;
  `tests/fixtures/` was regenerated to match.

### Documentation
- CHARMM36/CHARMM36m and pKAI+ citation metadata corrected and verified
  against source DOIs.
- README clarified: two/three independent conda environments (main,
  `py27`, `pdbfixer`), each with its own purpose, none nested in another.
- Naming unified to `pypkatool` throughout the codebase and docs.
