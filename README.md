# pypkatool

`pypkatool` is a command-line pipeline that takes a protein structure and a
target pH, runs [PyPKA](https://doi.org/10.1021/acs.jcim.0c00718)
(Poisson-Boltzmann + Monte Carlo pKa calculation), cross-validates every
predicted pKa against the independent machine-learning predictor
[pKAI+](https://doi.org/10.1021/acs.jctc.2c00308), and writes a
CHARMM-GUI-ready protonation-state table (which `RESI`/`PRES` to pick for
each titratable residue in PDB Reader).

Full documentation (installation, usage, force field parameters, tautomer
mapping, `fixstructure`, API reference) lives in [`docs/`](docs/), built with
Sphinx - see ["Building the docs"](#building-the-docs) below. This README is
a quickstart summary of the same content.

Output per run:

| File | Contents |
|---|---|
| `protonation_inputs.dat` | Residues that need a CHARMM-GUI action (RESI choice or PRES patch), plain text |
| `protonation_inputs.json` | Every titratable site, full detail, machine-readable |
| `detail.json` | Per-site tautomer population breakdown + the exact RTF blocks used |
| `crossvalidation_report.dat` | PyPKA vs. pKAI+ pKa and agreement, per site |
| `protocol.json` | Exact PyPKA parameters and package versions used (reproducibility) |
| `output_pypka/` | Raw PyPKA output: protonated PDB, full pH 0-14 titration curve |

## Installation

Requires [conda](https://docs.conda.io/) (miniconda or anaconda) on `PATH`.

This project uses **separate, independent conda environments** - none is
created "inside" another, and they serve different purposes:

| Environment | Created from | Do you activate it? | Required? |
|---|---|---|---|
| `pypkatool` | `environment.yml` | **Yes - every time** you run the tool: `conda activate pypkatool` | Always |
| `py27` | `environment-py27.yml` | **No, never.** `pypkatool` finds it on disk automatically at runtime. | Always |
| `pdbfixer` | `environment-pdbfixer.yml` | **No, never.** `pypkatool` finds it on disk automatically at runtime. | Only for the `fixstructure` command |

Why `py27` exists: PyPKA's compiled Poisson-Boltzmann backend (DelPhi4py)
needs `numpy<2` and `libgfortran4=7.5.0` in the main environment (the
compiled binary was built against NumPy's 1.x C API and links `GFORTRAN_7`,
which `libgfortran.so.5` / GFortran 10+ doesn't export) - that part lives in
`pypkatool`. Separately, PyPKA also shells out to a bare `python2.7`
interpreter internally (via `pdbmender`'s vendored `pdb2pqr.py`), which
cannot coexist with Python 3 in the same environment, hence the second,
minimal, interpreter-only environment.

Why `pdbfixer` exists: the optional `fixstructure` command (see
["Repairing fragmented structures"](#repairing-fragmented-structures)) uses
PDBFixer/OpenMM, which require `numpy>=2` - the opposite pin from `pypkatool`'s
own `numpy<2` requirement. The two cannot share an environment, so this is a
third, separate one. It is only needed if you plan to use `fixstructure`.

### One-time setup

```bash
# 1. Create the required environments (pdbfixer is optional - see below).
#    Order doesn't matter, and it doesn't matter whether any environment is
#    currently active: `conda env create` always builds a new, independent
#    environment from scratch.
conda env create -f environment.yml
conda env create -f environment-py27.yml
conda env create -f environment-pdbfixer.yml   # optional, only for fixstructure

# 2. Activate the MAIN environment (py27 is never activated) and install
#    pypkatool into it - pick ONE of the two:
conda activate pypkatool

# 2a. Editable install (recommended if you cloned this repo to modify or
#     update it): the command reads pypkatool/ from this checkout directly,
#     so `git pull` and local edits take effect immediately, no reinstall.
pip install -e .

# 2b. Regular install (recommended if you just want to use the CLI and
#     don't plan to touch the source): copies the package into the
#     environment's site-packages, same as any other pip package.
pip install .

# 3. Verify
pypkatool --help
```

`pypkatool` locates the Python 2.7 interpreter automatically by looking for a
conda environment literally named `py27` under your home directory
(`~/miniconda3/envs/py27/bin/python2.7` or `~/anaconda3/envs/py27/bin/python2.7`).
If you named it something else, either rename it to `py27` or prepend its
`bin/` to `PATH` yourself before running `pypkatool run`.

The `pdbfixer` environment (if created) is located the same way, by name,
under `~/miniconda3/envs/pdbfixer/bin/python` or
`~/anaconda3/envs/pdbfixer/bin/python`.

### Every time you want to use it

Only `pypkatool` needs activating - `py27` is "install once and forget":

```bash
conda activate pypkatool
pypkatool run my_protein.pdb --ph 7.0
```

## Usage

```bash
# Run the full pipeline on one PDB at pH 7.0
pypkatool run my_protein.pdb --ph 7.0

# Choose an output directory, CPU count, and protein interior dielectric
pypkatool run my_protein.pdb --ph 5.0 --outdir results/ --ncpus 8 --epsin 15

# Regenerate reports (e.g. at a different pH) from a previous run's raw
# PyPKA output, without rerunning the expensive PB+MC step
pypkatool reprocess results/ --ph 6.0 --pdb my_protein.pdb

# Repair a structurally incomplete PDB first, then feed the repaired
# structure into the normal pipeline (see "Repairing fragmented structures")
pypkatool fixstructure my_protein.pdb --outdir results/
pypkatool run results/my_protein_fixed.pdb --ph 7.0
```

Try it on the bundled example:

```bash
pypkatool run examples/denv2.pdb --ph 5.0
```

### `run` options

| Flag | Default | Meaning |
|---|---|---|
| `--ph` | *(required)* | Target pH for the most-probable-protonation-state report |
| `--outdir` | `pypkatool_<stem>_pH<ph>/` next to the input PDB | Output directory |
| `--ncpus` | all detected CPUs | Parallel PB solves in PyPKA |
| `--epsin` | `15` | Protein interior dielectric constant (PyPKA's literature-optimized default; RMSE 0.82 / MAE 0.57 on the PKAD benchmark) |
| `--charmm-input` | off | Set PyPKA's `ffinput=CHARMM` (see ["Force field parameters"](#force-field-parameters-used-for-the-pbmc-calculation)) |

### `fixstructure` options

| Flag | Default | Meaning |
|---|---|---|
| `pdb` | *(required)* | Input structure to repair |
| `--outdir` | input PDB's parent directory | Where to write `<stem>_fixed.pdb` |

## Force field parameters used for the PB+MC calculation

`pypkatool` always runs PyPKA with `ffID="CHARMM36m"`. This selects which
partial-charge/radius database (`DataBaseT.crg`/`.siz`) and which set of
tautomer definition files (`CHARMM36m/sts/*.st`, with their intrinsic
`pKmod` values) PyPKA uses for the electrostatics calculation itself. It is
not optional and cannot be overridden from the CLI: every label this tool
produces is a CHARMM36 RESI/PRES name, backed by the bundled
`top_all36_prot.rtf`, so the underlying pKa calculation has to be run with
CHARMM36m parameters for those labels to be self-consistent with the physics
that produced them. Setting `ffID` also determines whether `SER`/`THR` are
titrated at all: PyPKA disables `SER`/`THR` titration automatically whenever
`ffID` contains `"charmm36m"` (there is no `PRES` for a deprotonated `THR` in
standard CHARMM36, and `SER`/`THR` titration is not part of the CHARMM36m
parameterization PyPKA validates against).

A second, independent parameter, `ffinput`, controls how PyPKA interprets
residue names that are already in the *input* PDB before any calculation
happens - it has no effect on the electrostatics. `ffinput="CHARMM"`
(equivalently, `pypkatool run --charmm-input`) is needed only when the input
PDB already carries CHARMM protonation-state-specific residue names (`HSD`,
`HSE`, `HSP`, `ASPP`, `GLUP`, `CYSM`, ...) - for example, a structure
re-exported from a CHARMM-GUI PDB Reader step that has already had
protonation states assigned. Without it, PyPKA would not recognize those
residues as titratable sites at all. A standard PDB (RCSB, AlphaFold, or any
other source that spells residues canonically - `HIS`, `ASP`, `GLU`, ...)
is completely unaffected by this flag; it exists only for that CHARMM-GUI
round-trip case.

## How the tautomer -> CHARMM label mapping works

PyPKA reports, per titratable site, the most probable *tautomer* out of `N`
regular tautomers plus one *reference* tautomer (index `N+1`). Which physical
state that reference tautomer represents depends on site polarity - this is
easy to get backwards, so it's worth stating precisely (verified against the
installed PyPKA source, `titsite.py::Titsite.getRefProtState()`, and against
the signed atomic partial charges in the CHARMM36 `.st` tautomer files):

* **Cationic sites** (`HIS`, `LYS`, `NTR`): the reference tautomer is the
  **protonated**, positively charged state. It is also the CHARMM *default*
  RESI (`LYS`, `NTER`) - no patch needed. The regular tautomers are the
  neutral forms and need a patch (`LSN`, `NNEU`).
* **Anionic sites** (`ASP`, `GLU`, `CTR`, `CYS`, `TYR`, `SER`): the reference
  tautomer is the **deprotonated**, negatively charged state.
  * For `ASP`/`GLU`/`CTR`, that deprotonated/charged form *is* the CHARMM
    default RESI (`ASP`, `GLU`, `CTER`) - no patch needed; the protonated
    form needs a patch (`ASPP`, `GLUP`, `CNEU`).
  * For `CYS`/`TYR`/`SER` it's the opposite: the CHARMM default RESI is the
    *protonated*, neutral form (`CYS`, `TYR`, `SER`), so the deprotonated
    reference state is the one needing a patch (`CYSD`, `TYRD`, `SERD`).

`HIS` is handled separately since it toggles between two distinguishable
neutral tautomers rather than a single reference state: tautomer 1 ->
`HSD` (proton on ND1), tautomer 2 -> `HSE` (proton on NE2), reference ->
`HSP` (both nitrogens protonated).

Every label produced by the pipeline is additionally cross-checked against
the actual hydrogen atoms present in PyPKA's output protonated PDB (ground
truth) and against the bundled CHARMM36 RTF (`pypkatool/data/top_all36_prot.rtf`),
which is the authority for which labels exist and what atoms each implies. A
handful of deprotonated states have no patch in standard CHARMM36
(`THR` deprotonated has no `PRES THRD`; `TYR` deprotonated has no `PRES
TYRD`) - these are reported with a warning rather than silently mismapped.

## Repairing fragmented structures

PyPKA does not repair a structurally incomplete PDB - it requires every
present residue to have its expected atom set. A truncated terminus missing
its carboxyl oxygens, or a residue missing backbone/side-chain atoms, will
make the underlying preprocessing step fail (see ["Failure
modes"](#failure-modes)). For inputs like this, run `pypkatool fixstructure`
first to produce a repaired PDB, then feed that into `pypkatool run`.

`fixstructure` wraps [PDBFixer](https://github.com/openmm/pdbfixer) with a
fixed policy, chosen to repair damage without inventing biology that was
never part of the deposited structure:

1. **Missing heavy atoms are always filled in**, for any residue that is
   present in the structure - including residues at the very start or end of
   a chain (e.g. a C-terminal residue missing its `OXT`).
2. **Missing residues (a gap in the chain) are rebuilt only when the gap is
   internal** - i.e. there is at least one present residue on both sides of
   it. A gap at the very start or end of a chain is left untouched:
   `fixstructure` repairs holes inside a chain, it does not extend a chain
   to match a SEQRES record or a reference sequence.
3. **No hydrogens are added.** The output stays heavy-atom-only, the same
   convention as a standard deposited PDB (and the form PyPKA's own cleaning
   step expects as input).

| Input residue state | Chain position | Result |
|---|---|---|
| Present, missing some atoms (e.g. `OXT`, side-chain atoms) | Anywhere, including chain ends | Atoms filled in |
| Completely absent (gap in numbering) | Internal (residues present on both sides) | Residue rebuilt |
| Completely absent (gap in numbering) | At the start or end of a chain | Left untouched - chain is not extended |

```bash
pypkatool fixstructure my_protein.pdb --outdir results/
# -> results/my_protein_fixed.pdb
pypkatool run results/my_protein_fixed.pdb --ph 7.0
```

AlphaFold2/ColabFold predictions are heavy-atom-only (no hydrogens) and, for
the residue range you gave the model, have no internal gaps or missing
atoms by construction - `fixstructure` is normally unnecessary for them.
It is intended for experimentally determined structures with genuine
crystallographic disorder or truncated regions.

## Failure modes

If PyPKA's preprocessing step fails - most commonly because of a
structurally incomplete PDB - `pypkatool run` reports the underlying error
together with a suggestion to repair the structure first (see
["Repairing fragmented structures"](#repairing-fragmented-structures)),
rather than a raw Python traceback.

A gap that is silently *ignored* rather than causing a failure is also
possible: if an entire residue is absent from the input and `fixstructure`
was not run first, PyPKA simply omits that position - it will not appear
anywhere in the output tables, and no warning is issued. Compare the number
of sites reported against the residue count you expect for the protein if
this is a concern.

## Validation

Cross-validated against 12 benchmark proteins spanning pH 5 (flavivirus
envelope protein, endosomal pH sensing) and pH 7 (standard pKa benchmarks:
lysozyme, RNase A, barnase, BPTI, SNase, protein G, thioredoxin), matching
published experimental pKa values where available (e.g. lysozyme HIS15
pKa 6.33, experimental ~5.7) and reaching 96.7-100% sign agreement with the
independent pKAI+ predictor across ~1145 titratable sites (`SER`/`THR` are
not titrated under `ffID="CHARMM36m"` - see ["Force field
parameters"](#force-field-parameters-used-for-the-pbmc-calculation)). See
`tests/fixtures/` for the frozen per-protein output used by the regression
test suite, and `tests/data/` for the corresponding input PDBs.

Known limitation: for large multimeric/icosahedral assemblies (e.g. the
flavivirus envelope protein virion), pKa predictions computed on a
crystallographic dimer can miss inter-subunit electrostatic contacts present
in the full assembly, shifting a small number of buried/interfacial residues.
This is a physical limitation of static Poisson-Boltzmann on a partial
structure, not a mapping bug; the pipeline's cross-validation report is the
mechanism for catching such cases (see `crossvalidation_report.dat`).

## Tests

```bash
conda activate pypkatool
python tests/test_pypkatool.py -v
```

Three layers: unit tests on pure mapping/validation functions, regression
tests comparing against the frozen outputs in `tests/fixtures/`, and
adversarial tests that inject inconsistent inputs to confirm the RTF/PDB
cross-validation machinery (`reconcile_his_from_pdb`,
`reconcile_non_his_from_pdb`, `_validate_label_chain`) actually catches them.
No PyPKA rerun is needed to run this suite.

## Building the docs

The full documentation site (`docs/`) is built with [Sphinx](https://www.sphinx-doc.org/)
and does **not** require the `pypkatool`/`py27`/`pdbfixer` conda environments -
only a plain Python 3.10+ environment:

```bash
python -m venv .docs-venv && source .docs-venv/bin/activate   # or any Python 3.10+ env
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build
# open docs/_build/index.html
```

The API reference page uses `sphinx.ext.autodoc` on `pypkatool/core.py`'s
docstrings; `docs/conf.py` stubs out the `pkai` import-time check
(`_require_pkai()`) so this works without installing PyPKA/pKAI themselves.

## Repository layout

```
.
├── pypkatool/                 Installable package (CLI + pipeline)
│   ├── __init__.py            Exposes __version__ and main()
│   ├── core.py                All pipeline logic: PyPKA runner, tautomer->CHARMM
│   │                          mapping, PDB/RTF cross-validation, pKAI+
│   │                          cross-validation, report writers (see docstrings
│   │                          for full API detail)
│   └── data/                  Package data - bundled CHARMM36 protein topology and
│       │                      the fixstructure worker script
│       ├── top_all36_prot.rtf
│       └── fixstructure_worker.py   Runs inside the `pdbfixer` env (see _find_pdbfixer_python)
├── tests/
│   ├── test_pypkatool.py      Unit + regression + adversarial test suite
│   ├── data/                  Benchmark input PDBs (used for manual/example runs)
│   └── fixtures/               Frozen reference outputs (regression test fixtures,
│                               one directory per benchmark protein/pH)
├── examples/
│   └── denv2.pdb               Example input for a quick first run
├── docs/                      Sphinx documentation source (see "Building the docs")
│   ├── conf.py
│   ├── index.md                Landing page + toctree
│   ├── requirements.txt        Doc-build-only deps (independent of environment.yml)
│   └── ...                     installation/usage/force_field/tautomer_mapping/
│                               fixstructure/failure_modes/validation/api/citation/changelog
├── environment.yml            Main conda environment (Python 3.10 + PyPKA + pKAI)
├── environment-py27.yml       Python 2.7 helper environment (PyPKA internal dependency)
├── environment-pdbfixer.yml   PDBFixer/OpenMM helper environment (fixstructure command)
├── pyproject.toml             Package metadata + `pypkatool` console script
├── CITATION.cff               Citation metadata (GitHub "Cite this repository" button)
└── CHANGELOG.md                Notable changes, by version (Keep a Changelog format)
```

Note the two separate `data/` directories: `pypkatool/data/` ships *inside* the
installed package (the RTF topology file the pipeline needs at runtime),
while `tests/data/` is development-only input fixtures for the test suite
and is never imported by `pypkatool` itself.

## References

- PyPKA: Reis, P. B. P. S. et al. *J. Chem. Inf. Model.* 2020, 60, 4442-4448.
  [DOI: 10.1021/acs.jcim.0c00718](https://doi.org/10.1021/acs.jcim.0c00718)
- pKAI / pKAI+: Reis, P. B. P. S. et al. *J. Chem. Theory Comput.* 2022, 18, 3925-3935.
  [DOI: 10.1021/acs.jctc.2c00308](https://doi.org/10.1021/acs.jctc.2c00308)
- CHARMM36: Best, R. B.; Zhu, X.; Shim, J.; Lopes, P. E. M.; Mittal, J.;
  Feig, M.; MacKerell, A. D. Jr. *J. Chem. Theory Comput.* 2012, 8 (9),
  3257-3273. [DOI: 10.1021/ct300400x](https://doi.org/10.1021/ct300400x) -
  source of `top_all36_prot.rtf`'s RESI/PRES blocks for all protonation
  states used here.
- CHARMM36m: Huang, J.; Rauscher, S.; Nawrocki, G.; Ran, T.; Feig, M.;
  de Groot, B. L.; Grubmüller, H.; MacKerell, A. D. Jr. *Nat. Methods* 2017,
  14 (1), 71-73. [DOI: 10.1038/nmeth.4067](https://doi.org/10.1038/nmeth.4067) -
  the specific, more current CHARMM36 revision PyPKA's tautomer library
  (`CHARMM36m/sts/`) is built on.

## How to cite this repository

If you use `pypkatool` in published work, please cite it alongside the
methods it wraps (PyPKA, pKAI+, CHARMM36/CHARMM36m - see
["References"](#references) above), since those are what actually compute
the pKa values and define the protonation-state topology.

Citation metadata for `pypkatool` itself is kept in [`CITATION.cff`](CITATION.cff)
(GitHub reads this automatically and adds a "Cite this repository" button -
APA/BibTeX export - to the repo page). Manually:

```bibtex
@software{ropon_palacios_pypkatool,
  author  = {Ropón-Palacios, G.},
  title   = {pypkatool},
  version = {1.0.0},
  date    = {2026-07-08},
  url     = {https://github.com/groponp/PyPkaTool}
}
```

## Author

**Ropón-Palacios G.**
Department of Physics, UNESP.
[georcki.ropon@unesp.br](mailto:georcki.ropon@unesp.br)

## Disclaimer

This software is provided "as is", without warranty of any kind, express or
implied, including but not limited to the warranties of merchantability,
fitness for a particular purpose, and noninfringement - see the full text in
[LICENSE](LICENSE). pKa predictions and CHARMM protonation-state assignments
produced by this pipeline are computational estimates and must be checked
against experimental data and domain expertise before use in downstream
modeling; the authors assume no liability for outcomes derived from its use.

## License

MIT - see [LICENSE](LICENSE).
