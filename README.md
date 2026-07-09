# pypkatool

`pypkatool` is a command-line pipeline that takes a protein structure and a
target pH, runs [PyPKA](https://doi.org/10.1021/acs.jcim.0c00718)
(Poisson-Boltzmann + Monte Carlo pKa calculation), cross-validates every
predicted pKa against the independent machine-learning predictor
[pKAI+](https://doi.org/10.1021/acs.jctc.2c00308), and writes a
CHARMM-GUI-ready protonation-state table (which `RESI`/`PRES` to pick for
each titratable residue in PDB Reader).

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

This project uses **two separate, independent conda environments** - one is
not created "inside" the other, and they serve different purposes:

| Environment | Created from | Do you activate it? |
|---|---|---|
| `pypkatool` | `environment.yml` | **Yes - every time** you run the tool: `conda activate pypkatool` |
| `py27` | `environment-py27.yml` | **No, never.** `pypkatool` finds it on disk automatically at runtime. |

Why `py27` exists at all: PyPKA's compiled Poisson-Boltzmann backend
(DelPhi4py) needs `numpy<2` and `libgfortran4=7.5.0` in the main environment
(the compiled binary was built against NumPy's 1.x C API and links
`GFORTRAN_7`, which `libgfortran.so.5` / GFortran 10+ doesn't export) - that
part lives in `pypkatool`. Separately, PyPKA also shells out to a bare
`python2.7` interpreter internally (via `pdbmender`'s vendored `pdb2pqr.py`),
which cannot coexist with Python 3 in the same environment, hence the second,
minimal, interpreter-only environment.

### One-time setup

```bash
# 1. Create BOTH environments. Order doesn't matter, and it doesn't matter
#    whether any environment is currently active: `conda env create` always
#    builds a new, independent environment from scratch.
conda env create -f environment.yml
conda env create -f environment-py27.yml

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

## Validation

Cross-validated against 12 benchmark proteins spanning pH 5 (flavivirus
envelope protein, endosomal pH sensing) and pH 7 (standard pKa benchmarks:
lysozyme, RNase A, barnase, BPTI, SNase, protein G, thioredoxin), matching
published experimental pKa values where available (e.g. lysozyme HIS15
pKa 5.88, experimental ~5.7) and reaching 96.9-100% sign agreement with the
independent pKAI+ predictor across ~1684 titratable sites. See
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

## Repository layout

```
.
├── pypkatool/                 Installable package (CLI + pipeline)
│   ├── __init__.py            Exposes __version__ and main()
│   ├── core.py                All pipeline logic: PyPKA runner, tautomer->CHARMM
│   │                          mapping, PDB/RTF cross-validation, pKAI+
│   │                          cross-validation, report writers (see docstrings
│   │                          for full API detail)
│   └── data/                  Package data - bundled CHARMM36 protein topology
│       └── top_all36_prot.rtf
├── tests/
│   ├── test_pypkatool.py      Unit + regression + adversarial test suite
│   ├── data/                  Benchmark input PDBs (used for manual/example runs)
│   └── fixtures/               Frozen reference outputs (regression test fixtures,
│                               one directory per benchmark protein/pH)
├── examples/
│   └── denv2.pdb               Example input for a quick first run
├── environment.yml            Main conda environment (Python 3.10 + PyPKA + pKAI)
├── environment-py27.yml       Python 2.7 helper environment (PyPKA internal dependency)
└── pyproject.toml             Package metadata + `pypkatool` console script
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
