# Usage

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

## `run` options

| Flag | Default | Meaning |
|---|---|---|
| `--ph` | *(required)* | Target pH for the most-probable-protonation-state report |
| `--outdir` | `pypkatool_<stem>_pH<ph>/` next to the input PDB | Output directory |
| `--ncpus` | all detected CPUs | Parallel PB solves in PyPKA |
| `--epsin` | `15` | Protein interior dielectric constant (PyPKA's literature-optimized default; RMSE 0.82 / MAE 0.57 on the PKAD benchmark) |
| `--charmm-input` | off | Set PyPKA's `ffinput=CHARMM` (see {doc}`force_field`) |

## `reprocess` options

| Flag | Default | Meaning |
|---|---|---|
| `outdir` | *(required)* | A previous run's output directory (must contain the raw PyPKA titration data) |
| `--ph` | *(required)* | pH to regenerate reports at |
| `--pdb` | `<outdir>/<stem>.pdb` (stem taken from the protonated PDB filename) | Original input PDB, if it isn't at the default guessed location |
| `--epsin` | `15` | Only used for report metadata; the PB+MC step is not rerun |

## `fixstructure` options

| Flag | Default | Meaning |
|---|---|---|
| `pdb` | *(required)* | Input structure to repair |
| `--outdir` | input PDB's parent directory | Where to write `<stem>_fixed.pdb` |
