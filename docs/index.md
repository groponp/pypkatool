# pypkatool

`pypkatool` is a command-line pipeline that takes a protein structure and a
target pH, runs [PyPKA](https://doi.org/10.1021/acs.jcim.0c00718)
(Poisson-Boltzmann + Monte Carlo pKa calculation), cross-validates every
predicted pKa against the independent machine-learning predictor
[pKAI+](https://doi.org/10.1021/acs.jctc.2c00308), and writes a
CHARMM-GUI-ready protonation-state table (which `RESI`/`PRES` to pick for
each titratable residue in PDB Reader).

For a quick start, see the project [README](https://github.com/groponp/PyPkaTool#readme).
This site is the full manual and API reference.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
usage
force_field
tautomer_mapping
fixstructure
failure_modes
validation
api
citation
changelog
```

## Output per run

| File | Contents |
|---|---|
| `protonation_inputs.dat` | Residues that need a CHARMM-GUI action (RESI choice or PRES patch), plain text |
| `protonation_inputs.json` | Every titratable site, full detail, machine-readable |
| `detail.json` | Per-site tautomer population breakdown + the exact RTF blocks used |
| `crossvalidation_report.dat` | PyPKA vs. pKAI+ pKa and agreement, per site |
| `protocol.json` | Exact PyPKA parameters and package versions used (reproducibility) |
| `output_pypka/` | Raw PyPKA output: protonated PDB, full pH 0-14 titration curve |

## Indices

- {ref}`genindex`
- {ref}`modindex`
