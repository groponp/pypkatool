# Validation

Cross-validated against 12 benchmark proteins spanning pH 5 (flavivirus
envelope protein, endosomal pH sensing) and pH 7 (standard pKa benchmarks:
lysozyme, RNase A, barnase, BPTI, SNase, protein G, thioredoxin), matching
published experimental pKa values where available (e.g. lysozyme HIS15
pKa 6.33, experimental ~5.7) and reaching 96.7-100% sign agreement with the
independent pKAI+ predictor across ~1145 titratable sites (`SER`/`THR` are
not titrated under `ffID="CHARMM36m"` - see {doc}`force_field`). See
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
cross-validation machinery ({func}`pypkatool.core.reconcile_his_from_pdb`,
{func}`pypkatool.core.reconcile_non_his_from_pdb`, `_validate_label_chain`)
actually catches them. No PyPKA rerun is needed to run this suite.
