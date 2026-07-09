# Force field parameters used for the PB+MC calculation

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

## Where this is verified in the source

- `ffID` default and `SER`/`THR` gating: `pypka/config.py`,
  `Config.set_radii_charges_paths()`.
- `ffinput`/nomenclature call sites (all three, all gated on CHARMM-specific
  residue names): `pdbmender/utils.py::identify_tit_sites()`,
  `pypka/clean/cleaning.py` (two call sites).
- CHARMM-specific residue name set: `pdbmender/ffconverter.py::CHARMM_protomers`.

This was confirmed empirically as well as by source inspection: running all
12 benchmark proteins with `ffID="CHARMM36m"` alone versus
`ffID="CHARMM36m"` + `ffinput="CHARMM"` produced bit-identical pKa values,
tautomer populations, and CHARMM labels across all ~1145 titratable sites -
`ffinput` genuinely does nothing for standard PDB input.
