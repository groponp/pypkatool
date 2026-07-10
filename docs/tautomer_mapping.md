# How the tautomer -> CHARMM label mapping works

PyPKA reports, per titratable site, the most probable *tautomer* out of `N`
regular tautomers plus one *reference* tautomer (index `N+1`). "Reference" is
PyPKA's own internal term for this `N+1`-th slot (`titsite.py::getTautomers()`
docstring: "all tautomers instances except the tautomers of reference") - it
is not a CHARMM or RTF concept, and the RTF file itself never marks any RESI
or PRES block as a "reference" of anything. `pypkatool` is what connects the
two vocabularies: it reads which of PyPKA's tautomer indices won at the
target pH, and `_label()` maps that index to a concrete CHARMM label (e.g.
the reference tautomer maps to `HSP` for HIS - but `HSP` in the RTF is just
an ordinary `RESI` block like any other, with no special "reference" marker).

Which physical state the reference tautomer represents depends on site
polarity - this is easy to get backwards, so it's worth stating precisely
(verified against the installed PyPKA source,
`titsite.py::Titsite.getRefProtState()`, and against the signed atomic
partial charges in the CHARMM36 `.st` tautomer files):

* **Cationic sites** (`HIS`, `LYS`, `NTR`): the reference tautomer is the
  **protonated**, positively charged state. For `LYS` that is the CHARMM
  *default* `RESI` - no patch needed. `NTR`'s reference state maps to
  `NTER`, which - unlike `LYS` - is a `PRES` in the RTF, not a `RESI` (a
  chain's N-terminus is always applied as a patch on its first residue, it
  is never a standalone residue); it is, however, the *default* terminal
  patch CHARMM-GUI's PDB Reader applies automatically, so no explicit action
  is needed there either. The regular tautomers are the neutral forms and
  need an explicit patch (`LSN`, `NNEU`).
* **Anionic sites** (`ASP`, `GLU`, `CTR`, `CYS`, `TYR`, `SER`): the reference
  tautomer is the **deprotonated**, negatively charged state.
  * For `ASP`/`GLU`, that deprotonated/charged form *is* the CHARMM default
    `RESI` - no patch needed; the protonated form needs a patch (`ASPP`,
    `GLUP`). `CTR`'s reference state maps to `CTER`, which - like `NTER` -
    is a `PRES`, not a `RESI`, but is the default C-terminal patch applied
    automatically, so no explicit action is needed there either.
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

See {func}`pypkatool.core._label` for the exact mapping table implementation,
and {func}`pypkatool.core.reconcile_his_from_pdb` /
{func}`pypkatool.core.reconcile_non_his_from_pdb` for the PDB/RTF
cross-validation step.
