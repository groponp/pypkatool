# How the tautomer -> CHARMM label mapping works

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

See {func}`pypkatool.core._label` for the exact mapping table implementation,
and {func}`pypkatool.core.reconcile_his_from_pdb` /
{func}`pypkatool.core.reconcile_non_his_from_pdb` for the PDB/RTF
cross-validation step.
