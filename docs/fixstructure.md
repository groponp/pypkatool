# Repairing fragmented structures

PyPKA does not repair a structurally incomplete PDB - it requires every
present residue to have its expected atom set. A truncated terminus missing
its carboxyl oxygens, or a residue missing backbone/side-chain atoms, will
make the underlying preprocessing step fail (see {doc}`failure_modes`). For
inputs like this, run `pypkatool fixstructure` first to produce a repaired
PDB, then feed that into `pypkatool run`.

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
pypkatool run results/my_protein_fixed.pdb --pH 7.0
```

AlphaFold2/ColabFold predictions are heavy-atom-only (no hydrogens) and, for
the residue range you gave the model, have no internal gaps or missing
atoms by construction - `fixstructure` is normally unnecessary for them.
It is intended for experimentally determined structures with genuine
crystallographic disorder or truncated regions.

## Internal gaps require a reference sequence - `--pdbid`

Detecting an internal gap means comparing the chain's actual residues
against what *should* be there, and PDBFixer's `findMissingResidues()` gets
that reference sequence only from `SEQRES` records in the input PDB. A PDB
with no `SEQRES` (common for hand-edited or programmatically stripped test
files) gives PDBFixer nothing to compare against - confirmed directly:
`PDBFixer(filename=...).sequences` is `[]` on such a file, and
`findMissingResidues()` returns `{}` regardless of how many residues are
actually missing. **This means `fixstructure` silently reports 0 internal
gaps repaired on a `SEQRES`-less PDB, even when real gaps are present** - it
now also prints an explicit `WARNING` in that situation so this isn't
mistaken for "nothing was missing".

If the structure has a deposited RCSB entry, pass its 4-character code with
`--pdbid`: this fetches the deposited `SEQRES` from RCSB and uses it as the
reference sequence for gap detection, while the atoms/coordinates being
repaired still come entirely from your local file.

```bash
pypkatool fixstructure my_fragment.pdb --outdir results/ --pdbid 7A3S
```

Requires network access, and can take on the order of a minute for a large
entry (PDBFixer fetches the full deposited structure to derive the
sequence, not just a lightweight SEQRES query).

## Implementation

`fixstructure` runs `pypkatool/data/fixstructure_worker.py` as a subprocess
inside the separate `pdbfixer` conda environment (see
{func}`pypkatool.core._find_pdbfixer_python` and
{func}`pypkatool.core.fix_structure`). The internal/terminal-gap distinction
is implemented by inspecting the `(chain_index, position)` keys PDBFixer's
`findMissingResidues()` returns: a key with `position == 0` is before the
first present residue in its chain (N-terminal gap), and
`position == len(present_residues_in_chain)` is after the last one
(C-terminal gap) - both are filtered out before `addMissingAtoms()` is
called, so only strictly internal gaps (`0 < position < len(...)`) get
rebuilt.
